source("./custom/newmodel.R")
source("./custom/ingarch.string.R")
source("./custom/search.ingarch.R")

auto.ingarch <- function(y, 
                         max.p = 7,           # Maximum AR order
                         max.q = 7,           # Maximum MA order
                         start.p = 2,         # Starting AR order
                         start.q = 2,         # Starting MA order
                         distribution = c("poisson", "nbinom"),  # Distribution family
                         link = c("log", "identity"),        # Link function
                         xreg = NULL,         # External covariates
                         ic = c("aicc", "aic", "bic", "qic"),  # Information criterion
                         stepwise = TRUE,     # Whether to use stepwise selection
                         nmodels = 94,        # Maximum number of models to try
                         trace = FALSE,       # Whether to print model search progress
                         show_warnings = FALSE, # Whether to show warnings at the end
                         method = NULL,       # Estimation method
                         parallel = FALSE,    # Whether to use parallel processing
                         num.cores = 2,       # Number of cores for parallel processing
                         x = y,               # Data (same as y)
                         ...) {
  
  # Initial validation
  if (stepwise && parallel) {
    warning("Parallel computing is only implemented when stepwise=FALSE, the model will be fit in serial.")
    parallel <- FALSE
  }
  
  if (trace && parallel) {
    message("Tracing model searching in parallel is not supported.")
    trace <- FALSE
  }
  
  series <- deparse(substitute(y))
  x <- as.ts(x)
  
  if (!is.null(start.p) && start.p > max.p) {
    warning("start.p cannot be greater than max.p, setting start.p = max.p")
    start.p <- max.p
  }
  
  # Check for univariate series
  if (NCOL(x) > 1) {
    stop("auto.ingarch can only handle univariate time series")
  }
  
  # Check for count data
  if (!all(x[!is.na(x)] >= 0) || !all(x[!is.na(x)] == floor(x[!is.na(x)]))) {
    stop("INGARCH models require non-negative integer (count) data")
  }
  
  # Trim leading NAs and find length of non-missing data
  orig.x <- x
  missing <- is.na(x)
  firstnonmiss <- head(which(!missing), 1)
  lastnonmiss <- tail(which(!missing), 1)
  serieslength <- sum(!missing[firstnonmiss:lastnonmiss])
  
  # Handle xreg validation
  if(!is.null(xreg)) {
    if(!is.numeric(xreg))
      stop("xreg should be a numeric matrix or a numeric vector")
    
    xreg <- as.matrix(xreg)
    xreg <- xreg[firstnonmiss:NROW(xreg),,drop=FALSE]
    
    if(NROW(xreg) != length(x))
      stop("Number of rows in xreg must match length of time series")
    
    if(any(is.na(xreg)))
      stop("Missing values in external regressors are not allowed")
    
    if(is.null(colnames(xreg))) {
      colnames(xreg) <- paste0("xreg", 1:NCOL(xreg))
    }
    
    # Add check for non-negative xreg when using identity link
    link <- match.arg(link)
    if(link == "identity" && any(xreg < 0)) {
      stop("When using 'identity' link, all values in xreg must be non-negative")
    }
  }
  
  # Check for constant data
  if (all(x == x[1], na.rm = TRUE)) {
    if(all(is.na(x)))
      stop("All data are missing")
    
    const_value <- x[1]
    fit <- list(
      ts = orig.x,
      model = list(
        past_obs = 0,
        past_mean = 0,
        external = FALSE
      ),
      distribution = match.arg(distribution),
      link = match.arg(link),
      final_estimates = const_value,
      fitted.values = rep(const_value, length(x)),
      call = match.call()
    )
    
    if(fit$distribution == "nbinom") {
      fit$size <- Inf
    }
    
    class(fit) <- "tsglm"
    return(fit)
  }
  
  ic <- match.arg(ic)
  
  max.p <- min(max.p, floor(serieslength / 3))
  max.q <- min(max.q, floor(serieslength / 3))
  
  if (serieslength <= 3L) {
    ic <- "aic"
  }
  
  if (!stepwise) {
    # Use brute-force search instead of stepwise
    bestfit <- search.ingarch(
      x = x,
      max.p = max.p, 
      max.q = max.q,
      distribution = distribution,  # INGARCH-specific parameter
      link = link,                 # INGARCH-specific parameter
      xreg = xreg,
      ic = ic, 
      trace = trace,
      show_warnings = show_warnings, # Added to match warning system
      parallel = parallel, 
      num.cores = num.cores, 
      ...
    )
    
    # Set appropriate properties on bestfit
    bestfit$x <- orig.x
    bestfit$series <- series
    bestfit$ic <- NULL
    bestfit$call <- match.call()
    
    if (trace) {
      cat("\n\n Best model:", ingarch.string(bestfit, padding = TRUE), "\n\n")
    }
    
    if (!is.null(bestfit$coefficients)) {
      # Add minimum threshold for parameters to prevent convergence to zero
      bestfit$coefficients[bestfit$coefficients < 1e-4] <- 1e-4
    }
    
    return(bestfit)
  }
  
  # Starting model
  if (length(x) < 10L) {
    start.p <- min(start.p, 1L)
    start.q <- min(start.q, 1L)
  }
  
  p <- start.p <- min(start.p, max.p)
  q <- start.q <- min(start.q, max.q)
  
  # Initialize model warnings collection
  model_warnings <- list()
  
  # Initialize results matrix with NA values
  results <- matrix(NA, nrow = nmodels, ncol = 3)
  colnames(results) <- c("p", "q", "ic")
  
  # Initial model with intercept - capture warnings
  warnings_captured <- character(0)
  bestfit <- withCallingHandlers({
    myingarch(x, order = c(p, q), ic = ic, 
              trace = trace, xreg = xreg, 
              distr = distribution, link = link, ...)
  }, warning = function(w) {
    warnings_captured <<- c(warnings_captured, w$message)
    invokeRestart("muffleWarning")
  })
  
  # Store results and warnings
  results[1, ] <- c(p, q, bestfit$ic)
  model_warnings[[paste0("p", p, "_q", q)]] <- warnings_captured
  
  k <- 1
  
  # Fit (0,0) model - capture warnings
  k <- k + 1
  warnings_captured <- character(0)
  fit <- withCallingHandlers({
    myingarch(x, order = c(0, 0), ic = ic,
              trace = trace, xreg = xreg,
              distr = distribution, link = link, ...)
  }, warning = function(w) {
    warnings_captured <<- c(warnings_captured, w$message)
    invokeRestart("muffleWarning")
  })
  
  # Store results and warnings
  results[k, ] <- c(0, 0, fit$ic)
  model_warnings[[paste0("p", 0, "_q", 0)]] <- warnings_captured
  
  if (fit$ic < bestfit$ic) {
    bestfit <- fit
    p <- 0
    q <- 0
  }
  
  # Basic model with only past observations (p) - capture warnings
  if (max.p > 0 && (p != 1 || q != 0)) {
    k <- k + 1
    
    warnings_captured <- character(0)
    fit <- withCallingHandlers({
      myingarch(x, order = c(1, 0), ic = ic,
                trace = trace, xreg = xreg,
                distr = distribution, link = link, ...)
    }, warning = function(w) {
      warnings_captured <<- c(warnings_captured, w$message)
      invokeRestart("muffleWarning")
    })
    
    # Store results and warnings
    results[k, ] <- c(1, 0, fit$ic)
    model_warnings[[paste0("p", 1, "_q", 0)]] <- warnings_captured
    
    if (fit$ic < bestfit$ic) {
      bestfit <- fit
      p <- 1
      q <- 0
    }
  }
  
  # Basic model with only past means (q) - capture warnings
  if (max.q > 0 && (p != 0 || q != 1)) {
    k <- k + 1
    
    warnings_captured <- character(0)
    fit <- withCallingHandlers({
      myingarch(x, order = c(0, 1), ic = ic,
                trace = trace, xreg = xreg,
                distr = distribution, link = link, ...)
    }, warning = function(w) {
      warnings_captured <<- c(warnings_captured, w$message)
      invokeRestart("muffleWarning")
    })
    
    # Store results and warnings
    results[k, ] <- c(0, 1, fit$ic)
    model_warnings[[paste0("p", 0, "_q", 1)]] <- warnings_captured
    
    if (fit$ic < bestfit$ic) {
      bestfit <- fit
      p <- 0
      q <- 1
    }
  }
  
  startk <- 0
  
  if (trace) {
    cat("\nFitting through step-wise now...\n")
  }
  
  # Stepwise search loop
  while (startk < k && k < nmodels) {
    startk <- k
    
    # Try decreasing p - capture warnings
    if (p > 0 && newmodel(p - 1, q, results[1:k, ])) {
      k <- k + 1; 
      
      if(k > nmodels)
        break
      
      warnings_captured <- character(0)
      fit <- withCallingHandlers({
        myingarch(x, order = c(p - 1, q), ic = ic,
                  trace = trace, xreg = xreg,
                  distr = distribution, link = link, ...)
      }, warning = function(w) {
        warnings_captured <<- c(warnings_captured, w$message)
        invokeRestart("muffleWarning")
      })
      
      # Store results and warnings
      results[k, ] <- c(p - 1, q, fit$ic)
      model_warnings[[paste0("p", p-1, "_q", q)]] <- warnings_captured
      
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        p <- p - 1
        next
      }
    }
    
    # Try increasing q - capture warnings
    if (q < max.q && newmodel(p, q + 1, results[1:k, ])) {
      k <- k + 1;
      
      if(k > nmodels)
        break
      
      warnings_captured <- character(0)
      fit <- withCallingHandlers({
        myingarch(x, order = c(p, q + 1), ic = ic,
                  trace = trace, xreg = xreg,
                  distr = distribution, link = link, ...)
      }, warning = function(w) {
        warnings_captured <<- c(warnings_captured, w$message)
        invokeRestart("muffleWarning")
      })
      
      # Store results and warnings
      results[k, ] <- c(p, q + 1, fit$ic)
      model_warnings[[paste0("p", p, "_q", q+1)]] <- warnings_captured
      
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- q + 1
        next
      }
    }
    
    # Try decreasing both p and q - capture warnings
    if (q > 0 && p > 0 && newmodel(p - 1, q - 1, results[1:k, ])) {
      k <- k + 1;
      
      if(k > nmodels)
        break
      
      warnings_captured <- character(0)
      fit <- withCallingHandlers({
        myingarch(x, order = c(p - 1, q - 1), ic = ic,
                  trace = trace, xreg = xreg,
                  distr = distribution, link = link, ...)
      }, warning = function(w) {
        warnings_captured <<- c(warnings_captured, w$message)
        invokeRestart("muffleWarning")
      })
      
      # Store results and warnings
      results[k, ] <- c(p - 1, q - 1, fit$ic)
      model_warnings[[paste0("p", p-1, "_q", q-1)]] <- warnings_captured
      
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- q - 1
        p <- p - 1
        next
      }
    }
    
    # Try decreasing p and increasing q - capture warnings
    if (q < max.q && p > 0 && newmodel(p - 1, q + 1, results[1:k, ])) {
      k <- k + 1;
      
      if(k > nmodels)
        break
      
      warnings_captured <- character(0)
      fit <- withCallingHandlers({
        myingarch(x, order = c(p - 1, q + 1), ic = ic,
                  trace = trace, xreg = xreg,
                  distr = distribution, link = link, ...)
      }, warning = function(w) {
        warnings_captured <<- c(warnings_captured, w$message)
        invokeRestart("muffleWarning")
      })
      
      # Store results and warnings
      results[k, ] <- c(p - 1, q + 1, fit$ic)
      model_warnings[[paste0("p", p-1, "_q", q+1)]] <- warnings_captured
      
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- q + 1
        p <- p - 1
        next
      }
    }
    
    # Try increasing p and decreasing q - capture warnings
    if (q > 0 && p < max.p && newmodel(p + 1, q - 1, results[1:k, ])) {
      k <- k + 1;
      
      if(k > nmodels)
        break
      
      warnings_captured <- character(0)
      fit <- withCallingHandlers({
        myingarch(x, order = c(p + 1, q - 1), ic = ic,
                  trace = trace, xreg = xreg,
                  distr = distribution, link = link, ...)
      }, warning = function(w) {
        warnings_captured <<- c(warnings_captured, w$message)
        invokeRestart("muffleWarning")
      })
      
      # Store results and warnings
      results[k, ] <- c(p + 1, q - 1, fit$ic)
      model_warnings[[paste0("p", p+1, "_q", q-1)]] <- warnings_captured
      
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- q - 1
        p <- p + 1
        next
      }
    }
    
    # Try increasing both p and q - capture warnings
    if (q < max.q && p < max.p && newmodel(p + 1, q + 1, results[1:k, ])) {
      k <- k + 1;
      
      if(k > nmodels)
        break
      
      warnings_captured <- character(0)
      fit <- withCallingHandlers({
        myingarch(x, order = c(p + 1, q + 1), ic = ic,
                  trace = trace, xreg = xreg,
                  distr = distribution, link = link, ...)
      }, warning = function(w) {
        warnings_captured <<- c(warnings_captured, w$message)
        invokeRestart("muffleWarning")
      })
      
      # Store results and warnings
      results[k, ] <- c(p + 1, q + 1, fit$ic)
      model_warnings[[paste0("p", p+1, "_q", q+1)]] <- warnings_captured
      
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- q + 1
        p <- p + 1
        next
      }
    }
  }
  
  if(k >= nmodels){
    warning(sprintf("Stepwise search was stopped early due to reaching the model number limit: `nmodels = %i`", nmodels))
  }
  
  if (bestfit$ic == Inf) {
    if (trace) {
      cat("\n")
    }
    stop("No suitable INGARCH model found")
  }
  
  # Make sure these lines are still present before returning
  bestfit$x <- orig.x
  bestfit$series <- series
  bestfit$ic <- NULL
  bestfit$call <- call("auto.ingarch", y = as.name("ts_data"))
  
  # Store model search results and warnings
  bestfit$results <- results[1:k, , drop = FALSE]
  bestfit$model_warnings <- model_warnings
  
  # Find warnings specific to the best model
  best_p <- if(is.null(bestfit$model$past_obs)) 0 else length(bestfit$model$past_obs)
  best_q <- if(is.null(bestfit$model$past_mean)) 0 else length(bestfit$model$past_mean)
  bestfit$best_model_warnings <- model_warnings[[paste0("p", best_p, "_q", best_q)]]
  
  if (trace) {
    cat("\n\n Best model:", ingarch.string(bestfit, padding = TRUE), "\n\n")
  }
  
  # Display warnings for all models if requested
  if (show_warnings) {
    cat("\nWarnings for all models tested:\n")
    cat("------------------------------\n")
    
    # Get non-NA rows
    valid_results <- results[1:k, , drop = FALSE]
    valid_results <- valid_results[!is.na(valid_results[,1]), , drop = FALSE]
    
    # Order by IC value
    ordered_idx <- order(valid_results[,3])
    ordered_results <- valid_results[ordered_idx, , drop = FALSE]
    
    # Identify best model
    best_p <- if(is.null(bestfit$model$past_obs)) 0 else length(bestfit$model$past_obs)
    best_q <- if(is.null(bestfit$model$past_mean)) 0 else length(bestfit$model$past_mean)
    
    for (i in 1:nrow(ordered_results)) {
      p_val <- ordered_results[i, 1]
      q_val <- ordered_results[i, 2]
      ic_val <- ordered_results[i, 3]
      
      # Mark best model with a star
      is_best <- (p_val == best_p && q_val == best_q)
      model_label <- if(is_best) "★ BEST MODEL" else ""
      
      cat(sprintf("INGARCH(%d,%d) with IC=%.4f %s\n", 
                  p_val, q_val, ic_val, model_label))
      
      warnings_key <- paste0("p", p_val, "_q", q_val)
      if (warnings_key %in% names(model_warnings) && length(model_warnings[[warnings_key]]) > 0) {
        for (w in model_warnings[[warnings_key]]) {
          cat("  - ", w, "\n", sep="")
        }
      } else {
        cat("  No warnings\n")
      }
      cat("\n")
    }
  }
  
  if (!is.null(bestfit$coefficients)) {
    # Add minimum threshold for parameters to prevent convergence to zero
    bestfit$coefficients[bestfit$coefficients < 1e-4] <- 1e-4
  }
  
  return(bestfit)
}