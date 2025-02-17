source("./custom/myingarch.R")
source("./custom/ingarch.string.R")
source("./custom/newmodel.R")
library(tscount)
library(forecast)

auto.ingarch <- function(y, 
                           max.p = 5,           # Maximum AR order
                           max.q = 5,           # Maximum MA order
                           start.p = 2,         # Starting AR order
                           start.q = 2,         # Starting MA order
                           distribution = c("poisson", "nbinom"),  # Distribution family
                           link = c("log", "identity"),        # Link function
                           xreg = NULL, # External covariates
                           ic = c("aicc", "aic", "bic", "qic"),  # Information criterion
                           stepwise = TRUE,     # Whether to use stepwise selection
                           nmodels = 94,        # Maximum number of models to try
                           trace = FALSE,       # Whether to print model search progress
                           method = NULL,       # Estimation method
                           parallel = FALSE,    # Whether to use parallel processing
                           num.cores = 2,       # Number of cores for parallel processing
                           x = y,               # Data (same as y)
                           ...) {
  
  # Only non-stepwise parallel implemented so far
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
  
  # Trim initial missing values
  # x <- subset(x, start=firstnonmiss)
  
  if(!is.null(xreg)) {
    # Check if xreg is numeric
    if(!is.numeric(xreg))
      stop("xreg should be a numeric matrix or a numeric vector")
    
    # Convert to matrix format
    xreg <- as.matrix(xreg)
    
    # Trim to match the time series
    xreg <- xreg[firstnonmiss:NROW(xreg),,drop=FALSE]
    
    # Check dimensions match
    if(NROW(xreg) != length(x))
      stop("Number of rows in xreg must match length of time series")
    
    # Check for missing values
    if(any(is.na(xreg)))
      stop("Missing values in external regressors are not allowed")
    
    # Add column names if not present
    if(is.null(colnames(xreg))) {
      colnames(xreg) <- paste0("xreg", 1:NCOL(xreg))
    }
  }
  
  # Check for constant data
  if (all(x == x[1], na.rm = TRUE)) {
    if(all(is.na(x)))
      stop("All data are missing")
    
    # Get the constant value
    const_value <- x[1]
    
    # Create basic model structure following tscount conventions
    fit <- list(
      ts = orig.x,                        # Original time series
      model = list(
        past_obs = 0,                     # No past observations (p)
        past_mean = 0,                    # No past means (q)
        external = FALSE                  # No external regressors
      ),
      distribution = match.arg(distribution),
      link = match.arg(link),
      final_estimates = const_value,      # Parameter estimate
      fitted.values = rep(const_value, length(x)),
      call = match.call()
    )
    
    # Add distribution-specific parameters
    if(fit$distribution == "nbinom") {
      fit$size <- Inf  # For constant data, equivalent to Poisson
    }
    
    # Use tscount's class
    class(fit) <- "tsglm"
    return(fit)
  }
  
  ic <- match.arg(ic)
  
  max.p <- min(max.p, floor(serieslength / 3))
  max.q <- min(max.q, floor(serieslength / 3))
  
  # Use AIC if npar <= 3
  # AICc won't work for tiny samples.
  if (serieslength <= 3L) {
    ic <- "aic"
  }
  
  # Fixed differencing orders
  d <- D <- 0
  
  # Later we can set this in a way users can specify? Does it even make sense?
  constant <- TRUE
  
  # Need to change this to do a grid search for INGARCH models
  # For now we will ignore non-stepwise search
  # if (!stepwise) {
  #  bestfit <- search.arima(
  #    x, d, D, max.p, max.q, max.P, max.Q, max.order, stationary,
  #    ic, trace, approximation, method = method, xreg = xreg, offset = offset,
  #    allowdrift = allowdrift, allowmean = allowmean,
  #    parallel = parallel, num.cores = num.cores, ...
  #  )
  #  bestfit$call <- match.call()
  #  bestfit$call$x <- data.frame(x = x)
  #  bestfit$lambda <- lambda
  #  bestfit$x <- orig.x
  #  bestfit$series <- series
  #  bestfit$fitted <- forecast:::fitted.Arima(bestfit)
  #  if (trace) {
  #    cat("\n\n Best model:", forecast:::arima.string(bestfit, padding = TRUE), "\n\n")
  #  }
  #  return(bestfit)
  #}
  
  # Starting model
  if (length(x) < 10L) {
    start.p <- min(start.p, 1L)
    start.q <- min(start.q, 1L)
  }
  
  p <- start.p <- min(start.p, max.p)
  q <- start.q <- min(start.q, max.q)
  
  results <- matrix(NA, nrow = nmodels, ncol = 4)
  
  bestfit <- myingarch(x, order = c(p, q), constant = FALSE, ic = ic, 
                       trace = trace, xreg = xreg, 
                       distr = distribution, link = link, ...)
  
  results[1, ] <- c(p, q, constant, bestfit$ic)
  
  # Null model
  # Before we were passing constant equal true
  fit <- myingarch(x, order = c(0, 0), constant = FALSE, ic = ic, 
                   trace = trace, xreg = xreg,
                   distr = distribution, link = link, ...)
  
  results[2, ] <- c(0, 0, constant, fit$ic)
  
  if (fit$ic < bestfit$ic) {
    bestfit <- fit
    p <- q <- 0
  }
  
  k <- 2
  
  # Basic model with only past observations (p)
  if (max.p > 0) {
    fit <- myingarch(x, order = c(1, 0), constant = FALSE, ic = ic,
                     trace = trace, xreg = xreg,
                     distr = distribution, link = link, ...)
    results[k+1, ] <- c(max.p > 0, 0, constant, fit$ic)
    
    if (fit$ic < bestfit$ic) {
      bestfit <- fit
      p <- (max.p > 0)
      q <- 0
    }
    
    k <- k + 1
  }
  
  # Basic model with only past means (q)
  if (max.q > 0) {
    fit <- myingarch(x, order = c(0, 1), constant = FALSE, ic = ic,
                     trace = trace, xreg = xreg,
                     distr = distribution, link = link, ...)
    results[k+1, ] <- c(0, max.q > 0, constant, fit$ic)
    
    if (fit$ic < bestfit$ic) {
      bestfit <- fit
      p <- 0
      q <- (max.q > 0)
    }
    
    k <- k + 1
  }
  
  # Null model with no constant (This part is the same as the first model but with cosntant false)
  # Now we are creating a model (1,1)
  # Code below was inside the if constant
  # if (constant) {
  # }
  
  fit <- myingarch(x, order = c(1, 1), constant = FALSE, ic = ic,
                   trace = trace, xreg = xreg,
                   distr = distribution, link = link, ...)
  
  results[k+1, ] <- c(0, 0, 0, fit$ic)
  
  if (fit$ic < bestfit$ic) {
    bestfit <- fit
    p <- q <- 0
    # constant <- FALSE
  }
  
  k <- k + 1
  
  startk <- 0
  
  cat("\nFitting through step-wise now...\n")
  
  while (startk < k && k < nmodels) {
    startk <- k
    
    # Try decreasing p
    if (p > 0 && newmodel(p - 1, d, q, D, constant, results[1:k, ])) {
      k <- k + 1; if(k>nmodels) next
      fit <- myingarch(x, order = c(p - 1, q), constant = constant, ic = ic,
                       trace = trace, xreg = xreg,
                       distr = distribution, link = link, ...)
      results[k, ] <- c(p - 1, q, constant, fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        p <- (p - 1)
        next
      }
    }
    
    # Try decreasing q
    if (q > 0 && newmodel(p, d, q - 1, D, constant, results[1:k, ])) {
      k <- k + 1; if(k>nmodels) next
      fit <- myingarch(x, order = c(p, q - 1), constant = constant, ic = ic,
                       trace = trace, xreg = xreg,
                       distr = distribution, link = link, ...)
      results[k, ] <- c(p, q - 1, constant, fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- (q - 1)
        next
      }
    }
    
    # Try increasing p
    if (p < max.p && newmodel(p + 1, d, q, D, constant, results[1:k, ])) {
      k <- k + 1; if(k>nmodels) next
      fit <- myingarch(x, order = c(p + 1, q), constant = constant, ic = ic,
                       trace = trace, xreg = xreg,
                       distr = distribution, link = link, ...)
      results[k, ] <- c(p + 1, q, constant, fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        p <- (p + 1)
        next
      }
    }
    
    # Try increasing q
    if (q < max.q && newmodel(p, d, q + 1, D, constant, results[1:k, ])) {
      k <- k + 1; if(k>nmodels) next
      fit <- myingarch(x, order = c(p, q + 1), constant = constant, ic = ic,
                       trace = trace, xreg = xreg,
                       distr = distribution, link = link, ...)
      results[k, ] <- c(p, q + 1, constant, fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- (q + 1)
        next
      }
    }
    
    # Try decreasing both p and q
    if (q > 0 && p > 0 && newmodel(p - 1, d, q - 1, D, constant, results[1:k, ])) {
      k <- k + 1; if(k>nmodels) next
      fit <- myingarch(x, order = c(p - 1, q - 1), constant = constant, ic = ic,
                       trace = trace, xreg = xreg,
                       distr = distribution, link = link, ...)
      results[k, ] <- c(p - 1, q - 1, constant, fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- (q - 1)
        p <- (p - 1)
        next
      }
    }
    
    # Try decreasing p and increasing q
    if (q < max.q && p > 0 && newmodel(p - 1, d, q + 1, D, constant, results[1:k, ])) {
      k <- k + 1; if(k>nmodels) next
      fit <- myingarch(x, order = c(p - 1, q + 1), constant = constant, ic = ic,
                       trace = trace, xreg = xreg,
                       distr = distribution, link = link, ...)
      results[k, ] <- c(p - 1, q + 1, constant, fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- (q + 1)
        p <- (p - 1)
        next
      }
    }
    
    # Try increasing p and decreasing q
    if (q > 0 && p < max.p && newmodel(p + 1, D, q - 1, D, constant, results[1:k, ])) {
      k <- k + 1; if(k>nmodels) next
      fit <- myingarch(x, order = c(p + 1, q - 1), constant = constant, ic = ic,
                       trace = trace, xreg = xreg,
                       distr = distribution, link = link, ...)
      results[k, ] <- c(p + 1, q - 1, constant, fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- (q - 1)
        p <- (p + 1)
        next
      }
    }
    
    # Try increasing both p and q
    if (q < max.q && p < max.p && newmodel(p + 1, d, q + 1, D, constant, results[1:k, ])) {
      k <- k + 1; if(k>nmodels) next
      fit <- myingarch(x, order = c(p + 1, q + 1), constant = constant, ic = ic,
                       trace = trace, xreg = xreg,
                       distr = distribution, link = link, ...)
      results[k, ] <- c(p + 1, q + 1, constant, fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- (q + 1)
        p <- (p + 1)
        next
      }
    }
    
    # Try toggling constant term
    if (!is.null(constant)) {
      if (newmodel(p, d, q, D, !constant, results[1:k, ])) {
        k <- k + 1; if(k>nmodels) next
        fit <- myingarch(x, order = c(p, q), constant = !constant, ic = ic,
                         trace = trace, xreg = xreg,
                         distr = distribution, link = link, ...)
        results[k, ] <- c(p, q, !constant, fit$ic)
        if (fit$ic < bestfit$ic) {
          bestfit <- fit
          constant <- !constant
        }
      }
    }
  }
  
  if(k > nmodels){
    warning(sprintf("Stepwise search was stopped early due to reaching the model number limit: `nmodels = %i`", nmodels))
  }
  
  
  # Nothing fitted
  if (bestfit$ic == Inf) {
    if (trace) {
      cat("\n")
    }
    stop("No suitable INGARCH model found")
  }
  
  # Return best fit
  bestfit$x <- orig.x
  bestfit$series <- series
  bestfit$ic <- NULL
  bestfit$call <- match.call()
  bestfit$call$x <- data.frame(x = x)
  
  if (trace) {
    cat("\n\n Best model:", ingarch.string(bestfit, padding = TRUE), "\n\n")
  }
  
  return(bestfit)
}