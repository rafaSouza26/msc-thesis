source("./custom/newmodel.R")
source("./custom/ingarch.string.R")

auto.ingarch <- function(y, 
                         max.p = 5,           # Maximum AR order
                         max.q = 5,           # Maximum MA order
                         start.p = 2,         # Starting AR order
                         start.q = 2,         # Starting MA order
                         distribution = c("poisson", "nbinom"),  # Distribution family
                         link = c("log", "identity"),        # Link function
                         xreg = NULL,         # External covariates
                         ic = c("aicc", "aic", "bic", "qic"),  # Information criterion
                         stepwise = TRUE,     # Whether to use stepwise selection
                         nmodels = 94,        # Maximum number of models to try
                         trace = FALSE,       # Whether to print model search progress
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
  
  # Starting model
  if (length(x) < 10L) {
    start.p <- min(start.p, 1L)
    start.q <- min(start.q, 1L)
  }
  
  p <- start.p <- min(start.p, max.p)
  q <- start.q <- min(start.q, max.q)
  
  # Initialize results matrix - only store model configurations and ICs
  results <- matrix(NA, nrow = nmodels, ncol = 3)
  colnames(results) <- c("p", "q", "ic")
  
  # Initial model with intercept
  bestfit <- myingarch(x, order = c(p, q), ic = ic, trace = trace, 
                       xreg = xreg, distr = distribution, link = link, ...)
  
  # Store results
  results[1, ] <- c(p, q, bestfit$ic)
  
  k <- 1
  
  # Fit (0,0) model - explicitly passing order=c(0,0)
  k <- k + 1
  fit <- myingarch(x, order = c(0, 0), ic = ic, trace = trace, 
                   xreg = xreg, distr = distribution, link = link, ...)
  
  # Store results
  results[k, ] <- c(0, 0, fit$ic)
  
  if (fit$ic < bestfit$ic) {
    bestfit <- fit
    p <- 0
    q <- 0
  }
  
  # Basic model with only past observations (p)
  if (max.p > 0 && (p != 1 || q != 0)) {
    k <- k + 1
    
    fit <- myingarch(x, order = c(1, 0), ic = ic, trace = trace, 
                     xreg = xreg, distr = distribution, link = link, ...)
    
    # Store results
    results[k, ] <- c(1, 0, fit$ic)
    
    if (fit$ic < bestfit$ic) {
      bestfit <- fit
      p <- 1
      q <- 0
    }
  }
  
  # Basic model with only past means (q)
  if (max.q > 0 && (p != 0 || q != 1)) {
    k <- k + 1
    
    fit <- myingarch(x, order = c(0, 1), ic = ic, trace = trace, 
                     xreg = xreg, distr = distribution, link = link, ...)
    
    # Store results
    results[k, ] <- c(0, 1, fit$ic)
    
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
    
    # Try decreasing p
    if (p > 0 && newmodel(p - 1, q, results[1:k, ])) {
      k <- k + 1; 
      
      if(k > nmodels)
        break
      
      fit <- myingarch(x, order = c(p - 1, q), ic = ic, trace = trace, 
                       xreg = xreg, distr = distribution, link = link, ...)
      
      # Store results
      results[k, ] <- c(p - 1, q, fit$ic)
      
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        p <- p - 1
        next
      }
    }
    
    # Try increasing q
    if (q < max.q && newmodel(p, q + 1, results[1:k, ])) {
      k <- k + 1;
      
      if(k > nmodels)
        break
      
      fit <- myingarch(x, order = c(p, q + 1), ic = ic, trace = trace, 
                       xreg = xreg, distr = distribution, link = link, ...)
      
      # Store results
      results[k, ] <- c(p, q + 1, fit$ic)
      
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- (q + 1)
        next
      }
    }
    
    # Try decreasing both p and q
    if (q > 0 && p > 0 && newmodel(p - 1, q - 1, results[1:k, ])) {
      k <- k + 1;
      
      if(k > nmodels)
        break
      
      fit <- myingarch(x, order = c(p - 1, q - 1), ic = ic, trace = trace, 
                       xreg = xreg, distr = distribution, link = link, ...)
      
      # Store results
      results[k, ] <- c(p - 1, q - 1, fit$ic)
      
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- (q - 1)
        p <- (p - 1)
        next
      }
    }
    
    # Try decreasing p and increasing q
    if (q < max.q && p > 0 && newmodel(p - 1, q + 1, results[1:k, ])) {
      k <- k + 1;
      
      if(k > nmodels)
        break
      
      fit <- myingarch(x, order = c(p - 1, q + 1), ic = ic, trace = trace, 
                       xreg = xreg, distr = distribution, link = link, ...)
      
      # Store results
      results[k, ] <- c(p - 1, q + 1, fit$ic)
      
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- (q + 1)
        p <- (p - 1)
        next
      }
    }
    
    # Try increasing p and decreasing q
    if (q > 0 && p < max.p && newmodel(p + 1, q - 1, results[1:k, ])) {
      k <- k + 1;
      
      if(k > nmodels)
        break
      
      fit <- myingarch(x, order = c(p + 1, q - 1), ic = ic, trace = trace, 
                       xreg = xreg, distr = distribution, link = link, ...)
      
      # Store results
      results[k, ] <- c(p + 1, q - 1, fit$ic)
      
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- (q - 1)
        p <- (p + 1)
        next
      }
    }
    
    # Try increasing both p and q
    if (q < max.q && p < max.p && newmodel(p + 1, q + 1, results[1:k, ])) {
      k <- k + 1;
      
      if(k > nmodels)
        break
      
      fit <- myingarch(x, order = c(p + 1, q + 1), ic = ic, trace = trace, 
                       xreg = xreg, distr = distribution, link = link, ...)
      
      # Store results
      results[k, ] <- c(p + 1, q + 1, fit$ic)
      
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        q <- (q + 1)
        p <- (p + 1)
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
  
  # Store model search results without warnings
  bestfit$results <- results[1:k, , drop = FALSE]
  
  if (trace) {
    cat("\n\n Best model:", ingarch.string(bestfit, padding = TRUE), "\n\n")
  }
  
  if (!is.null(bestfit$coefficients)) {
    # Add minimum threshold for parameters to prevent convergence to zero
    bestfit$coefficients[bestfit$coefficients < 1e-4] <- 1e-4
  }
  
  return(bestfit)
}