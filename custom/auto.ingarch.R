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
                         constant = TRUE,     # Whether to include constant term
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
  
  results <- matrix(NA, nrow = nmodels, ncol = 4)
  
  # Initial model with user-specified constant
  bestfit <- myingarch(x, order = c(p, q), constant = constant, ic = ic, 
                       trace = trace, xreg = xreg, 
                       distr = distribution, link = link, ...)
  
  results[1, ] <- c(p, q, constant, bestfit$ic)
  
  
  k <- 1
  
  # Fit (0,0) model
  fit <- myingarch(x, constant = constant, ic = ic,
                   trace = trace, xreg = xreg,
                   distr = distribution, link = link,...)
  results[k, ] <- c(0, 0, constant, fit$ic)
  
  if (fit$ic < bestfit$ic) {
    bestfit <- fit
    # This is symbolical!! Not passing any past_obs or past_mean is the real (0, 0)
    # Do not pass (0, 0)
    p <- 0
    q <- 0
  }
  
  k <- k + 1
  
  # Basic model with only past observations (p)
  if (max.p > 0 && (p != 1 || q != 0)) {
    k <- k + 1
    fit <- myingarch(x, order = c(1, NULL), constant = constant, ic = ic,
                     trace = trace, xreg = xreg,
                     distr = distribution, link = link, ...)
    results[k, ] <- c(1, 0, constant, fit$ic)
    
    if (fit$ic < bestfit$ic) {
      bestfit <- fit
      # This is symbolical!! Not passing any past_obs or past_mean is the real (0, 0)
      p <- 1
      q <- 0
    }
  }
  
  # Basic model with only past means (q)
  if (max.q > 0 && (p != 0 || q != 1)) {
    k <- k + 1
    fit <- myingarch(x, order = c(NULL, 1), constant = constant, ic = ic,
                     trace = trace, xreg = xreg,
                     distr = distribution, link = link, ...)
    results[k, ] <- c(0, 1, constant, fit$ic)
    
    if (fit$ic < bestfit$ic) {
      bestfit <- fit
      # This is symbolical!! Not passing any past_obs or past_mean is the real (0, 0)
      p <- 0
      q <- 1
    }
  }
  
  # Fit (1,1) model
  if (p != 1 || q != 1) {
    k <- k + 1
    fit <- myingarch(x, order = c(1, 1), constant = constant, ic = ic,
                     trace = trace, xreg = xreg,
                     distr = distribution, link = link,...)
    results[k, ] <- c(1, 1, constant, fit$ic)
    
    if (fit$ic < bestfit$ic) {
      bestfit <- fit
      p <- 1
      q <- 1
    }
  }
  
  startk <- 0
  
  cat("\nFitting through step-wise now...\n")
  
  # Stepwise search loop
  while (startk < k && k < nmodels) {
    startk <- k
    
    # Try decreasing p
    if (p > 0 && newmodel(p - 1, q, constant, results[1:k, ])) {
      k <- k + 1; if(k>nmodels) next
      fit <- myingarch(x, order = c(p - 1, q), constant = constant, ic = ic,
                       trace = trace, xreg = xreg,
                       distr = distribution, link = link, ...)
      results[k, ] <- c(p - 1, q, constant, fit$ic)
      if (fit$ic < bestfit$ic) {
        bestfit <- fit
        p <- p - 1
        next
      }
    }
    
    # Try increasing q
    if (q < max.q && newmodel(p, q + 1, constant, results[1:k, ])) {
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
    if (q > 0 && p > 0 && newmodel(p - 1, q - 1, constant, results[1:k, ])) {
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
    if (q < max.q && p > 0 && newmodel(p - 1, q + 1, constant, results[1:k, ])) {
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
    if (q > 0 && p < max.p && newmodel(p + 1, q - 1, constant, results[1:k, ])) {
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
    if (q < max.q && p < max.p && newmodel(p + 1, q + 1, constant, results[1:k, ])) {
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
    
    # Try toggling constant at each step if it improves the model
    if (newmodel(p, q, !constant, results[1:k, ])) {
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
  
  if(k > nmodels){
    warning(sprintf("Stepwise search was stopped early due to reaching the model number limit: `nmodels = %i`", nmodels))
  }
  
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
  
  if (!is.null(bestfit$coefficients)) {
    # Add minimum threshold for parameters to prevent convergence to zero
    bestfit$coefficients[bestfit$coefficients < 1e-4] <- 1e-4
  }
  
  return(bestfit)
}