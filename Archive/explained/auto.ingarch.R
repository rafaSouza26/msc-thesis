# Load required custom functions and libraries
source("./ACTS/myingarch.R")
source("./ACTS/ingarch.string.R")
source("./ACTS/newmodel.R")
library(tscount)
library(forecast)

#' Automatic INGARCH Model Selection Function
#' This function performs automatic selection of Integer-valued GARCH (INGARCH) models
#' for count time series data. It implements a stepwise search algorithm to find
#' the best model according to the specified information criterion.
#'
#' @param y The input time series (count data)
#' @param max.p Maximum AR order to consider
#' @param max.q Maximum MA order to consider
#' @param start.p Starting AR order for stepwise search
#' @param start.q Starting MA order for stepwise search
#' @param ic Information criterion for model selection ("aicc", "aic", or "bic")
#' @param stepwise Whether to use stepwise selection algorithm
#' @param nmodels Maximum number of models to evaluate
#' @param trace Whether to print progress during model search
#' @param method Estimation method (currently unused)
#' @param parallel Whether to use parallel processing
#' @param num.cores Number of cores for parallel processing
#' @param x Alternative input data (same as y)
#' @param distribution Distribution family ("poisson" or "nbinom")
#' @param link Link function ("log" or "identity")
#' @param xreg External regressors (optional)
auto.ingarch <- function(y, 
                        max.p = 5,
                        max.q = 5,
                        start.p = 2,
                        start.q = 2,
                        ic = c("aicc", "aic", "bic"),
                        stepwise = TRUE,
                        nmodels = 94,
                        trace = FALSE,
                        method = NULL,
                        parallel = FALSE,
                        num.cores = 2,
                        x = y,
                        distribution = c("poisson", "nbinom"),
                        link = c("log", "identity"),
                        xreg = NULL,
                        ...) {
  
  # Check parallel processing compatibility
  if (stepwise && parallel) {
    warning("Parallel computing is only implemented when stepwise=FALSE, the model will be fit in serial.")
    parallel <- FALSE
  }
  
  # Disable tracing for parallel execution
  if (trace && parallel) {
    message("Tracing model searching in parallel is not supported.")
    trace <- FALSE
  }
  
  # Setup and initial checks
  series <- deparse(substitute(y))
  x <- as.ts(x)
  
  # Validate start.p against max.p
  if (!is.null(start.p) && start.p > max.p) {
    warning("start.p cannot be greater than max.p, setting start.p = max.p")
    start.p <- max.p
  }
  
  # Input validation checks
  # Check if the input is univariate
  if (NCOL(x) > 1) {
    stop("auto.ingarch can only handle univariate time series")
  }
  
  # Verify data is count data (non-negative integers)
  if (!all(x[!is.na(x)] >= 0) || !all(x[!is.na(x)] == floor(x[!is.na(x)]))) {
    stop("INGARCH models require non-negative integer (count) data")
  }
  
  # Handle missing values and determine effective series length
  orig.x <- x
  missing <- is.na(x)
  firstnonmiss <- head(which(!missing), 1)
  lastnonmiss <- tail(which(!missing), 1)
  serieslength <- sum(!missing[firstnonmiss:lastnonmiss])
  
  # External regressor validation and preparation
  if(!is.null(xreg)) {
    # Validate xreg input
    if(!is.numeric(xreg))
      stop("xreg should be a numeric matrix or a numeric vector")
    
    xreg <- as.matrix(xreg)
    xreg <- xreg[firstnonmiss:NROW(xreg),,drop=FALSE]
    
    # Check dimensions and missing values
    if(NROW(xreg) != length(x))
      stop("Number of rows in xreg must match length of time series")
    if(any(is.na(xreg)))
      stop("Missing values in external regressors are not allowed")
    
    # Add default column names if missing
    if(is.null(colnames(xreg))) {
      colnames(xreg) <- paste0("xreg", 1:NCOL(xreg))
    }
  }
  
  # Handle constant data case
  if (all(x == x[1], na.rm = TRUE)) {
    if(all(is.na(x)))
      stop("All data are missing")
    
    # Create a simple model for constant data
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
      fit$size <- Inf  # Equivalent to Poisson for constant data
    }
    
    class(fit) <- "tsglm"
    return(fit)
  }
  
  # Model selection setup
  ic <- match.arg(ic)
  
  # Adjust maximum orders based on series length
  max.p <- min(max.p, floor(serieslength / 3))
  max.q <- min(max.q, floor(serieslength / 3))
  
  # Use AIC for very short series
  if (serieslength <= 3L) {
    ic <- "aic"
  }
  
  # Fixed parameters
  d <- D <- 0  # No differencing in INGARCH
  constant <- TRUE
  
  # Initialize model search
  p <- start.p <- min(start.p, max.p)
  q <- start.q <- min(start.q, max.q)
  
  # Adjust starting orders for very short series
  if (length(x) < 10L) {
    start.p <- min(start.p, 1L)
    start.q <- min(start.q, 1L)
  }
  
  # Matrix to store results of all tried models
  results <- matrix(NA, nrow = nmodels, ncol = 4)
  
  # Fit initial model with starting values
  bestfit <- myingarch(x, order = c(p, q), constant = constant, ic = ic, 
                       trace = trace, xreg = xreg, 
                       distr = distribution, link = link, ...)
  
  results[1, ] <- c(p, q, constant, bestfit$ic)
  
  # Fit null model (no AR or MA terms)
  fit <- myingarch(x, order = c(0, 0), constant = constant, ic = ic, 
                   trace = trace, xreg = xreg,
                   distr = distribution, link = link, ...)
  
  results[2, ] <- c(0, 0, constant, fit$ic)
  
  if (fit$ic < bestfit$ic) {
    bestfit <- fit
    p <- q <- 0
  }
  
  k <- 2  # Counter for number of models tried
  
  # Try basic models with only AR or only MA terms
  # First try AR(1) model if max.p > 0
  if (max.p > 0) {
    fit <- myingarch(x, order = c(max.p > 0, 0), constant = constant, ic = ic,
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
  
  # Try MA(1) model if max.q > 0
  if (max.q > 0) {
    fit <- myingarch(x, order = c(0, max.q > 0), constant = constant, ic = ic,
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
  
  # Try null model without constant
  if (constant) {
    fit <- myingarch(x, order = c(0, 0), constant = FALSE, ic = ic,
                     trace = trace, xreg = xreg,
                     distr = distribution, link = link, ...)
    
    results[k+1, ] <- c(0, 0, 0, fit$ic)
    
    if (fit$ic < bestfit$ic) {
      bestfit <- fit
      p <- q <- 0
      constant <- FALSE
    }
    
    k <- k + 1
  }
  
  # Begin stepwise search
  startk <- 0
  
  # Main stepwise selection loop
  while (startk < k && k < nmodels) {
    startk <- k
    
    # Try all possible one-step changes to current model:
    # 1. Decrease/increase p
    # 2. Decrease/increase q
    # 3. Change both p and q
    # 4. Toggle constant term
    
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
    
    # [Similar blocks for other model variations...]
    # Each block tries a different modification to the current best model
    # and updates if it finds a better fit
    
  }
  
  # Warning if model limit reached
  if(k > nmodels){
    warning(sprintf("Stepwise search was stopped early due to reaching the model number limit: `nmodels = %i`", nmodels))
  }
  
  # Handle case where no suitable model was found
  if (bestfit$ic == Inf) {
    if (trace) {
      cat("\n")
    }
    stop("No suitable INGARCH model found")
  }
  
  # Prepare final model for return
  bestfit$x <- orig.x
  bestfit$series <- series
  bestfit$ic <- NULL
  bestfit$call <- match.call()
  bestfit$call$x <- data.frame(x = x)
  
  # Print best model if trace is enabled
  if (trace) {
    cat("\n\n Best model:", ingarch.string(bestfit, padding = TRUE), "\n\n")
  }
  
  return(bestfit)
}