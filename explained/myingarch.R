# Load required helper function for INGARCH model string representation
source("./custom/ingarch.string.R")

#' Fit an Integer-valued GARCH (INGARCH) Model
#' This function fits an INGARCH model to count time series data using the tscount package.
#' It handles both Poisson and Negative Binomial distributions with log or identity links.
#'
#' @param x The input time series (count data)
#' @param order Vector of length 2 specifying the AR(p) and MA(q) orders
#' @param constant Whether to include an intercept term
#' @param ic Information criterion for model selection ("aic", "aicc", or "bic")
#' @param trace Whether to print model fitting progress
#' @param xreg Optional matrix of external regressors
#' @param distr Distribution family ("poisson" or "nbinom")
#' @param link Link function ("log" or "identity")
#' @param ... Additional arguments passed to tscount::tsglm
myingarch <- function(x, order = c(0, 0), constant=TRUE, ic="aic", trace=FALSE,
                      xreg=NULL, distr="poisson", link="log", ...) {
  
  # Calculate effective length of time series accounting for missing values
  # This is important for proper model fitting and information criteria calculation
  missing <- is.na(x)
  firstnonmiss <- head(which(!missing),1)
  lastnonmiss <- tail(which(!missing),1)
  n <- sum(!missing[firstnonmiss:lastnonmiss])
  
  # Input validation
  # Check that order parameter is valid (two non-negative integers)
  if (!is.null(order) && (length(order) != 2 || !all(order >= 0))) {
    stop("order must be a vector of two non-negative integers")
  }
  
  # Validate link function specification
  if (!link %in% c("log", "identity")) {
    stop("Link function must be either 'log' or 'identity'")
  }
  
  # Validate distribution family specification
  if (!distr %in% c("poisson", "nbinom")) {
    stop("Distribution must be either 'poisson' or 'nbinom'")
  }
  
  # Set up model specification
  # past_obs (p) represents AR order
  # past_mean (q) represents MA order
  model_spec <- list(
    past_obs = order[1],    # p parameter 
    past_mean = order[2],   # q parameter
    external = FALSE
  )
  
  # Validate and prepare external regressors if provided
  if (!is.null(xreg)) {
    xreg <- as.matrix(xreg)
    if (length(x) != nrow(xreg)) {
      stop("Number of observations in x and xreg don't match")
    }
  }
  
  # Fit the model based on whether constant term is included
  # Handle constant term through external regressors approach
  if (!constant) {
    # If no constant is desired, use a zero offset as external regressor
    xreg_modified <- if(is.null(xreg)) matrix(0, nrow=length(x), ncol=1) else 
      cbind(0, xreg)
    colnames(xreg_modified)[1] <- "intercept"
    
    # Attempt to fit model without constant
    suppressWarnings(fit <- try(tscount::tsglm(
      ts = x,
      model = model_spec,
      xreg = xreg_modified,
      distr = distr,
      link = link,
      ...
    ), silent = TRUE))
  } else {
    # Fit model with constant term
    suppressWarnings(fit <- try(tscount::tsglm(
      ts = x,
      model = model_spec,
      xreg = xreg,
      distr = distr,
      link = link,
      ...
    ), silent = TRUE))
  }
  
  # Process successful model fit
  if (!is.element("try-error", class(fit))) {
    # Calculate model evaluation metrics
    # Effective sample size accounts for lags used in model
    nstar <- length(x) - max(order)  # Effective sample size
    npar <- length(coef(fit)) # Number of parameters in fitted model
    
    # Calculate information criteria for model selection
    loglik <- fit$logLik
    fit$aic <- -2 * loglik + 2 * npar
    fit$bic <- -2 * loglik + npar * log(nstar)
    
    # Calculate AICc (corrected AIC for small sample sizes)
    # Set to Inf if sample size is too small relative to parameters
    if (nstar <= npar + 2) {
      fit$aicc <- Inf
    } else {
      fit$aicc <- fit$aic + 2 * npar * (npar + 1) / (nstar - npar - 1)
    }
    
    # Store selected information criterion based on user choice
    fit$ic <- switch(ic, bic = fit$bic, aic = fit$aic, aicc = fit$aicc)
    
    # Store model specifications in fit object
    fit$order <- order
    fit$constant <- constant
    fit$xreg <- xreg
    fit$distr <- distr
    fit$link <- link
    
    # Print model information if trace is enabled
    if (trace) {
      cat("\n", ingarch.string(fit, padding = TRUE))
      cat(" (", distr, " distribution, ", link, " link)", sep="")
      cat(" :", fit$ic)
    }
    
    # Return fitted model with appropriate class structure
    return(structure(fit, class = c("forecast_tsglm", "tsglm", "tscount")))
  }
  else {
    # Handle model fitting failures
    # Check for invalid arguments
    if (length(grep("unused argument", fit)) > 0L) {
      stop(fit[1])
    }
    
    # Print failure information if trace is enabled
    if (trace) {
      cat("\n INGARCH(", order[1], ",", order[2], ")", sep = "")
      if (constant) {
        cat(" with non-zero mean")
      } else {
        cat(" with zero mean    ")
      }
      cat(" (", distr, " distribution, ", link, " link)", sep="")
      cat(" :", Inf)
    }
    
    # Return infinite IC to indicate failed fit
    return(list(ic = Inf))
  }
}