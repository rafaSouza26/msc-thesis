myingarch <- function(x, order = c(0, 0), constant=TRUE, ic="aic", trace=FALSE,
                      xreg=NULL, distr="poisson", link="log", ...) {
  # Ensure orders are non-negative integers
  if (!all(order >= 0) || !all(order == floor(order))) {
    stop("Model orders must be non-negative integers")
  }
  
  # Prepare model specification
  model_spec <- list(
    past_obs = if(order[1] > 0) 1:order[1] else NULL,    # p parameter
    past_mean = if(order[2] > 0) 1:order[2] else NULL,   # q parameter
    external = !is.null(xreg)
  )
  
  # Handle external regressors if provided
  if (!is.null(xreg)) {
    xreg <- as.matrix(xreg)
    if (length(x) != nrow(xreg)) {
      stop("Number of observations in x and xreg don't match")
    }
  }
  
  # Fit model
  if (!constant) {
    # If no constant is desired, add zero offset
    xreg_modified <- if(is.null(xreg)) matrix(0, nrow=length(x), ncol=1) else 
      cbind(0, xreg)
    colnames(xreg_modified)[1] <- "intercept"
    
    fit <- try(tscount::tsglm(
      ts = x,
      model = model_spec,
      xreg = xreg_modified,
      distr = distr,
      link = link,
      ...
    ), silent = TRUE)
  } else {
    fit <- try(tscount::tsglm(
      ts = x,
      model = model_spec,
      xreg = xreg,
      distr = distr,
      link = link,
      ...
    ), silent = TRUE)
  }
  
  if (!inherits(fit, "try-error")) {
    # Calculate effective sample size and number of parameters
    nstar <- length(x) - max(order)  # Effective sample size
    npar <- length(coef(fit)) # Get number of parameters from fitted model
    
    # Calculate information criteria
    loglik <- fit$logLik
    fit$aic <- -2 * loglik + 2 * npar
    fit$bic <- -2 * loglik + npar * log(nstar)
    
    if (nstar <= npar + 2) {
      fit$aicc <- Inf
    } else {
      fit$aicc <- fit$aic + 2 * npar * (npar + 1) / (nstar - npar - 1)
    }
    
    fit$ic <- switch(ic, bic = fit$bic, aic = fit$aic, aicc = fit$aicc)
    
    # Store model parameters
    fit$order <- order
    fit$constant <- constant
    fit$xreg <- xreg
    fit$distr <- distr
    fit$link <- link
    
    if (trace) {
      cat("\n", ingarch.string(fit, padding = TRUE))
      cat(" (", distr, " distribution, ", link, " link)", sep="")
      cat(" :", fit$ic)
    }
    
    return(structure(fit, class = c("forecast_tsglm", "tsglm", "tscount")))
  } else {
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
    return(list(ic = Inf))
  }
}