myingarch <- function(x, order = c(NULL, NULL), constant = TRUE, ic = "aic", trace = FALSE,
                      xreg = NULL, distr = "poisson", link = "log", ...) {
  # Input validation
  if (!is.numeric(x)) {
    stop("Input time series 'x' must be numeric")
  }
  
  if (!is.logical(constant)) {
    stop("constant must be logical (TRUE/FALSE)")
  }
  
  # Match arguments
  distr <- match.arg(distr, choices = c("poisson", "nbinom"))
  link <- match.arg(link, choices = c("log", "identity"))
  ic <- match.arg(ic, choices = c("aic", "aicc", "bic", "qic"))
  
  # Convert x to time series if it isn't already
  x <- as.ts(x)
  
  # Properly handle order parameters for p and q with safer conditional checks
  # Make sure we check both NULL and NA values
  p_order <- if(!is.null(order) && length(order) >= 1 && !is.null(order[1]) && !is.na(order[1]) && order[1] > 0) 1:order[1] else NULL
  q_order <- if(!is.null(order) && length(order) >= 2 && !is.null(order[2]) && !is.na(order[2]) && order[2] > 0) 1:order[2] else NULL
  
  model_spec <- list(
    past_obs = p_order,    # p parameter
    past_mean = q_order,   # q parameter
    external = !is.null(xreg)
  )
  
  # Handle external regressors
  if (!is.null(xreg)) {
    if (!is.numeric(xreg)) {
      stop("xreg should be a numeric matrix or a numeric vector")
    }
    xreg <- as.matrix(xreg)
    if (length(x) != nrow(xreg)) {
      stop("Number of observations in x and xreg don't match")
    }
  }
  
  # Get the maximum lag for effective sample size calculation
  max_order <- max(
    if(is.null(p_order)) 0 else max(p_order),
    if(is.null(q_order)) 0 else max(q_order)
  )
  
  # Try fitting model
  tryCatch({
    # Create the model formula
    fit <- tscount::tsglm(
      ts = x,
      model = model_spec,
      xreg = xreg,
      distr = distr,
      link = link,
      ...
    )
    
    # Set class
    class(fit) <- c("forecast_tsglm", "tsglm", "tscount")
    
    # Ensure all required elements are present
    fit$ts <- x
    fit$model <- model_spec
    fit$xreg <- xreg
    fit$distr <- distr
    fit$link <- link
    fit$constant <- constant
    
    # Add required elements that might be missing
    if (is.null(fit$coefficients)) fit$coefficients <- numeric(0)
    if (is.null(fit$residuals)) {
      fit$residuals <- x - fit$fitted.values
    }
    
    # Calculate effective sample size and number of parameters
    nstar <- length(x) - max_order
    if (nstar <= 0) {
      stop("Insufficient data for the specified model order")
    }
    
    # Count parameters
    npar <- length(coef(fit))
    if (!constant) {
      npar <- npar - 1  # Subtract the constrained intercept
    }
    
    # Calculate information criteria
    loglik <- fit$logLik
    if (!is.finite(loglik)) {
      warning("Non-finite log-likelihood obtained")
      return(list(ic = Inf))
    }
    
    # Calculate AIC
    fit$aic <- -2 * loglik + 2 * npar
    
    # Calculate BIC
    fit$bic <- -2 * loglik + npar * log(nstar)
    
    # Calculate AICc
    if (nstar <= npar + 2) {
      fit$aicc <- Inf
      warning("Sample size too small for AICc calculation")
    } else {
      fit$aicc <- fit$aic + 2 * npar * (npar + 1) / (nstar - npar - 1)
    }
    
    # Set information criterion based on specified ic
    fit$ic <- switch(ic,
                     bic = fit$bic,
                     aic = fit$aic,
                     aicc = fit$aicc,
                     qic = fit$aic)  # QIC defaults to AIC for now
    
    if (trace) {
      cat("\n", ingarch.string(fit, padding = TRUE))
      cat(" (", distr, " distribution, ", link, " link)", sep="")
      cat(" :", fit$ic)
    }
    
    return(fit)
    
  }, error = function(e) {
    if (trace) {
      p_val <- if(is.null(order) || length(order) < 1 || is.null(order[1]) || is.na(order[1])) 0 else order[1]
      q_val <- if(is.null(order) || length(order) < 2 || is.null(order[2]) || is.na(order[2])) 0 else order[2]
      cat("\n INGARCH(", p_val, ",", q_val, ")", sep = "")
      cat(if(constant) " with non-zero mean" else " with zero mean    ")
      cat(" (", distr, " distribution, ", link, " link)", sep="")
      cat(" :", Inf)
      cat("\nError: ", conditionMessage(e), "\n")
    }
    return(list(
      ic = Inf,
      error = conditionMessage(e),
      order = order,
      constant = constant,
      distr = distr,
      link = link
    ))
  })
}