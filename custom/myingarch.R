myingarch <- function(x, order = c(NULL, NULL), ic = "aic", trace = FALSE,
                      xreg = NULL, distr = "poisson", link = "log", ...) {
  
  # Input validation
  if (!is.numeric(x)) {
    stop("Input time series 'x' must be numeric")
  }
  
  # Match arguments
  distr <- match.arg(distr, choices = c("poisson", "nbinom"))
  link <- match.arg(link, choices = c("log", "identity"))
  ic <- match.arg(ic, choices = c("aic", "aicc", "bic", "qic"))
  
  # Convert x to time series if it isn't already
  x <- as.ts(x)
  
  # Store original order values for trace output
  # Interpret 0 as NULL for p and q
  orig_p <- if(!is.null(order) && length(order) >= 1) {
    if(is.null(order[1]) || is.na(order[1]))
      0 
    else
      order[1]
  } else 0
  
  orig_q <- if(!is.null(order) && length(order) >= 2) {
    if(is.null(order[2]) || is.na(order[2]))
      0 
    else
      order[2]
  } else 0
  
  # Handle 0 values: convert them to NULL for tscount compatibility
  # For non-zero values, create proper sequence
  p_order <- if(!is.null(order) && length(order) >= 1) {
    if(is.null(order[1]) || is.na(order[1]) || order[1] == 0) NULL 
    else 1:order[1]
  } else NULL
  
  q_order <- if(!is.null(order) && length(order) >= 2) {
    if(is.null(order[2]) || is.na(order[2]) || order[2] == 0) NULL 
    else 1:order[2]
  } else NULL
  
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
    fit$constant <- TRUE  # Always set to TRUE since we're always including the intercept
    
    # Store original order values for reference
    fit$orig_order <- c(orig_p, orig_q)
    
    # Add required elements that might be missing
    if (is.null(fit$coefficients))
      fit$coefficients <- numeric(0)
    
    if (is.null(fit$residuals)) {
      fit$residuals <- x - fit$fitted.values
    }
    
    # Calculate effective sample size and number of parameters
    nstar <- length(x) - max_order
    if (nstar <= 0) {
      stop("Insufficient data for the specified model order")
    }
    
    # Count parameters - always include intercept
    npar <- length(coef(fit))
    
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
      # Use the original order values for trace output
      cat("\n INGARCH(", orig_p, ",", orig_q, ")", sep = "")
      cat(" with non-zero mean")  # Always show "with non-zero mean"
      cat(" (", distr, " distribution, ", link, " link)", sep="")
      cat(" :", fit$ic)
    }
    
    return(fit)
    
  }, error = function(e) {
    if (trace) {
      cat("\n INGARCH(", orig_p, ",", orig_q, ")", sep = "")
      cat(" with non-zero mean")  # Always show "with non-zero mean"
      cat(" (", distr, " distribution, ", link, " link)", sep="")
      cat(" :", Inf)
      cat("\nError: ", conditionMessage(e), "\n")
    }
    return(list(
      ic = Inf,
      error = conditionMessage(e),
      order = c(orig_p, orig_q),
      constant = TRUE,  # Always set to TRUE
      distr = distr,
      link = link
    ))
  })
}