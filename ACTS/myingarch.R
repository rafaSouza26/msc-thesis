myingarch <- function(x, order = c(NULL, NULL), ic = "aic", trace = FALSE,
                      xreg = NULL, external = FALSE, # 'external' parameter added, default FALSE
                      distr = "poisson", link = "log", ...) {
  
  # Input validation
  if (!is.numeric(x)) {
    stop("Input time series 'x' must be numeric")
  }
  
  # Match arguments
  distr <- match.arg(distr, choices = c("poisson", "nbinom"))
  link <- match.arg(link, choices = c("log", "identity"))
  ic <- match.arg(ic, choices = c("aic", "aicc", "bic", "qic"))
  
  # MODIFIED: Validation for the 'external' parameter to allow NULL
  if (!is.null(external) && !is.logical(external)) {
    stop("Argument 'external' must be NULL, logical TRUE/FALSE, or a vector of logicals.")
  }
  # This warning correctly handles external = NULL because any(NULL) is FALSE.
  if (any(external) && is.null(xreg)) { 
    warning("'external' suggests active external regressors (e.g., is TRUE or a vector containing TRUE), but 'xreg' is NULL. 'tsglm' might interpret this as effectively no external regressors being active.")
  }
  
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
    external = external    # Use the 'external' function argument (can be TRUE, FALSE, vector, or NULL)
  )
  
  # Handle external regressors (data validation)
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
      model = model_spec, # model_spec now includes the 'external' argument value
      xreg = xreg,
      distr = distr,
      link = link,
      ...
    )
    
    # Set class but maintain the original tsglm class to ensure compatibility
    class(fit) <- c("forecast_tsglm", "tsglm", "tscount")
    
    # Ensure all required elements are present
    required_names <- c("coefficients", "start", "residuals", "fitted.values", 
                        "linear.predictors", "response", "logLik", "score", 
                        "info.matrix", "info.matrix_corrected", "call", "n_obs", 
                        "n_eff", "ts", "model", "xreg", "distr", "distrcoefs", "sigmasq")
    
    # Check which required elements are missing
    missing_elements <- setdiff(required_names, names(fit))
    
    # Add any missing elements
    if (length(missing_elements) > 0) {
      for (element in missing_elements) {
        # Add default values for missing elements
        fit[[element]] <- switch(element,
                                 "coefficients" = numeric(0),
                                 "start" = c(1, frequency(x)),
                                 "residuals" = x - fit$fitted.values, # Simplified, ensure fitted.values exists
                                 "fitted.values" = if(is.null(fit$fitted.values)) rep(mean(x, na.rm = TRUE), length(x)) else fit$fitted.values,
                                 "linear.predictors" = if(is.null(fit$linear.predictors)) rep(0, length(x)) else fit$linear.predictors,
                                 "response" = if(is.null(fit$response)) x else fit$response,
                                 "logLik" = if(is.null(fit$logLik)) -Inf else fit$logLik,
                                 "score" = if(is.null(fit$score)) matrix(0, nrow = 1, ncol = 1) else fit$score,
                                 "info.matrix" = if(is.null(fit$info.matrix)) matrix(0, nrow = 1, ncol = 1) else fit$info.matrix,
                                 "info.matrix_corrected" = if(is.null(fit$info.matrix_corrected)) matrix(0, nrow = 1, ncol = 1) else fit$info.matrix_corrected,
                                 "call" = match.call(),
                                 "n_obs" = length(x),
                                 "n_eff" = length(x) - max_order,
                                 "ts" = x,
                                 "model" = model_spec, # model_spec here will reflect the new 'external' value
                                 "xreg" = xreg,
                                 "distr" = distr,
                                 "distrcoefs" = if(distr == "nbinom" && is.null(fit$distrcoefs)) list(size = 10) else fit$distrcoefs, # Default size param
                                 "sigmasq" = if(is.null(fit$sigmasq)) 1 else fit$sigmasq,
                                 NULL)  # Default case
      }
    }
    
    # Store original order values for reference
    fit$orig_order <- c(orig_p, orig_q)
    
    # Calculate effective sample size and number of parameters
    nstar <- length(x) - max_order
    if (nstar <= 0) {
      stop("Insufficient data for the specified model order")
    }
    
    # Count parameters - always include intercept
    npar <- length(stats::coef(fit)) # Use stats::coef for safety
    
    # Calculate information criteria
    loglik <- fit$logLik
    if (!is.finite(loglik)) {
      warning("Non-finite log-likelihood obtained")
      # Create a minimal list with ic = Inf for consistency in return type for auto.ingarch
      # This part was previously returning list(ic=Inf), let's ensure it has enough structure
      # if auto.ingarch relies on more fields from this error path.
      # For now, keeping it simple as per original structure, but if auto.ingarch needs full dummy_fit, this might need adjustment.
      # The original dummy_fit on error is more complete. This path is for non-finite logLik *after* a successful fit.
      fit$ic <- Inf
      fit$aic <- Inf
      fit$aicc <- Inf
      fit$bic <- Inf
      return(fit) # Return the fit object with Inf ICs
    }
    
    # Calculate AIC (summary.tsglm is robust)
    fit$aic <- summary(fit)$AIC
    
    # Calculate BIC
    fit$bic <- summary(fit)$BIC
    
    # Calculate AICc
    if (nstar <= npar + 2) { # npar should be count of estimated parameters
      fit$aicc <- Inf
      warning("Sample size too small for AICc calculation (nstar <= npar + 2)")
    } else {
      fit$aicc <- fit$aic + 2 * npar * (npar + 1) / (nstar - npar - 1)
    }
    
    # Set information criterion based on specified ic
    fit$ic <- switch(ic,
                     bic = fit$bic,
                     aic = fit$aic,
                     aicc = fit$aicc,
                     qic = fit$aic)  # QIC defaults to AIC for now (as per original)
    
    if (trace) {
      cat("\n INGARCH(", orig_p, ",", orig_q, ")", sep = "")
      cat(" (", distr, " distribution, ", link, " link)", sep="")
      # Updated trace for external regressors
      if (!is.null(xreg) && !is.null(model_spec$external) && any(model_spec$external)) {
        cat(" with external regressors")
      }
      cat(" :", fit$ic)
    }
    
    return(fit)
    
  }, error = function(e) {
    if (trace) {
      cat("\n INGARCH(", orig_p, ",", orig_q, ")", sep = "")
      cat(" (", distr, " distribution, ", link, " link)", sep="")
      # Updated trace for external regressors in error case
      if (!is.null(xreg) && !is.null(model_spec$external) && any(model_spec$external)) {
        cat(" with external regressors")
      }
      cat(" :", Inf)
      cat("\nError: ", conditionMessage(e), "\n")
    }
    
    # Ensure max_order is defined even if error occurred before its calculation for dummy_fit$n_eff
    # This might happen if x is not numeric, for example.
    # However, max_order calculation is placed after basic x validation.
    # If an error occurs really early, match.call() might be the only safe thing.
    # The dummy_fit creation is robust as is for most errors post initial checks.
    
    dummy_fit <- list(
      ic = Inf,
      error = conditionMessage(e),
      # order = c(orig_p, orig_q), # orig_order is what's added to successful fit.
      orig_order = c(orig_p, orig_q), # Match what's added to successful fit
      coefficients = numeric(0),
      start = c(1, if(is.ts(x)) frequency(x) else 1), # frequency(x) only if x is ts
      residuals = rep(NA, length(x)), # length(x) might error if x is not vector-like
      fitted.values = rep(NA, length(x)),
      linear.predictors = rep(NA, length(x)),
      response = if(is.numeric(x)) x else numeric(0), # Ensure response is numeric or empty
      logLik = -Inf,
      score = matrix(0, nrow = 1, ncol = 1),
      info.matrix = matrix(0, nrow = 1, ncol = 1),
      info.matrix_corrected = matrix(0, nrow = 1, ncol = 1),
      call = match.call(), # This is good
      n_obs = if(is.numeric(x) || is.ts(x)) length(x) else 0,
      n_eff = if(is.numeric(x) || is.ts(x)) length(x) - max_order else 0 - max_order, # max_order should be defined
      ts = if(is.numeric(x)) x else numeric(0),
      model = model_spec, # model_spec here will reflect the new 'external' value
      xreg = xreg, # This might be NULL
      distr = distr, # This will be matched
      distrcoefs = if(distr == "nbinom") list(size = 10) else NULL, # Default size, consistent
      sigmasq = 1
    )
    
    class(dummy_fit) <- c("forecast_tsglm", "tsglm", "tscount") # Consistent class
    
    return(dummy_fit)
  })
}