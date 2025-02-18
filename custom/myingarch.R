myingarch <- function(x, order = c(NULL, NULL), constant = TRUE, ic = "aic", trace = FALSE,
                      xreg = NULL, distr = "poisson", link = "log", ...) {
  # Input validation
  if (!is.numeric(x)) {
    stop("Input time series 'x' must be numeric")
  }
  
  if (!is.logical(constant)) {
    stop("constant must be logical (TRUE/FALSE)")
  }
  
  # if (!all(order >= 0) || !all(is.null(order))) {
  #   stop("Model orders must be non-negative integers")
  # }
  
  # Match arguments
  distr <- match.arg(distr, choices = c("poisson", "nbinom"))
  link <- match.arg(link, choices = c("log", "identity"))
  ic <- match.arg(ic, choices = c("aic", "aicc", "bic", "qic"))
  
  # Convert x to time series if it isn't already
  x <- as.ts(x)
  
  # Prepare model specification
  cat(if(!is.null(order[1]) && order[1] > 0) 1:order[1] else NULL)
  
  model_spec <- list(
    past_obs = if(!is.null(order[1]) && order[1] > 0) 1:order[1] else NULL,    # p parameter
    past_mean = if(!is.null(order[2])) 1:order[2] else NULL,   # q parameter
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
  
  # Try fitting model
  tryCatch({
    # Prepare xreg for constant handling
    final_xreg <- xreg
    
    # if (!constant) {
    #   final_xreg <- if(is.null(xreg)) {
    #     matrix(0, nrow = length(x), ncol = 1)
    #   } else {
    #     cbind(0, xreg)
    #   }
    #   colnames(final_xreg)[1] <- "intercept"
    # }
    
    # Create the model formula
    fit <- tscount::tsglm(
      ts = x,
      model = model_spec,
      xreg = final_xreg,
      distr = distr,
      link = link,
      ...
    )
    
    # Set class
    class(fit) <- c("forecast_tsglm", "tsglm", "tscount")
    
    # Ensure all required elements are present
    fit$ts <- x
    fit$model <- model_spec
    fit$xreg <- final_xreg
    fit$distr <- distr
    fit$link <- link
    
    # Add required elements that might be missing
    if (is.null(fit$coefficients)) fit$coefficients <- fit$coefficients
    if (is.null(fit$start)) fit$start <- fit$start
    if (is.null(fit$n.iter)) fit$n.iter <- 1
    if (is.null(fit$runtime)) fit$runtime <- 0
    if (is.null(fit$score)) fit$score <- numeric(0)
    if (is.null(fit$info.matrix)) fit$info.matrix <- matrix(0, 0, 0)
    if (is.null(fit$info.matrix_corrected)) fit$info.matrix_corrected <- matrix(0, 0, 0)
    if (is.null(fit$inter_par)) fit$inter_par <- list()
    if (is.null(fit$residuals)) {
      fit$residuals <- x - fit$fitted.values
    }
    if (is.null(fit$response)) fit$response <- x
    if (is.null(fit$linear.predictors)) {
      fit$linear.predictors <- fit$fitted.values
    }
    
    # Calculate effective sample size and number of parameters
    nstar <- length(x) - max(order)
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
    
    # Store additional model parameters
    # fit$order <- order
    # fit$constant <- constant
    # fit$n_obs <- length(x)
    # fit$n_eff <- nstar
    # fit$call <- match.call()
    # fit$distrcoefs <- if(distr == "nbinom") list(size = fit$distrcoefs$size) else NULL
    # fit$sigmasq <- if(distr == "nbinom") fit$sigmasq else 0
    
    if (trace) {
      cat("\n", ingarch.string(fit, padding = TRUE))
      cat(" (", distr, " distribution, ", link, " link)", sep="")
      cat(" :", fit$ic)
    }
    
    return(fit)
    
  }, error = function(e) {
    if (trace) {
      cat("\n INGARCH(", order[1], ",", order[2], ")", sep = "")
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