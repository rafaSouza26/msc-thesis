myingarch <- function(x, order = c(0, 0), constant = TRUE, ic = "aic", trace = FALSE,
                      xreg = NULL, distr = "poisson", link = "log", ...) {
  # Input validation
  if (!is.numeric(x)) {
    stop("Input time series 'x' must be numeric")
  }
  
  if (!is.logical(constant)) {
    stop("constant must be logical (TRUE/FALSE)")
  }
  
  if (!all(order >= 0) || !all(order == floor(order))) {
    stop("Model orders must be non-negative integers")
  }
  
  # Match arguments
  distr <- match.arg(distr, choices = c("poisson", "nbinom"))
  link <- match.arg(link, choices = c("log", "identity"))
  ic <- match.arg(ic, choices = c("aic", "aicc", "bic", "qic"))
  
  # Prepare model specification
  model_spec <- list(
    past_obs = if(order[1] > 0) 1:order[1] else NULL,    # p parameter
    past_mean = if(order[2] > 0) 1:order[2] else NULL,   # q parameter
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
    if (!constant) {
      xreg_modified <- if(is.null(xreg)) {
        matrix(0, nrow = length(x), ncol = 1)
      } else {
        cbind(0, xreg)
      }
      colnames(xreg_modified)[1] <- "intercept"
      
      fit <- tscount::tsglm(
        ts = x,
        model = model_spec,
        xreg = xreg_modified,
        distr = distr,
        link = link,
        ...
      )
    } else {
      fit <- tscount::tsglm(
        ts = x,
        model = model_spec,
        xreg = xreg,
        distr = distr,
        link = link,
        ...
      )
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
    
    # Add required elements for tsglm class
    fit$ts <- x
    fit$model <- model_spec
    fit$distribution <- distr
    fit$link <- link
    fit$start <- fit$start
    fit$n.iter <- fit$n.iter
    fit$runtime <- fit$runtime
    fit$score <- if(!is.null(fit$score)) fit$score else NULL
    fit$info <- if(!is.null(fit$info)) fit$info else NULL
    fit$inter_par <- NULL
    fit$residuals <- fit$residuals
    fit$response <- x
    
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
    
    # Store model parameters
    fit$order <- order
    fit$constant <- constant
    fit$xreg <- if (!constant && !is.null(xreg)) xreg[, -1, drop = FALSE] else xreg
    
    if (trace) {
      cat("\n", ingarch.string(fit, padding = TRUE))
      cat(" (", distr, " distribution, ", link, " link)", sep="")
      cat(" :", fit$ic)
    }
    
    # Set appropriate class and return
    class(fit) <- c("forecast_tsglm", "tsglm", "tscount")
    return(fit)
    
  }, error = function(e) {
    if (trace) {
      cat("\n INGARCH(", order[1], ",", order[2], ")", sep = "")
      cat(if(constant) " with non-zero mean" else " with zero mean    ")
      cat(" (", distr, " distribution, ", link, " link)", sep="")
      cat(" :", Inf)
      cat("\nError: ", conditionMessage(e), "\n")
    }
    return(list(ic = Inf))
  })
}