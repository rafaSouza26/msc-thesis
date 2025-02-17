source("./custom/ingarch.string.R")

myingarch <- function(x, order = c(0, 0), constant=TRUE, ic="aic", trace=FALSE,
                      xreg=NULL, distr="poisson", link="log", ...) {
  # Length of non-missing interior
  missing <- is.na(x)
  firstnonmiss <- head(which(!missing),1)
  lastnonmiss <- tail(which(!missing),1)
  n <- sum(!missing[firstnonmiss:lastnonmiss])
  
  if (!is.null(order) && (length(order) != 2 || !all(order >= 0))) {
    stop("order must be a vector of two non-negative integers")
  }
  
  # Check valid link function
  if (!link %in% c("log", "identity")) {
    stop("Link function must be either 'log' or 'identity'")
  }
  
  # Check valid distribution
  if (!distr %in% c("poisson", "nbinom")) {
    stop("Distribution must be either 'poisson' or 'nbinom'")
  }
  
  
  # Prepare model fit
  model_spec <- list(
    past_obs = if(order[1] > 0) 1:order[1] else NULL,    # p parameter 
    past_mean = if(order[2] > 0) 1:order[2] else NULL,   # q parameter
    external = FALSE
  )
  
  # Handle external regressors if provided
  if (!is.null(xreg)) {
    xreg <- as.matrix(xreg)
    if (length(x) != nrow(xreg)) {
      stop("Number of observations in x and xreg don't match")
    }
  }
  
  # Fit model based on constant parameter
  # Instead of using include.intercept, we'll modify the model specification
  # A partida o modelo vai ter sempre constante, por enquanto não precisamos de ter essa diferenciação
  # Comentar
  if (!constant) {
    # If no constant is desired, we use an offset that's always zero
    xreg_modified <- if(is.null(xreg)) matrix(0, nrow=length(x), ncol=1) else 
      cbind(0, xreg)
    colnames(xreg_modified)[1] <- "intercept"
    
    # No INGARCH temos de usar [1:13] em vez de 13 para calcular os coeficientes
    # senão ele calcula apenas o coef 13
    suppressWarnings(fit <- try(tscount::tsglm(
      ts = x,
      model = model_spec,
      xreg = xreg_modified,
      distr = distr,
      link = link,
      ...
    ), silent = TRUE))
  } else {
    # No INGARCH temos de usar [1:13] em vez de 13 para calcular os coeficientes
    # senão ele calcula apenas o coef 13
    suppressWarnings(fit <- try(tscount::tsglm(
      ts = x,
      model = model_spec,
      xreg = xreg,
      distr = distr,
      link = link,
      ...
    ), silent = TRUE))
  }
  
  if (!is.element("try-error", class(fit))) {
    # Calculate effective sample size and number of parameters
    nstar <- length(x) - max(order)  # Effective sample size
    npar <- length(coef(fit)) # Get number of parameters from fitted model
    
    # Calculate information criteria
    loglik <- fit$logLik
    fit$aic <- -2 * loglik + 2 * npar
    fit$bic <- -2 * loglik + npar * log(nstar)
    
    # Adicionar aqui o AICc e QIC
    
    # Add check for aicc calculation
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
  }
  else {
    if (length(grep("unused argument", fit)) > 0L) {
      stop(fit[1])
    }
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