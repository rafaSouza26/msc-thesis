ingarch.sim <- function(n, param = list(intercept = 1, past_obs = NULL, past_mean = NULL, 
                                        xreg = NULL), 
                        model = list(past_obs = NULL, past_mean = NULL, external = FALSE), 
                        xreg = NULL,
                        link = c("identity", "log"),
                        distr = c("poisson", "nbinom"),
                        distrcoefs = NULL,
                        fit = NULL, 
                        n_start = 50) {
  
  # Enhanced input validation with specific error messages
  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n != floor(n)) {
    stop("n must be a positive integer")
  }
  
  if (n_start < 0 || n_start != floor(n_start)) {
    stop("n_start must be a non-negative integer")
  }
  
  # Validate fit object if provided
  if (!is.null(fit)) {
    if (!inherits(fit, "tsglm")) {
      stop("fit must be an object of class 'tsglm'")
    }
    
    tryCatch({
      model <- fit$model
      link <- fit$link
      distr <- fit$distr
      distrcoefs <- fit$distrcoefs
      param <- list(
        intercept = fit$coefficients[1],
        past_obs = fit$coefficients[2:(1+length(model$past_obs))],
        past_mean = fit$coefficients[(2+length(model$past_obs)):length(fit$coefficients)]
      )
    }, error = function(e) {
      stop("Failed to extract parameters from fit object: ", conditionMessage(e))
    })
  }
  
  # Validate link and distribution
  link <- match.arg(link)
  distr <- match.arg(distr)
  
  # Validate parameters
  if (!is.numeric(param$intercept) || length(param$intercept) != 1) {
    stop("intercept must be a single numeric value")
  }
  
  if (!is.null(param$past_obs)) {
    if (!is.numeric(param$past_obs)) {
      stop("past_obs coefficients must be numeric")
    }
    if (any(param$past_obs < 0)) {
      stop("past_obs coefficients must be non-negative")
    }
  }
  
  if (!is.null(param$past_mean)) {
    if (!is.numeric(param$past_mean)) {
      stop("past_mean coefficients must be numeric")
    }
    if (any(param$past_mean < 0)) {
      stop("past_mean coefficients must be non-negative")
    }
  }
  
  # Validate distribution parameters
  if (distr == "nbinom") {
    if (is.null(distrcoefs) || is.null(distrcoefs$size)) {
      stop("size parameter must be provided for negative binomial distribution")
    }
    if (!is.numeric(distrcoefs$size) || distrcoefs$size <= 0) {
      stop("size parameter must be positive for negative binomial distribution")
    }
  }
  
  # Validate external regressors
  if (!is.null(xreg)) {
    if (!is.numeric(xreg)) {
      stop("xreg must be numeric")
    }
    xreg <- as.matrix(xreg)
    if (nrow(xreg) != n) {
      stop(sprintf("Number of rows in xreg (%d) must match n (%d)", nrow(xreg), n))
    }
  }
  
  # Set up tsglm.sim parameters
  ts_param <- list(
    intercept = param$intercept,
    past_obs = param$past_obs,
    past_mean = param$past_mean
  )
  
  ts_model <- list(
    past_obs = if(is.null(model$past_obs)) NULL else seq_along(model$past_obs),
    past_mean = if(is.null(model$past_mean)) NULL else seq_along(model$past_mean),
    external = model$external
  )
  
  # Try simulation with error handling
  tryCatch({
    result <- tsglm.sim(
      n = n,
      param = ts_param,
      model = ts_model,
      xreg = xreg,
      link = link,
      distr = distr,
      distrcoefs = distrcoefs,
      n_start = n_start
    )
    
    # Store additional information
    result$parameters <- param
    result$model <- model
    result$link <- link
    result$distr <- distr
    
    # Set class and return
    class(result) <- c("tsglm.sim", "tsglm", "tscount")
    return(result)
    
  }, error = function(e) {
    stop("Simulation failed: ", conditionMessage(e))
  })
}