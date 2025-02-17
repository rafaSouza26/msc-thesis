ingarch.sim <- function(n, param = list(intercept = 1, past_obs = NULL, past_mean = NULL, 
                                        xreg = NULL), 
                        model = list(past_obs = NULL, past_mean = NULL, external = FALSE), 
                        xreg = NULL,
                        link = c("identity", "log"),
                        distr = c("poisson", "nbinom"),
                        distrcoefs = NULL,
                        fit = NULL, 
                        n_start = 50) {
  
  # Input validation
  if(!is.numeric(n) || n <= 0)
    stop("n must be a positive integer")
  
  # If fit object is provided, use its parameters
  if(!is.null(fit)) {
    if(!inherits(fit, "tsglm"))
      stop("fit must be an object of class 'tsglm'")
    
    model <- fit$model
    link <- fit$link
    distr <- fit$distr
    distrcoefs <- fit$distrcoefs
    param <- list(
      intercept = fit$coefficients[1],
      past_obs = fit$coefficients[2:(1+length(model$past_obs))],
      past_mean = fit$coefficients[(2+length(model$past_obs)):length(fit$coefficients)]
    )
  }
  
  # Parameter validation
  link <- match.arg(link)
  distr <- match.arg(distr)
  
  if(!is.null(param$past_obs) && (!is.numeric(param$past_obs) || any(param$past_obs < 0)))
    stop("past_obs coefficients must be non-negative")
  if(!is.null(param$past_mean) && (!is.numeric(param$past_mean) || any(param$past_mean < 0)))
    stop("past_mean coefficients must be non-negative")
  if(distr == "nbinom" && (is.null(distrcoefs$size) || distrcoefs$size <= 0))
    stop("size parameter must be positive for negative binomial distribution")
  
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
  
  # Call tsglm.sim
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
  
  # Store the true parameters in the result
  result$parameters <- param
  result$model <- model
  result$link <- link
  result$distr <- distr
  
  # Return result with class to match specification
  class(result) <- c("tsglm.sim", "tsglm", "tscount")
  
  return(result)
}