# First function: ingarch.sim
library(tscount)

#' Simulate Count Time Series from an INGARCH Model
#' This function generates simulated count time series data following an INGARCH model specification.
#' It can either use parameters directly specified or extract them from a fitted model object.
#' 
#' @param n Integer. Number of observations to simulate
#' @param param List containing model parameters:
#'        - intercept: Base level
#'        - past_obs: AR coefficients
#'        - past_mean: MA coefficients
#'        - xreg: External regressors
#' @param model List specifying model structure
#' @param xreg External regressors matrix
#' @param link Link function type
#' @param distr Distribution family
#' @param distrcoefs Distribution-specific parameters
#' @param fit Optional fitted model object to use for simulation
#' @param n_start Burn-in period length
ingarch.sim <- function(n, param = list(intercept = 1, past_obs = NULL, past_mean = NULL, 
                                        xreg = NULL), 
                        model = list(past_obs = NULL, past_mean = NULL, external = FALSE), 
                        xreg = NULL,
                        link = c("identity", "log"),
                        distr = c("poisson", "nbinom"),
                        distrcoefs = NULL,
                        fit = NULL, 
                        n_start = 50) {
  
  # Basic input validation
  if(!is.numeric(n) || n <= 0)
    stop("n must be a positive integer")
  
  # Extract parameters from fitted model if provided
  if(!is.null(fit)) {
    # Validate fit object type
    if(!inherits(fit, "tsglm"))
      stop("fit must be an object of class 'tsglm'")
    
    # Extract model components from fit object
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
  
  # Parameter validation and processing
  link <- match.arg(link)
  distr <- match.arg(distr)
  
  # Validate coefficient constraints
  if(!is.null(param$past_obs) && (!is.numeric(param$past_obs) || any(param$past_obs < 0)))
    stop("past_obs coefficients must be non-negative")
  if(!is.null(param$past_mean) && (!is.numeric(param$past_mean) || any(param$past_mean < 0)))
    stop("past_mean coefficients must be non-negative")
  if(distr == "nbinom" && (is.null(distrcoefs$size) || distrcoefs$size <= 0))
    stop("size parameter must be positive for negative binomial distribution")
  
  # Prepare parameters for tsglm.sim function
  ts_param <- list(
    intercept = param$intercept,
    past_obs = param$past_obs,
    past_mean = param$past_mean
  )
  
  # Convert model specification to sequence format required by tsglm.sim
  ts_model <- list(
    past_obs = if(is.null(model$past_obs)) NULL else seq_along(model$past_obs),
    past_mean = if(is.null(model$past_mean)) NULL else seq_along(model$past_mean),
    external = model$external
  )
  
  # Generate simulated data using tsglm.sim
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
  
  # Add required components to output
  result$ts <- result$ts
  result$linear.predictors <- result$linear.predictors
  result$xreg.effects <- result$xreg.effects
  
  # Set appropriate class for the result
  class(result) <- c("tsglm.sim", "tsglm", "tscount")
  
  return(result)
}
