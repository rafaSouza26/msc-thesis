ingarch.string <- function(object, padding=FALSE) {
  # Extract model components from object
  past_obs <- object$model$past_obs    # p parameter
  past_mean <- object$model$past_mean  # q parameter
  
  # Get lengths of past_obs and past_mean
  p <- if(is.null(past_obs)) 0 else length(past_obs)
  q <- if(is.null(past_mean)) 0 else length(past_mean)
  
  # Basic INGARCH string
  result <- paste("INGARCH(", p, ",", q, ")", sep = "")
  
  # Handle regression components
  if (!is.null(object$xreg)) {
    if (NCOL(object$xreg) == 1 && is.element("drift", colnames(object$xreg))) {
      result <- paste(result, "with drift        ")
    } else {
      result <- paste("Regression with", result, "errors")
    }
  } else {
    # Check for intercept
    if (!is.null(object$intercept) && object$intercept) {
      result <- paste(result, "with non-zero mean")
    } else {
      result <- paste(result, "                  ")
    }
  }
  
  if (!padding) {
    # Strip trailing spaces
    result <- gsub("[ ]*$", "", result)
  }
  
  return(result)
}