ingarch.string <- function(object, padding=FALSE) {
  # Extract model components from object
  past_obs <- object$model$past_obs    # p parameter
  past_mean <- object$model$past_mean  # q parameter
  
  # Get lengths of past_obs and past_mean
  # Properly handle NULL values to report as 0
  p <- if(is.null(past_obs)) 0 else length(past_obs)
  q <- if(is.null(past_mean)) 0 else length(past_mean)
  
  # Basic INGARCH string
  base_string <- sprintf("INGARCH(%d,%d)", p, q)
  
  # Handle regression components and constant term
  if (!is.null(object$xreg)) {
    result <- paste("Regression with", base_string, "errors")
  } else {
    # Check if constant is in the object - handle cases where it might not be
    if (("constant" %in% names(object) && object$constant) || 
        (!("constant" %in% names(object)) && !is.null(object$coefficients) && length(object$coefficients) > 0)) {
      result <- paste(base_string, "with non-zero mean")
    } else {
      result <- paste(base_string, "with zero mean    ")
    }
  }
  
  if (!padding) {
    # Strip trailing spaces
    result <- gsub("[ ]*$", "", result)
  }
  
  return(result)
}