ingarch.string <- function(object, padding=FALSE) {
  # Check if we have stored original order values and use them if available
  if (!is.null(object$orig_order)) {
    p <- object$orig_order[1]
    q <- object$orig_order[2]
  } else {
    # Extract model components from object
    past_obs <- object$model$past_obs    # p parameter
    past_mean <- object$model$past_mean  # q parameter
    
    # Get lengths of past_obs and past_mean
    # Properly handle NULL values to report as 0
    p <- if(is.null(past_obs)) 0 else length(past_obs)
    q <- if(is.null(past_mean)) 0 else length(past_mean)
  }
  
  # Basic INGARCH string
  base_string <- sprintf("INGARCH(%d,%d)", p, q)
  
  # Handle regression components
  if (!is.null(object$xreg) && ncol(object$xreg) > 0) {
    result <- paste("Regression with", base_string, "errors")
  } else {
    result <- base_string
  }
  
  if (!padding) {
    # Strip trailing spaces
    result <- gsub("[ ]*$", "", result)
  }
  
  return(result)
}