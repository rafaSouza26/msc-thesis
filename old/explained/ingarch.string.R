# Second function: ingarch.string
#' Generate String Representation of INGARCH Model
#' This function creates a formatted string description of an INGARCH model,
#' including its order, mean structure, distribution, and link function.
#' 
#' @param object The INGARCH model object to describe
#' @param padding Whether to include padding spaces for alignment
#' @return A string describing the model
ingarch.string <- function(object, padding=FALSE) {
  # Extract model order parameters
  past_obs <- object$model$past_obs    # p parameter (AR order)
  past_mean <- object$model$past_mean  # q parameter (MA order)
  
  # Create base model string with order
  result <- paste("INGARCH(", past_obs, ",", past_mean, ")", sep = "")
  
  # Add information about regression components
  if (!is.null(object$xreg)) {
    # Check if it's a drift model
    if (NCOL(object$xreg) == 1 && is.element("drift", colnames(object$xreg))) {
      result <- paste(result, "with drift        ")
    } else {
      # General regression case
      result <- paste("Regression with", result, "errors")
    }
  } else {
    # Handle mean structure
    if (!is.null(object$intercept) && object$intercept) {
      result <- paste(result, "with non-zero mean")
    } else {
      result <- paste(result, "                  ")
    }
  }
  
  # Add distribution family information
  result <- paste(result, "(", object$distr, " distribution)", sep="")
  
  # Add link function information
  result <- paste(result, "(", object$link, " link)", sep="")
  
  # Remove trailing spaces if padding is not requested
  if (!padding) {
    result <- gsub("[ ]*$", "", result)
  }
  
  return(result)
}
