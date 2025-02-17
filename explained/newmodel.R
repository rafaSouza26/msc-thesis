# Third function: newmodel
#' Check if Model Configuration is New
#' This helper function checks if a particular combination of model parameters
#' has already been evaluated during the model selection process.
#' 
#' @param p AR order
#' @param d Difference order
#' @param q MA order
#' @param D Seasonal difference order
#' @param constant Whether to include constant term
#' @param results Matrix of previously evaluated models
#' @return Logical indicating if this is a new model configuration
newmodel <- function(p, d, q, D, constant, results) {
  # Get number of previously evaluated models
  n <- nrow(results)
  
  # Check each previous result
  for (i in 1:n) {
    # Skip incomplete results
    if(!all(is.na(results[i, seq(7)]))) {
      # Compare current parameters with stored results
      if (all(c(p, d, q, D, constant) == results[i, 1:7])) {
        return(FALSE)  # Model configuration already exists
      }
    }
  }
  
  # If we get here, this is a new configuration
  return(TRUE)
}