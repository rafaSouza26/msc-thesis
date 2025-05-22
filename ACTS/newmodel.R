newmodel <- function(p, q, results) {
  # Function to check if a specific model configuration has already been tested
  # p: AR order
  # q: MA order
  # results: matrix with p, q, and IC values of already tested models
  
  # Extract only the rows where actual values exist for p and q
  valid_rows <- which(!is.na(results[,1]) & !is.na(results[,2]))
  
  if (length(valid_rows) == 0) {
    return(TRUE)  # No valid models tested yet
  }
  
  # Extract valid results
  valid_results <- results[valid_rows, 1:2, drop = FALSE]
  
  # Check if this model configuration has already been tested
  for (i in 1:nrow(valid_results)) {
    # Compare p and q values with small tolerance for floating point
    if (abs(valid_results[i, 1] - p) < 1e-10 && 
        abs(valid_results[i, 2] - q) < 1e-10) {
      return(FALSE)  # Model has already been tested
    }
  }
  
  return(TRUE)  # Model has not been tested yet
}