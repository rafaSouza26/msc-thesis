newmodel <- function(p, d, q, D, constant, results) {
  n <- nrow(results)
  for (i in 1:n) {
    if(!all(is.na(results[i, 1:3]))) {  # Check only p, q, constant
      if (all(c(p, q, constant) == results[i, 1:3])) {  # Compare only stored values
        return(FALSE)
      }
    }
  }
  return(TRUE)
}