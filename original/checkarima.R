# Check that Arima object has positive coefficient variances without returning warnings
checkarima <- function(object) {
  suppressWarnings(test <- any(is.nan(sqrt(diag(object$var.coef)))))
  return(test)
}