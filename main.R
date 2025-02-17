# Load required packages
library(tscount)

# Load our custom functions
source("./custom/ingarch.sim.R")
source("./custom/myingarch.R")
source("./custom/auto.ingarch.R")

# Set a random seed for reproducibility
set.seed(123)

# First, let's generate some data with known parameters
# We'll use a simple INGARCH(1,1) model where today's count depends on:
#  - yesterday's actual count (past_obs = 0.3)
#  - yesterday's expected count (past_mean = 0.4)
#  - a baseline level (intercept = 1.5)
simulated_data <- ingarch.sim(
  n = 1000,           # Generate 200 observations
  param = list(
    intercept = 1.5, # Baseline level
    past_obs = 0.3,  # Effect of yesterday's count
    past_mean = 0.4  # Effect of yesterday's expectation
  ),
  model = list(
    past_obs = 1,    # Use lag-1 for past observations
    past_mean = 1,   # Use lag-1 for past means
    external = FALSE # No external regressors
  ),
  link = "identity", # Use identity link function
  distr = "poisson"  # Use Poisson distribution
)

# Convert to time series object for model fitting
ts_data <- ts(simulated_data$ts)

# Print the true parameters we used to generate the data
cat("\nTrue parameters used in simulation:\n")
cat("Intercept:", simulated_data$parameters$intercept, "\n")
cat("Past observations coefficient (p=1):", simulated_data$parameters$past_obs, "\n")
cat("Past means coefficient (q=1):", simulated_data$parameters$past_mean, "\n")

# Now let's try to recover these parameters using auto.ingarch
# We'll set maximum orders slightly higher than the true values to see
# if the function can identify the correct model order
cat("\nFitting model with auto.ingarch...\n")
fitted_model <- try({
  auto.ingarch(
    y = ts_data,
    max.p = 2,        # Maximum AR order (true is 1)
    max.q = 2,        # Maximum MA order (true is 1)
    distribution = "poisson",
    link = "identity",
    trace = TRUE,
    nmodels = 100
  )
})

# Check if model fitting was successful and display results
if (!inherits(fitted_model, "try-error")) {
  cat("\nFitted model parameters:\n")
  print(summary(fitted_model))
  
  # Create plot comparing actual vs fitted values
  par(mfrow=c(2,1))
  
  # Plot the time series
  plot(ts_data, 
       type = "l",
       main = "Actual vs Fitted Values",
       ylab = "Count",
       xlab = "Time")
  lines(fitted(fitted_model), col = "red")
  legend("topright", 
         legend = c("Actual", "Fitted"),
         col = c("black", "red"),
         lty = 1)
  
  # Plot residuals
  plot(residuals(fitted_model),
       type = "h",
       main = "Model Residuals",
       ylab = "Residual",
       xlab = "Time")
  abline(h = 0, col = "red", lty = 2)
} else {
  cat("\nError in model fitting:", fitted_model, "\n")
}