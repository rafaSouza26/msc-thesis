# Load necessary libraries
library(tscount)  # Make sure this is installed

# Set seed for reproducibility
set.seed(456)

# Number of realizations and observations
n_realizations <- 1000
n_observations <- 1000
n_start <- 50  # Burn-in period

# Load covariates from the provided data
# Assuming d_1_data is available or can be loaded from a file
# d_1_data <- readRDS("d_1_data.rds")  # Uncomment if needed

# Extract the necessary covariates (first 1000 observations)
covariates <- d_1_data[1:1000, c("mean.Temp.Lag5", "mean.DTemp.Lag0", "PM25.Lag0", "mean.NO2.Lag2")]

# Parameters from the Aveiro Full model
p <- 2
q <- 6
intercept <- 0.07797
past_obs_coef <- c(0.216433, 0.149956)  # beta_1, beta_2
past_mean_coef <- c(-0.26131, 0.290165, 0.043232, -0.26497, 0.086309, 0.716672)  # alpha_1 to alpha_6
xreg_coef <- c(0.095653, -0.0073, 0.007343, 0.002856, -0.00503)  # Temp, DTemp, PM, NO2
sigmasq <- 0.02256  # Overdispersion parameter
size <- 1/sigmasq  # NB size parameter

# Storage for results
selected_p <- numeric(n_realizations)
selected_q <- numeric(n_realizations)
auto_time <- numeric(n_realizations)
grid_time <- numeric(n_realizations)
auto_models <- numeric(n_realizations)
grid_models <- numeric(n_realizations)

# Create model and parameter lists for simulation
model <- list(past_obs = c(1, 2), past_mean = c(1, 2, 3, 4, 5, 6), external = FALSE)
param <- list(
  intercept = intercept,
  past_obs = past_obs_coef,
  past_mean = past_mean_coef,
  xreg = xreg_coef
)

# Run simulation study
for (i in 1:n_realizations) {
  cat(sprintf("Running realization %d of %d\n", i, n_realizations))
  
  # Simulate INGARCH time series with covariates
  sim_data <- ingarch.sim(
    n = n_observations,
    param = param,
    model = model,
    xreg = as.matrix(covariates),
    link = "log",
    distr = "nbinom",
    size = size,
    n_start = n_start
  )
  
  # Model selection using auto.ingarch
  t1 <- Sys.time()
  auto_result <- auto.ingarch(
    y = sim_data$ts,
    max.p = 7,
    max.q = 7,
    distribution = "nbinom",
    link = "log",
    xreg = as.matrix(covariates),
    ic = "aic",
    stepwise = TRUE
  )
  t2 <- Sys.time()
  auto_time[i] <- as.numeric(difftime(t2, t1, units = "secs"))
  
  # Get selected orders from auto.ingarch
  selected_p[i] <- length(auto_result$model$past_obs)
  selected_q[i] <- length(auto_result$model$past_mean)
  auto_models[i] <- nrow(auto_result$results)
  
  # Model selection using grid search
  t1 <- Sys.time()
  grid_result <- auto.ingarch(
    y = sim_data$ts,
    max.p = 7,
    max.q = 7,
    distribution = "nbinom",
    link = "log",
    xreg = as.matrix(covariates),
    ic = "aic",
    stepwise = FALSE
  )
  t2 <- Sys.time()
  grid_time[i] <- as.numeric(difftime(t2, t1, units = "secs"))
  
  # Count models tested in grid search
  grid_models[i] <- sum(!is.na(grid_result$results[,1]))
  
  # Save intermediate results every 100 realizations
  if (i %% 100 == 0) {
    results_df <- data.frame(
      realization = 1:i,
      selected_p = selected_p[1:i],
      selected_q = selected_q[1:i],
      auto_time = auto_time[1:i],
      grid_time = grid_time[1:i],
      auto_models = auto_models[1:i],
      grid_models = grid_models[1:i]
    )
    saveRDS(results_df, file = "ingarch_covariates_simulation_results_interim.rds")
  }
}

# Summarize results
results_df <- data.frame(
  realization = 1:n_realizations,
  selected_p = selected_p,
  selected_q = selected_q,
  auto_time = auto_time,
  grid_time = grid_time,
  auto_models = auto_models,
  grid_models = grid_models
)

# Save full results
saveRDS(results_df, file = "ingarch_covariates_simulation_results.rds")

# Print summary statistics
cat("\n\nSummary of Results (INGARCH with covariates):\n")
cat("=========================================\n\n")

cat("Orders Selected (auto.ingarch):\n")
table_auto <- table(selected_p, selected_q)
print(table_auto)
cat("\n")

cat("Mean Computing Time (seconds):\n")
cat(sprintf("auto.ingarch: %.2f\n", mean(auto_time)))
cat(sprintf("grid.search: %.2f\n", mean(grid_time)))
cat("\n")

cat("Mean Number of Models Tested:\n")
cat(sprintf("auto.ingarch: %.2f\n", mean(auto_models)))
cat(sprintf("grid.search: %.2f\n", mean(grid_models)))

# Create plots
pdf("ingarch_covariates_simulation_results.pdf", width = 10, height = 8)

# Distribution of selected orders
par(mfrow = c(2, 2))
hist(selected_p, breaks = seq(0, 7, by = 1), main = "Distribution of selected p",
     xlab = "Order p", xlim = c(0, 7), col = "lightblue")
abline(v = p, col = "red", lwd = 2)
legend("topright", legend = c("True p"), col = "red", lwd = 2)

hist(selected_q, breaks = seq(0, 7, by = 1), main = "Distribution of selected q",
     xlab = "Order q", xlim = c(0, 7), col = "lightblue")
abline(v = q, col = "red", lwd = 2)
legend("topright", legend = c("True q"), col = "red", lwd = 2)

# Computing time comparison
boxplot(auto_time, grid_time, names = c("auto.ingarch", "grid.search"),
        main = "Computing Time Comparison", ylab = "Seconds", col = c("lightblue", "lightgreen"))

# Number of models tested
boxplot(auto_models, grid_models, names = c("auto.ingarch", "grid.search"),
        main = "Number of Models Tested", ylab = "Count", col = c("lightblue", "lightgreen"))

dev.off()