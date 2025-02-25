# Script to test if auto.ingarch correctly identifies the order of INGARCH models
# Load required packages
library(tscount)

# Load our custom functions
source("./custom/ingarch.sim.R")
source("./custom/myingarch.R")
source("./custom/auto.ingarch.R")
source("./custom/ingarch.string.R")
source("./custom/newmodel.R")

# Set a random seed for reproducibility
set.seed(456)

# Function to run a single test with specified parameters
test_model_identification <- function(true_p, true_q, intercept, p_coef, q_coef, 
                                      n_obs = 500, link = "identity", distr = "poisson") {
  # Set up the model parameters
  if (true_p == 0) {
    past_obs_param <- NULL
    past_obs_model <- NULL
  } else {
    past_obs_param <- p_coef
    past_obs_model <- 1:true_p
  }
  
  if (true_q == 0) {
    past_mean_param <- NULL
    past_mean_model <- NULL
  } else {
    past_mean_param <- q_coef
    past_mean_model <- 1:true_q
  }
  
  # Generate simulated data
  simulated_data <- tryCatch({
    ingarch.sim(
      n = n_obs,
      param = list(
        intercept = intercept,
        past_obs = past_obs_param,
        past_mean = past_mean_param
      ),
      model = list(
        past_obs = past_obs_model,
        past_mean = past_mean_model,
        external = FALSE
      ),
      link = link,
      distr = distr
    )
  }, error = function(e) {
    cat("Error in simulation:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (is.null(simulated_data)) {
    return(NULL)
  }
  
  # Convert to time series object for model fitting
  ts_data <- ts(simulated_data$ts)
  
  cat("\nTesting INGARCH(", true_p, ",", true_q, ") model identification...\n", sep="")
  cat("True parameters:\n")
  cat("  Intercept:", intercept, "\n")
  if (true_p > 0) {
    for (i in 1:length(p_coef)) {
      cat("  Past observation (p=", i, "):", p_coef[i], "\n", sep="")
    }
  }
  if (true_q > 0) {
    for (i in 1:length(q_coef)) {
      cat("  Past mean (q=", i, "):", q_coef[i], "\n", sep="")
    }
  }
  
  # Set max_p and max_q higher than the true values
  max_p <- max(true_p + 2, 3)
  max_q <- max(true_q + 2, 3)
  
  # Attempt to fit the model
  start_time <- Sys.time()
  fitted_model <- tryCatch({
    auto.ingarch(
      y = ts_data,
      max.p = max_p,
      max.q = max_q,
      distribution = distr,
      link = link,
      trace = FALSE
    )
  }, error = function(e) {
    cat("Error in model fitting:", conditionMessage(e), "\n")
    return(NULL)
  })
  end_time <- Sys.time()
  
  if (is.null(fitted_model)) {
    return(NULL)
  }
  
  # Extract estimated p and q
  est_p <- if(is.null(fitted_model$model$past_obs)) 0 else length(fitted_model$model$past_obs)
  est_q <- if(is.null(fitted_model$model$past_mean)) 0 else length(fitted_model$model$past_mean)
  
  # Check if orders were correctly identified
  order_correct <- (est_p == true_p) && (est_q == true_q)
  
  cat("\nResults:\n")
  cat("  Estimated model:", ingarch.string(fitted_model), "\n")
  cat("  Order correctly identified:", ifelse(order_correct, "YES", "NO"), "\n")
  
  if (order_correct) {
    cat("  Estimated parameters:\n")
    cat("    Intercept:", fitted_model$coefficients[1], "\n")
    
    if (est_p > 0) {
      past_obs_idx <- 2:(est_p+1)
      for (i in 1:est_p) {
        cat("    Past observation (p=", i, "):", 
            fitted_model$coefficients[past_obs_idx[i]], "\n", sep="")
      }
    }
    
    if (est_q > 0) {
      past_mean_idx <- (est_p+2):(est_p+est_q+1)
      for (i in 1:est_q) {
        cat("    Past mean (q=", i, "):", 
            fitted_model$coefficients[past_mean_idx[i]], "\n", sep="")
      }
    }
  }
  
  cat("  Time taken:", format(end_time - start_time, digits=2), "\n")
  
  return(list(
    true_p = true_p,
    true_q = true_q,
    est_p = est_p,
    est_q = est_q,
    order_correct = order_correct,
    fitted_model = fitted_model,
    time_taken = end_time - start_time
  ))
}

# Run a series of tests with different model orders
cat("\n====== TESTING INGARCH MODEL ORDER IDENTIFICATION ======\n")

# Test cases to run
test_cases <- list(
  # INGARCH(0,0) - constant mean model
  list(true_p = 0, true_q = 0, intercept = 5, p_coef = NULL, q_coef = NULL),
  
  # INGARCH(1,0) - only AR effect
  list(true_p = 1, true_q = 0, intercept = 1.2, p_coef = 0.6, q_coef = NULL),
  
  # INGARCH(0,1) - only MA effect
  list(true_p = 0, true_q = 1, intercept = 1.5, p_coef = NULL, q_coef = 0.7),
  
  # INGARCH(1,1) - both AR and MA effects
  list(true_p = 1, true_q = 1, intercept = 1.2, p_coef = 0.5, q_coef = 0.4),
  
  # INGARCH(2,1) - higher order AR
  list(true_p = 2, true_q = 1, intercept = 1.0, p_coef = c(0.4, 0.3), q_coef = 0.2),
  
  # INGARCH(1,2) - higher order MA
  list(true_p = 1, true_q = 2, intercept = 1.0, p_coef = 0.3, q_coef = c(0.3, 0.2))
)

# Run all test cases
results <- list()
for (i in 1:length(test_cases)) {
  with(test_cases[[i]], {
    cat("\n----------------------------------------------------------------\n")
    cat("Test case", i, "of", length(test_cases), "\n")
    results[[i]] <- test_model_identification(true_p, true_q, intercept, p_coef, q_coef)
  })
}

# Get valid results (not NULL)
valid_results <- results[!sapply(results, is.null)]

# Summarize overall results
cat("\n====== SUMMARY OF RESULTS ======\n")
successful_tests <- length(valid_results)
cat("Successfully ran", successful_tests, "of", length(test_cases), "tests\n")

correct_identifications <- sum(sapply(valid_results, function(x) x$order_correct))
cat("Correctly identified model order in", correct_identifications, "of", successful_tests, "tests\n")

# Print detailed summary table
cat("\nDetailed results:\n")
cat("----------------------------------------------------------------------\n")
cat(sprintf("%-12s %-12s %-12s %-12s %s\n", 
            "True Model", "Est. Model", "Correct?", "Time (s)", "Notes"))
cat("----------------------------------------------------------------------\n")

for (i in 1:length(test_cases)) {
  if (!is.null(results[[i]])) {
    r <- results[[i]]
    true_model <- sprintf("INGARCH(%d,%d)", r$true_p, r$true_q)
    est_model <- sprintf("INGARCH(%d,%d)", r$est_p, r$est_q)
    correct <- ifelse(r$order_correct, "YES", "NO")
    time_taken <- round(as.numeric(r$time_taken, units="secs"), 2)
    
    cat(sprintf("%-12s %-12s %-12s %-12.2f %s\n", 
                true_model, est_model, correct, time_taken, ""))
  } else {
    true_p <- test_cases[[i]]$true_p
    true_q <- test_cases[[i]]$true_q
    true_model <- sprintf("INGARCH(%d,%d)", true_p, true_q)
    
    cat(sprintf("%-12s %-12s %-12s %-12s %s\n", 
                true_model, "ERROR", "N/A", "N/A", "Fitting failed"))
  }
}

# Find a successful model to plot (if any)
plot_index <- NA
for (i in 1:length(results)) {
  if (!is.null(results[[i]]) && results[[i]]$order_correct) {
    plot_index <- i
    break
  }
}

# Create a plot comparing one of the successful models (if any)
if (!is.na(plot_index)) {
  r <- results[[plot_index]]
  model <- r$fitted_model
  ts_data <- model$ts
  
  pdf("ingarch_model_fit.pdf", width=10, height=8)
  par(mfrow=c(2,1))
  
  # Plot the time series
  plot(ts_data, 
       type = "l",
       main = paste("Actual vs Fitted Values -", ingarch.string(model)),
       ylab = "Count",
       xlab = "Time")
  lines(fitted(model), col = "red")
  legend("topright", 
         legend = c("Actual", "Fitted"),
         col = c("black", "red"),
         lty = 1)
  
  # Plot residuals
  plot(residuals(model),
       type = "h",
       main = "Model Residuals",
       ylab = "Residual",
       xlab = "Time")
  abline(h = 0, col = "red", lty = 2)
  dev.off()
  
  cat("\nPlot saved to 'ingarch_model_fit.pdf'\n")
}