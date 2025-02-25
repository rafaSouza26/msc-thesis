#!/usr/bin/env Rscript
# Extensive testing of auto.ingarch order identification accuracy
# This script runs multiple simulations for each model specification
# and calculates the percentage of times the algorithm correctly
# identifies the true model order

# Load required packages
library(tscount)

# Load custom functions
source("./custom/ingarch.sim.R")
source("./custom/myingarch.R")
source("./custom/auto.ingarch.R")
source("./custom/ingarch.string.R")
source("./custom/newmodel.R")

# Function to simulate data and test model identification
simulate_and_identify <- function(true_p, true_q, intercept, p_coef, q_coef, 
                                  n_obs = 500, link = "identity", distr = "poisson",
                                  max_p = 3, max_q = 3, verbose = FALSE) {
  
  # Set up the model parameters
  if (true_p == 0) {
    past_obs_param <- NULL
  } else {
    past_obs_param <- p_coef
  }
  
  if (true_q == 0) {
    past_mean_param <- NULL
  } else {
    past_mean_param <- q_coef
  }
  
  # Generate simulated data
  tryCatch({
    simulated_data <- ingarch.sim(
      n = n_obs,
      param = list(
        intercept = intercept,
        past_obs = past_obs_param,
        past_mean = past_mean_param
      ),
      model = list(
        past_obs = if(true_p > 0) 1:true_p else NULL,
        past_mean = if(true_q > 0) 1:true_q else NULL,
        external = FALSE
      ),
      link = link,
      distr = distr
    )
    
    # Convert to time series object
    ts_data <- ts(simulated_data$ts)
    
    # Fit model with auto.ingarch
    if (verbose) {
      cat("\nTesting INGARCH(", true_p, ",", true_q, ") model...\n", sep="")
    }
    
    # Attempt to fit the model
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
      if (verbose) {
        cat("Error fitting model:", conditionMessage(e), "\n")
      }
      return(NULL)
    })
    
    # Process results
    if (!is.null(fitted_model)) {
      # Extract estimated p and q
      est_p <- if(is.null(fitted_model$model$past_obs)) 0 else length(fitted_model$model$past_obs)
      est_q <- if(is.null(fitted_model$model$past_mean)) 0 else length(fitted_model$model$past_mean)
      
      # Check if orders were correctly identified
      order_correct <- (est_p == true_p) && (est_q == true_q)
      
      if (verbose) {
        cat("  Estimated model: INGARCH(", est_p, ",", est_q, ")\n", sep="")
        cat("  Order correctly identified:", ifelse(order_correct, "YES", "NO"), "\n")
      }
      
      # Return result
      return(list(
        success = TRUE,
        order_correct = order_correct,
        est_p = est_p,
        est_q = est_q
      ))
    } else {
      return(list(
        success = FALSE,
        order_correct = FALSE,
        est_p = NA,
        est_q = NA
      ))
    }
  }, error = function(e) {
    if (verbose) {
      cat("Error in simulation:", conditionMessage(e), "\n")
    }
    return(list(
      success = FALSE,
      order_correct = FALSE,
      est_p = NA,
      est_q = NA
    ))
  })
}

# Function to run multiple simulations for a given model configuration
run_monte_carlo <- function(true_p, true_q, intercept_range, p_coef_range, q_coef_range,
                            n_sims = 100, n_obs = 500, link = "identity", distr = "poisson",
                            max_p = 3, max_q = 3, seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Store results
  results <- data.frame(
    sim = 1:n_sims,
    true_p = true_p,
    true_q = true_q,
    est_p = NA,
    est_q = NA,
    success = FALSE,
    order_correct = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Progress bar
  pb <- txtProgressBar(min = 0, max = n_sims, style = 3)
  
  for (i in 1:n_sims) {
    # Generate random parameters within ranges
    intercept <- runif(1, intercept_range[1], intercept_range[2])
    
    if (true_p > 0) {
      p_min <- p_coef_range[1]
      p_max <- p_coef_range[2]
      p_coef <- runif(true_p, p_min, p_max)
    } else {
      p_coef <- NULL
    }
    
    if (true_q > 0) {
      q_min <- q_coef_range[1]
      q_max <- q_coef_range[2]
      q_coef <- runif(true_q, q_min, q_max)
    } else {
      q_coef <- NULL
    }
    
    # Run simulation
    sim_result <- simulate_and_identify(
      true_p = true_p,
      true_q = true_q,
      intercept = intercept,
      p_coef = p_coef,
      q_coef = q_coef,
      n_obs = n_obs,
      link = link,
      distr = distr,
      max_p = max_p,
      max_q = max_q,
      verbose = FALSE
    )
    
    # Store results
    results$success[i] <- sim_result$success
    results$order_correct[i] <- sim_result$order_correct
    results$est_p[i] <- sim_result$est_p
    results$est_q[i] <- sim_result$est_q
    
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  # Calculate summary statistics
  successful_runs <- sum(results$success)
  correct_identifications <- sum(results$order_correct)
  accuracy <- correct_identifications / n_sims
  
  # Print summary
  cat("\nINGARCH(", true_p, ",", true_q, ") Summary:\n", sep="")
  cat("  Total simulations:", n_sims, "\n")
  cat("  Successful runs:", successful_runs, "(", round(successful_runs/n_sims*100, 1), "%)\n", sep="")
  cat("  Correct identifications:", correct_identifications, "(", round(correct_identifications/n_sims*100, 1), "%)\n", sep="")
  
  # Calculate distribution of estimated orders
  est_orders <- table(results$est_p, results$est_q, useNA = "ifany")
  dimnames(est_orders) <- list(p = rownames(est_orders), q = colnames(est_orders))
  
  return(list(
    model = paste0("INGARCH(", true_p, ",", true_q, ")"),
    results = results,
    successful_runs = successful_runs,
    correct_identifications = correct_identifications,
    accuracy = accuracy,
    est_orders = est_orders
  ))
}

# Main testing function
extensive_ingarch_testing <- function(n_sims = 100, n_obs = 500, max_order = 3, seed = 42) {
  # Set seed for reproducibility
  set.seed(seed)
  
  # Define test cases
  test_cases <- list(
    # INGARCH(0,0) - constant mean model
    list(true_p = 0, true_q = 0, 
         intercept_range = c(2, 10), 
         p_coef_range = NULL, 
         q_coef_range = NULL),
    
    # INGARCH(1,0) - only AR effect
    list(true_p = 1, true_q = 0, 
         intercept_range = c(0.5, 3), 
         p_coef_range = c(0.3, 0.7), 
         q_coef_range = NULL),
    
    # INGARCH(0,1) - only MA effect
    list(true_p = 0, true_q = 1, 
         intercept_range = c(0.5, 3), 
         p_coef_range = NULL, 
         q_coef_range = c(0.3, 0.7)),
    
    # INGARCH(1,1) - both AR and MA effects
    list(true_p = 1, true_q = 1, 
         intercept_range = c(0.5, 2), 
         p_coef_range = c(0.2, 0.6), 
         q_coef_range = c(0.2, 0.6)),
    
    # INGARCH(2,0) - higher order AR
    list(true_p = 2, true_q = 0, 
         intercept_range = c(0.5, 2), 
         p_coef_range = c(0.2, 0.5), 
         q_coef_range = NULL),
    
    # INGARCH(0,2) - higher order MA
    list(true_p = 0, true_q = 2, 
         intercept_range = c(0.5, 2), 
         p_coef_range = NULL, 
         q_coef_range = c(0.2, 0.5)),
    
    # INGARCH(2,1) - higher order mixed
    list(true_p = 2, true_q = 1, 
         intercept_range = c(0.5, 2), 
         p_coef_range = c(0.2, 0.4), 
         q_coef_range = c(0.2, 0.4)),
    
    # INGARCH(1,2) - higher order mixed
    list(true_p = 1, true_q = 2, 
         intercept_range = c(0.5, 2), 
         p_coef_range = c(0.2, 0.4), 
         q_coef_range = c(0.2, 0.4))
  )
  
  # Start time
  start_time <- Sys.time()
  
  # Run all test cases
  cat(paste0("Starting extensive INGARCH testing with ", n_sims, " simulations per model\n"))
  cat("=====================================================================\n")
  
  # Create directory for results if it doesn't exist
  dir.create("ingarch_results", showWarnings = FALSE)
  
  # Run monte carlo simulations for each test case
  all_results <- list()
  for (i in 1:length(test_cases)) {
    cat(paste0("\nTest case ", i, " of ", length(test_cases), "\n"))
    
    tc <- test_cases[[i]]
    result <- run_monte_carlo(
      true_p = tc$true_p,
      true_q = tc$true_q,
      intercept_range = tc$intercept_range,
      p_coef_range = tc$p_coef_range,
      q_coef_range = tc$q_coef_range,
      n_sims = n_sims,
      n_obs = n_obs,
      max_p = max_order,
      max_q = max_order,
      seed = seed + i
    )
    
    all_results[[i]] <- result
    
    # Save individual test case results
    model_name <- paste0("INGARCH(", tc$true_p, ",", tc$true_q, ")")
    file_name <- file.path("ingarch_results", paste0("model_", tc$true_p, "_", tc$true_q, ".rds"))
    saveRDS(result, file = file_name)
    
    cat("\nModel:", model_name, "- Results saved to", file_name, "\n")
    cat("Estimated order distribution:\n")
    print(result$est_orders)
    cat("\n")
  }
  
  # End time
  end_time <- Sys.time()
  time_taken <- end_time - start_time
  
  # Summarize overall results
  accuracy_vector <- sapply(all_results, function(x) x$accuracy)
  models <- sapply(all_results, function(x) x$model)
  
  summary_df <- data.frame(
    model = models,
    accuracy = accuracy_vector,
    success_rate = sapply(all_results, function(x) x$successful_runs / n_sims),
    stringsAsFactors = FALSE
  )
  
  # Print overall summary
  cat("\n=====================================================================\n")
  cat("OVERALL SUMMARY - INGARCH ORDER IDENTIFICATION PERFORMANCE\n")
  cat("=====================================================================\n")
  cat("Total time taken:", format(time_taken), "\n")
  cat("Average accuracy across all models:", round(mean(accuracy_vector) * 100, 2), "%\n\n")
  
  # Print model-specific accuracy
  cat("Accuracy by model:\n")
  for (i in 1:length(all_results)) {
    cat(sprintf("  %-12s: %6.2f%%\n", summary_df$model[i], summary_df$accuracy[i] * 100))
  }
  
  # Generate plot of accuracy by model
  pdf(file.path("ingarch_results", "accuracy_by_model.pdf"), width = 10, height = 6)
  par(mar = c(8, 4, 4, 2) + 0.1)
  barplot(summary_df$accuracy * 100, names.arg = summary_df$model, 
          main = "INGARCH Order Identification Accuracy by Model",
          ylab = "Accuracy (%)", col = "skyblue",
          ylim = c(0, 100), las = 2)
  abline(h = seq(0, 100, by = 10), lty = 3, col = "gray")
  abline(h = mean(accuracy_vector) * 100, lwd = 2, col = "red")
  legend("topright", legend = "Average", lwd = 2, col = "red")
  dev.off()
  
  # Generate confusion matrix-style visualization for each model
  for (i in 1:length(all_results)) {
    model_name <- all_results[[i]]$model
    est_orders <- all_results[[i]]$est_orders
    
    # Convert to dataframe for easier plotting
    est_matrix <- as.data.frame.table(est_orders)
    colnames(est_matrix) <- c("est_p", "est_q", "frequency")
    
    # Convert to numeric
    est_matrix$est_p <- as.numeric(as.character(est_matrix$est_p))
    est_matrix$est_q <- as.numeric(as.character(est_matrix$est_q))
    
    # Calculate percentages
    est_matrix$percentage <- est_matrix$frequency / sum(est_matrix$frequency) * 100
    
    # Create heatmap
    pdf(file.path("ingarch_results", paste0("heatmap_", gsub("[()]", "_", model_name), ".pdf")), 
        width = 7, height = 6)
    
    # Get dimensions for the plot
    p_values <- sort(unique(est_matrix$est_p))
    q_values <- sort(unique(est_matrix$est_q))
    
    # Create a matrix of percentages
    heat_matrix <- matrix(0, nrow = length(p_values), ncol = length(q_values))
    rownames(heat_matrix) <- p_values
    colnames(heat_matrix) <- q_values
    
    for (j in 1:nrow(est_matrix)) {
      p_idx <- which(p_values == est_matrix$est_p[j])
      q_idx <- which(q_values == est_matrix$est_q[j])
      heat_matrix[p_idx, q_idx] <- est_matrix$percentage[j]
    }
    
    # Plot heatmap
    layout(matrix(c(1,2), ncol=2), widths=c(4,1))
    par(mar=c(5,4,4,1))
    
    # Extract true p and q
    true_p <- as.numeric(gsub("INGARCH\\((\\d+),(\\d+)\\)", "\\1", model_name))
    true_q <- as.numeric(gsub("INGARCH\\((\\d+),(\\d+)\\)", "\\2", model_name))
    
    image(1:length(p_values), 1:length(q_values), t(heat_matrix), 
          col = colorRampPalette(c("white", "orange", "red"))(100),
          xlab = "Estimated p", ylab = "Estimated q", axes = FALSE,
          main = paste(model_name, "Order Estimation"))
    
    axis(1, at = 1:length(p_values), labels = p_values)
    axis(2, at = 1:length(q_values), labels = q_values)
    
    # Add grid lines
    grid(length(p_values), length(q_values), lty = 1, col = "gray")
    
    # Highlight the true order
    true_p_idx <- which(p_values == true_p)
    true_q_idx <- which(q_values == true_q)
    
    if (length(true_p_idx) > 0 && length(true_q_idx) > 0) {
      points(true_p_idx, true_q_idx, pch = 16, cex = 2)
      points(true_p_idx, true_q_idx, pch = 1, cex = 2, lwd = 2)
    }
    
    # Add percentages as text
    for (p_idx in 1:length(p_values)) {
      for (q_idx in 1:length(q_values)) {
        if (heat_matrix[p_idx, q_idx] > 0) {
          text(p_idx, q_idx, sprintf("%.1f%%", heat_matrix[p_idx, q_idx]),
               col = ifelse(heat_matrix[p_idx, q_idx] > 50, "white", "black"))
        }
      }
    }
    
    # Add colorbar
    par(mar=c(5,1,4,3))
    colorbar <- as.matrix(seq(0, 100, length.out = 100))
    image(1, seq(0, 100, length.out = 100), colorbar, col = colorRampPalette(c("white", "orange", "red"))(100),
          axes = FALSE, xlab = "", ylab = "")
    axis(4, at = seq(0, 100, length.out = 6), labels = sprintf("%.0f%%", seq(0, 100, length.out = 6)))
    box()
    
    dev.off()
  }
  
  # Save all results
  saveRDS(all_results, file = file.path("ingarch_results", "all_results.rds"))
  saveRDS(summary_df, file = file.path("ingarch_results", "summary.rds"))
  
  # Return summary dataframe
  return(summary_df)
}

# Function to load and analyze previously saved results
analyze_results <- function(results_dir = "ingarch_results") {
  # Check if directory exists
  if (!dir.exists(results_dir)) {
    stop("Results directory does not exist")
  }
  
  # Load summary results
  summary_file <- file.path(results_dir, "summary.rds")
  if (file.exists(summary_file)) {
    summary_df <- readRDS(summary_file)
    
    # Print summary
    cat("\n=====================================================================\n")
    cat("INGARCH ORDER IDENTIFICATION PERFORMANCE SUMMARY\n")
    cat("=====================================================================\n")
    cat("Average accuracy across all models:", round(mean(summary_df$accuracy) * 100, 2), "%\n\n")
    
    # Print model-specific accuracy
    cat("Accuracy by model:\n")
    for (i in 1:nrow(summary_df)) {
      cat(sprintf("  %-12s: %6.2f%%\n", summary_df$model[i], summary_df$accuracy[i] * 100))
    }
    
    return(summary_df)
  } else {
    stop("Summary file does not exist")
  }
}

# Run tests if this script is executed directly
if (!interactive()) {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  n_sims <- 100  # Default number of simulations per model
  n_obs <- 500   # Default number of observations per simulation
  max_order <- 3 # Default maximum order to consider
  seed <- 42     # Default random seed
  
  if (length(args) >= 1) n_sims <- as.numeric(args[1])
  if (length(args) >= 2) n_obs <- as.numeric(args[2])
  if (length(args) >= 3) max_order <- as.numeric(args[3])
  if (length(args) >= 4) seed <- as.numeric(args[4])
  
  # Run the tests
  results <- extensive_ingarch_testing(n_sims, n_obs, max_order, seed)
} else {
  # If being sourced, provide instructions
  cat("Functions available:\n")
  cat("  extensive_ingarch_testing(n_sims = 100, n_obs = 500, max_order = 3, seed = 42)\n")
  cat("  analyze_results(results_dir = 'ingarch_results')\n")
  cat("\nExample usage:\n")
  cat("  results <- extensive_ingarch_testing(n_sims = 50)\n")
  cat("  # Later, analyze the results:\n")
  cat("  summary <- analyze_results()\n")
}