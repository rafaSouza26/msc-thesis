# INGARCH Model Simulation Study - Without Covariates
# Comparing model selection methods: auto.ingarch with stepwise vs grid search

# Load required packages
library(tscount)    # For INGARCH models
library(ggplot2)    # For visualization
library(dplyr)      # For data manipulation
library(readxl)     # For reading Excel files

# Source custom functions
source("./custom/auto.ingarch.R")
source("./custom/ingarch.sim.R")
source("./custom/newmodel.R")
source("./custom/ingarch.string.R")
source("./custom/search.ingarch.R")

# Function to extract model parameters from Excel file
extract_model_params <- function(data_path, model_row = 1) {
  # Load the Excel file
  excel_data <- readxl::read_excel(data_path)
  
  # Print Excel file contents for verification
  cat("Excel file columns:", colnames(excel_data), "\n")
  cat("Excel file rows:", nrow(excel_data), "\n")
  
  # Extract parameters from the specified row (default is first row - Modelo Ingarch_Aveiro)
  model_data <- excel_data[model_row, ]
  
  # Extract parameters
  p <- as.numeric(model_data$p)
  q <- as.numeric(model_data$q)
  sigmasq <- as.numeric(model_data$sigmasq)
  intercept <- as.numeric(model_data$Intercept)
  
  # Extract beta parameters (past observations coefficients)
  beta_cols <- grep("^beta", colnames(excel_data))
  betas <- as.numeric(unlist(model_data[, beta_cols]))
  
  # Extract alpha parameters (past means coefficients)
  alpha_cols <- grep("^alpha", colnames(excel_data))
  alphas <- as.numeric(unlist(model_data[, alpha_cols]))
  
  # Return the extracted parameters
  return(list(
    p = p,
    q = q,
    sigmasq = sigmasq,
    intercept = intercept,
    betas = betas,
    alphas = alphas
  ))
}

# Main simulation study - Without Covariates
run_simulation_study_no_covariates <- function() {
  set.seed(12345)  # For reproducibility
  
  # 1. Extract model parameters from Excel file
  cat("Extracting model parameters from Excel file...\n")
  params <- extract_model_params("./data/modelosAveiro.xlsx", model_row = 1)
  
  # Print model parameters for verification
  cat("Model parameters:\n")
  cat("p =", params$p, "\n")
  cat("q =", params$q, "\n")
  cat("sigmasq =", params$sigmasq, "\n")
  cat("intercept =", params$intercept, "\n")
  cat("beta coefficients:", params$betas, "\n")
  cat("alpha coefficients:", params$alphas, "\n")
  
  # Set up parameters for ingarch.sim
  ingarch_params <- list(
    intercept = params$intercept,
    past_obs = params$betas[1:params$p],
    past_mean = params$alphas[1:params$q]
  )
  
  # Set up model for ingarch.sim
  ingarch_model <- list(
    past_obs = 1:params$p,
    past_mean = 1:params$q,
    external = FALSE
  )
  
  # Calculate size parameter for negative binomial (size = 1/sigmasq)
  size_param <- 1/params$sigmasq
  
  # 2. Simulate INGARCH without covariates
  cat("Simulating 1000 INGARCH realizations without covariates...\n")
  sims <- vector("list", 1000)
  
  for(i in 1:1000) {
    if(i %% 100 == 0) cat("Generating simulation", i, "of 1000\n")
    
    # Get the simulation result
    sim_result <- ingarch.sim(
      n = 1000,                  # 1000 observations
      param = ingarch_params,    # Parameters from Excel file
      model = ingarch_model,     # Model specification
      link = "log",              # Log link
      distr = "nbinom",          # Negative Binomial distribution
      size = size_param,         # Size parameter
      n_start = 100              # Burn-in period
    )
    
    # Extract only the time series values from the simulation result
    # Determine how to extract based on the class and structure of sim_result
    if(i == 1) {
      cat("Class of simulation result:", class(sim_result), "\n")
      if(is.list(sim_result)) {
        cat("Names of sim_result components:", names(sim_result), "\n")
      }
    }
    
    # Store only the time series part, not the entire object
    if(is.numeric(sim_result)) {
      # If it's already a numeric vector, use it directly
      sims[[i]] <- as.numeric(sim_result)
    } else if(is.list(sim_result)) {
      # If it's a list, try to extract the time series component
      if("ts" %in% names(sim_result)) {
        sims[[i]] <- as.numeric(sim_result$ts)
      } else if("series" %in% names(sim_result)) {
        sims[[i]] <- as.numeric(sim_result$series)
      } else if("values" %in% names(sim_result)) {
        sims[[i]] <- as.numeric(sim_result$values)
      } else if(length(sim_result) > 0 && is.numeric(sim_result[[1]])) {
        # If the first element is numeric, use that
        sims[[i]] <- as.numeric(sim_result[[1]])
      } else {
        stop(paste("Cannot extract numeric time series from ingarch.sim output in simulation", i))
      }
    } else {
      stop(paste("Unexpected output type from ingarch.sim in simulation", i))
    }
    
    # Verify the extracted data
    if(i == 1) {
      cat("Class of extracted simulation data:", class(sims[[i]]), "\n")
      cat("Length of extracted simulation data:", length(sims[[i]]), "\n")
      cat("First few values:", head(sims[[i]]), "\n")
    }
  }
  
  # 3. Run model selection for each simulation
  cat("Running model selection for INGARCH without covariates...\n")
  
  # Initialize results storage
  results <- list(
    stepwise = vector("list", length(sims)),
    grid_search = vector("list", length(sims))
  )
  
  # Perform model selection
  for(i in 1:length(sims)) {
    if(i %% 10 == 0) cat("Processing simulation", i, "of", length(sims), "\n")
    
    # Get simulated data (now a numeric vector, not a list)
    sim_data <- sims[[i]]
    
    # Run auto.ingarch with stepwise=TRUE
    start_time_stepwise <- Sys.time()
    fit_stepwise <- tryCatch({
      auto.ingarch(
        y = sim_data,
        max.p = 7,
        max.q = 7,
        distribution = "nbinom",
        link = "log",
        ic = "aicc",
        stepwise = TRUE,
        trace = FALSE,
        show_warnings = FALSE
      )
    }, error = function(e) {
      cat("Error in stepwise fit for simulation", i, ":", conditionMessage(e), "\n")
      return(NULL)
    })
    end_time_stepwise <- Sys.time()
    
    # Record results for stepwise
    if(!is.null(fit_stepwise)) {
      p_stepwise <- if(is.null(fit_stepwise$model$past_obs)) 0 else length(fit_stepwise$model$past_obs)
      q_stepwise <- if(is.null(fit_stepwise$model$past_mean)) 0 else length(fit_stepwise$model$past_mean)
      n_models_stepwise <- if(is.null(fit_stepwise$results)) 0 else nrow(fit_stepwise$results)
      
      results$stepwise[[i]] <- list(
        p = p_stepwise,
        q = q_stepwise,
        time = as.numeric(difftime(end_time_stepwise, start_time_stepwise, units = "secs")),
        n_models = n_models_stepwise,
        aic = tryCatch(AIC(fit_stepwise), error = function(e) NA),
        bic = tryCatch(BIC(fit_stepwise), error = function(e) NA)
      )
    } else {
      results$stepwise[[i]] <- list(
        p = NA,
        q = NA,
        time = as.numeric(difftime(end_time_stepwise, start_time_stepwise, units = "secs")),
        n_models = NA,
        aic = NA,
        bic = NA
      )
    }
    
    # Run auto.ingarch with stepwise=FALSE (grid search)
    start_time_grid <- Sys.time()
    fit_grid <- tryCatch({
      auto.ingarch(
        y = sim_data,
        max.p = 7,
        max.q = 7,
        distribution = "nbinom",
        link = "log",
        ic = "aicc",
        stepwise = FALSE,
        trace = FALSE,
        show_warnings = FALSE
      )
    }, error = function(e) {
      cat("Error in grid search fit for simulation", i, ":", conditionMessage(e), "\n")
      return(NULL)
    })
    end_time_grid <- Sys.time()
    
    # Record results for grid search
    if(!is.null(fit_grid)) {
      p_grid <- if(is.null(fit_grid$model$past_obs)) 0 else length(fit_grid$model$past_obs)
      q_grid <- if(is.null(fit_grid$model$past_mean)) 0 else length(fit_grid$model$past_mean)
      
      # For grid search, number of models tested is always (max.p+1)*(max.q+1)
      n_models_grid <- (7+1) * (7+1)
      
      results$grid_search[[i]] <- list(
        p = p_grid,
        q = q_grid,
        time = as.numeric(difftime(end_time_grid, start_time_grid, units = "secs")),
        n_models = n_models_grid,
        aic = tryCatch(AIC(fit_grid), error = function(e) NA),
        bic = tryCatch(BIC(fit_grid), error = function(e) NA)
      )
    } else {
      results$grid_search[[i]] <- list(
        p = NA,
        q = NA,
        time = as.numeric(difftime(end_time_grid, start_time_grid, units = "secs")),
        n_models = NA,
        aic = NA,
        bic = NA
      )
    }
  }
  
  # 4. Summarize results
  cat("Summarizing results for INGARCH without covariates...\n")
  
  # Convert lists to data frames, handling possible NAs
  stepwise_results_df <- do.call(rbind, lapply(1:length(results$stepwise), function(i) {
    x <- results$stepwise[[i]]
    if(is.null(x)) return(NULL)
    data.frame(
      method = "stepwise",
      sim_id = i,
      p = x$p,
      q = x$q,
      time = x$time,
      n_models = x$n_models,
      aic = x$aic,
      bic = x$bic
    )
  }))
  
  grid_results_df <- do.call(rbind, lapply(1:length(results$grid_search), function(i) {
    x <- results$grid_search[[i]]
    if(is.null(x)) return(NULL)
    data.frame(
      method = "grid_search",
      sim_id = i,
      p = x$p,
      q = x$q,
      time = x$time,
      n_models = x$n_models,
      aic = x$aic,
      bic = x$bic
    )
  }))
  
  results_df <- rbind(stepwise_results_df, grid_results_df)
  
  # Save results
  saveRDS(sims, "ingarch_no_covariates_simulations.rds")
  saveRDS(results_df, "ingarch_no_covariates_results.rds")
  
  # Calculate and print summary statistics, handling NAs
  summary_stats <- results_df %>%
    group_by(method) %>%
    summarize(
      mean_p = mean(p, na.rm = TRUE),
      median_p = median(p, na.rm = TRUE),
      sd_p = sd(p, na.rm = TRUE),
      mean_q = mean(q, na.rm = TRUE),
      median_q = median(q, na.rm = TRUE),
      sd_q = sd(q, na.rm = TRUE),
      mean_time = mean(time, na.rm = TRUE),
      median_time = median(time, na.rm = TRUE),
      mean_models = mean(n_models, na.rm = TRUE),
      median_models = median(n_models, na.rm = TRUE),
      na_count = sum(is.na(p) | is.na(q))  # Count of failures
    )
  
  # Order selection frequencies
  order_freq <- results_df %>%
    filter(!is.na(p) & !is.na(q)) %>%  # Remove NAs
    group_by(method, p, q) %>%
    summarize(count = n()) %>%
    group_by(method) %>%
    mutate(freq = count / sum(count) * 100)
  
  # Save summary statistics
  write.csv(summary_stats, "ingarch_no_covariates_summary.csv", row.names = FALSE)
  write.csv(order_freq, "ingarch_no_covariates_order_freq.csv", row.names = FALSE)
  
  # Print summary
  cat("\n===== Summary Statistics =====\n")
  print(summary_stats)
  
  cat("\n===== Most Frequent (p,q) Orders =====\n")
  top_orders <- order_freq %>%
    arrange(method, desc(freq)) %>%
    group_by(method) %>%
    slice_head(n = 5)
  print(top_orders)
  
  cat("\nResults saved to: ingarch_no_covariates_results.rds\n")
  
  return(results_df)
}

# Run the study
main <- function() {
  results <- run_simulation_study_no_covariates()
  cat("Simulation study completed!\n")
}

# Execute main function
main()