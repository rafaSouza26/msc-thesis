# INGARCH Model Simulation Study - With Covariates
# Comparing model selection methods: auto.ingarch with stepwise vs grid search

# Load required packages
# install.packages(c("tscount", "ggplot2", "dplyr", "readxl", "doParallel", "foreach", "progressr")) # Uncomment to install if needed
library(tscount)   # For INGARCH models
library(ggplot2)   # For visualization
library(dplyr)     # For data manipulation
library(readxl)    # For reading Excel files
library(doParallel)# For parallel processing
library(foreach)   # For parallel loops
library(progressr) # For progress reporting with foreach

# Source custom functions
# Ensure these paths are correct relative to your script's location
# Make sure these functions are available in the global environment
# or adjust the .export argument in the foreach loop accordingly.
# --- IMPORTANT: Ensure ALL functions called by auto.ingarch are sourced here ---
source("./custom/auto.ingarch.R")
source("./custom/ingarch.sim.R") # Ensure this uses xreg=param$external internally now
source("./custom/newmodel.R")     # Assuming myingarch might be defined here or in search.ingarch?
source("./custom/ingarch.string.R")
source("./custom/search.ingarch.R") # Assuming myingarch might be defined here or in newmodel?

# Function to extract model parameters (including covariates) from Excel file
# Minimal printing inside this function now
extract_model_params_with_covariates <- function(data_path, model_row = 2,
                                                 covariate_cols = c("Temp", "DTemp", "PM", "NO2")) {
  # Load the Excel file
  excel_data <- readxl::read_excel(data_path)
  
  # Print Excel file contents for verification (matches original script)
  cat("Excel file columns:", paste(colnames(excel_data), collapse=", "), "\n")
  cat("Excel file rows:", nrow(excel_data), "\n")
  
  # Extract parameters from the specified row
  if (model_row > nrow(excel_data)) {
    stop("model_row index out of bounds for the loaded Excel file.")
  }
  model_data <- excel_data[model_row, ]
  # Removed 'Extracting parameters from row X' print
  
  # Basic parameters
  p <- as.numeric(model_data$p)
  q <- as.numeric(model_data$q)
  sigmasq <- as.numeric(model_data$sigmasq)
  intercept <- as.numeric(model_data$Intercept)
  
  # Extract beta parameters
  beta_cols <- grep("^beta", colnames(excel_data), value = TRUE)
  if(length(beta_cols) >= p && p > 0) {
    betas <- as.numeric(unlist(model_data[, beta_cols[1:p]]))
  } else if (p == 0) {
    betas <- NULL
  } else {
    warning(paste("Found fewer beta columns (", length(beta_cols), ") than specified by p=", p, ". Check Excel file."))
    betas <- as.numeric(unlist(model_data[, beta_cols])) # Take what's available
  }
  
  # Extract alpha parameters
  alpha_cols <- grep("^alpha", colnames(excel_data), value = TRUE)
  if(length(alpha_cols) >= q && q > 0) {
    alphas <- as.numeric(unlist(model_data[, alpha_cols[1:q]]))
  } else if (q == 0) {
    alphas <- NULL
  } else {
    warning(paste("Found fewer alpha columns (", length(alpha_cols), ") than specified by q=", q, ". Check Excel file."))
    alphas <- as.numeric(unlist(model_data[, alpha_cols])) # Take what's available
  }
  
  # Extract external covariate coefficients
  external_coefs <- numeric(length(covariate_cols))
  names(external_coefs) <- covariate_cols
  missing_cov_cols <- c()
  for(col_name in covariate_cols) {
    if(col_name %in% colnames(model_data)) {
      external_coefs[col_name] <- as.numeric(model_data[[col_name]])
    } else {
      warning(paste("Covariate coefficient column '", col_name, "' not found in Excel file.", sep=""))
      external_coefs[col_name] <- NA
      missing_cov_cols <- c(missing_cov_cols, col_name)
    }
  }
  if(any(is.na(external_coefs))) {
    stop("Could not extract all required covariate coefficients: ", paste(missing_cov_cols, collapse=", "))
  }
  
  # Return the extracted parameters (no printing here)
  return(list(
    p = p, q = q, sigmasq = sigmasq, intercept = intercept,
    betas = betas, alphas = alphas, external_coefs = external_coefs
  ))
}

# Main simulation study - With Covariates and Parallel Model Selection
run_simulation_study_with_covariates_parallel <- function() {
  set.seed(54321) # Use a different seed for reproducibility
  
  # --- 1. Setup Covariates ---
  # Minimal printing during setup
  cat("Loading covariate data...\n") # Simplified message
  covariate_data_path <- "./data/count_covariates_data.RData"
  if (!file.exists(covariate_data_path)) {
    stop("Covariate data file not found: ", covariate_data_path)
  }
  load(covariate_data_path)
  if (!exists("d_1_data")) {
    stop("Data frame object 'd_1_data' not found within the loaded RData file: ", covariate_data_path)
  }
  if (!is.data.frame(d_1_data)) {
    stop("Object 'd_1_data' loaded from RData file is not a data frame.")
  }
  # Removed 'loaded successfully' print
  
  cov_names_r <- c("mean.Temp.Lag5", "mean.DTemp.Lag0", "PM25.Lag0", "mean.NO2.Lag2")
  if (!all(cov_names_r %in% colnames(d_1_data))) {
    missing_cols <- cov_names_r[!cov_names_r %in% colnames(d_1_data)]
    stop("Missing required covariate columns in d_1_data: ", paste(missing_cols, collapse=", "))
  }
  n_obs_sim <- 10
  if (nrow(d_1_data) < n_obs_sim) {
    stop(paste("Covariate data frame 'd_1_data' has only", nrow(d_1_data), "rows, but", n_obs_sim, "are required."))
  }
  xreg_matrix <- as.matrix(d_1_data[1:n_obs_sim, cov_names_r])
  # Removed 'Covariate matrix dimensions' print
  if(anyNA(xreg_matrix)) {
    warning("NA values found in the selected covariate matrix (first ", n_obs_sim, " rows). Check data or handle NAs.", immediate. = TRUE)
  }
  
  # --- 2. Extract Model Parameters ---
  cat("Extracting model parameters from Excel file...\n") # Matches original script
  excel_cov_names <- c("Temp", "DTemp", "PM", "NO2")
  excel_file_path <- "./data/modelosAveiro.xlsx"
  if (!file.exists(excel_file_path)) {
    stop("Excel parameter file not found: ", excel_file_path)
  }
  params <- extract_model_params_with_covariates(
    data_path = excel_file_path,
    model_row = 2,
    covariate_cols = excel_cov_names
  )
  
  # Print model parameters for verification (matches original format + external)
  cat("Model parameters:\n")
  cat("p =", params$p, "\n")
  cat("q =", params$q, "\n")
  cat("sigmasq =", params$sigmasq, "\n")
  cat("intercept =", params$intercept, "\n")
  cat("beta coefficients:", if(is.null(params$betas)) "None" else paste(params$betas, collapse=" "), "\n")
  cat("alpha coefficients:", if(is.null(params$alphas)) "None" else paste(params$alphas, collapse=" "), "\n")
  cat("external coefficients:", paste(names(params$external_coefs), params$external_coefs, sep="=", collapse=", "), "\n") # Keep this extra info
  
  # Setup simulation parameters
  ingarch_params <- list(
    intercept = params$intercept, past_obs = params$betas,
    past_mean = params$alphas, external = params$external_coefs # Name 'external' here is fine
  )
  ingarch_model <- list(
    past_obs = if(params$p > 0) 1:params$p else NULL,
    past_mean = if(params$q > 0) 1:params$q else NULL,
    external = TRUE
  )
  if (params$sigmasq <= 0) {
    stop("Invalid sigmasq value extracted from Excel (must be > 0): ", params$sigmasq)
  }
  size_param <- 1/params$sigmasq
  # Removed 'Negative Binomial size parameter' print
  
  # --- Ensure names consistency for external regressors (Optional but recommended) ---
  if (!is.null(ingarch_params$external) && !is.null(colnames(xreg_matrix))) {
    if (!identical(names(ingarch_params$external), colnames(xreg_matrix))) {
      warning("Names of extracted external coefficients do not match column names of xreg_matrix. Forcing names based on order.", immediate. = TRUE)
      names(ingarch_params$external) <- colnames(xreg_matrix)
      # Removed print statement showing corrected names
    }
  } else if (!is.null(ingarch_params$external) && is.null(colnames(xreg_matrix))) {
    warning("External coefficients provided but xreg_matrix has no column names.", immediate. = TRUE)
  }
  
  
  # --- 3. Simulate INGARCH with covariates ---
  n_sims <- 1 # Set desired number of simulations (e.g., 100 for testing, 1000 for full run)
  cat(paste("Simulating", n_sims, "INGARCH realizations with covariates...\n")) # Adapted message
  sims <- vector("list", n_sims)
  
  for(i in 1:n_sims) {
    # Simulation progress message (matches original frequency)
    if(i %% 100 == 0) cat("Generating simulation", i, "of", n_sims, "\n")
    
    sim_result <- tryCatch({
      # Call ingarch.sim - make sure it internally passes 'xreg=param$external' to tsglm.sim
      ingarch.sim(
        n = n_obs_sim, param = ingarch_params, model = ingarch_model,
        xreg = xreg_matrix, link = "log", distr = "nbinom",
        size = size_param, n_start = 100
      )
    }, error = function(e) {
      # Minimal error print during simulation generation
      cat("Error during ingarch.sim for simulation", i, ":", conditionMessage(e), "\n")
      return(NULL)
    })
    
    if(is.null(sim_result)) {
      sims[[i]] <- NULL
      next
    }
    
    # Print details only for the first simulation (matches original script)
    if(i == 1) {
      cat("Class of simulation result:", class(sim_result), "\n")
      if(is.list(sim_result)) {
        cat("Names of sim_result components:", names(sim_result), "\n")
      }
    }
    
    # Extract time series (using robust logic)
    extracted_series <- NULL
    if(is.numeric(sim_result)) {
      extracted_series <- as.numeric(sim_result)
    } else if(inherits(sim_result, "tsglm.sim") || is.list(sim_result)) {
      # Prioritize 'ts' if available, then 'response', etc.
      if("ts" %in% names(sim_result)) extracted_series <- as.numeric(sim_result$ts)
      else if("response" %in% names(sim_result)) extracted_series <- as.numeric(sim_result$response)
      else if("series" %in% names(sim_result)) extracted_series <- as.numeric(sim_result$series)
      else if("values" %in% names(sim_result)) extracted_series <- as.numeric(sim_result$values)
      else if(length(sim_result) > 0 && is.numeric(sim_result[[1]])) extracted_series <- as.numeric(sim_result[[1]])
    }
    
    if(is.null(extracted_series) || length(extracted_series) != n_obs_sim) {
      warning(paste("Issue extracting valid series for simulation", i, "- Storing NULL."), immediate. = TRUE)
      sims[[i]] <- NULL
    } else {
      sims[[i]] <- extracted_series
    }
    
    # Print details only for the first *successful* extraction (matches original script)
    if(i == 1 && !is.null(sims[[i]])) {
      cat("Class of extracted simulation data:", class(sims[[i]]), "\n")
      cat("Length of extracted simulation data:", length(sims[[i]]), "\n")
      cat("First few values:", head(sims[[i]]), "\n")
    }
  } # End simulation generation loop
  
  
  # --- 4. Run model selection for each simulation (PARALLELIZED) ---
  cat("Running PARALLEL model selection for INGARCH with covariates...\n") # Adapted message
  
  # Define maximum orders
  max_order_p <- 7
  max_order_q <- 7
  
  # --- Setup Parallel Backend ---
  num_cores <- detectCores(logical = FALSE)
  print(paste("Detected cores:", num_cores))
  cores_to_use <- max(1, num_cores - 1) # Use n-1 cores
  print(paste("Registering parallel backend with", cores_to_use, "cores."))
  cl <- makeCluster(cores_to_use)
  registerDoParallel(cl)
  print(paste("Workers registered:", getDoParWorkers()))
  
  # --- Setup Progress Reporting ---
  handlers(global = TRUE) # Enable progressr handlers (e.g., text bar)
  # Define the progressor based on the number of simulations
  p <- progressr::progressor(along = 1:length(sims))
  
  # --- Parallel Model Selection Loop ---
  start_parallel_time <- Sys.time()
  
  # Define custom functions to export if they are not in loaded packages
  # Ensure these functions exist in the current environment before the loop
  custom_funcs_to_export <- c("auto.ingarch", "ingarch.sim", "newmodel",
                              "ingarch.string", "search.ingarch", "myingarch")
  # Filter out any that might not exist to avoid errors during export
  custom_funcs_to_export <- custom_funcs_to_export[sapply(custom_funcs_to_export, exists, envir = .GlobalEnv)]
  if(!("myingarch" %in% custom_funcs_to_export)) {
    warning("Function 'myingarch' was not found in the global environment for export!", immediate. = TRUE)
  }
  if(length(custom_funcs_to_export) < 6) { # Adjust count if more functions expected
    warning("One or more expected custom functions were not found for export.", immediate. = TRUE)
  }
  
  # Wrap the foreach call in with_progress for the progress bar
  with_progress({
    parallel_results_list <- foreach(
      i = 1:length(sims),
      .packages = c("tscount", "stats"), # Add stats explicitly for AIC/BIC, although often loaded
      .export = c("xreg_matrix", "max_order_p", "max_order_q", custom_funcs_to_export), # Export necessary variables and custom functions
      .errorhandling = 'pass' # Return errors instead of stopping
    ) %dopar% {
      
      # Get the data for the current simulation
      sim_data <- sims[[i]]
      
      # Initialize results for this iteration
      iter_stepwise_result <- list(p = NA, q = NA, time = 0, n_models = NA, aic = NA, bic = NA, status = "Not Run")
      iter_grid_result <- list(p = NA, q = NA, time = 0, n_models = NA, aic = NA, bic = NA, status = "Not Run")
      
      # --- Pre-computation Checks ---
      if (is.null(sim_data)) {
        iter_stepwise_result$status <- "Input Sim Failed"
        iter_grid_result$status <- "Input Sim Failed"
        # --- Signal progress update even if skipping ---
        p(sprintf("i=%d (skipped)", i))
        return(list(stepwise = iter_stepwise_result, grid_search = iter_grid_result))
      }
      
      if (length(sim_data) != nrow(xreg_matrix)) {
        iter_stepwise_result$status <- "Input Length Mismatch"
        iter_grid_result$status <- "Input Length Mismatch"
        # --- Signal progress update even if skipping ---
        p(sprintf("i=%d (skipped)", i))
        return(list(stepwise = iter_stepwise_result, grid_search = iter_grid_result))
      }
      
      # --- Stepwise ---
      fit_status_stepwise <- "Fit Failed" # Default status
      start_time_stepwise <- Sys.time()
      fit_stepwise <- tryCatch({
        model <- auto.ingarch(y = sim_data, xreg = xreg_matrix, max.p = max_order_p, max.q = max_order_q,
                              distribution = "nbinom", link = "log", ic = "aicc",
                              stepwise = TRUE, trace = FALSE, show_warnings = FALSE)
        fit_status_stepwise <- "Success"
        model # Return the model on success
      }, error = function(e) {
        error_msg <- conditionMessage(e)
        fit_status_stepwise <<- paste("Error:", error_msg)
        return(NULL) # Return NULL on error
      })
      end_time_stepwise <- Sys.time()
      time_stepwise <- as.numeric(difftime(end_time_stepwise, start_time_stepwise, units = "secs"))
      
      # Record stepwise results for this iteration
      iter_stepwise_result$time <- time_stepwise
      iter_stepwise_result$status <- fit_status_stepwise # Assign status (Success or Error: message)
      if (!is.null(fit_stepwise) && inherits(fit_stepwise, "tsglm")) {
        p_stepwise <- if (is.null(fit_stepwise$model$past_obs)) 0 else length(fit_stepwise$model$past_obs)
        q_stepwise <- if (is.null(fit_stepwise$model$past_mean)) 0 else length(fit_stepwise$model$past_mean)
        n_models_stepwise <- if (is.null(fit_stepwise$results)) NA else nrow(fit_stepwise$results)
        iter_stepwise_result$p <- p_stepwise
        iter_stepwise_result$q <- q_stepwise
        iter_stepwise_result$n_models <- n_models_stepwise
        iter_stepwise_result$aic <- tryCatch(stats::AIC(fit_stepwise), error = function(e) NA)
        iter_stepwise_result$bic <- tryCatch(stats::BIC(fit_stepwise), error = function(e) NA)
      }
      
      # --- Grid Search ---
      fit_status_grid <- "Fit Failed" # Default status
      start_time_grid <- Sys.time()
      fit_grid <- tryCatch({
        model <- auto.ingarch(y = sim_data, xreg = xreg_matrix, max.p = max_order_p, max.q = max_order_q,
                              distribution = "nbinom", link = "log", ic = "aicc",
                              stepwise = FALSE, trace = FALSE, show_warnings = FALSE)
        fit_status_grid <- "Success"
        model # Return the model on success
      }, error = function(e) {
        error_msg <- conditionMessage(e)
        fit_status_grid <<- paste("Error:", error_msg)
        return(NULL) # Return NULL on error
      })
      end_time_grid <- Sys.time()
      time_grid <- as.numeric(difftime(end_time_grid, start_time_grid, units = "secs"))
      
      # Record grid search results for this iteration
      n_models_grid_expected <- (max_order_p + 1) * (max_order_q + 1)
      iter_grid_result$time <- time_grid
      iter_grid_result$status <- fit_status_grid # Assign status (Success or Error: message)
      if (!is.null(fit_grid) && inherits(fit_grid, "tsglm")) {
        p_grid <- if (is.null(fit_grid$model$past_obs)) 0 else length(fit_grid$model$past_obs)
        q_grid <- if (is.null(fit_grid$model$past_mean)) 0 else length(fit_grid$model$past_mean)
        iter_grid_result$p <- p_grid
        iter_grid_result$q <- q_grid
        iter_grid_result$n_models <- n_models_grid_expected # Grid search tests all models
        iter_grid_result$aic <- tryCatch(stats::AIC(fit_grid), error = function(e) NA)
        iter_grid_result$bic <- tryCatch(stats::BIC(fit_grid), error = function(e) NA)
      } else {
        iter_grid_result$n_models <- NA
      }
      
      # --- Signal progress update for this iteration ---
      p(sprintf("i=%d", i)) # Update progress bar for this iteration
      
      # Return a list containing both results for this iteration 'i'
      return(list(stepwise = iter_stepwise_result, grid_search = iter_grid_result))
      
    } # End foreach loop
  }) # End with_progress block
  
  end_parallel_time <- Sys.time()
  parallel_duration <- end_parallel_time - start_parallel_time
  cat("Parallel model selection finished.\n")
  print(paste("Total parallel execution time:", format(parallel_duration)))
  
  # --- Stop the Cluster ---
  # Important to release resources
  stopCluster(cl)
  # registerDoSEQ() # Optional: Register sequential backend if needed later
  cat("Parallel cluster stopped.\n")
  
  
  # --- Process Parallel Results ---
  results <- list(
    stepwise = vector("list", length(sims)),
    grid_search = vector("list", length(sims))
  )
  
  for (i in 1:length(parallel_results_list)) {
    # Check if the result for iteration 'i' is an error object itself
    if (inherits(parallel_results_list[[i]], "error")) {
      error_msg <- conditionMessage(parallel_results_list[[i]])
      results$stepwise[[i]] <- list(p = NA, q = NA, time = NA, n_models = NA, aic = NA, bic = NA, status = paste("Outer Error:", error_msg))
      results$grid_search[[i]] <- list(p = NA, q = NA, time = NA, n_models = NA, aic = NA, bic = NA, status = paste("Outer Error:", error_msg))
    } else if (is.list(parallel_results_list[[i]]) && !is.null(parallel_results_list[[i]]$stepwise) && !is.null(parallel_results_list[[i]]$grid_search)) {
      # If it's the expected list structure
      results$stepwise[[i]] <- parallel_results_list[[i]]$stepwise
      results$grid_search[[i]] <- parallel_results_list[[i]]$grid_search
    } else {
      # Handle unexpected structure or NULL results
      results$stepwise[[i]] <- list(p = NA, q = NA, time = NA, n_models = NA, aic = NA, bic = NA, status = "Unexpected Result Structure")
      results$grid_search[[i]] <- list(p = NA, q = NA, time = NA, n_models = NA, aic = NA, bic = NA, status = "Unexpected Result Structure")
    }
  }
  cat("Results restructured.\n")
  
  
  # --- 5. Summarize results ---
  cat("Summarizing results for INGARCH with covariates...\n") # Adapted message
  
  # Convert lists to data frames (using robust function from previous version)
  results_to_df <- function(results_list, method_name) {
    df_list <- lapply(1:length(results_list), function(i) {
      res <- results_list[[i]]
      # Add default values if res is NULL or missing components
      if (is.null(res)) {
        res <- list(p=NA, q=NA, time=NA, n_models=NA, aic=NA, bic=NA, status="Result List Null")
      } else {
        # Ensure all expected fields exist, assign NA if not
        res$p <- ifelse("p" %in% names(res), res$p, NA)
        res$q <- ifelse("q" %in% names(res), res$q, NA)
        res$time <- ifelse("time" %in% names(res), res$time, NA)
        res$n_models <- ifelse("n_models" %in% names(res), res$n_models, NA)
        res$aic <- ifelse("aic" %in% names(res), res$aic, NA)
        res$bic <- ifelse("bic" %in% names(res), res$bic, NA)
        res$status <- ifelse("status" %in% names(res), as.character(res$status), "Status Missing") # Ensure status is character
      }
      data.frame( method = method_name, sim_id = i,
                  p = res$p, q = res$q, time = res$time, n_models = res$n_models,
                  aic = res$aic, bic = res$bic, status = res$status ) # Ensure status is included
    })
    do.call(rbind, df_list) # Base R approach
  }
  
  # Make sure dplyr is loaded if using bind_rows
  # library(dplyr)
  stepwise_results_df <- results_to_df(results$stepwise, "stepwise")
  grid_results_df <- results_to_df(results$grid_search, "grid_search")
  results_df <- rbind(stepwise_results_df, grid_results_df)
  
  # output_dir
  output_dir <- "./simulation_with_covariates_output" 
  if (!dir.exists(output_dir)) { dir.create(output_dir) }
  
  # Define filenames
  sims_filename <- file.path(output_dir, "ingarch_with_covariates_simulations.rds")
  results_filename_csv <- file.path(output_dir, "ingarch_with_covariates_results_parallel.csv") # Changed to .csv
  summary_filename <- file.path(output_dir, "ingarch_with_covariates_summary_parallel.csv")
  order_freq_filename <- file.path(output_dir, "ingarch_with_covariates_order_freq_parallel.csv")
  
  # Save simulations (still as RDS as per original structure for this file)
  saveRDS(sims, sims_filename) 
  
  # Save the results_df data frame as CSV only
  write.csv(results_df, results_filename_csv, row.names = FALSE) 
  
  # Calculate summary statistics
  # Make sure status column exists and is character before filtering
  results_df$status <- as.character(results_df$status)
  
  summary_stats <- results_df %>%
    filter(!is.na(method)) %>% # Basic filter
    group_by(method) %>%
    summarize(
      total_sims = n(),
      successful_fits = sum(status == "Success", na.rm = TRUE),
      mean_p = mean(p[status == "Success"], na.rm = TRUE),
      median_p = median(p[status == "Success"], na.rm = TRUE),
      sd_p = sd(p[status == "Success"], na.rm = TRUE),
      mean_q = mean(q[status == "Success"], na.rm = TRUE),
      median_q = median(q[status == "Success"], na.rm = TRUE),
      sd_q = sd(q[status == "Success"], na.rm = TRUE),
      mean_time_secs = mean(time, na.rm = TRUE),
      median_time_secs = median(time, na.rm = TRUE),
      sd_time_secs = sd(time, na.rm = TRUE),
      mean_models_tested = mean(n_models[status == "Success"], na.rm = TRUE), # Only average models tested for successful fits
      median_models_tested = median(n_models[status == "Success"], na.rm = TRUE),
      sd_models_tested = sd(n_models[status == "Success"], na.rm = TRUE),
      failure_count = sum(status != "Success", na.rm = TRUE) # Count non-success statuses
    )
  
  # Order selection frequencies
  order_freq <- results_df %>%
    filter(status == "Success") %>% # Only consider successful fits
    group_by(method, p, q) %>%
    summarize(count = n(), .groups = 'drop') %>%
    group_by(method) %>%
    mutate(freq = count / sum(count) * 100) %>%
    arrange(method, desc(freq))
  
  # Save summaries
  write.csv(summary_stats, summary_filename, row.names = FALSE)
  write.csv(order_freq, order_freq_filename, row.names = FALSE)
  
  # Print summaries (matches original script format)
  cat("\n===== Summary Statistics (Parallel Run) =====\n") # Use generic title
  print(summary_stats)
  cat("\n===== Most Frequent (p,q) Orders (Parallel Run) =====\n") # Use generic title
  if(nrow(order_freq) > 0) {
    top_orders <- order_freq %>% group_by(method) %>% slice_head(n = 5)
    print(top_orders)
  } else {
    cat("No successful fits to determine frequent orders.\n")
  }
  
  
  # Optional: Print failure summary (useful for debugging)
  cat("\n===== Failure/Status Summary (Parallel Run) =====\n")
  # Summarize the actual statuses reported
  failure_summary <- results_df %>%
    group_by(method, status) %>%
    summarise(count = n(), .groups = 'drop') %>%
    arrange(method, desc(count))
  print(failure_summary)
  
  # Final message (matches original script format, adapted path)
  cat(paste0("\nResults saved to directory: ", output_dir, "\n"))
  
  # Ask user to check warnings
  cat("\n--- Please run warnings() to see any specific warnings generated during execution. ---\n")
  
  return(results_df)
}

# --- Execute Main Function ---
# Wrap the execution
main_parallel <- function() {
  results <- run_simulation_study_with_covariates_parallel()
  cat("Parallel simulation study completed!\n") # Matches original script
}

# Call the main function to run the study
main_parallel()

