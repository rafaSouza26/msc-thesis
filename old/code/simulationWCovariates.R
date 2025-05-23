# INGARCH Model Simulation Study - With Covariates
# Comparing model selection methods: auto.ingarch with stepwise vs grid search

# Load required packages
# install.packages(c("tscount", "ggplot2", "dplyr", "readxl")) # Uncomment to install if needed
library(tscount)   # For INGARCH models
library(ggplot2)   # For visualization
library(dplyr)     # For data manipulation
library(readxl)    # For reading Excel files

# Source custom functions
# Ensure these paths are correct relative to your script's location
source("./ACTS/auto.ingarch.R")
source("./ACTS/ingarch.sim.R") # Ensure this uses xreg=param$external internally now
source("./ACTS/newmodel.R")
source("./ACTS/ingarch.string.R")
source("./ACTS/search.ingarch.R")

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

# Main simulation study - With Covariates
run_simulation_study_with_covariates <- function() {
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
  n_obs_sim <- 1000
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
  cat("beta coefficients:", if(is.null(params$betas)) "None" else params$betas, "\n")
  cat("alpha coefficients:", if(is.null(params$alphas)) "None" else params$alphas, "\n")
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
  n_sims <- 1000 # Set desired number of simulations
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
  
  
  # --- 4. Run model selection for each simulation ---
  cat("Running model selection for INGARCH with covariates...\n") # Adapted message
  
  results <- list(stepwise = vector("list", length(sims)),
                  grid_search = vector("list", length(sims)))
  max_order_p <- 7
  max_order_q <- 7
  
  for(i in 1:length(sims)) {
    # Model selection progress (matches original frequency)
    if(i %% 10 == 0) cat("Processing simulation", i, "of", length(sims), "\n")
    
    sim_data <- sims[[i]]
    if (is.null(sim_data)) {
      # Minimal message for skipped simulation
      # cat("Skipping model selection for failed simulation", i, "\n") # Optionally remove this too
      results$stepwise[[i]] <- list(p = NA, q = NA, time = 0, n_models = NA, aic = NA, bic = NA, status = "Simulation Failed")
      results$grid_search[[i]] <- list(p = NA, q = NA, time = 0, n_models = NA, aic = NA, bic = NA, status = "Simulation Failed")
      next
    }
    if (length(sim_data) != nrow(xreg_matrix)) {
      warning(paste("Length mismatch for simulation", i, "before fitting. Skipping."), immediate. = TRUE)
      results$stepwise[[i]] <- list(p = NA, q = NA, time = 0, n_models = NA, aic = NA, bic = NA, status = "Length Mismatch")
      results$grid_search[[i]] <- list(p = NA, q = NA, time = 0, n_models = NA, aic = NA, bic = NA, status = "Length Mismatch")
      next
    }
    
    # --- Stepwise ---
    fit_status_stepwise <- "Fit Failed"
    start_time_stepwise <- Sys.time()
    fit_stepwise <- tryCatch({
      model <- auto.ingarch(y = sim_data, xreg = xreg_matrix, max.p = max_order_p, max.q = max_order_q,
                            distribution = "nbinom", link = "log", ic = "aicc",
                            stepwise = TRUE, trace = FALSE, show_warnings = FALSE)
      fit_status_stepwise <- "Success"
      model
    }, error = function(e) {
      # Keep error message print (matches original script)
      cat("Error in stepwise fit for simulation", i, ":", conditionMessage(e), "\n")
      fit_status_stepwise <<- paste("Error:", substr(conditionMessage(e), 1, 50)) # Store brief error
      return(NULL)
    })
    end_time_stepwise <- Sys.time()
    time_stepwise <- as.numeric(difftime(end_time_stepwise, start_time_stepwise, units = "secs"))
    
    # Record stepwise results
    if(!is.null(fit_stepwise) && inherits(fit_stepwise, "tsglm")) {
      p_stepwise <- if(is.null(fit_stepwise$model$past_obs)) 0 else length(fit_stepwise$model$past_obs)
      q_stepwise <- if(is.null(fit_stepwise$model$past_mean)) 0 else length(fit_stepwise$model$past_mean)
      n_models_stepwise <- if(is.null(fit_stepwise$results)) NA else nrow(fit_stepwise$results)
      results$stepwise[[i]] <- list(p = p_stepwise, q = q_stepwise, time = time_stepwise, n_models = n_models_stepwise, status = fit_status_stepwise,
                                    aic = tryCatch(AIC(fit_stepwise), error = function(e) NA), bic = tryCatch(BIC(fit_stepwise), error = function(e) NA))
    } else {
      results$stepwise[[i]] <- list(p = NA, q = NA, time = time_stepwise, n_models = NA, status = fit_status_stepwise, aic = NA, bic = NA)
    }
    
    # --- Grid Search ---
    fit_status_grid <- "Fit Failed"
    start_time_grid <- Sys.time()
    fit_grid <- tryCatch({
      model <- auto.ingarch(y = sim_data, xreg = xreg_matrix, max.p = max_order_p, max.q = max_order_q,
                            distribution = "nbinom", link = "log", ic = "aicc",
                            stepwise = FALSE, trace = FALSE, show_warnings = FALSE)
      fit_status_grid <- "Success"
      model
    }, error = function(e) {
      # Keep error message print (matches original script)
      cat("Error in grid search fit for simulation", i, ":", conditionMessage(e), "\n")
      fit_status_grid <<- paste("Error:", substr(conditionMessage(e), 1, 50)) # Store brief error
      return(NULL)
    })
    end_time_grid <- Sys.time()
    time_grid <- as.numeric(difftime(end_time_grid, start_time_grid, units = "secs"))
    
    # Record grid search results
    n_models_grid_expected <- (max_order_p + 1) * (max_order_q + 1)
    if(!is.null(fit_grid) && inherits(fit_grid, "tsglm")) {
      p_grid <- if(is.null(fit_grid$model$past_obs)) 0 else length(fit_grid$model$past_obs)
      q_grid <- if(is.null(fit_grid$model$past_mean)) 0 else length(fit_grid$model$past_mean)
      results$grid_search[[i]] <- list(p = p_grid, q = q_grid, time = time_grid, n_models = n_models_grid_expected, status = fit_status_grid,
                                       aic = tryCatch(AIC(fit_grid), error = function(e) NA), bic = tryCatch(BIC(fit_grid), error = function(e) NA))
    } else {
      results$grid_search[[i]] <- list(p = NA, q = NA, time = time_grid, n_models = NA, status = fit_status_grid, aic = NA, bic = NA)
    }
  } # End model selection loop
  
  
  # --- 5. Summarize results ---
  cat("Summarizing results for INGARCH with covariates...\n") # Adapted message
  
  # Convert lists to data frames (using robust function from previous version)
  results_to_df <- function(results_list, method_name) {
    df_list <- lapply(1:length(results_list), function(i) {
      res <- results_list[[i]]
      if (is.null(res)) { res <- list(p=NA, q=NA, time=NA, n_models=NA, aic=NA, bic=NA, status="Result List Null") }
      data.frame( method = method_name, sim_id = i,
                  p = ifelse("p" %in% names(res), res$p, NA), q = ifelse("q" %in% names(res), res$q, NA),
                  time = ifelse("time" %in% names(res), res$time, NA), n_models = ifelse("n_models" %in% names(res), res$n_models, NA),
                  aic = ifelse("aic" %in% names(res), res$aic, NA), bic = ifelse("bic" %in% names(res), res$bic, NA),
                  status = ifelse("status" %in% names(res), res$status, NA) )
    })
    do.call(rbind, df_list)
  }
  stepwise_results_df <- results_to_df(results$stepwise, "stepwise")
  grid_results_df <- results_to_df(results$grid_search, "grid_search")
  results_df <- rbind(stepwise_results_df, grid_results_df)
  
  # Save results (no printing here, final message later)
  output_dir <- "./simulation_output"
  if (!dir.exists(output_dir)) { dir.create(output_dir) }
  sims_filename <- file.path(output_dir, "ingarch_with_covariates_simulations.rds")
  results_filename <- file.path(output_dir, "ingarch_with_covariates_results.rds")
  summary_filename <- file.path(output_dir, "ingarch_with_covariates_summary.csv")
  order_freq_filename <- file.path(output_dir, "ingarch_with_covariates_order_freq.csv")
  saveRDS(sims, sims_filename)
  saveRDS(results_df, results_filename)
  
  # Calculate summary statistics
  summary_stats <- results_df %>% filter(!is.na(method)) %>% group_by(method) %>%
    summarize(
      total_sims = n(), successful_fits = sum(status == "Success", na.rm = TRUE),
      mean_p = mean(p[status == "Success"], na.rm = TRUE), median_p = median(p[status == "Success"], na.rm = TRUE), sd_p = sd(p[status == "Success"], na.rm = TRUE),
      mean_q = mean(q[status == "Success"], na.rm = TRUE), median_q = median(q[status == "Success"], na.rm = TRUE), sd_q = sd(q[status == "Success"], na.rm = TRUE),
      mean_time_secs = mean(time, na.rm = TRUE), median_time_secs = median(time, na.rm = TRUE), sd_time_secs = sd(time, na.rm = TRUE),
      mean_models_tested = mean(n_models, na.rm = TRUE), median_models_tested = median(n_models, na.rm = TRUE), sd_models_tested = sd(n_models, na.rm = TRUE),
      failure_count = sum(status != "Success", na.rm = TRUE)
    )
  
  # Order selection frequencies
  order_freq <- results_df %>% filter(status == "Success") %>% group_by(method, p, q) %>%
    summarize(count = n(), .groups = 'drop') %>% group_by(method) %>%
    mutate(freq = count / sum(count) * 100) %>% arrange(method, desc(freq))
  
  # Save summaries
  write.csv(summary_stats, summary_filename, row.names = FALSE)
  write.csv(order_freq, order_freq_filename, row.names = FALSE)
  
  # Print summaries (matches original script format)
  cat("\n===== Summary Statistics =====\n") # Use generic title
  print(summary_stats)
  cat("\n===== Most Frequent (p,q) Orders =====\n") # Use generic title
  top_orders <- order_freq %>% group_by(method) %>% slice_head(n = 5)
  print(top_orders)
  
  # Optional: Print failure summary (useful for debugging)
  # cat("\n===== Failure/Status Summary =====\n")
  # failure_summary <- results_df %>% group_by(method, status) %>%
  #      summarise(count = n(), .groups = 'drop') %>% arrange(method, desc(count))
  # print(failure_summary)
  
  # Final message (matches original script format, adapted path)
  cat(paste0("\nResults saved to: ", results_filename, "\n"))
  # Removed separate prints for each file
  
  return(results_df)
}

# Run the study
main_with_covariates <- function() {
  results <- run_simulation_study_with_covariates()
  cat("Simulation study completed!\n") # Matches original script
}

# Execute main function
main_with_covariates()