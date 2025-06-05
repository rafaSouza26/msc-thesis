# --------------------------------------------------------------------------
# R SCRIPT TO FIT INGARCH MODELS (Combined CSV, Parallel, Coeffs, RDS input, External XREG)
# --------------------------------------------------------------------------

# --- 0. Load necessary libraries ---
# install.packages(c("tscount", "foreach", "doParallel")) # Uncomment if not installed
library(tscount)
library(foreach)
library(doParallel)

# --- 1. Define Model Parameters ---
# Model specifications (INGARCH)
# tsglm will fit phi_1, phi_2 (or beta_1, beta_2 if xreg) and theta_2, ..., theta_6
model_spec <- list(past_obs = 1:2, past_mean = 2:6) 
link_function <- "log"
distribution_function <- "nbinom"

# --- 1b. Load Common Covariate Matrix (xreg) from external RData file ---
common_xreg_data <- NULL # Initialize
covariates_rdata_file <- "./Data/count_covariates_data.RData"
covariates_object_name <- "d_1_data"
xreg_column_names_from_file <- c("mean.Temp.Lag5", "mean.DTemp.Lag0", "PM25.Lag0", "mean.NO2.Lag2")
num_rows_for_xreg <- 1000

cat("Loading common covariates from:", covariates_rdata_file, "\n")
loaded_cov_env <- new.env()
tryCatch({
  load(covariates_rdata_file, envir = loaded_cov_env)
  if (exists(covariates_object_name, envir = loaded_cov_env)) {
    all_covariates_df <- loaded_cov_env[[covariates_object_name]]
    if (!is.data.frame(all_covariates_df)) {
      stop(paste0("Object '", covariates_object_name, "' in '", covariates_rdata_file, "' is not a data frame."))
    }
    if (!all(xreg_column_names_from_file %in% colnames(all_covariates_df))) {
      missing_cols <- xreg_column_names_from_file[!xreg_column_names_from_file %in% colnames(all_covariates_df)]
      stop(paste0("Missing required covariate columns in '", covariates_object_name, "': ", paste(missing_cols, collapse=", ")))
    }
    if (nrow(all_covariates_df) < num_rows_for_xreg) {
      stop(paste0("Object '", covariates_object_name, "' has only ", nrow(all_covariates_df), 
                  " rows, but ", num_rows_for_xreg, " are required for xreg."))
    }
    common_xreg_data <- all_covariates_df[1:num_rows_for_xreg, xreg_column_names_from_file, drop = FALSE]
    if(!is.data.frame(common_xreg_data) && !is.matrix(common_xreg_data)){
      common_xreg_data <- as.data.frame(common_xreg_data) 
    }
    cat("Successfully loaded and processed 'common_xreg_data'.\n")
    cat("Dimensions of common_xreg_data:", dim(common_xreg_data), "\n")
    cat("Column names of common_xreg_data:", colnames(common_xreg_data), "\n")
  } else {
    stop(paste0("Object '", covariates_object_name, "' not found in '", covariates_rdata_file, "'."))
  }
}, error = function(e) {
  warning(paste("Error loading or processing common covariates:", e$message))
  warning("common_xreg_data will be NULL. Models requiring xreg may fail or produce unexpected results.")
  common_xreg_data <- NULL 
})

# --- 1c. Define Input RDS Files and Scenario Names ---
# Each element: input_rds_path, scenario_name, needs_xreg (TRUE/FALSE)
file_processing_list <- list(
  list(
    input_rds_path = "./Results/Simulation/no_covariates/Model1/ingarch_no_covariates_simulations.rds",
    scenario_name = "M1", # Scenario name
    needs_xreg = FALSE
  ),
  list(
    input_rds_path = "./Results/Simulation/no_covariates/Model2/ingarch_no_covariates_simulations.rds",
    scenario_name = "M2", # Scenario name
    needs_xreg = FALSE
  ),
  list(
    input_rds_path = "./Results/Simulation/with_covariates/Original/ingarch_with_covariates_simulations.rds",
    scenario_name = "with_covariates", # Scenario name
    needs_xreg = TRUE 
  )
)
# Define the name for the single, combined output CSV file
final_output_csv_filename <- "Models_INGARCH_Fitted_All_Scenarios.csv"


# --- 1d. Setup Parallel Backend ---
if (is.null(common_xreg_data) && any(sapply(file_processing_list, function(x) x$needs_xreg))) {
  stop("FATAL ERROR: common_xreg_data could not be loaded, but some files require it. Stopping execution.")
}
num_cores <- detectCores()
cat("Number of available cores:", num_cores, "\n")
cores_to_use <- max(1, num_cores - 1) 
cat("Registering parallel backend with", cores_to_use, "cores.\n")
cl <- makeCluster(cores_to_use)
registerDoParallel(cl)


# --- Function to process a single input file and RETURN its results (PARALLELIZED) ---
process_simulation_file <- function(input_rds_path, 
                                    scenario_name_param, # New: scenario name
                                    model_spec_to_use, 
                                    link_func, 
                                    distr_func, 
                                    xreg_data_to_pass = NULL) { 
  cat("Processing input file:", input_rds_path, "for Scenario:", scenario_name_param, "\n")
  simulations <- NULL
  
  tryCatch({
    simulations <- readRDS(input_rds_path)
  }, error = function(e) {
    warning(paste("Error reading RDS file:", input_rds_path, "-", e$message, ". Skipping this file."))
    return(NULL) 
  })
  
  if (is.null(simulations) || !is.list(simulations)) {
    warning(paste("Content of", input_rds_path, "is not a list or could not be loaded. Skipping this file."))
    return(NULL)
  }
  
  num_sims <- length(simulations)
  if (num_sims == 0) {
    cat("No simulations found in", input_rds_path, ". Skipping this file.\n")
    return(NULL)
  }
  
  cat("  Starting parallel processing for", num_sims, "simulations from", basename(input_rds_path), "...\n")
  
  current_input_path_for_warnings <- input_rds_path 
  
  file_results_df <- foreach(
    i = 1:num_sims, 
    .combine = 'rbind', 
    .export = c("model_spec_to_use", "link_func", "distr_func", 
                "xreg_data_to_pass", "simulations", "current_input_path_for_warnings"), 
    .packages = "tscount" 
  ) %dopar% {
    
    ts_data <- simulations[[i]] 
    
    result_row <- data.frame(
      # Scenario column will be added after foreach loop
      SimulationIndex = i, AIC = NA_real_, selected_p = NA_integer_, selected_q = NA_integer_,
      alfa1 = NA_real_, alfa2 = NA_real_,
      beta1 = NA_real_, beta2 = NA_real_, beta3 = NA_real_, beta4 = NA_real_, beta5 = NA_real_, beta6 = NA_real_
    )
    
    if (!is.numeric(ts_data) || length(ts_data) != 1000) {
      return(result_row)
    }
    
    fit_model <- NULL 
    fit_model <- tryCatch({
      tsglm(
        ts = ts_data,
        model = model_spec_to_use,
        link = link_func,
        distr = distr_func,
        xreg = xreg_data_to_pass 
      )
    }, error = function(e) {
      return(NULL) 
    })
    
    if (!is.null(fit_model)) {
      coeffs <- coef(fit_model)
      uses_xreg_in_this_fit <- !is.null(xreg_data_to_pass) 
      
      selected_p_order <- if (length(model_spec_to_use$past_obs) > 0) max(model_spec_to_use$past_obs) else 0
      selected_q_order <- if (length(model_spec_to_use$past_mean) > 0) max(model_spec_to_use$past_mean) else 0
      
      result_row$AIC <- AIC(fit_model)
      result_row$selected_p <- selected_p_order
      result_row$selected_q <- selected_q_order
      
      if (uses_xreg_in_this_fit) { 
        if ("beta_1" %in% names(coeffs)) result_row$alfa1 <- coeffs["beta_1"]
        if ("beta_2" %in% names(coeffs)) result_row$alfa2 <- coeffs["beta_2"]
      } else { 
        if ("phi_1" %in% names(coeffs)) result_row$alfa1 <- coeffs["phi_1"]
        if ("phi_2" %in% names(coeffs)) result_row$alfa2 <- coeffs["phi_2"]
      }
      
      if ("theta_2" %in% names(coeffs)) result_row$beta2 <- coeffs["theta_2"]
      if ("theta_3" %in% names(coeffs)) result_row$beta3 <- coeffs["theta_3"]
      if ("theta_4" %in% names(coeffs)) result_row$beta4 <- coeffs["theta_4"]
      if ("theta_5" %in% names(coeffs)) result_row$beta5 <- coeffs["theta_5"]
      if ("theta_6" %in% names(coeffs)) result_row$beta6 <- coeffs["theta_6"]
    }
    return(result_row) 
  } 
  
  cat("  Finished parallel processing for", basename(input_rds_path), ".\n")
  
  if (!is.null(file_results_df) && nrow(file_results_df) > 0) {
    file_results_df$Scenario <- scenario_name_param # Add Scenario column
    
    cat("Results for", basename(input_rds_path), "(Scenario:", scenario_name_param, ") processed.\n")
    cat("Summary of AICs for this scenario:\n")
    if(any(!is.na(file_results_df$AIC))) {
      print(summary(file_results_df$AIC[!is.na(file_results_df$AIC)]))
    } else {
      cat("  No valid AICs to summarize for this scenario.\n")
    }
    return(file_results_df) # Return the data frame with results
  } else {
    cat("No results dataframe generated for file:", input_rds_path, "(Scenario:", scenario_name_param, ").\n")
    return(NULL) # Return NULL if no results
  }
}


# --- Main Processing Loop ---
cat("Starting main processing loop...\n")
all_results_list <- list() # Initialize a list to store all data frames

for (file_config in file_processing_list) {
  current_xreg_data <- NULL
  if (file_config$needs_xreg) {
    if (!is.null(common_xreg_data)) { 
      current_xreg_data <- common_xreg_data 
    } else {
      warning(paste("Skipping xreg for", file_config$input_rds_path, 
                    "because common_xreg_data failed to load. Model will be fit without xreg."))
      # If xreg is absolutely mandatory for this config, consider stopping or skipping this file_config
    }
  }
  
  # Call process_simulation_file and store its result
  single_file_df <- process_simulation_file(
    input_rds_path = file_config$input_rds_path,
    scenario_name_param = file_config$scenario_name, # Pass scenario name
    model_spec_to_use = model_spec, 
    link_func = link_function, 
    distr_func = distribution_function, 
    xreg_data_to_pass = current_xreg_data
  )
  
  if (!is.null(single_file_df) && nrow(single_file_df) > 0) {
    all_results_list[[length(all_results_list) + 1]] <- single_file_df
  }
  cat("Finished processing for scenario:", file_config$scenario_name, "\n\n")
}

# --- Combine and Save All Results to a Single CSV File ---
cat("Combining all results...\n")
if (length(all_results_list) > 0) {
  final_combined_df <- do.call(rbind, all_results_list)
  
  # Define desired column order, matching typical output structure
  desired_column_order <- c("Scenario", "SimulationIndex", "AIC", 
                            "alfa1", "alfa2", 
                            "beta1", "beta2", "beta3", "beta4", "beta5", "beta6",
                            "selected_p", "selected_q")
  
  # Ensure all desired columns are present before trying to reorder
  # Note: rbind will create columns if they exist in any df, filling with NA if not in others.
  # So, all columns from the result_row definition should be present.
  
  # Get current columns and reorder, placing desired ones first
  current_cols <- colnames(final_combined_df)
  ordered_cols <- character(0)
  
  for(col in desired_column_order){
    if(col %in% current_cols){
      ordered_cols <- c(ordered_cols, col)
    }
  }
  # Add any other columns that might exist but weren't in desired_column_order (e.g. xreg coeffs if added later)
  other_existing_cols <- current_cols[!current_cols %in% ordered_cols]
  final_ordered_cols <- c(ordered_cols, other_existing_cols)
  
  final_combined_df <- final_combined_df[, final_ordered_cols, drop = FALSE]
  
  tryCatch({
    write.csv(final_combined_df, file = final_output_csv_filename, row.names = FALSE, na = "") # Write NA as empty string
    cat("\nAll combined results successfully saved to '", final_output_csv_filename, "'\n")
  }, error = function(e) {
    warning(paste("Error saving combined results to CSV '", final_output_csv_filename, "':", e$message))
  })
} else {
  cat("\nNo results were generated from any file to combine and save.\n")
}


# --- Shutdown Parallel Backend ---
cat("Stopping parallel backend...\n")
if(exists("cl") && !is.null(cl)) { 
  tryCatch({
    stopCluster(cl)
    cat("Parallel backend stopped.\n")
  }, error = function(e){
    warning(paste("Error stopping cluster:", e$message))
  })
  cl <- NULL # Set cl to NULL after stopping
} else {
  cat("No parallel backend (cl) to stop or it was already stopped/not initialized.\n")
}

cat("Script execution finished.\n")
