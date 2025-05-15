# --- Script: auto.ingarch for Districts (stepwise = FALSE) ---
# Stores Alphas and Betas as comma-separated strings

# --- 0. Setup ---
# Load necessary libraries
if (!require(tscount)) {
  install.packages("tscount")
  library(tscount)
}

# Suppress warnings globally (use with caution, inspect if unexpected results)
options(warn = -1)
# trace = FALSE and show_warnings = FALSE will be set in auto.ingarch call

# Source custom functions from the 'custom' subdirectory
tryCatch({
  source("./custom/auto.ingarch.R") 
  cat("Custom functions (auto.ingarch and its dependencies) loaded successfully.\n")
}, error = function(e) {
  stop("Error loading custom functions: ", e$message, 
       "\nPlease ensure 'auto.ingarch.R' and its dependencies are accessible, typically in a './custom/' subdirectory relative to this script.")
})


# --- 1. Define File Paths and Parameters ---
# !!! USER: MODIFY THESE PATHS AND SETTINGS AS NEEDED !!!
input_rdata_file <- "./data/count_covariates_data.RData" # Adjusted default input file name and path
output_csv_file_stepwise_false <- "district_ingarch_results_stepwise_FALSE_stringcoeffs.csv" # Output name

# Number of districts to process from the loaded list (e.g., 18 for all, or a smaller number for testing)
NUM_DISTRICTS_TO_PROCESS <- 1 # SET THIS TO THE DESIRED NUMBER 

# auto.ingarch parameters
MAX_P <- 7
MAX_Q <- 7
DISTRIBUTION <- "nbinom"
LINK_FUNCTION <- "log"
IC_CRITERION <- "aic"
STEPWISE_METHOD <- FALSE # << KEY CHANGE FOR THIS SCRIPT
TRACE_OUTPUT <- FALSE 
SHOW_MODEL_WARNINGS <- FALSE 

NUM_COVARIATES <- 4 # Assuming 4 covariates as per description
COEFF_ROUNDING_DIGITS <- 6 # For rounding coefficients in string

# --- 2. Load District Data ---
cat("Loading district data from:", input_rdata_file, "\n")
if (!file.exists(input_rdata_file)) {
  stop("Input RData file not found: ", input_rdata_file, 
       "\nPlease ensure the path is correct and the file 'count_covariates_data.RData' exists (e.g., in a './data/' subdirectory).")
}
load_env <- new.env() # Environment to load RData into
load(input_rdata_file, envir = load_env) 

district_data_list <- NULL
loaded_objects <- ls(envir = load_env)

# Logic to correctly get district_data_list from loaded RData
if (length(loaded_objects) == 1 && is.list(load_env[[loaded_objects[1]]])) {
  # Case 1: RData file contains a single list object (e.g. a list of 18 district dataframes)
  district_data_list <- load_env[[loaded_objects[1]]]
  cat("Loaded data as a single list object:", loaded_objects[1], "\n")
} else {
  # Case 2: RData file contains multiple d_X_data objects (or similar pattern like d_1_data, d_2_data)
  # This attempts to find objects matching the pattern "d_\\d+_data" (e.g., d_1_data)
  district_df_names <- grep("^d_\\d+_data$", loaded_objects, value = TRUE)
  if (length(district_df_names) > 0) {
    district_data_list <- mget(district_df_names, envir = load_env)
    cat("Collated district data from objects:", paste(district_df_names, collapse=", "), "\n")
  } else if ("district_data_list" %in% loaded_objects && is.list(load_env[["district_data_list"]])) {
    # Fallback: if an object specifically named "district_data_list" was loaded
    district_data_list <- load_env[["district_data_list"]]
    cat("Loaded data from an object named 'district_data_list'.\n")
  }
}

if (is.null(district_data_list) || !is.list(district_data_list) || length(district_data_list) == 0) {
  stop("Failed to load district data as a list from '", input_rdata_file ,"'. 
       Please check its structure. It should contain either:
       1. A single list object where each element is a district's data frame.
       2. Multiple data frame objects named like 'd_1_data', 'd_2_data', etc.
       Currently loaded objects: ", paste(loaded_objects, collapse=", "))
}

total_districts_available <- length(district_data_list)
cat("Total districts available in loaded data:", total_districts_available, "\n")

num_districts_to_run <- total_districts_available
if (!is.null(NUM_DISTRICTS_TO_PROCESS) && is.numeric(NUM_DISTRICTS_TO_PROCESS) && NUM_DISTRICTS_TO_PROCESS > 0) {
  if (NUM_DISTRICTS_TO_PROCESS > total_districts_available) {
    cat("Warning: Requested to process", NUM_DISTRICTS_TO_PROCESS, "districts, but only", total_districts_available, "are available. Processing all available.\n")
  } else {
    num_districts_to_run <- NUM_DISTRICTS_TO_PROCESS
  }
}
cat("Will process the first", num_districts_to_run, "districts.\n")

district_names_all <- names(district_data_list)
if (is.null(district_names_all) || length(district_names_all) != total_districts_available) {
  district_names_all <- paste0("District_", 1:total_districts_available) # Fallback names
}
district_names_to_run <- district_names_all[1:num_districts_to_run]
district_data_to_run <- district_data_list[1:num_districts_to_run]


# --- 3. Initialize Results Storage ---
results_collector_list <- list()

# --- 4. Main Loop: Process Each District ---
cat("Starting auto.ingarch modeling for each district (stepwise = FALSE)...\n") # << Log Message Updated

for (i in 1:num_districts_to_run) {
  district_name <- district_names_to_run[i]
  cat("Processing:", district_name, "(", i, "/", num_districts_to_run, ")...\n")
  
  current_data <- district_data_to_run[[i]]
  
  temp_xreg_colnames_for_df <- paste0("Cov", 1:NUM_COVARIATES, "_coeff") 
  if (is.data.frame(current_data) && ncol(current_data) >= (2 + NUM_COVARIATES)) {
    current_xreg_original_names <- colnames(current_data[, 3:(2 + NUM_COVARIATES), drop = FALSE])
    if(!is.null(current_xreg_original_names) && length(current_xreg_original_names) == NUM_COVARIATES && !any(duplicated(current_xreg_original_names)) && all(nchar(current_xreg_original_names)>0) ) {
      temp_xreg_colnames_for_df <- paste0(make.names(current_xreg_original_names), "_coeff")
    }
  }
  
  all_result_colnames <- c("DistrictName", "p_order", "q_order", "AIC", "BIC", "AICc", 
                           "NumModelsTested", "Intercept", "Betas_String", "Alphas_String", 
                           temp_xreg_colnames_for_df, "Error")
  
  error_row <- as.data.frame(matrix(NA, nrow = 1, ncol = length(all_result_colnames)))
  colnames(error_row) <- all_result_colnames
  error_row$DistrictName <- district_name
  
  if (!is.data.frame(current_data) || ncol(current_data) < (2 + NUM_COVARIATES)) {
    cat("  Skipping", district_name, "- data is not a data frame or has insufficient columns.\n")
    error_row$Error <- "Invalid data structure or insufficient columns"
    results_collector_list[[district_name]] <- error_row
    next
  }
  
  count_data <- current_data[[2]] 
  xreg_data <- as.matrix(current_data[, 3:(2 + NUM_COVARIATES), drop = FALSE]) 
  
  original_xreg_colnames <- colnames(xreg_data)
  if (is.null(original_xreg_colnames) || length(original_xreg_colnames) != NUM_COVARIATES || any(duplicated(original_xreg_colnames)) || !all(nchar(original_xreg_colnames)>0)) {
    original_xreg_colnames <- paste0("Cov", 1:NUM_COVARIATES) 
    colnames(xreg_data) <- original_xreg_colnames
  }
  covariate_coeff_colnames_for_df <- paste0(make.names(original_xreg_colnames), "_coeff")
  
  fit_model <- NULL
  error_message <- NA_character_
  
  tryCatch({
    fit_model <- auto.ingarch(
      y = count_data,
      xreg = xreg_data,
      max.p = MAX_P,
      max.q = MAX_Q,
      distribution = DISTRIBUTION,
      link = LINK_FUNCTION,
      ic = IC_CRITERION,
      stepwise = STEPWISE_METHOD, # This is FALSE
      trace = TRACE_OUTPUT,
      show_warnings = SHOW_MODEL_WARNINGS
    )
  }, error = function(e) {
    cat("  Error fitting model for", district_name, ":", conditionMessage(e), "\n")
    error_message <<- conditionMessage(e) 
  })
  
  if (!is.null(fit_model) && inherits(fit_model, "tsglm")) {
    p_order <- if (is.null(fit_model$model$past_obs)) 0 else length(fit_model$model$past_obs)
    q_order <- if (is.null(fit_model$model$past_mean)) 0 else length(fit_model$model$past_mean)
    
    aic_val <- tryCatch(fit_model$aic, error = function(e) NA)
    bic_val <- tryCatch(fit_model$bic, error = function(e) NA)
    aicc_val <- tryCatch(fit_model$aicc, error = function(e) NA)
    
    num_models <- NA
    # For stepwise=FALSE, auto.ingarch calls search.ingarch.
    # The 'results' element of the returned model is a matrix of tested models (p, q, ic).
    if (!is.null(fit_model$results) && is.matrix(fit_model$results)) {
      num_models <- nrow(fit_model$results) 
    } else { 
      # Fallback if results matrix isn't populated as expected, calculate from max orders
      num_models <- (MAX_P + 1) * (MAX_Q + 1) 
    }
    
    coeffs <- fit_model$coefficients
    intercept_val <- ifelse("intercept" %in% names(coeffs), coeffs["intercept"], NA_real_)
    
    # Betas String (same logic as stepwise=TRUE version)
    betas_string <- NA_character_
    if (p_order > 0 && !is.null(coeffs) && length(coeffs) > 0) {
      actual_betas_coeffs <- coeffs[grep("^beta", names(coeffs))]
      if (length(actual_betas_coeffs) > 0) {
        beta_names <- names(actual_betas_coeffs)
        beta_indices <- sapply(beta_names, function(name) {
          if (name == "beta") return(1L)
          num_part <- sub("beta", "", name)
          if (grepl("^[0-9]+$", num_part)) return(as.integer(num_part))
          return(NA_integer_)
        }, USE.NAMES = FALSE)
        valid_indices <- !is.na(beta_indices)
        actual_betas_coeffs <- actual_betas_coeffs[valid_indices]
        beta_indices <- beta_indices[valid_indices]
        if(length(actual_betas_coeffs) > 0){
          ordered_betas <- actual_betas_coeffs[order(beta_indices)]
          betas_string <- paste(round(ordered_betas, COEFF_ROUNDING_DIGITS), collapse = ", ")
        } else { betas_string <- "" }
      } else { betas_string <- "" }
    } else if (p_order == 0) { betas_string <- "" }
    
    # Alphas String (same logic as stepwise=TRUE version)
    alphas_string <- NA_character_
    if (q_order > 0 && !is.null(coeffs) && length(coeffs) > 0) {
      actual_alphas_coeffs <- coeffs[grep("^alpha", names(coeffs))]
      if (length(actual_alphas_coeffs) > 0) {
        alpha_names <- names(actual_alphas_coeffs)
        alpha_indices <- sapply(alpha_names, function(name) {
          if (name == "alpha") return(1L)
          num_part <- sub("alpha", "", name)
          if (grepl("^[0-9]+$", num_part)) return(as.integer(num_part))
          return(NA_integer_)
        }, USE.NAMES = FALSE)
        valid_indices <- !is.na(alpha_indices)
        actual_alphas_coeffs <- actual_alphas_coeffs[valid_indices]
        alpha_indices <- alpha_indices[valid_indices]
        if(length(actual_alphas_coeffs) > 0){
          ordered_alphas <- actual_alphas_coeffs[order(alpha_indices)]
          alphas_string <- paste(round(ordered_alphas, COEFF_ROUNDING_DIGITS), collapse = ", ")
        } else { alphas_string <- "" }
      } else { alphas_string <- "" }
    } else if (q_order == 0) { alphas_string <- "" }
    
    cov_coeffs_padded <- rep(NA_real_, NUM_COVARIATES)
    if(!is.null(coeffs) && length(coeffs) > 0){
      for(k_idx in 1:NUM_COVARIATES){
        original_cov_name <- original_xreg_colnames[k_idx]
        tsglm_cov_name_pattern <- paste0("xreg", k_idx) 
        if (original_cov_name %in% names(coeffs)) {
          cov_coeffs_padded[k_idx] <- coeffs[original_cov_name]
        } else if (tsglm_cov_name_pattern %in% names(coeffs)) { 
          cov_coeffs_padded[k_idx] <- coeffs[tsglm_cov_name_pattern]
        }
      }
    }
    names(cov_coeffs_padded) <- covariate_coeff_colnames_for_df
    
    results_df_row_list <- list(
      DistrictName = district_name,
      p_order = p_order,
      q_order = q_order,
      AIC = aic_val,
      BIC = bic_val,
      AICc = aicc_val,
      NumModelsTested = num_models,
      Intercept = intercept_val,
      Betas_String = betas_string,
      Alphas_String = alphas_string
    )
    for(k_idx in 1:NUM_COVARIATES){
      results_df_row_list[[covariate_coeff_colnames_for_df[k_idx]]] <- cov_coeffs_padded[k_idx]
    }
    results_df_row_list[["Error"]] <- error_message
    
    current_results_df <- as.data.frame(results_df_row_list, stringsAsFactors = FALSE)
    results_collector_list[[district_name]] <- current_results_df
    
  } else {
    cat("  Failed to fit model for", district_name, "or model object is not tsglm.\n")
    error_row$Error <- ifelse(is.na(error_message), "Fit failed or model not tsglm", error_message)
    results_collector_list[[district_name]] <- error_row
  }
}

# --- 5. Combine and Save Results ---
cat("Combining results...\n")
if (length(results_collector_list) > 0) {
  base_cols <- c("DistrictName", "p_order", "q_order", "AIC", "BIC", "AICc", 
                 "NumModelsTested", "Intercept", "Betas_String", "Alphas_String")
  all_dynamic_cov_cols <- unique(unlist(lapply(results_collector_list, function(df) {
    grep("_coeff$", names(df), value = TRUE)
  })))
  final_ordered_colnames <- c(base_cols, sort(all_dynamic_cov_cols), "Error")
  
  standardized_results_list <- lapply(results_collector_list, function(df) {
    new_df <- data.frame(matrix(NA, nrow = 1, ncol = length(final_ordered_colnames)))
    colnames(new_df) <- final_ordered_colnames
    for (col_name in names(df)) {
      if (col_name %in% final_ordered_colnames) {
        new_df[1, col_name] <- df[1, col_name]
      }
    }
    new_df$DistrictName <- as.character(df$DistrictName)
    return(new_df)
  })
  
  final_results_df <- do.call(rbind, standardized_results_list)
  rownames(final_results_df) <- NULL 
  
  cat("Saving results to CSV:", output_csv_file_stepwise_false, "\n")
  tryCatch({
    write.csv(final_results_df, output_csv_file_stepwise_false, row.names = FALSE, na = "")
    cat("Successfully saved results.\n")
  }, error = function(e) {
    cat("Error saving CSV file:", conditionMessage(e), "\n")
  })
} else {
  cat("No results to save.\n")
}

# Restore warning settings
options(warn = 0)

cat("--- Script finished for stepwise = FALSE (string coeffs) ---\n") # << Log Message Updated