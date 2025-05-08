# --- Apply auto.ingarch to District Data with Covariates by Index (Sequential) ---
# Version: 2025-05-07 (Simplified with index-based covariates) # User's version
# Modified on: 2025-05-08 (Prioritize n_total_models, Save all tested models)

# Record script start time
script_start_time <- Sys.time()

# --- 1. Load required packages ---
library(tscount)
library(dplyr)
library(readxl) # Note: readxl is loaded but main data uses .RData

# --- 2. Source custom functions ---
cat("Sourcing custom functions...\n")
tryCatch(source("./custom/auto.ingarch.R"), error = function(e) stop("Failed to source auto.ingarch.R: ", e))
tryCatch(source("./custom/newmodel.R"), error = function(e) stop("Failed to source newmodel.R: ", e))
tryCatch(source("./custom/ingarch.string.R"), error = function(e) stop("Failed to source ingarch.string.R: ", e))
tryCatch(source("./custom/search.ingarch.R"), error = function(e) stop("Failed to source search.ingarch.R: ", e))

# --- 3. Define Parameters and Setup ---
cat("Defining parameters...\n")
MAX_P <- 7
MAX_Q <- 7
DISTRIBUTION <- "nbinom"
LINK <- "log"
IC <- "aicc" # Information Criterion used for model selection
STEPWISE <- TRUE  # Always use stepwise approach

DATA_PATH <- "./data/count_covariates_data.RData"
COUNT_VAR_NAME <- "Count"
# Using indices instead of names for data columns
COUNT_VAR_INDEX <- 2  # Count is in column 2
COVARIATE_INDICES <- 3:6  # Covariates are in columns 3, 4, 5, 6
MIN_EXPECTED_COLUMNS <- max(c(COUNT_VAR_INDEX, COVARIATE_INDICES))

OUTPUT_DIR <- "./district_fitting_output_simplified"
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}
# Files for BEST model summaries
RESULTS_RDS_FILE <- file.path(OUTPUT_DIR, "district_fitting_results_simplified.rds")
RESULTS_CSV_FILE <- file.path(OUTPUT_DIR, "district_fitting_summary_simplified.csv")
# NEW Files for ALL TESTED models
ALL_TESTED_MODELS_RDS_FILE <- file.path(OUTPUT_DIR, "all_district_tested_models_details.rds")
ALL_TESTED_MODELS_CSV_FILE <- file.path(OUTPUT_DIR, "all_district_tested_models_details.csv")


# --- Configuration for Districts to Process ---
# Can easily adjust how many districts to process
NUM_DISTRICTS_TO_PROCESS <- 18  # Change this number as needed (1-18)
ALL_DISTRICT_CODES_AVAILABLE <- 1:18

# Ensure NUM_DISTRICTS_TO_PROCESS is valid
if (NUM_DISTRICTS_TO_PROCESS > length(ALL_DISTRICT_CODES_AVAILABLE)) {
  warning(paste("NUM_DISTRICTS_TO_PROCESS (", NUM_DISTRICTS_TO_PROCESS,
                ") is greater than available districts (", length(ALL_DISTRICT_CODES_AVAILABLE),
                "). Processing all available districts instead."))
  NUM_DISTRICTS_TO_PROCESS <- length(ALL_DISTRICT_CODES_AVAILABLE)
}
if (NUM_DISTRICTS_TO_PROCESS < 1) {
  warning("NUM_DISTRICTS_TO_PROCESS is less than 1. Setting to 1.")
  NUM_DISTRICTS_TO_PROCESS <- 1
}

DISTRICT_CODES_TO_PROCESS <- ALL_DISTRICT_CODES_AVAILABLE[1:NUM_DISTRICTS_TO_PROCESS]

cat("Processing", NUM_DISTRICTS_TO_PROCESS, "district(s). Codes:", paste(DISTRICT_CODES_TO_PROCESS, collapse=", "), "\n")

# District names lookup table
DISTRICT_NAMES <- c("Aveiro", "Beja", "Braga", "Bragança", "Castelo Branco",
                    "Coimbra", "Évora", "Faro", "Guarda", "Leiria", "Lisboa",
                    "Portalegre", "Porto", "Santarém", "Setúbal",
                    "Viana do Castelo", "Vila Real", "Viseu")
district_lookup <- setNames(DISTRICT_NAMES, ALL_DISTRICT_CODES_AVAILABLE)

# --- 4. Load Data ---
cat("Loading district data from:", DATA_PATH, "\n")
if (!file.exists(DATA_PATH)) {
  stop("Data file not found: ", DATA_PATH)
}
load(DATA_PATH, envir = .GlobalEnv)
cat("Data loaded into the Global Environment.\n")

# --- 5. Verify Data Objects & Basic Structure Check ---
cat("Verifying data objects and basic structure for selected districts...\n")
missing_data_objs <- c()
# Create names for only the districts we will process
all_df_names_to_process <- paste0("d_", DISTRICT_CODES_TO_PROCESS, "_data")
all_struct_ok <- TRUE

for (df_name in all_df_names_to_process) {
  if (!exists(df_name, envir = .GlobalEnv)) {
    missing_data_objs <- c(missing_data_objs, df_name)
    all_struct_ok <- FALSE
  } else {
    df_temp <- get(df_name, envir = .GlobalEnv)
    if (!is.data.frame(df_temp)){
      cat("ERROR: Object '", df_name, "' is not a data frame.\n")
      all_struct_ok <- FALSE
    } else if (ncol(df_temp) < MIN_EXPECTED_COLUMNS) {
      cat("ERROR: Data frame '", df_name, "' has ", ncol(df_temp), " columns, but expected at least ", MIN_EXPECTED_COLUMNS, ".\n")
      all_struct_ok <- FALSE
    } else if (names(df_temp)[COUNT_VAR_INDEX] != COUNT_VAR_NAME) {
      cat("ERROR: Data frame '", df_name, "' column ", COUNT_VAR_INDEX, " is named '", names(df_temp)[COUNT_VAR_INDEX], "' but expected '", COUNT_VAR_NAME, "'.\n")
      all_struct_ok <- FALSE
    }
  }
}

if (length(missing_data_objs) > 0) {
  stop("Stopping script. Missing required data frame objects for selected districts: ", paste(missing_data_objs, collapse=", "))
}
if (!all_struct_ok) {
  stop("Stopping script due to data structure issues in selected districts. Please check error messages above and your data file/indices.")
} else {
  cat("Basic structure check passed for selected districts.\n")
}

# --- 6. Sequential Model Fitting ---
cat("Starting sequential model fitting for", length(DISTRICT_CODES_TO_PROCESS), "districts...\n")

# Initialize results list for BEST models
results_list <- list()
# NEW: Initialize list for ALL tested models
all_tested_models_list <- list()

# Sequential loop over each district
for (i in DISTRICT_CODES_TO_PROCESS) {
  cat("\nProcessing District", i, "(", district_lookup[as.character(i)], ")...\n")
  
  # Initialize result for this district's BEST model
  iter_result <- list(
    district_code = i, p = NA_integer_, q = NA_integer_,
    aic = NA_real_, bic = NA_real_, aicc = NA_real_, # Note: aicc here is for the best model
    time_secs = NA_real_, n_models_stepwise = 0,
    betas_str = NA_character_, alphas_str = NA_character_,
    covariates_used_names = NA_character_,
    status = "Not Run"
  )
  
  df_name <- paste0("d_", i, "_data")
  xreg <- NULL
  cov_names_for_this_district <- character(0)
  
  current_data <- tryCatch({
    df <- get(df_name, envir = .GlobalEnv)
    if (!is.data.frame(df)) stop("Object is not a data frame.")
    df
  }, error = function(e) {
    iter_result$status <<- paste("Data Error:", conditionMessage(e))
    return(NULL)
  })
  
  if (!is.null(current_data)) {
    if (ncol(current_data) < MIN_EXPECTED_COLUMNS) {
      iter_result$status <- "Insufficient columns in data"
    } else {
      y <- current_data[, COUNT_VAR_INDEX]
      
      if(max(COVARIATE_INDICES) <= ncol(current_data)) {
        xreg <- as.matrix(current_data[, COVARIATE_INDICES])
        cov_names_for_this_district <- colnames(current_data)[COVARIATE_INDICES]
        iter_result$covariates_used_names <- paste(cov_names_for_this_district, collapse=", ")
      } else {
        iter_result$status <- "Covariate indices out of bounds for this dataframe."
        xreg <- NULL 
      }
      
      if (anyNA(y)) {
        iter_result$status <- "NA values in response (y)"
      } else if (!is.null(xreg) && anyNA(xreg)) {
        iter_result$status <- "NA values in covariates (xreg)"
      } else if (iter_result$status == "Not Run") {
        fit_status <- "Fit Failed"
        start_time_fit <- Sys.time()
        
        fit_model <- tryCatch({
          cat("  Fitting auto.ingarch model with covariates...\n")
          model <- auto.ingarch(
            y = y, xreg = xreg, max.p = MAX_P, max.q = MAX_Q,
            distribution = DISTRIBUTION, link = LINK, ic = IC,
            stepwise = STEPWISE, trace = FALSE, show_warnings = FALSE
          )
          if (!inherits(model, "tsglm")) {
            stop("auto.ingarch did not return a valid model object")
          }
          fit_status <<- "Success"
          model
        }, error = function(e) {
          fit_status <<- paste("Fitting Error:", conditionMessage(e))
          cat("  Model fitting error:", conditionMessage(e), "\n")
          return(NULL)
        })
        
        end_time_fit <- Sys.time()
        time_fit <- as.numeric(difftime(end_time_fit, start_time_fit, units = "secs"))
        
        iter_result$time_secs <- time_fit
        iter_result$status <- fit_status
        
        if (fit_status == "Success" && !is.null(fit_model) && inherits(fit_model, "tsglm")) {
          iter_result$p <- if (is.null(fit_model$model$past_obs)) 0 else length(fit_model$model$past_obs)
          iter_result$q <- if (is.null(fit_model$model$past_mean)) 0 else length(fit_model$model$past_mean)
          
          iter_result$aic <- tryCatch(stats::AIC(fit_model), error = function(e) NA_real_)
          iter_result$bic <- tryCatch(stats::BIC(fit_model), error = function(e) NA_real_)
          
          if (!is.null(fit_model[[IC]])) { # IC value for the BEST model
            iter_result[[IC]] <- fit_model[[IC]] # Store it in the correctly named column (e.g. iter_result$aicc)
          } else if (!is.null(fit_model$results) && IC %in% names(fit_model$results)) {
            # This path for IC of best model might be redundant if auto.ingarch returns IC directly
            best_model_row <- which.min(fit_model$results[[IC]]) 
            if(length(best_model_row) == 1) iter_result[[IC]] <- fit_model$results[[IC]][best_model_row]
          }
          
          # Count number of models evaluated (prioritizing n_total_models)
          if (!is.null(fit_model$n_total_models)) {
            iter_result$n_models_stepwise <- fit_model$n_total_models
          } else if (!is.null(fit_model$results) && (is.matrix(fit_model$results) || is.data.frame(fit_model$results))) {
            iter_result$n_models_stepwise <- nrow(as.data.frame(fit_model$results))
          } else {
            iter_result$n_models_stepwise <- 1 
          }
          
          # --- NEW: Store ALL models tested for this district ---
          if (!is.null(fit_model$results) && (is.matrix(fit_model$results) || is.data.frame(fit_model$results))) {
            if (nrow(as.data.frame(fit_model$results)) > 0) { # Check if there are any tested models logged
              district_tested_models_df <- as.data.frame(fit_model$results)
              # Rename columns: auto.ingarch $results has "p", "q", "ic"
              # The 'ic' column contains values of the criterion defined by the 'IC' variable (e.g., "aicc")
              colnames(district_tested_models_df) <- c("tested_p", "tested_q", paste0(IC, "_value"))
              
              district_tested_models_df$district_code <- i
              district_tested_models_df$district_name <- district_lookup[as.character(i)]
              
              # Reorder for clarity
              district_tested_models_df <- district_tested_models_df[, c("district_code", "district_name", 
                                                                         "tested_p", "tested_q", 
                                                                         paste0(IC, "_value"))]
              all_tested_models_list[[length(all_tested_models_list) + 1]] <- district_tested_models_df
            }
          }
          # --- END NEW SECTION for ALL tested models ---
          
          all_coeffs <- stats::coef(fit_model)
          if (iter_result$p > 0) {
            beta_coef_names <- paste0("beta_", 1:iter_result$p)
            actual_beta_coeffs <- all_coeffs[names(all_coeffs) %in% beta_coef_names]
            iter_result$betas_str <- if (length(actual_beta_coeffs) > 0) paste(round(actual_beta_coeffs, 5), collapse=",") else ""
          } else { iter_result$betas_str <- "" }
          
          if (iter_result$q > 0) {
            alpha_coef_names <- paste0("alpha_", 1:iter_result$q)
            actual_alpha_coeffs <- all_coeffs[names(all_coeffs) %in% alpha_coef_names]
            iter_result$alphas_str <- if (length(actual_alpha_coeffs) > 0) paste(round(actual_alpha_coeffs, 5), collapse=",") else ""
          } else { iter_result$alphas_str <- "" }
          
          cat("  Model fitted successfully: p=", iter_result$p, ", q=", iter_result$q,  
              ", ", IC, "=", round(iter_result[[IC]], 2), # Use the IC variable for printing
              ", Models evaluated=", iter_result$n_models_stepwise, "\n", sep="")
        } else {
          cat("  Model fitting failed with status:", fit_status, "\n")
        }
      }
    }
  }
  
  results_list[[length(results_list) + 1]] <- iter_result # For BEST model
  cat("Completed District", i, "with status:", iter_result$status, "\n")
}

# --- 7. Process and Combine Results (BEST Models) ---
cat("\nProcessing results for BEST models...\n")
# (This section remains largely the same, processing 'results_list' for best models)
if (length(results_list) > 0) {
  results_df <- bind_rows(lapply(results_list, function(res) {
    if ("n_models_stepwise" %in% names(res) && is.na(res$n_models_stepwise)) {
      res$n_models_stepwise <- 0
    }
    as.data.frame(res, stringsAsFactors = FALSE)
  }))
  
  results_df <- results_df %>%
    mutate(district_name = ifelse(!is.na(district_code), district_lookup[as.character(district_code)], NA_character_)) %>%
    select(district_code, district_name, p, q, aic, bic, aicc, time_secs, # Ensure 'aicc' or general IC column is here
           n_models_stepwise, betas_str, alphas_str, covariates_used_names, status)
} else {
  results_df <- data.frame(
    district_code = integer(0), district_name = character(0), p = integer(0), q = integer(0),
    aic = numeric(0), bic = numeric(0), aicc = numeric(0), # Ensure 'aicc' or general IC column
    time_secs = numeric(0), n_models_stepwise = integer(0),
    betas_str = character(0), alphas_str = character(0),
    covariates_used_names = character(0), status = character(0),
    stringsAsFactors = FALSE
  )
  cat("No BEST model results to process or save.\n")
}

# --- 7.1. Process and Combine ALL Tested Model Results --- (NEW SECTION)
cat("\nProcessing ALL tested model results...\n")
all_tested_models_df <- data.frame() 

if (length(all_tested_models_list) > 0) {
  all_tested_models_df <- bind_rows(all_tested_models_list)
  cat(nrow(all_tested_models_df), "total models documented as tested across all processed districts.\n")
} else {
  cat("No detailed tested model results to process or save.\n")
  # Create empty dataframe with expected structure for consistency
  empty_df_cols <- list(
    district_code = integer(0),
    district_name = character(0),
    tested_p = integer(0),
    tested_q = integer(0)
  )
  ic_col_name <- paste0(IC, "_value") # IC is defined in Section 3 (e.g., "aicc")
  empty_df_cols[[ic_col_name]] <- numeric(0)
  all_tested_models_df <- as.data.frame(empty_df_cols, stringsAsFactors = FALSE)
}

# --- 8. Save Results (BEST Models) ---
# (This section remains largely the same, saving 'results_df')
if(nrow(results_df) > 0) {
  cat("Saving BEST model results to:", OUTPUT_DIR, "\n")
  tryCatch({  
    saveRDS(results_df, file = RESULTS_RDS_FILE)
    cat("Detailed BEST model results saved to:", RESULTS_RDS_FILE, "\n")
  }, error = function(e) { cat("ERROR saving BEST models RDS file:", conditionMessage(e), "\n") })
  
  tryCatch({  
    results_df$n_models_stepwise <- as.numeric(results_df$n_models_stepwise)
    results_df$n_models_stepwise[is.na(results_df$n_models_stepwise)] <- 0
    
    write.csv(results_df, file = RESULTS_CSV_FILE, row.names = FALSE, na = "")
    cat("Summary BEST model results saved to:", RESULTS_CSV_FILE, "\n")
  }, error = function(e) { cat("ERROR saving BEST models CSV file:", conditionMessage(e), "\n") })
}

# --- 8.1. Save ALL Tested Model Results --- (NEW SECTION)
if(nrow(all_tested_models_df) > 0) {
  cat("Saving ALL tested model details to:", OUTPUT_DIR, "\n")
  tryCatch({  
    saveRDS(all_tested_models_df, file = ALL_TESTED_MODELS_RDS_FILE)
    cat("All tested model details saved to:", ALL_TESTED_MODELS_RDS_FILE, "\n")
  }, error = function(e) { cat("ERROR saving all tested models RDS file:", conditionMessage(e), "\n") })
  
  tryCatch({  
    write.csv(all_tested_models_df, file = ALL_TESTED_MODELS_CSV_FILE, row.names = FALSE, na = "")
    cat("All tested model details saved to:", ALL_TESTED_MODELS_CSV_FILE, "\n")
  }, error = function(e) { cat("ERROR saving all tested models CSV file:", conditionMessage(e), "\n") })
} else if (length(all_tested_models_list) == 0) { # Only print if list was empty
  cat("No data for all tested models to save.\n")
}


# --- 9. Display Summary (BEST Models) ---
cat("\n===== Fitting Summary (BEST Models - Simplified, Index-based Covariates) =====\n")
# (This section remains largely the same, displaying 'results_df')
if(nrow(results_df) > 0) {
  # Ensure the IC column used (e.g. aicc) is part of the display
  cols_to_display <- c("district_code", "district_name", "p", "q", IC, 
                       "time_secs", "n_models_stepwise", 
                       "betas_str", "alphas_str", "covariates_used_names", "status")
  # Make sure IC column actually exists in results_df
  if (!IC %in% colnames(results_df) && "aicc" %in% colnames(results_df) && IC == "aicc") {
    # Default to aicc if IC specific col not found but aicc is there and IC is aicc
  } else if (!IC %in% colnames(results_df)) {
    warning(paste("Specified IC column '", IC, "' not found in results_df for display. Defaulting to 'aicc' if available or skipping."))
    if ("aicc" %in% colnames(results_df)) {
      cols_to_display <- unique(c("district_code", "district_name", "p", "q", "aicc", 
                                  "time_secs", "n_models_stepwise", 
                                  "betas_str", "alphas_str", "covariates_used_names", "status"))
    } else {
      cols_to_display <- unique(c("district_code", "district_name", "p", "q",
                                  "time_secs", "n_models_stepwise",
                                  "betas_str", "alphas_str", "covariates_used_names", "status"))
    }
  }
  
  print(as.data.frame(results_df[,intersect(cols_to_display, colnames(results_df))]),
        row.names = FALSE, max = nrow(results_df) + 5)
} else {
  cat("No BEST model results to display.\n")
}

cat("\n===== Status Summary (BEST Models) =====\n")
# (This section remains largely the same, summarizing 'results_df$status')
if(nrow(results_df) > 0) {
  results_df$status <- as.character(results_df$status) 
  status_summary <- results_df %>%
    group_by(status) %>%
    summarise(count = n(), .groups = 'drop') %>%
    arrange(desc(count))
  print(status_summary)
} else {
  cat("No status summary to display.\n")
}

cat("\n--- Script finished --- \n")

# Calculate and print total run time
script_end_time <- Sys.time()
total_run_time <- difftime(script_end_time, script_start_time, units = "auto")
cat("Total run time:", format(total_run_time), attr(total_run_time, "units"), "\n")