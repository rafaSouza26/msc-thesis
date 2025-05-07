# --- Apply auto.ingarch to District Data with Covariates by Index (Sequential) ---
# Version: 2025-05-07 (Simplified with index-based covariates)

# Record script start time
script_start_time <- Sys.time()

# --- 1. Load required packages ---
library(tscount)
library(dplyr)
library(readxl)

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
IC <- "aicc"
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
RESULTS_RDS_FILE <- file.path(OUTPUT_DIR, "district_fitting_results_simplified.rds")
RESULTS_CSV_FILE <- file.path(OUTPUT_DIR, "district_fitting_summary_simplified.csv")

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

# Initialize results list
results_list <- list()

# Sequential loop over each district
for (i in DISTRICT_CODES_TO_PROCESS) {
  cat("\nProcessing District", i, "(", district_lookup[as.character(i)], ")...\n")
  
  # Initialize result for this district
  iter_result <- list(
    district_code = i, p = NA_integer_, q = NA_integer_,
    aic = NA_real_, bic = NA_real_, aicc = NA_real_,
    time_secs = NA_real_, n_models_stepwise = 0,  # Initialize to 0 instead of NA
    betas_str = NA_character_, alphas_str = NA_character_,
    covariates_used_names = NA_character_,
    status = "Not Run"
  )
  
  df_name <- paste0("d_", i, "_data")
  xreg <- NULL
  cov_names_for_this_district <- character(0)
  
  # Get the current district's data
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
      # Extract response variable (count)
      y <- current_data[, COUNT_VAR_INDEX]
      
      # Extract covariates using indices
      if(max(COVARIATE_INDICES) <= ncol(current_data)) {
        xreg <- as.matrix(current_data[, COVARIATE_INDICES])
        cov_names_for_this_district <- colnames(current_data)[COVARIATE_INDICES]
        iter_result$covariates_used_names <- paste(cov_names_for_this_district, collapse=", ")
      } else {
        iter_result$status <- "Covariate indices out of bounds for this dataframe."
        xreg <- NULL # Cannot form xreg
      }
      
      # Check for NA values
      if (anyNA(y)) {
        iter_result$status <- "NA values in response (y)"
      } else if (!is.null(xreg) && anyNA(xreg)) {
        iter_result$status <- "NA values in covariates (xreg)"
      } else if (iter_result$status == "Not Run") {
        # Fit the model
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
        
        # Extract results from the model if successful
        if (fit_status == "Success" && !is.null(fit_model) && inherits(fit_model, "tsglm")) {
          # Extract p and q
          iter_result$p <- if (is.null(fit_model$model$past_obs)) 0 else length(fit_model$model$past_obs)
          iter_result$q <- if (is.null(fit_model$model$past_mean)) 0 else length(fit_model$model$past_mean)
          
          # Extract information criteria
          iter_result$aic <- tryCatch(stats::AIC(fit_model), error = function(e) NA_real_)
          iter_result$bic <- tryCatch(stats::BIC(fit_model), error = function(e) NA_real_)
          
          if (!is.null(fit_model[[IC]])) {
            iter_result$aicc <- fit_model[[IC]]
          } else if (!is.null(fit_model$results) && IC %in% names(fit_model$results)) {
            best_model_row <- which.min(fit_model$results[[IC]])
            if(length(best_model_row) == 1) iter_result$aicc <- fit_model$results[[IC]][best_model_row]
          }
          
          # Count number of models evaluated
          if(!is.null(fit_model$results) && inherits(fit_model$results, "data.frame")) {
            iter_result$n_models_stepwise <- nrow(fit_model$results)
          } else if(!is.null(fit_model$n_total_models)) {
            # Backup method to get number of models evaluated
            iter_result$n_models_stepwise <- fit_model$n_total_models
          } else {
            # Set a default value instead of NA
            iter_result$n_models_stepwise <- 1
          }
          
          # Extract beta and alpha coefficients
          all_coeffs <- stats::coef(fit_model)
          
          # Beta coefficients (past observations)
          if (iter_result$p > 0) {
            beta_coef_names <- paste0("beta_", 1:iter_result$p)
            actual_beta_coeffs <- all_coeffs[names(all_coeffs) %in% beta_coef_names]
            iter_result$betas_str <- if (length(actual_beta_coeffs) > 0) 
              paste(round(actual_beta_coeffs, 5), collapse=",") else ""
          } else { 
            iter_result$betas_str <- "" 
          }
          
          # Alpha coefficients (past means)
          if (iter_result$q > 0) {
            alpha_coef_names <- paste0("alpha_", 1:iter_result$q)
            actual_alpha_coeffs <- all_coeffs[names(all_coeffs) %in% alpha_coef_names]
            iter_result$alphas_str <- if (length(actual_alpha_coeffs) > 0) 
              paste(round(actual_alpha_coeffs, 5), collapse=",") else ""
          } else { 
            iter_result$alphas_str <- "" 
          }
          
          cat("  Model fitted successfully: p=", iter_result$p, ", q=", iter_result$q, 
              ", AICc=", round(iter_result$aicc, 2), 
              ", Models evaluated=", iter_result$n_models_stepwise, "\n", sep="")
        } else {
          cat("  Model fitting failed with status:", fit_status, "\n")
        }
      }
    }
  }
  
  # Add this district's results to the overall list
  results_list[[length(results_list) + 1]] <- iter_result
  cat("Completed District", i, "with status:", iter_result$status, "\n")
}

# --- 7. Process and Combine Results ---
cat("\nProcessing results...\n")

if (length(results_list) > 0) {
  # Convert list of results to data frame
  results_df <- bind_rows(lapply(results_list, function(res) {
    # Ensure n_models_stepwise is numeric before binding
    if ("n_models_stepwise" %in% names(res) && is.na(res$n_models_stepwise)) {
      res$n_models_stepwise <- 0
    }
    as.data.frame(res, stringsAsFactors = FALSE)
  }))
  
  # Add district names
  results_df <- results_df %>%
    mutate(district_name = ifelse(!is.na(district_code), district_lookup[as.character(district_code)], NA_character_)) %>%
    select(district_code, district_name, p, q, aic, bic, aicc, time_secs,
           n_models_stepwise, betas_str, alphas_str, covariates_used_names, status)
} else {
  # Create empty structure if no results
  results_df <- data.frame(
    district_code = integer(0), district_name = character(0), p = integer(0), q = integer(0),
    aic = numeric(0), bic = numeric(0), aicc = numeric(0),
    time_secs = numeric(0), n_models_stepwise = integer(0),
    betas_str = character(0), alphas_str = character(0),
    covariates_used_names = character(0), status = character(0),
    stringsAsFactors = FALSE  # Ensure strings don't get converted to factors
  )
  cat("No results to process or save.\n")
}

# --- 8. Save Results ---
if(nrow(results_df) > 0) {
  cat("Saving results to:", OUTPUT_DIR, "\n")
  tryCatch({ 
    saveRDS(results_df, file = RESULTS_RDS_FILE)
    cat("Detailed results saved to:", RESULTS_RDS_FILE, "\n")
  }, error = function(e) { cat("ERROR saving RDS file:", conditionMessage(e), "\n") })
  
  tryCatch({ 
    # Ensure n_models_stepwise is numeric before writing to CSV
    results_df$n_models_stepwise <- as.numeric(results_df$n_models_stepwise)
    # Replace NA with 0 for n_models_stepwise
    results_df$n_models_stepwise[is.na(results_df$n_models_stepwise)] <- 0
    
    write.csv(results_df, file = RESULTS_CSV_FILE, row.names = FALSE, na = "")
    cat("Summary results saved to:", RESULTS_CSV_FILE, "\n")
  }, error = function(e) { cat("ERROR saving CSV file:", conditionMessage(e), "\n") })
}

# --- 9. Display Summary ---
cat("\n===== Fitting Summary (Simplified, Index-based Covariates) =====\n")
if(nrow(results_df) > 0) {
  print(as.data.frame(results_df[, c("district_code", "district_name", "p", "q", "aicc",
                                     "time_secs", "betas_str", "alphas_str", "covariates_used_names", "status")]),
        row.names = FALSE, max = nrow(results_df) + 5)
} else {
  cat("No results to display.\n")
}

cat("\n===== Status Summary =====\n")
if(nrow(results_df) > 0) {
  results_df$status <- as.character(results_df$status) # Ensure character for grouping
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