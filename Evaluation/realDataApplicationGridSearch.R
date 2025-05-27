# --- Apply auto.ingarch to District Data with Covariates by Index (Sequential) ---
# Version: 2025-05-16 (Grid Search, Single Detailed CSV output, AIC as standard IC)
# Version: 2025-05-23 (Grid Search, Compact CSV with Base Covariate Names)

# Record script start time
script_start_time <- Sys.time()

# --- 1. Load required packages ---
library(tscount)
library(dplyr)
library(readxl) # Kept if needed by custom functions or future use

# --- 2. Source custom functions ---
cat("Sourcing custom functions...\n")
tryCatch(source("./ACTS/auto.ingarch.R"), error = function(e) stop("Failed to source auto.ingarch.R: ", e))
tryCatch(source("./ACTS/newmodel.R"), error = function(e) stop("Failed to source newmodel.R: ", e))
tryCatch(source("./ACTS/ingarch.string.R"), error = function(e) stop("Failed to source ingarch.string.R: ", e))
tryCatch(source("./ACTS/search.ingarch.R"), error = function(e) stop("Failed to source search.ingarch.R: ", e))
tryCatch(source("./ACTS/myingarch.R"), error = function(e) stop("Failed to source myingarch.R: ", e))


# --- 3. Define Parameters and Setup ---
cat("Defining parameters...\n")
MAX_P <- 7
MAX_Q <- 7
DISTRIBUTION <- "nbinom"
LINK <- "log"
IC <- "aic" # Ensure AIC is used
STEPWISE <- FALSE # Changed for Grid Search
MAX_ORDER_GRID_SEARCH <- MAX_P + MAX_Q # Max sum of p and q for grid search
TRACE <- FALSE
SHOW_WARNINGS <- FALSE

DATA_PATH <- "./data/count_covariates_data.RData"
COUNT_VAR_NAME <- "Count"
COUNT_VAR_INDEX <- 2
COVARIATE_INDICES <- 3:6
MIN_EXPECTED_COLUMNS <- max(c(COUNT_VAR_INDEX, COVARIATE_INDICES))

# <<< NEW: Define Base Covariate Names for fixed columns >>>
BASE_COVARIATE_NAMES <- c("mean.Temp", "mean.DTemp", "PM25", "mean.NO2")

OUTPUT_DIR <- "./real_data_application_gs" # Changed for Grid Search output
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}
DETAILED_RESULTS_CSV_FILE <- file.path(OUTPUT_DIR, "district_fitting_detailed_results_gridsearch.csv")


NUM_DISTRICTS_TO_PROCESS <- 18
ALL_DISTRICT_CODES_AVAILABLE <- 1:18

if (NUM_DISTRICTS_TO_PROCESS > length(ALL_DISTRICT_CODES_AVAILABLE)) {
  warning(paste("NUM_DISTRICTS_TO_PROCESS (", NUM_DISTRICTS_TO_PROCESS,
                ") is greater than available districts. Processing all."))
  NUM_DISTRICTS_TO_PROCESS <- length(ALL_DISTRICT_CODES_AVAILABLE)
}
if (NUM_DISTRICTS_TO_PROCESS < 1) {
  warning("NUM_DISTRICTS_TO_PROCESS < 1. Setting to 1.")
  NUM_DISTRICTS_TO_PROCESS <- 1
}
DISTRICT_CODES_TO_PROCESS <- ALL_DISTRICT_CODES_AVAILABLE[1:NUM_DISTRICTS_TO_PROCESS]
cat("Processing", NUM_DISTRICTS_TO_PROCESS, "district(s) using GRID SEARCH. Codes:", paste(DISTRICT_CODES_TO_PROCESS, collapse=", "), "\n")

DISTRICT_NAMES <- c("Aveiro", "Beja", "Braga", "Bragança", "Castelo Branco",
                    "Coimbra", "Évora", "Faro", "Guarda", "Leiria", "Lisboa",
                    "Portalegre", "Porto", "Santarém", "Setúbal",
                    "Viana do Castelo", "Vila Real", "Viseu")
district_lookup <- setNames(DISTRICT_NAMES, ALL_DISTRICT_CODES_AVAILABLE)

# <<< NEW: Helper function to get base covariate name >>>
get_base_covariate_name <- function(full_cov_name) {
  # Attempts to extract base name like "mean.Temp" from "mean.Temp.Lag5"
  base_name_candidate <- sub("\\.Lag[0-9A-Za-z]*$", "", full_cov_name)
  if (base_name_candidate %in% BASE_COVARIATE_NAMES) {
    return(base_name_candidate)
  } else {
    # Fallback for more complex names or if primary sub fails
    for(bcn in BASE_COVARIATE_NAMES){
      if(startsWith(full_cov_name, bcn)){ # Check if full name starts with a known base name
        return(bcn)
      }
    }
    warning(paste("Could not map '", full_cov_name, "' to a known base covariate name. Returning original.", sep=""))
    return(full_cov_name) # Return original if no mapping found
  }
}

cat("Loading district data from:", DATA_PATH, "\n")
if (!file.exists(DATA_PATH)) stop("Data file not found: ", DATA_PATH)
load(DATA_PATH, envir = .GlobalEnv)
cat("Data loaded.\n")

cat("Verifying data objects...\n")
# ... (Data verification logic remains largely the same) ...
missing_data_objs <- c()
all_df_names_to_process <- paste0("d_", DISTRICT_CODES_TO_PROCESS, "_data")
all_struct_ok <- TRUE
for (df_name in all_df_names_to_process) {
  if (!exists(df_name, envir = .GlobalEnv)) {
    missing_data_objs <- c(missing_data_objs, df_name); all_struct_ok <- FALSE
  } else {
    df_temp <- get(df_name, envir = .GlobalEnv)
    if (!is.data.frame(df_temp)){ cat("ERROR: '", df_name, "' not a DF.\n"); all_struct_ok <- FALSE
    } else if (ncol(df_temp) < MIN_EXPECTED_COLUMNS) { cat("ERROR: '", df_name, "' has ", ncol(df_temp), " cols, need >=", MIN_EXPECTED_COLUMNS, ".\n"); all_struct_ok <- FALSE
    } else if (names(df_temp)[COUNT_VAR_INDEX] != COUNT_VAR_NAME) { cat("ERROR: '", df_name, "' col ", COUNT_VAR_INDEX, " is '", names(df_temp)[COUNT_VAR_INDEX], "' not '", COUNT_VAR_NAME, "'.\n"); all_struct_ok <- FALSE
    }
  }
}
if (length(missing_data_objs) > 0) stop("Missing DFs: ", paste(missing_data_objs, collapse=", "))
if (!all_struct_ok) stop("Data structure issues. Check messages.")
cat("Data structure check passed.\n")


cat("Starting sequential model fitting (Grid Search) for", length(DISTRICT_CODES_TO_PROCESS), "districts...\n")
results_list <- list()

base_iter_result_cols <- c(
  "district_code", "p", "q", "aic_value", "bic_value", "time_secs",
  "n_models_tested", "intercept", "sigmasq", "betas_str", "alphas_str",
  "covariates_used_names", "status", "error_message"
)

for (i in DISTRICT_CODES_TO_PROCESS) {
  cat("\nProcessing District", i, "(", district_lookup[as.character(i)], ")...\n")
  
  # <<< MODIFIED: Initialize iter_result to include BASE_COVARIATE_NAMES >>>
  iter_result <- vector("list", length(base_iter_result_cols) + length(BASE_COVARIATE_NAMES))
  names(iter_result) <- c(base_iter_result_cols, BASE_COVARIATE_NAMES)
  
  iter_result[sapply(iter_result, is.null)] <- NA # General NA init
  
  # Specific NA types for base_iter_result_cols
  iter_result$district_code <- i # Set district code early
  iter_result$p <- NA_integer_
  iter_result$q <- NA_integer_
  iter_result$aic_value <- NA_real_
  iter_result$bic_value <- NA_real_
  iter_result$time_secs <- NA_real_
  iter_result$n_models_tested <- 0L 
  iter_result$intercept <- NA_real_
  iter_result$sigmasq <- NA_real_
  iter_result$betas_str <- ""      # NA_character_ typically
  iter_result$alphas_str <- ""     # NA_character_
  iter_result$covariates_used_names <- NA_character_
  iter_result$status <- "Not Run"
  iter_result$error_message <- ""  # NA_character_
  
  # Initialize BASE_COVARIATE_NAMES columns to NA_real_
  for(base_cov_name in BASE_COVARIATE_NAMES) {
    iter_result[[base_cov_name]] <- NA_real_
  }
  # End specific NA initializations
  
  df_name <- paste0("d_", i, "_data")
  y <- NULL
  xreg <- NULL
  cov_names_for_this_district <- character(0) # Full original names with lags
  
  current_data <- tryCatch(get(df_name, envir = .GlobalEnv), error = function(e) {
    iter_result$status <<- "Data Load Error"; iter_result$error_message <<- conditionMessage(e); NULL
  })
  
  if (is.null(current_data) || !is.data.frame(current_data)) {
    if(iter_result$status == "Not Run") iter_result$status <- "Data Not Found/Invalid"
    results_list[[length(results_list) + 1]] <- iter_result
    cat("  Skipping district", i, "due to data loading issue.\n")
    next
  }
  
  if (ncol(current_data) < MIN_EXPECTED_COLUMNS) {
    iter_result$status <- "Insufficient columns"
    iter_result$error_message <- "Data frame has too few columns."
  } else {
    y <- current_data[, COUNT_VAR_INDEX]
    
    if(max(COVARIATE_INDICES) <= ncol(current_data) && length(COVARIATE_INDICES) > 0) {
      xreg_candidate <- current_data[, COVARIATE_INDICES, drop = FALSE]
      if(all(sapply(xreg_candidate, is.numeric))) {
        xreg <- as.matrix(xreg_candidate)
        cov_names_for_this_district <- colnames(xreg) # These are full names like "mean.Temp.Lag5"
        if(is.null(cov_names_for_this_district) || any(cov_names_for_this_district == "" | is.na(cov_names_for_this_district))) {
          warning(paste("District", i, ": Covariate names missing/invalid. Using generic names."))
          cov_names_for_this_district <- paste0("xreg_coef", 1:ncol(xreg))
          colnames(xreg) <- cov_names_for_this_district
        }
        iter_result$covariates_used_names <- paste(cov_names_for_this_district, collapse=", ")
        # NO dynamic iter_result[[cn]] addition here for full names
      } else {
        iter_result$status <- "Non-numeric Covariates"; iter_result$error_message <- "Covariates non-numeric."; xreg <- NULL
      }
    } else if (length(COVARIATE_INDICES) == 0) {
      xreg <- NULL; iter_result$covariates_used_names <- ""
    } else {
      iter_result$status <- "Covariate Indices OOB"; iter_result$error_message <- "Covariate indices OOB."; xreg <- NULL
    }
    
    if (iter_result$status == "Not Run") {
      if (anyNA(y)) {
        iter_result$status <- "NA in y"; iter_result$error_message <- "NA values in response."
      } else if (!is.null(xreg) && anyNA(xreg)) {
        iter_result$status <- "NA in xreg"; iter_result$error_message <- "NA values in covariates."
      }
    }
  }
  
  if (iter_result$status == "Not Run") {
    fit_status_local <- "Fit Failed"; error_msg_fit_local <- ""
    start_time_fit <- Sys.time()
    
    fit_model <- tryCatch({
      cat("  Fitting auto.ingarch (stepwise=", STEPWISE, ", ic=", IC, ", max.order=", MAX_ORDER_GRID_SEARCH, ")...\n")
      model <- auto.ingarch(
        y = y, xreg = xreg, max.p = MAX_P, max.q = MAX_Q,
        distribution = DISTRIBUTION, link = LINK, ic = IC,
        stepwise = STEPWISE, # This is FALSE for grid search
        max.order = MAX_ORDER_GRID_SEARCH,
        trace = TRACE, show_warnings = SHOW_WARNINGS
      )
      if (!inherits(model, "tsglm")) stop("auto.ingarch did not return a valid tsglm model object")
      fit_status_local <- "Success"
      model
    }, error = function(e) {
      error_msg_fit_local <<- conditionMessage(e)
      cat("  Model fitting error for district", i, ":", error_msg_fit_local, "\n")
      NULL
    })
    
    end_time_fit <- Sys.time()
    iter_result$time_secs <- as.numeric(difftime(end_time_fit, start_time_fit, units = "secs"))
    iter_result$status <- if(is.null(fit_model)) "Fitting Error" else fit_status_local # Simplified
    iter_result$error_message <- if(is.null(fit_model)) substr(error_msg_fit_local,1,250) else ""
    
    if (iter_result$status == "Success" && !is.null(fit_model)) {
      iter_result$p <- if(is.null(fit_model$model$past_obs)) 0L else length(fit_model$model$past_obs)
      iter_result$q <- if(is.null(fit_model$model$past_mean)) 0L else length(fit_model$model$past_mean)
      iter_result$aic_value <- tryCatch(stats::AIC(fit_model), error = function(e) NA_real_)
      iter_result$bic_value <- tryCatch(stats::BIC(fit_model), error = function(e) NA_real_)
      
      iter_result$n_models_tested <- if(!is.null(fit_model$n_models_evaluated)) {
        fit_model$n_models_evaluated
      } else {
        (MAX_P + 1) * (MAX_Q + 1) # Default for full grid if not specified
      }
      
      all_coeffs <- tryCatch(stats::coef(fit_model), error = function(e) NULL)
      if(!is.null(all_coeffs)){
        iter_result$intercept <- if("intercept" %in% names(all_coeffs)) all_coeffs[["intercept"]] else if("(Intercept)" %in% names(all_coeffs)) all_coeffs[["(Intercept)"]] else NA_real_
        # ... (beta_str and alpha_str logic as before)
        if (iter_result$p > 0) {
          beta_coefs <- sapply(1:iter_result$p, function(k) {
            b_name1 <- paste0("beta_", k); b_name2 <- paste0("past_obs",k)
            if(b_name1 %in% names(all_coeffs)) all_coeffs[[b_name1]] else if(b_name2 %in% names(all_coeffs)) all_coeffs[[b_name2]] else NA_real_
          })
          iter_result$betas_str <- paste(round(na.omit(beta_coefs), 5), collapse=",")
        } # else it remains ""
        if (iter_result$q > 0) {
          alpha_coefs <- sapply(1:iter_result$q, function(k) {
            a_name1 <- paste0("alpha_", k); a_name2 <- paste0("past_mean",k)
            if(a_name1 %in% names(all_coeffs)) all_coeffs[[a_name1]] else if(a_name2 %in% names(all_coeffs)) all_coeffs[[a_name2]] else NA_real_
          })
          iter_result$alphas_str <- paste(round(na.omit(alpha_coefs), 5), collapse=",")
        } # else it remains ""
        
        # <<< MODIFIED: Assign covariate coefficients to base name columns >>>
        if (length(cov_names_for_this_district) > 0 && !is.null(xreg)) {
          for (k_xreg_name_original in cov_names_for_this_district) { # e.g., "mean.Temp.Lag5"
            if (k_xreg_name_original %in% names(all_coeffs)) {
              coeff_value <- all_coeffs[[k_xreg_name_original]]
              base_name <- get_base_covariate_name(k_xreg_name_original) # e.g., "mean.Temp"
              if (base_name %in% BASE_COVARIATE_NAMES) {
                iter_result[[base_name]] <- coeff_value
              } else {
                # This case is handled by get_base_covariate_name warning if mapping fails to a known base name.
                # If base_name was returned as k_xreg_name_original and is not in BASE_COVARIATE_NAMES:
                if (k_xreg_name_original == base_name) { 
                  cat("  Info: Coef for '", k_xreg_name_original, "' not mapped to a standard base column. Value not stored in fixed columns.\n", sep="")
                }
              }
            }
          }
        }
      }
      
      if(fit_model$distr == "nbinom"){
        iter_result$sigmasq <- if(!is.null(fit_model$sigmasq)) fit_model$sigmasq else if (!is.null(fit_model$distrcoefs$size) && fit_model$distrcoefs$size > 0) 1/fit_model$distrcoefs$size else NA_real_
      } # else it remains NA_real_ (or whatever it was initialized to if not nbinom)
      cat("  Model fitted successfully: p=", iter_result$p, ", q=", iter_result$q,
          ", Selected IC (", IC, ")=", round(iter_result$aic_value, 2),
          ", Models evaluated=", iter_result$n_models_tested, "\n", sep="")
    } else {
      cat("  Model fitting failed for district", i, "with status:", iter_result$status, "\n")
    }
  }
  results_list[[length(results_list) + 1]] <- iter_result
  cat("Completed District", i, "with status:", iter_result$status, "\n")
}

# --- Process Results ---
cat("\nProcessing results...\n")

# <<< MODIFIED: Define final column order using BASE_COVARIATE_NAMES >>>
# all_possible_xreg_cols is now effectively BASE_COVARIATE_NAMES for fixed columns
final_col_order <- c("district_code", "district_name", "p", "q", "aic_value", "bic_value", "time_secs",
                     "n_models_tested", "intercept", "sigmasq", "betas_str", "alphas_str",
                     BASE_COVARIATE_NAMES, # Use the predefined base names
                     "covariates_used_names", "status", "error_message")
final_col_order <- unique(final_col_order)

if (length(results_list) > 0) {
  results_df <- bind_rows(lapply(results_list, function(res_item) {
    res_item_as_list <- as.list(res_item)
    for(col_n in final_col_order){
      if(!col_n %in% names(res_item_as_list) && col_n != "district_name") {
        if(col_n %in% c("p","q","n_models_tested","district_code")) {
          res_item_as_list[[col_n]] <- NA_integer_
        } else if(col_n %in% c("aic_value","bic_value","time_secs","intercept","sigmasq") || col_n %in% BASE_COVARIATE_NAMES) { # Check against base names
          res_item_as_list[[col_n]] <- NA_real_
        } else {
          res_item_as_list[[col_n]] <- NA_character_
        }
      }
    }
    as.data.frame(res_item_as_list, stringsAsFactors = FALSE)
  }))
  
  results_df <- results_df %>%
    mutate(district_name = ifelse(!is.na(district_code), district_lookup[as.character(district_code)], NA_character_))
  
  for(col_to_ensure in final_col_order){
    if(!col_to_ensure %in% names(results_df)){
      cat("Notice: Column '", col_to_ensure, "' from final_col_order was not in results_df. Adding it now.\n", sep="")
      if(col_to_ensure %in% c("p","q","n_models_tested","district_code")) results_df[[col_to_ensure]] <- NA_integer_
      else if(col_to_ensure %in% c("aic_value","bic_value","time_secs","intercept","sigmasq") || col_to_ensure %in% BASE_COVARIATE_NAMES) results_df[[col_to_ensure]] <- NA_real_
      else results_df[[col_to_ensure]] <- NA_character_
    }
  }
  # Ensure correct column order, selecting only existing columns from final_col_order
  actual_cols_to_select <- intersect(final_col_order, names(results_df))
  results_df <- results_df[, actual_cols_to_select, drop = FALSE]
  
} else {
  # Create an empty data frame with the correct column structure
  df_cols <- lapply(final_col_order, function(col_name) {
    if(col_name %in% c("p","q","n_models_tested","district_code")) return(integer(0))
    if(col_name %in% c("aic_value","bic_value","time_secs","intercept","sigmasq") || col_name %in% BASE_COVARIATE_NAMES) return(numeric(0))
    return(character(0))
  })
  names(df_cols) <- final_col_order
  results_df <- as.data.frame(df_cols, stringsAsFactors = FALSE)
  cat("No results to process or save.\n")
}

# --- Save Detailed Results CSV ---
if(nrow(results_df) > 0) {
  cat("Saving detailed results to:", OUTPUT_DIR, "\n")
  tryCatch({
    write.csv(results_df, file = DETAILED_RESULTS_CSV_FILE, row.names = FALSE, na = "")
    cat("Detailed results saved to CSV:", DETAILED_RESULTS_CSV_FILE, "\n")
  }, error = function(e) { cat("ERROR saving detailed CSV file:", conditionMessage(e), "\n") })
}

# --- Display Summary ---
cat("\n===== Fitting Summary (Grid Search, IC=", IC, ") =====\n")
if(nrow(results_df) > 0) {
  cols_to_print <- c("district_code", "district_name", "p", "q", "aic_value", "time_secs", "status")
  # <<< MODIFIED: Use BASE_COVARIATE_NAMES for printing summary >>>
  existing_xreg_cols_in_df <- intersect(BASE_COVARIATE_NAMES, names(results_df))
  for(dyn_col in existing_xreg_cols_in_df){ # dyn_col is now a base name
    if(any(!is.na(results_df[[dyn_col]]))) {
      cols_to_print <- c(cols_to_print, dyn_col)
    }
  }
  cols_to_print_existing <- intersect(cols_to_print, names(results_df))
  if(length(cols_to_print_existing) > 0){
    print(results_df[, cols_to_print_existing, drop = FALSE], row.names = FALSE, max = nrow(results_df) + 5)
  } else {
    cat("No columns selected for printing summary or results_df is empty.\n")
  }
} else {
  cat("No results to display.\n")
}

cat("\n===== Status Summary =====\n")
if(nrow(results_df) > 0 && "status" %in% names(results_df)) { # Check status column exists
  results_df$status <- as.character(results_df$status)
  status_summary <- results_df %>%
    group_by(status) %>%
    summarise(count = n(), .groups = 'drop') %>%
    arrange(desc(count))
  print(status_summary)
} else {
  cat("No status summary to display (results_df empty or 'status' column missing).\n")
}

cat("\n--- Script finished --- \n")
script_end_time <- Sys.time()
total_run_time <- difftime(script_end_time, script_start_time, units = "auto")
cat("Total run time:", format(total_run_time), attr(total_run_time, "units"), "\n")