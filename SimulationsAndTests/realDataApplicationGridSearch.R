# --- Apply auto.ingarch to District Data with Covariates by Index (Sequential) ---
# Version: 2025-05-16 (Grid Search, Single Detailed CSV output, AIC as standard IC)

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
tryCatch(source("./custom/myingarch.R"), error = function(e) stop("Failed to source myingarch.R: ", e))


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
SHOW_WARNINGS <- FALSE # Set to TRUE to see warnings from auto.ingarch fits

DATA_PATH <- "./data/count_covariates_data.RData"
COUNT_VAR_NAME <- "Count"
COUNT_VAR_INDEX <- 2
COVARIATE_INDICES <- 3:6
MIN_EXPECTED_COLUMNS <- max(c(COUNT_VAR_INDEX, COVARIATE_INDICES))

OUTPUT_DIR <- "./real_data_application_gs" # Changed for Grid Search output
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}
# Only one detailed CSV file will be saved
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

cat("Loading district data from:", DATA_PATH, "\n")
if (!file.exists(DATA_PATH)) stop("Data file not found: ", DATA_PATH)
load(DATA_PATH, envir = .GlobalEnv)
cat("Data loaded.\n")

cat("Verifying data objects...\n")
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
  
  iter_result <- vector("list", length(base_iter_result_cols))
  names(iter_result) <- base_iter_result_cols
  iter_result[sapply(iter_result, is.null)] <- NA 
  iter_result$district_code <- i
  iter_result$p <- NA_integer_
  iter_result$q <- NA_integer_
  iter_result$n_models_tested <- 0L # Will be (MAX_P+1)*(MAX_Q+1) for successful grid search
  iter_result$status <- "Not Run"
  iter_result$error_message <- ""
  iter_result$betas_str <- ""
  iter_result$alphas_str <- ""
  
  df_name <- paste0("d_", i, "_data")
  y <- NULL 
  xreg <- NULL 
  cov_names_for_this_district <- character(0)
  
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
    
    if(max(COVARIATE_INDICES) <= ncol(current_data)) {
      xreg <- as.matrix(current_data[, COVARIATE_INDICES, drop = FALSE])
      cov_names_for_this_district <- colnames(xreg)
      if(is.null(cov_names_for_this_district) || any(cov_names_for_this_district == "")) { 
        cov_names_for_this_district <- paste0("xreg_coef", 1:ncol(xreg))
        colnames(xreg) <- cov_names_for_this_district 
      }
      iter_result$covariates_used_names <- paste(cov_names_for_this_district, collapse=", ")
      for(cn in make.names(cov_names_for_this_district, unique=TRUE)) iter_result[[cn]] <- NA_real_
    } else {
      iter_result$status <- "Covariate Indices OOB"
      iter_result$error_message <- "Covariate indices out of bounds, xreg not formed."
      xreg <- NULL 
    }
    
    if (iter_result$status == "Not Run") {
      if (anyNA(y)) {
        iter_result$status <- "NA in y"
        iter_result$error_message <- "NA values found in response variable."
      } else if (!is.null(xreg) && anyNA(xreg)) { 
        iter_result$status <- "NA in xreg"
        iter_result$error_message <- "NA values found in covariates."
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
        max.order = MAX_ORDER_GRID_SEARCH, # Pass max.order for grid search
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
    iter_result$status <- if(is.null(fit_model)) paste("Fitting Error") else fit_status_local
    iter_result$error_message <- if(is.null(fit_model)) substr(error_msg_fit_local,1,250) else ""
    
    if (iter_result$status == "Success" && !is.null(fit_model)) {
      iter_result$p <- if(is.null(fit_model$model$past_obs)) 0L else length(fit_model$model$past_obs)
      iter_result$q <- if(is.null(fit_model$model$past_mean)) 0L else length(fit_model$model$past_mean)
      
      iter_result$aic_value <- tryCatch(stats::AIC(fit_model), error = function(e) NA_real_)
      iter_result$bic_value <- tryCatch(stats::BIC(fit_model), error = function(e) NA_real_)
      
      # For grid search, n_models_tested is (MAX_P+1)*(MAX_Q+1) if search.ingarch completed
      # auto.ingarch with stepwise=FALSE calls search.ingarch, which should evaluate all these.
      # The n_models_evaluated from search.ingarch might be stored in fit_model$n_models_evaluated
      iter_result$n_models_tested <- if(!is.null(fit_model$n_models_evaluated)) {
        fit_model$n_models_evaluated
      } else {
        (MAX_P + 1) * (MAX_Q + 1) # Theoretical max for full grid if not otherwise specified
      }
      
      
      all_coeffs <- tryCatch(stats::coef(fit_model), error = function(e) NULL)
      if(!is.null(all_coeffs)){
        iter_result$intercept <- if("intercept" %in% names(all_coeffs)) all_coeffs[["intercept"]] else if("(Intercept)" %in% names(all_coeffs)) all_coeffs[["(Intercept)"]] else NA_real_
        
        if (iter_result$p > 0) {
          beta_coefs <- sapply(1:iter_result$p, function(k) {
            b_name1 <- paste0("beta_", k); b_name2 <- paste0("past_obs",k)
            if(b_name1 %in% names(all_coeffs)) all_coeffs[[b_name1]]
            else if(b_name2 %in% names(all_coeffs)) all_coeffs[[b_name2]]
            else NA_real_
          })
          iter_result$betas_str <- paste(round(na.omit(beta_coefs), 5), collapse=",")
        } else iter_result$betas_str <- ""
        
        if (iter_result$q > 0) {
          alpha_coefs <- sapply(1:iter_result$q, function(k) {
            a_name1 <- paste0("alpha_", k); a_name2 <- paste0("past_mean",k)
            if(a_name1 %in% names(all_coeffs)) all_coeffs[[a_name1]]
            else if(a_name2 %in% names(all_coeffs)) all_coeffs[[a_name2]]
            else NA_real_
          })
          iter_result$alphas_str <- paste(round(na.omit(alpha_coefs), 5), collapse=",")
        } else iter_result$alphas_str <- ""
        
        if (length(cov_names_for_this_district) > 0 && !is.null(xreg)) { 
          for (k_xreg_name in cov_names_for_this_district) {
            valid_col_name <- make.names(k_xreg_name, unique=TRUE) 
            iter_result[[valid_col_name]] <- if(k_xreg_name %in% names(all_coeffs)) all_coeffs[[k_xreg_name]] else NA_real_
          }
        }
      } 
      
      if(fit_model$distr == "nbinom"){
        iter_result$sigmasq <- if(!is.null(fit_model$sigmasq)) fit_model$sigmasq else if (!is.null(fit_model$distrcoefs$size) && fit_model$distrcoefs$size > 0) 1/fit_model$distrcoefs$size else NA_real_
      } else {
        iter_result$sigmasq <- NA_real_
      }
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

cat("\nProcessing results...\n")

all_possible_xreg_cols <- character(0)
if(exists("d_1_data")){ 
  if(max(COVARIATE_INDICES) <= ncol(d_1_data)) {
    all_possible_xreg_cols <- make.names(colnames(d_1_data[,COVARIATE_INDICES, drop=FALSE]), unique=TRUE)
  }
}
final_col_order <- c("district_code", "district_name", "p", "q", "aic_value", "bic_value", "time_secs", 
                     "n_models_tested", "intercept", "sigmasq", "betas_str", "alphas_str", 
                     all_possible_xreg_cols, 
                     "covariates_used_names", "status", "error_message")
final_col_order <- unique(final_col_order)

if (length(results_list) > 0) {
  results_df <- bind_rows(lapply(results_list, function(res_item) {
    for(col_n in final_col_order){
      if(!col_n %in% names(res_item) && col_n != "district_name") { 
        res_item[[col_n]] <- NA 
        if(col_n %in% c("p","q","n_models_tested","district_code")) res_item[[col_n]] <- NA_integer_
        else if(col_n %in% c("aic_value","bic_value","time_secs","intercept","sigmasq") || col_n %in% all_possible_xreg_cols) res_item[[col_n]] <- NA_real_
        else res_item[[col_n]] <- NA_character_ 
      }
    }
    as.data.frame(res_item, stringsAsFactors = FALSE)
  }))
  
  results_df <- results_df %>%
    mutate(district_name = ifelse(!is.na(district_code), district_lookup[as.character(district_code)], NA_character_))
  
  for(col_to_ensure in final_col_order){
    if(!col_to_ensure %in% names(results_df)){
      if(col_to_ensure %in% c("p","q","n_models_tested","district_code")) results_df[[col_to_ensure]] <- NA_integer_
      else if(col_to_ensure %in% c("aic_value","bic_value","time_secs","intercept","sigmasq") || col_to_ensure %in% all_possible_xreg_cols) results_df[[col_to_ensure]] <- NA_real_
      else results_df[[col_to_ensure]] <- NA_character_
    }
  }
  results_df <- results_df[, final_col_order, drop = FALSE]
  
} else {
  df_cols <- lapply(final_col_order, function(col_name) {
    if(col_name %in% c("p","q","n_models_tested","district_code")) return(integer(0))
    if(col_name %in% c("aic_value","bic_value","time_secs","intercept","sigmasq") || col_name %in% all_possible_xreg_cols) return(numeric(0))
    return(character(0))
  })
  names(df_cols) <- final_col_order
  results_df <- as.data.frame(df_cols, stringsAsFactors = FALSE)
  cat("No results to process or save.\n")
}

if(nrow(results_df) > 0) {
  cat("Saving detailed results to:", OUTPUT_DIR, "\n")
  tryCatch({ 
    write.csv(results_df, file = DETAILED_RESULTS_CSV_FILE, row.names = FALSE, na = "")
    cat("Detailed results saved to CSV:", DETAILED_RESULTS_CSV_FILE, "\n")
  }, error = function(e) { cat("ERROR saving detailed CSV file:", conditionMessage(e), "\n") })
  
  # Removed saving of separate summary CSV file
}

cat("\n===== Fitting Summary (Grid Search, IC=", IC, ") =====\n")
if(nrow(results_df) > 0) {
  cols_to_print <- c("district_code", "district_name", "p", "q", "aic_value", "time_secs", "status")
  existing_xreg_cols_in_df <- intersect(all_possible_xreg_cols, names(results_df))
  for(dyn_col in existing_xreg_cols_in_df){
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
script_end_time <- Sys.time()
total_run_time <- difftime(script_end_time, script_start_time, units = "auto")
cat("Total run time:", format(total_run_time), attr(total_run_time, "units"), "\n")

