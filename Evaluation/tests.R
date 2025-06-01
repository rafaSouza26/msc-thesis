# INGARCH Model Simulation Study - Without Covariates (Sequential, DGP 1,3)

# Ensure all necessary libraries are loaded.
# dplyr and readxl might only be needed if you revert to using extract_model_params
library(tscount)
library(dplyr)
library(readxl) # Kept in case extract_model_params is used in other contexts

# Ensure these paths are correct for your environment
source("./ACTS/auto.ingarch.R")
source("./ACTS/ingarch.sim.R")
source("./ACTS/newmodel.R")
source("./ACTS/ingarch.string.R")
source("./ACTS/search.ingarch.R")
source("./ACTS/myingarch.R")

# extract_model_params function remains available if needed for other scenarios,
# but will be bypassed for the (1,3) DGP in this script.
extract_model_params <- function(data_path, model_row = 1) {
  excel_data <- readxl::read_excel(data_path)
  model_data <- excel_data[model_row, ]
  
  p <- as.numeric(model_data$p)
  q <- as.numeric(model_data$q)
  sigmasq <- as.numeric(model_data$sigmasq)
  intercept <- as.numeric(model_data$Intercept)
  
  beta_cols <- grep("^beta", colnames(excel_data), ignore.case = TRUE)
  betas <- if(length(beta_cols) > 0 && ncol(model_data[, beta_cols, drop=FALSE]) > 0) as.numeric(unlist(model_data[, beta_cols])) else numeric(0)
  
  alpha_cols <- grep("^alpha", colnames(excel_data), ignore.case = TRUE)
  alphas <- if(length(alpha_cols) > 0 && ncol(model_data[, alpha_cols, drop=FALSE]) > 0) as.numeric(unlist(model_data[, alpha_cols])) else numeric(0)
  
  return(list(
    p = p, q = q, sigmasq = sigmasq, intercept = intercept,
    betas = betas, alphas = alphas
  ))
}

run_simulation_study_sequential_1_3 <- function() {
  set.seed(12345) # For reproducibility
  
  num_simulations <- 100   # Number of simulation runs
  sim_length <- 500     # Length of each simulated time series
  progress_print_frequency <- max(1, floor(num_simulations / 10))
  
  cat("Starting INGARCH simulation study (without covariates) - SEQUENTIAL EXECUTION.\n")
  cat("Data Generating Process: INGARCH(1,3)\n")
  cat("Configuration: num_simulations =", num_simulations, ", sim_length =", sim_length, "\n")
  
  # --- Define Fixed DGP Parameters for INGARCH(1,3) ---
  true_p_dgp <- 1
  true_q_dgp <- 3
  dgp_intercept <- 0.5          # Example intercept
  dgp_betas <- c(0.2)           # Example beta_1 for p=1
  dgp_alphas <- c(0.15, 0.1, 0.05) # Example alpha_1, alpha_2, alpha_3 for q=3
  dgp_sigmasq <- 0.1            # Example sigmasq (for NB dispersion: size = 1/sigmasq)
  
  cat("\nStep 1: Using fixed Data Generating Process (DGP) parameters.\n")
  cat("  True p_order=", true_p_dgp, ", True q_order=", true_q_dgp, ", True sigmasq=", dgp_sigmasq, "\n")
  cat("  DGP Intercept=", dgp_intercept, "\n")
  cat("  DGP Betas (past_obs):", paste(dgp_betas, collapse=", "), "\n")
  cat("  DGP Alphas (past_mean):", paste(dgp_alphas, collapse=", "), "\n")
  
  ingarch_params_dgp <- list(
    intercept = dgp_intercept,
    past_obs = dgp_betas,    # Should match true_p_dgp
    past_mean = dgp_alphas   # Should match true_q_dgp
  )
  ingarch_model_dgp <- list(
    past_obs = if(true_p_dgp > 0) 1:true_p_dgp else NULL,
    past_mean = if(true_q_dgp > 0) 1:true_q_dgp else NULL,
    external = FALSE
  )
  size_param_dgp <- if(!is.na(dgp_sigmasq) && dgp_sigmasq > 0) 1/dgp_sigmasq else {
    warning("Defined dgp_sigmasq is NA, zero, or non-positive; cannot calculate 'size' for nbinom.")
    NA
  }
  if(is.na(size_param_dgp)) stop("size_param_dgp for nbinom is NA. Halting simulation.")
  cat("  Calculated 'size' parameter for nbinom simulation (1/sigmasq):", size_param_dgp, "\n")
  
  cat(paste("\nStep 2: Simulating", num_simulations, "INGARCH(1,3) realizations of length", sim_length, "...\n"))
  sims <- vector("list", num_simulations)
  for(i in 1:num_simulations) {
    if(i == 1 || i == num_simulations || (i %% progress_print_frequency == 0 && num_simulations > 10 && progress_print_frequency > 0)) {
      cat("  Generating simulation", i, "of", num_simulations, "...\n")
    }
    sim_result <- ingarch.sim(n = sim_length,
                              param = ingarch_params_dgp,
                              model = ingarch_model_dgp,
                              link = "log", distr = "nbinom", # Assuming log link and nbinom for simulation
                              size = size_param_dgp, n_start = 100)
    if (!is.null(sim_result) && is.list(sim_result) && "ts" %in% names(sim_result)) {
      sims[[i]] <- as.numeric(sim_result$ts)
    } else {
      warning(paste("Cannot extract numeric time series from ingarch.sim output in simulation", i, ". Skipping this simulation."))
      sims[[i]] <- NULL
    }
  }
  sims <- sims[!sapply(sims, is.null)] # Remove any NULL simulations
  actual_num_simulations <- length(sims)
  if (actual_num_simulations == 0) {
    stop("No valid simulations were generated. Halting.")
  }
  if(actual_num_simulations < num_simulations){
    cat("Warning: Only", actual_num_simulations, "valid simulations were generated out of", num_simulations, "requested.\n")
  }
  cat("Finished generating", actual_num_simulations, "valid simulations.\n")
  
  cat("\nStep 3: Running model selection SEQUENTIALLY...\n")
  
  # Parameters for auto.ingarch fitting (can be kept as they were or adjusted)
  task_max.p <- 4
  task_max.q <- 4
  task_max.order_stepwise <- 5  # Max sum p+q for initial models in stepwise
  task_max.order_grid <- 4     # Max sum p+q for non-stepwise (search.ingarch)
  task_distribution <- "nbinom" # Distribution for fitting
  task_link <- "log"            # Link function for fitting
  task_ic <- "aic"              # Information criterion for auto.ingarch selection
  
  define_results_template <- function(max_p_cols, max_q_cols) {
    na_betas_list <- list(); if(max_p_cols > 0) { na_betas_list <- as.list(rep(NA_real_, max_p_cols)); names(na_betas_list) <- paste0("beta", 1:max_p_cols) }
    na_alphas_list <- list(); if(max_q_cols > 0) { na_alphas_list <- as.list(rep(NA_real_, max_q_cols)); names(na_alphas_list) <- paste0("alpha", 1:max_q_cols) }
    template <- c(list(p_order=NA_integer_, q_order=NA_integer_, time=NA_real_, n_models_tested=NA_integer_, aic=NA_real_, bic=NA_real_, intercept=NA_real_, sigmasq=NA_real_), na_betas_list, na_alphas_list, list(status="Not Run", error_message=""))
    template[sapply(template, is.null)] <- NA
    for (name in c("p_order", "q_order", "n_models_tested")) template[[name]] <- NA_integer_
    for (name in c("time", "aic", "bic", "intercept", "sigmasq")) template[[name]] <- NA_real_
    if(max_p_cols > 0) for(k in 1:max_p_cols) template[[paste0("beta",k)]] <- NA_real_
    if(max_q_cols > 0) for(k in 1:max_q_cols) template[[paste0("alpha",k)]] <- NA_real_
    template$status <- "Not Run"; template$error_message <- ""
    return(template)
  }
  results_template <- define_results_template(task_max.p, task_max.q)
  
  populate_results_from_fit <- function(fit_object, template_data,
                                        max_p_cols_in_template, max_q_cols_in_template) {
    populated_data <- template_data
    if (!is.null(fit_object) && is.null(fit_object$error) && inherits(fit_object, "tsglm")) {
      populated_data$status <- "Success"; populated_data$error_message <- ""
      model_summary_obj <- tryCatch(summary(fit_object), error = function(e) NULL)
      coeffs_vector <- NULL
      if (!is.null(model_summary_obj) && !is.null(model_summary_obj$coefficients) &&
          is.matrix(model_summary_obj$coefficients) && "Estimate" %in% colnames(model_summary_obj$coefficients)) {
        coeffs_vector <- model_summary_obj$coefficients[, "Estimate", drop=TRUE]
        names(coeffs_vector) <- rownames(model_summary_obj$coefficients)
      } else { coeffs_vector <- tryCatch(stats::coef(fit_object), error = function(e) NULL) }
      
      current_p_order <- 0L; current_q_order <- 0L
      if (!is.null(fit_object$model)) { # $model should contain $past_obs and $past_mean as vectors of lags
        p_lags <- fit_object$model$past_obs
        q_lags <- fit_object$model$past_mean
        current_p_order <- if(is.null(p_lags) || length(p_lags) == 0) 0L else as.integer(max(p_lags))[1]
        current_q_order <- if(is.null(q_lags) || length(q_lags) == 0) 0L else as.integer(max(q_lags))[1]
      }
      populated_data$p_order <- current_p_order
      populated_data$q_order <- current_q_order
      populated_data$n_models_tested <- if(!is.null(fit_object$n_total_models)) as.integer(fit_object$n_total_models)[1] else NA_integer_
      populated_data$aic <- as.numeric(tryCatch(stats::AIC(fit_object), error = function(e) NA_real_))[1]
      populated_data$bic <- as.numeric(tryCatch(stats::BIC(fit_object), error = function(e) NA_real_))[1]
      
      if (!is.null(coeffs_vector)) {
        if ("(Intercept)" %in% names(coeffs_vector)) { populated_data$intercept <- as.numeric(coeffs_vector["(Intercept)"])[1]
        } else if ("intercept" %in% names(coeffs_vector)) { populated_data$intercept <- as.numeric(coeffs_vector["intercept"])[1] }
        
        # Match coefficient names: tsglm uses beta_LAG, alpha_LAG or past_obsLAG, past_meanLAG
        # The template uses betaK, alphaK.
        # We need to map from model coefficient names to template names.
        # This part assumes fit_object$model$past_obs = c(1,2,...) if p > 0
        
        # Populate Betas (past_obs)
        if (max_p_cols_in_template > 0 && current_p_order > 0 && !is.null(fit_object$model$past_obs)) {
          for (lag_val in fit_object$model$past_obs) { # Iterate through actual lags used
            if (lag_val > max_p_cols_in_template) next # Only store up to template capacity
            template_beta_col <- paste0("beta", lag_val)
            # Standard names from summary.tsglm are like 'past_obs1', 'past_obs2'
            # Or if your myingarch/auto.ingarch renames them to beta_1, beta_2
            model_coef_name_past_obs <- paste0("past_obs", lag_val)
            model_coef_name_beta <- paste0("beta_", lag_val) # Common in some ARIMA outputs
            
            if (model_coef_name_past_obs %in% names(coeffs_vector)) {
              populated_data[[template_beta_col]] <- as.numeric(coeffs_vector[model_coef_name_past_obs])[1]
            } else if (model_coef_name_beta %in% names(coeffs_vector)) {
              populated_data[[template_beta_col]] <- as.numeric(coeffs_vector[model_coef_name_beta])[1]
            }
          }
        }
        # Populate Alphas (past_mean)
        if (max_q_cols_in_template > 0 && current_q_order > 0 && !is.null(fit_object$model$past_mean)) {
          for (lag_val in fit_object$model$past_mean) { # Iterate through actual lags used
            if (lag_val > max_q_cols_in_template) next
            template_alpha_col <- paste0("alpha", lag_val)
            model_coef_name_past_mean <- paste0("past_mean", lag_val)
            model_coef_name_alpha <- paste0("alpha_", lag_val)
            
            if (model_coef_name_past_mean %in% names(coeffs_vector)) {
              populated_data[[template_alpha_col]] <- as.numeric(coeffs_vector[model_coef_name_past_mean])[1]
            } else if (model_coef_name_alpha %in% names(coeffs_vector)) {
              populated_data[[template_alpha_col]] <- as.numeric(coeffs_vector[model_coef_name_alpha])[1]
            }
          }
        }
      } # end if !is.null(coeffs_vector)
      
      if (!is.null(fit_object$distr) && fit_object$distr == "nbinom") {
        if(!is.null(fit_object$sigmasq)) { # If myingarch/auto.ingarch added it
          populated_data$sigmasq <- as.numeric(fit_object$sigmasq)[1]
        } else if (!is.null(model_summary_obj) && !is.null(model_summary_obj$dispersion) && "estimate" %in% names(model_summary_obj$dispersion)){
          alpha_disp <- as.numeric(model_summary_obj$dispersion["estimate"])[1] # This is the overdispersion alpha for var=mu+alpha*mu^2
          if(!is.na(alpha_disp)) populated_data$sigmasq <- alpha_disp # Storing alpha directly as sigmasq based on original script's naming
        } else if (!is.null(fit_object$distrcoefs) && "size" %in% names(fit_object$distrcoefs) && is.numeric(fit_object$distrcoefs$size) && fit_object$distrcoefs$size > 0) {
          # Fallback if summary doesn't have dispersion but fit object has size
          populated_data$sigmasq <- 1 / as.numeric(fit_object$distrcoefs$size)[1] # sigmasq as 1/size
        }
      } else { populated_data$sigmasq <- NA_real_ }
    } else { # Fit object error or not tsglm
      populated_data$status <- if(!is.null(fit_object$status_message)) fit_object$status_message else "Fit Error/Null"
      populated_data$error_message <- if(!is.null(fit_object$message)) as.character(fit_object$message)[1] else "Fit object error or NULL"
    }
    # Ensure all template fields exist, even if NULL from fit_object processing
    for(name_iter in names(template_data)){
      if(!name_iter %in% names(populated_data) || is.null(populated_data[[name_iter]])){
        populated_data[[name_iter]] <- template_data[[name_iter]] # Use template's NA default
      }
      if(name_iter %in% c("error_message", "status") && (is.na(populated_data[[name_iter]]) || is.null(populated_data[[name_iter]]))) {
        populated_data[[name_iter]] <- if(name_iter == "error_message") "" else "Undefined"
      }
    }
    return(populated_data)
  }
  
  # --- Sequential Loop for Model Fitting ---
  sequential_results_list <- vector("list", actual_num_simulations)
  
  for (i in 1:actual_num_simulations) {
    cat("  Processing simulation", i, "of", actual_num_simulations, "sequentially...\n")
    sim_data_ts <- ts(sims[[i]]) # Ensure it's a ts object for auto.ingarch
    
    # Define template for this iteration's results
    # (worker_results_template was used in parallel, can define results_template once outside loop too)
    # results_template is already defined outside
    
    # --- Stepwise Search ---
    current_stepwise_template <- results_template # Start with a fresh template
    stepwise_fit_object <- NULL
    cat("    Fitting stepwise model for simulation", i, "...\n")
    stepwise_time_taken <- system.time({
      stepwise_fit_object <- tryCatch({
        auto.ingarch(y = sim_data_ts, max.p = task_max.p, max.q = task_max.q,
                     max.order = task_max.order_stepwise,
                     distribution = task_distribution, link = task_link, ic = task_ic,
                     stepwise = TRUE, trace = FALSE, show_warnings = FALSE, parallel = FALSE)
      }, error = function(e) {
        cat("      ERROR in stepwise auto.ingarch for sim", i, ":", conditionMessage(e), "\n")
        return(list(error = TRUE, message = conditionMessage(e), status_message = "Stepwise Fit Error"))
      })
    })["elapsed"]
    current_stepwise_template$time <- as.numeric(stepwise_time_taken)
    iter_stepwise_result <- populate_results_from_fit(stepwise_fit_object, current_stepwise_template, task_max.p, task_max.q)
    # populate_results_from_fit might overwrite time, so re-assign if necessary or ensure it doesn't
    iter_stepwise_result$time <- current_stepwise_template$time
    
    
    # --- Grid Search (Non-stepwise) ---
    current_grid_template <- results_template # Start with a fresh template
    grid_fit_object <- NULL
    cat("    Fitting grid search model for simulation", i, "...\n")
    grid_time_taken <- system.time({
      grid_fit_object <- tryCatch({
        auto.ingarch(y = sim_data_ts, max.p = task_max.p, max.q = task_max.q,
                     max.order = task_max.order_grid, # This max.order is for search.ingarch
                     distribution = task_distribution, link = task_link, ic = task_ic,
                     stepwise = FALSE, trace = FALSE, show_warnings = FALSE, parallel = FALSE)
      }, error = function(e) {
        cat("      ERROR in grid auto.ingarch for sim", i, ":", conditionMessage(e), "\n")
        return(list(error = TRUE, message = conditionMessage(e), status_message = "Grid Fit Error"))
      })
    })["elapsed"]
    current_grid_template$time <- as.numeric(grid_time_taken)
    iter_grid_result <- populate_results_from_fit(grid_fit_object, current_grid_template, task_max.p, task_max.q)
    iter_grid_result$time <- current_grid_template$time # Re-assign
    
    # This logic for n_models_tested for grid search was in the parallel version, might still be relevant
    # if auto.ingarch/search.ingarch doesn't always populate n_total_models.
    # However, search.ingarch itself has n_models_evaluated.
    # Let's rely on n_total_models from the fit object if populated by auto.ingarch.
    # if(iter_grid_result$status == "Success" && (is.na(iter_grid_result$n_models_tested) || iter_grid_result$n_models_tested == 0)){
    #   # This is a rough count for a full grid if max.order wasn't restrictive
    #   # iter_grid_result$n_models_tested <- sum(outer(0:task_max.p, 0:task_max.q, "+") <= task_max.order_grid)
    # }
    
    sequential_results_list[[i]] <- list(stepwise = iter_stepwise_result, grid = iter_grid_result)
  } # End of sequential for loop
  
  cat("Sequential model fitting finished.\n")
  
  # --- Process and Save Results (similar to parallel version) ---
  results_list_stepwise <- vector("list", actual_num_simulations)
  results_list_grid <- vector("list", actual_num_simulations)
  
  for(i in 1:actual_num_simulations) {
    # No need to check for 'error' class as tryCatch now returns a list with $error
    current_iter_result <- sequential_results_list[[i]]
    if (is.list(current_iter_result) && !is.null(current_iter_result$stepwise) && !is.null(current_iter_result$grid)) {
      results_list_stepwise[[i]] <- current_iter_result$stepwise
      results_list_grid[[i]] <- current_iter_result$grid
    } else {
      cat("Unexpected result structure from sequential iteration", i, "\n")
      error_template <- define_results_template(task_max.p, task_max.q)
      error_template$status <- "Iteration Result Invalid Structure"
      error_template$error_message <- "Iteration did not return the expected list structure."
      results_list_stepwise[[i]] <- error_template
      results_list_grid[[i]] <- error_template
    }
  }
  
  cat("\nStep 4: Summarizing and saving results...\n")
  target_column_names <- c("method", "sim_id", names(results_template))
  
  convert_to_df_with_bind_rows <- function(list_of_result_lists_input, method_name_str, template_for_cols, target_cols_config) {
    # This function seems robust enough, keeping as is from original script.
    # Ensure it correctly handles the structure of list_of_result_lists_input.
    list_of_1row_dfs <- lapply(1:length(list_of_result_lists_input), function(j) {
      current_row_data_list <- list_of_result_lists_input[[j]]
      final_df_row_list <- vector("list", length(target_cols_config))
      names(final_df_row_list) <- target_cols_config
      final_df_row_list[["method"]] <- method_name_str
      final_df_row_list[["sim_id"]] <- as.integer(j)
      for(col_name_template in names(template_for_cols)){
        if(col_name_template %in% names(current_row_data_list) && !is.null(current_row_data_list[[col_name_template]])){
          val_to_assign <- current_row_data_list[[col_name_template]]
          if(length(val_to_assign) > 1 && !is.list(val_to_assign)) val_to_assign <- val_to_assign[1] # Avoid issues with vectors in cells
          final_df_row_list[[col_name_template]] <- val_to_assign
        } else { final_df_row_list[[col_name_template]] <- template_for_cols[[col_name_template]] }
        if(col_name_template == "error_message" && (is.na(final_df_row_list[[col_name_template]]) || is.null(final_df_row_list[[col_name_template]]))) {
          final_df_row_list[[col_name_template]] <- "" }
        if(col_name_template == "status" && (is.na(final_df_row_list[[col_name_template]]) || is.null(final_df_row_list[[col_name_template]]))) {
          final_df_row_list[[col_name_template]] <- "Undefined" }
      }
      for(target_col in target_cols_config){ # Ensure all target columns are present
        if(!target_col %in% names(final_df_row_list) || is.null(final_df_row_list[[target_col]])){
          if(target_col %in% names(template_for_cols)) final_df_row_list[[target_col]] <- template_for_cols[[target_col]]
          else final_df_row_list[[target_col]] <- NA }
      }
      # Ensure no list columns which bind_rows handles poorly unless all are lists
      is_list_col <- sapply(final_df_row_list, is.list)
      if(any(is_list_col)){
        # This indicates an issue, as results_template mostly defines atomic types.
        # For simplicity, convert problematic list columns to character.
        warning(paste("List column(s) found for sim_id",j,"method",method_name_str, ":", names(final_df_row_list)[is_list_col]))
        final_df_row_list[is_list_col] <- lapply(final_df_row_list[is_list_col], function(x) paste(unlist(x), collapse=";"))
      }
      
      return(as.data.frame(final_df_row_list, stringsAsFactors = FALSE))
    })
    list_of_1row_dfs_filtered <- list_of_1row_dfs[!sapply(list_of_1row_dfs, is.null)]
    if (length(list_of_1row_dfs_filtered) > 0) {
      # Align columns before bind_rows for robustness
      all_names <- unique(unlist(lapply(list_of_1row_dfs_filtered, names)))
      list_of_1row_dfs_aligned <- lapply(list_of_1row_dfs_filtered, function(df_row) {
        missing_cols <- setdiff(all_names, names(df_row))
        for (col in missing_cols) {
          if(col %in% names(template_for_cols)) df_row[[col]] <- template_for_cols[[col]]
          else df_row[[col]] <- NA # Default NA for any unexpected new columns
        }
        # Ensure correct order of columns
        return(df_row[, all_names, drop = FALSE])
      })
      # Check for non-atomic columns again before bind_rows
      # It might be better that populate_results_from_fit guarantees atomic results for each cell
      # For now, this should mostly work.
      return(suppressWarnings(dplyr::bind_rows(list_of_1row_dfs_aligned))) # suppress warnings from type coercions if any
    } else { # Return empty df with correct structure
      empty_df_structure <- template_for_cols; empty_df_structure$method <- character(0); empty_df_structure$sim_id <- integer(0)
      empty_df_final_structure <- empty_df_structure[target_cols_config]
      return(as.data.frame(empty_df_final_structure))
    }
  }
  
  stepwise_results_df <- convert_to_df_with_bind_rows(results_list_stepwise, "stepwise", results_template, target_column_names)
  grid_results_df     <- convert_to_df_with_bind_rows(results_list_grid, "grid_search", results_template, target_column_names)
  
  results_df <- dplyr::bind_rows(stepwise_results_df, grid_results_df)
  
  output_dir_name <- "ingarch_dgp13_sequential_results" # Changed output directory name
  if (!dir.exists(output_dir_name)) {
    dir.create(output_dir_name, recursive = TRUE)
  }
  
  # Define file paths using the new directory
  results_csv_path <- file.path(output_dir_name, "ingarch_dgp13_results.csv")
  sims_rds_path <- file.path(output_dir_name, "ingarch_dgp13_simulations.rds")
  summary_stats_csv_path <- file.path(output_dir_name, "ingarch_dgp13_summary_stats.csv")
  order_freq_csv_path <- file.path(output_dir_name, "ingarch_dgp13_order_frequencies.csv")
  
  if (nrow(results_df) > 0) {
    if("status" %in% names(results_df) && is.factor(results_df$status)) results_df$status <- as.character(results_df$status)
    if("error_message" %in% names(results_df) && is.factor(results_df$error_message)) results_df$error_message <- as.character(results_df$error_message)
    
    write.csv(results_df, results_csv_path, row.names = FALSE, na = "")
    cat("  Detailed results saved to:", results_csv_path, "\n")
  } else { cat("  No results to save to CSV.\n") }
  
  saveRDS(sims, sims_rds_path)
  cat("  Simulated datasets saved to:", sims_rds_path, "\n")
  
  # --- Summarize Results (logic largely unchanged) ---
  if (nrow(results_df) > 0) {
    numeric_cols_for_summary <- c("p_order", "q_order", "time", "n_models_tested", "aic", "bic", "intercept", "sigmasq")
    # Ensure beta/alpha columns for summary match what results_template can hold (up to task_max.p/q)
    if(task_max.p > 0) numeric_cols_for_summary <- c(numeric_cols_for_summary, paste0("beta", 1:min(task_max.p, ncol(results_df[startsWith(names(results_df), "beta")]))))
    if(task_max.q > 0) numeric_cols_for_summary <- c(numeric_cols_for_summary, paste0("alpha", 1:min(task_max.q, ncol(results_df[startsWith(names(results_df), "alpha")]))))
    numeric_cols_for_summary <- intersect(numeric_cols_for_summary, names(results_df)) # Only existing columns
    
    for(col_s in numeric_cols_for_summary){
      if(col_s %in% names(results_df)){ # Should always be true due to intersect
        results_df[[col_s]] <- suppressWarnings(as.numeric(as.character(results_df[[col_s]])))
      }
    }
    results_df$status <- as.character(results_df$status) # Ensure character
    
    summary_stats <- results_df %>%
      filter(!is.na(method)) %>%
      group_by(method) %>%
      summarize(
        total_sims = n(),
        successful_fits = sum(status == "Success", na.rm = TRUE),
        mean_p_order = mean(p_order[status == "Success"], na.rm = TRUE),
        median_p_order = median(p_order[status == "Success"], na.rm = TRUE),
        sd_p_order = sd(p_order[status == "Success"], na.rm = TRUE),
        mean_q_order = mean(q_order[status == "Success"], na.rm = TRUE),
        median_q_order = median(q_order[status == "Success"], na.rm = TRUE),
        sd_q_order = sd(q_order[status == "Success"], na.rm = TRUE),
        mean_time_secs = mean(time, na.rm = TRUE),
        median_time_secs = median(time, na.rm = TRUE),
        sd_time_secs = sd(time, na.rm = TRUE),
        mean_fitted_models = mean(n_models_tested[status == "Success"], na.rm = TRUE),
        median_fitted_models = median(n_models_tested[status == "Success"], na.rm = TRUE),
        sd_fitted_models = sd(n_models_tested[status == "Success"], na.rm = TRUE),
        failure_count = sum(status != "Success" & !(is.na(status) | status == ""), na.rm = TRUE),
        .groups = 'drop'
      )
    
    order_freq <- results_df %>%
      filter(status == "Success", !is.na(p_order), !is.na(q_order)) %>%
      group_by(method, p_order, q_order) %>%
      summarize(count = n(), .groups = 'drop') %>%
      group_by(method) %>%
      mutate(freq = count / sum(count) * 100) %>%
      arrange(method, desc(freq))
    
    write.csv(summary_stats, summary_stats_csv_path, row.names = FALSE, na = "")
    write.csv(order_freq, order_freq_csv_path, row.names = FALSE, na = "")
    cat("  Summary statistics and order frequencies saved to CSV files in '", output_dir_name, "'.\n")
    
    cat("\n===== Summary Statistics (Sequential, DGP 1,3) =====\n"); print(summary_stats)
    cat("\n===== Most Frequent (p_order,q_order) Orders (Top 5 per method - Sequential, DGP 1,3) =====\n")
    if(nrow(order_freq) > 0) {
      top_orders <- order_freq %>% group_by(method) %>% slice_head(n = 5)
      print(as.data.frame(top_orders)) # Print as data frame for better console view
    } else { cat("No successful fits to determine frequent orders.\n") }
    cat("\n===== Status Summary (Sequential, DGP 1,3) =====\n")
    print(as.data.frame(results_df %>% group_by(method, status) %>% summarise(count = n(), .groups = 'drop') %>% arrange(method, desc(count))))
    
  } else { cat("  No results available to summarize.\n") }
  
  return(results_df)
}

# Main function to call the sequential simulation study
main_sequential_1_3_no_cov <- function() {
  # Clean up any old flags if they exist (from original script)
  if(exists(".coef_names_printed_flag_v4", envir = .GlobalEnv)){
    rm(".coef_names_printed_flag_v4", envir = .GlobalEnv)
  }
  
  overall_start_time <- Sys.time()
  results_data_frame <- tryCatch(run_simulation_study_sequential_1_3(),
                                 error = function(e) {
                                   cat("\nERROR during SEQUENTIAL simulation study run (DGP 1,3, no covariates):\n")
                                   print(e)
                                   # No cluster to stop in sequential version
                                   return(NULL)
                                 })
  overall_end_time <- Sys.time()
  overall_duration <- overall_end_time - overall_start_time
  
  if(!is.null(results_data_frame)){
    cat("\nSequential simulation study (DGP 1,3, no covariates) completed!\n")
  } else {
    cat("\nSequential simulation study (DGP 1,3, no covariates) encountered an error and did not complete successfully.\n")
  }
  cat("Total duration of the study:", format(overall_duration), "\n")
}

# Execute the main function
main_sequential_1_3_no_cov()