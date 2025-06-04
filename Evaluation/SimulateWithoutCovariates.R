# INGARCH Model Simulation Study - Without Covariates

library(tscount)
library(dplyr)
library(readxl)
library(doParallel)
library(foreach)


source("./ACTS/auto.ingarch.R")
source("./ACTS/ingarch.sim.R")
source("./ACTS/newmodel.R")
source("./ACTS/ingarch.string.R")
source("./ACTS/search.ingarch.R")
source("./ACTS/myingarch.R")

extract_model_params <- function(data_path, model_row = 1) {
  excel_data <- readxl::read_excel(data_path)
  model_data <- excel_data[model_row, ]
  
  p <- as.numeric(model_data$p)
  q <- as.numeric(model_data$q)
  cat("\nOrders: (", p, ",", q, ")")
  sigmasq <- as.numeric(model_data$sigmasq) # This is used as 'alpha' in var=mu+alpha*mu^2 for NB
  cat("\nSigmasq: ", sigmasq)
  intercept <- as.numeric(model_data$Intercept)
  cat("\nIntercept: ", intercept)
  
  beta_cols <- grep("^beta", colnames(excel_data), ignore.case = TRUE)
  betas <- if(length(beta_cols) > 0 && ncol(model_data[, beta_cols, drop=FALSE]) > 0) as.numeric(unlist(model_data[, beta_cols])) else numeric(0)
  cat("\nBetas: ")
  print(betas)
  alpha_cols <- grep("^alpha", colnames(excel_data), ignore.case = TRUE)
  alphas <- if(length(alpha_cols) > 0 && ncol(model_data[, alpha_cols, drop=FALSE]) > 0) as.numeric(unlist(model_data[, alpha_cols])) else numeric(0)
  cat("\nAlphas: ")
  print(alphas)
  
  return(list(
    p = p, q = q, sigmasq = sigmasq, intercept = intercept,
    betas = betas, alphas = alphas
  ))
}

run_simulation_study_no_covariates_parallel <- function() {
  set.seed(12345)
  
  num_simulations <- 1000
  sim_length <- 1000
  progress_print_frequency <- max(1, floor(num_simulations / 10))
  
  cat("Starting INGARCH simulation study (without covariates) - PARALLEL EXECUTION.\n")
  cat("Configuration: num_simulations =", num_simulations, ", sim_length =", sim_length, "\n")
  
  cat("Step 1: Extracting model parameters from Excel file...\n")
  params <- extract_model_params("./data/modelosAveiro.xlsx", model_row = 2) # Ensure this path is correct
  cat("  True p_order=", params$p, ", True q_order=", params$q, ", True sigmasq (alpha for NB variance form mu+alpha*mu^2)=", params$sigmasq, "\n") # MODIFIED comment
  
  ingarch_params <- list(
    intercept = params$intercept,
    past_obs = if(params$p > 0 && length(params$betas) >= params$p) params$betas[1:params$p] else NULL,
    past_mean = if(params$q > 0 && length(params$alphas) >= params$q) params$alphas[1:params$q] else NULL
    # No 'xreg' element as this simulation is without covariates
  )
  ingarch_model <- list(
    past_obs = if(params$p > 0) 1:params$p else NULL,
    past_mean = if(params$q > 0) 1:params$q else NULL,
    external = NULL
  )
  
  distrcoefs_param <- if(!is.na(params$sigmasq) && params$sigmasq > 0) {
    1 / params$sigmasq 
  } else {
    warning("sigmasq from Excel is NA, zero, or non-positive; cannot calculate 'distrcoefs' (size) for nbinom.")
    NA
  }
  if(is.na(distrcoefs_param)) stop("Calculated 'distrcoefs' (size) for nbinom is NA. Halting simulation.")
  cat("  Calculated 'distrcoefs' (size = 1/sigmasq) for nbinom simulation:", distrcoefs_param, "\n") # MODIFIED message
  
  cat(paste("\nStep 2: Simulating", num_simulations, "INGARCH realizations of length", sim_length, "...\n"))
  sims <- vector("list", num_simulations)
  for(i in 1:num_simulations) {
    if(i == 1 || i == num_simulations || (i %% progress_print_frequency == 0 && num_simulations > 10 && progress_print_frequency > 0 )) { # Added check for progress_print_frequency > 0
      cat("  Generating simulation", i, "of", num_simulations, "...\n")
    }

    sim_result <- ingarch.sim(n = sim_length, 
                              param = ingarch_params, 
                              model = ingarch_model, 
                              link = "log", 
                              distr = "nbinom", 
                              distrcoefs = distrcoefs_param, # Use the calculated distrcoefs_param
                              n_start = 100)
    if (!is.null(sim_result) && is.list(sim_result) && "ts" %in% names(sim_result)) {
      sims[[i]] <- as.numeric(sim_result$ts)
    } else {
      warning(paste("Cannot extract numeric time series from ingarch.sim output in simulation", i, ". Skipping this simulation."))
      sims[[i]] <- NULL
    }
  }
  sims <- sims[!sapply(sims, is.null)]
  actual_num_simulations <- length(sims)
  if (actual_num_simulations == 0) {
    stop("No valid simulations were generated. Halting.")
  }
  if(actual_num_simulations < num_simulations){ # Added this warning from your other script
    cat("Warning: Only", actual_num_simulations, "valid simulations were generated out of", num_simulations, "requested.\n")
  }
  cat("Finished generating", actual_num_simulations, "valid simulations.\n")
  
  cat("\nStep 3: Running model selection in PARALLEL...\n")
  
  task_max.p <- 7
  task_max.q <- 7
  task_max.order_stepwise <- 5
  task_max.order_grid <- 14 # Max order for grid search (non-stepwise)
  task_distribution <- "nbinom"
  task_link <- "log"
  task_ic <- "aic"
  
  define_results_template <- function(max_p_cols, max_q_cols) {
    na_betas_list <- list(); if(max_p_cols > 0) { na_betas_list <- as.list(rep(NA_real_, max_p_cols)); names(na_betas_list) <- paste0("beta", 1:max_p_cols) }
    na_alphas_list <- list(); if(max_q_cols > 0) { na_alphas_list <- as.list(rep(NA_real_, max_q_cols)); names(na_alphas_list) <- paste0("alpha", 1:max_q_cols) }
    template <- c(list(p_order=NA_integer_, q_order=NA_integer_, time=NA_real_, n_models_tested=NA_integer_, aic=NA_real_, bic=NA_real_, intercept=NA_real_, sigmasq=NA_real_), na_betas_list, na_alphas_list, list(status="Not Run", error_message=""))
    template[sapply(template, is.null)] <- NA # Should not be needed if list elements are NA_ type
    for (name in c("p_order", "q_order", "n_models_tested")) template[[name]] <- NA_integer_
    for (name in c("time", "aic", "bic", "intercept", "sigmasq")) template[[name]] <- NA_real_
    if(max_p_cols > 0) for(k in 1:max_p_cols) template[[paste0("beta",k)]] <- NA_real_
    if(max_q_cols > 0) for(k in 1:max_q_cols) template[[paste0("alpha",k)]] <- NA_real_
    template$status <- "Not Run"
    template$error_message <- ""
    return(template)
  }
  results_template <- define_results_template(task_max.p, task_max.q)
  
  populate_results_from_fit <- function(fit_object, template_data,
                                        max_p_cols_in_template, max_q_cols_in_template) {
    populated_data <- template_data # Start with the template
    if (!is.null(fit_object) && is.null(fit_object$error) && inherits(fit_object, "tsglm")) {
      populated_data$status <- "Success"
      populated_data$error_message <- ""
      model_summary_obj <- tryCatch(summary(fit_object), error = function(e) NULL)
      coeffs_vector <- NULL
      if (!is.null(model_summary_obj) && !is.null(model_summary_obj$coefficients) &&
          is.matrix(model_summary_obj$coefficients) && "Estimate" %in% colnames(model_summary_obj$coefficients)) {
        coeffs_vector <- model_summary_obj$coefficients[, "Estimate", drop=TRUE]
        names(coeffs_vector) <- rownames(model_summary_obj$coefficients)
      } else {
        coeffs_vector <- tryCatch(stats::coef(fit_object), error = function(e) NULL) # Fallback
      }
      
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
        
        # Populate Betas (past_obs)
        if (max_p_cols_in_template > 0 && current_p_order > 0 && !is.null(fit_object$model$past_obs)) {
          for (lag_val in fit_object$model$past_obs) { # Iterate through actual lags used
            if (lag_val > max_p_cols_in_template) next # Only store up to template capacity
            template_beta_col <- paste0("beta", lag_val)
            model_coef_name_past_obs <- paste0("past_obs", lag_val) # tsglm default name
            model_coef_name_beta <- paste0("beta_", lag_val) # Alternative if auto.ingarch renames
            
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
            model_coef_name_past_mean <- paste0("past_mean", lag_val) # tsglm default name
            model_coef_name_alpha <- paste0("alpha_", lag_val) # Alternative
            
            if (model_coef_name_past_mean %in% names(coeffs_vector)) {
              populated_data[[template_alpha_col]] <- as.numeric(coeffs_vector[model_coef_name_past_mean])[1]
            } else if (model_coef_name_alpha %in% names(coeffs_vector)) {
              populated_data[[template_alpha_col]] <- as.numeric(coeffs_vector[model_coef_name_alpha])[1]
            }
          }
        }
      } # end if !is.null(coeffs_vector)
      
      # sigmasq for NB (this is alpha in var=mu+alpha*mu^2)
      if (!is.null(fit_object$distr) && fit_object$distr == "nbinom") {
        if(!is.null(fit_object$sigmasq)) { # If myingarch/auto.ingarch added it directly (as alpha)
          populated_data$sigmasq <- as.numeric(fit_object$sigmasq)[1]
        } else if (!is.null(model_summary_obj) && !is.null(model_summary_obj$dispersion) && "estimate" %in% names(model_summary_obj$dispersion)){
          alpha_disp <- as.numeric(model_summary_obj$dispersion["estimate"])[1] # This is the overdispersion alpha
          if(!is.na(alpha_disp)) populated_data$sigmasq <- alpha_disp 
        } else if (!is.null(fit_object$distrcoefs)) { # Check if new ingarch.sim output has distrcoefs
          if (is.numeric(fit_object$distrcoefs) && fit_object$distrcoefs > 0) { # Assuming distrcoefs is the size
            populated_data$sigmasq <- 1 / as.numeric(fit_object$distrcoefs)[1] # sigmasq as 1/size
          } else if (is.list(fit_object$distrcoefs) && "size" %in% names(fit_object$distrcoefs) && is.numeric(fit_object$distrcoefs$size) && fit_object$distrcoefs$size > 0) {
            populated_data$sigmasq <- 1 / as.numeric(fit_object$distrcoefs$size)[1] # sigmasq as 1/size
          }
        }
      } else { # Not nbinom
        populated_data$sigmasq <- NA_real_
      }
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
  
  # Parallel setup
  num_cores_detected <- detectCores(logical = FALSE) # Physical cores
  cores_to_use <- max(1, num_cores_detected - 1)
  if(actual_num_simulations < cores_to_use) cores_to_use <- actual_num_simulations # Don't use more cores than tasks
  
  cl <- makeCluster(cores_to_use)
  registerDoParallel(cl)
  cat(paste("Registered parallel backend with", getDoParWorkers(), "cores for", actual_num_simulations, "simulations.\n"))
  
  # Export necessary functions and objects to workers
  custom_funcs_to_export <- c("auto.ingarch", "ingarch.sim", "newmodel",
                              "ingarch.string", "search.ingarch", "myingarch",
                              "populate_results_from_fit", "define_results_template")
  # Filter out functions that might not exist in the global environment for some reason
  custom_funcs_to_export <- custom_funcs_to_export[sapply(custom_funcs_to_export, exists, envir = .GlobalEnv, mode = "function")]
  
  parallel_results_list <- foreach(
    i = 1:actual_num_simulations,
    .packages = c("tscount", "stats", "dplyr"), # Ensure workers have necessary packages
    .export = c(custom_funcs_to_export, # Exporting filtered custom functions
                "sims", "task_max.p", "task_max.q", "task_max.order_stepwise", "task_max.order_grid",
                "task_distribution", "task_link", "task_ic"), # Exporting task parameters
    .errorhandling = 'pass' # Pass errors as results
  ) %dopar% {
    # This code runs on each worker
    worker_results_template <- define_results_template(task_max.p, task_max.q) # Define template inside worker
    sim_data_ts <- ts(sims[[i]]) # Ensure it's a ts object
    
    # --- Stepwise Search on Worker ---
    current_stepwise_template <- worker_results_template # Fresh template
    stepwise_fit_object <- NULL
    stepwise_time_taken <- system.time({
      stepwise_fit_object <- tryCatch({
        auto.ingarch(y = sim_data_ts, max.p = task_max.p, max.q = task_max.q, max.order = task_max.order_stepwise,
                     distribution = task_distribution, link = task_link, ic = task_ic,
                     stepwise = TRUE, trace = FALSE, show_warnings = FALSE, parallel = FALSE) # parallel=FALSE inside worker
      }, error = function(e) {
        # Return a list indicating error, message, and status
        return(list(error = TRUE, message = conditionMessage(e), status_message = "Stepwise Fit Error"))
      })
    })["elapsed"]
    current_stepwise_template$time <- as.numeric(stepwise_time_taken)
    iter_stepwise_result <- populate_results_from_fit(stepwise_fit_object, current_stepwise_template, task_max.p, task_max.q)
    iter_stepwise_result$time <- current_stepwise_template$time # Ensure time is correctly set
    
    # --- Grid Search (Non-stepwise) on Worker ---
    current_grid_template <- worker_results_template # Fresh template
    grid_fit_object <- NULL
    grid_time_taken <- system.time({
      grid_fit_object <- tryCatch({
        auto.ingarch(y = sim_data_ts, max.p = task_max.p, max.q = task_max.q, max.order = task_max.order_grid,
                     distribution = task_distribution, link = task_link, ic = task_ic,
                     stepwise = FALSE, trace = FALSE, show_warnings = FALSE, parallel = FALSE) # parallel=FALSE inside worker
      }, error = function(e) {
        return(list(error = TRUE, message = conditionMessage(e), status_message = "Grid Fit Error"))
      })
    })["elapsed"]
    current_grid_template$time <- as.numeric(grid_time_taken)
    iter_grid_result <- populate_results_from_fit(grid_fit_object, current_grid_template, task_max.p, task_max.q)
    iter_grid_result$time <- current_grid_template$time # Ensure time is correctly set
    
    # Estimate n_models_tested for grid search if not provided by auto.ingarch/search.ingarch
    if(iter_grid_result$status == "Success" && (is.na(iter_grid_result$n_models_tested) || iter_grid_result$n_models_tested == 0)){
      # This is a simplified count; search.ingarch might be more precise if it returns this.
      # It counts models where p+q <= max.order_grid.
      count = 0
      for(p_count in 0:task_max.p){
        for(q_count in 0:task_max.q){
          if(p_count + q_count <= task_max.order_grid) {
            count = count + 1
          }
        }
      }
      iter_grid_result$n_models_tested <- count
    }
    
    return(list(stepwise = iter_stepwise_result, grid = iter_grid_result))
  } # End of foreach parallel loop
  
  stopCluster(cl)
  cat("Parallel processing finished.\n")
  
  # --- Process and Save Results ---
  results_list_stepwise <- vector("list", actual_num_simulations)
  results_list_grid <- vector("list", actual_num_simulations)
  
  for(i in 1:actual_num_simulations) {
    if (inherits(parallel_results_list[[i]], "error")) {
      cat("Error in parallel iteration", i, ":", conditionMessage(parallel_results_list[[i]]), "\n")
      error_template <- define_results_template(task_max.p, task_max.q) # Use the globally defined template
      error_template$status <- "Parallel Iteration Error"
      error_template$error_message <- conditionMessage(parallel_results_list[[i]])
      results_list_stepwise[[i]] <- error_template
      results_list_grid[[i]] <- error_template
    } else if (is.list(parallel_results_list[[i]]) && !is.null(parallel_results_list[[i]]$stepwise) && !is.null(parallel_results_list[[i]]$grid)) {
      results_list_stepwise[[i]] <- parallel_results_list[[i]]$stepwise
      results_list_grid[[i]] <- parallel_results_list[[i]]$grid
    } else { # Should not happen if foreach errorhandling='pass' returns list with error element
      cat("Unexpected result structure from parallel iteration", i, "\n")
      error_template <- define_results_template(task_max.p, task_max.q)
      error_template$status <- "Worker Result Invalid Structure"
      error_template$error_message <- "Worker did not return the expected list structure."
      results_list_stepwise[[i]] <- error_template
      results_list_grid[[i]] <- error_template
    }
  }
  
  cat("\nStep 4: Summarizing and saving results...\n")
  target_column_names <- c("method", "sim_id", names(results_template)) # results_template is defined outside loops
  
  convert_to_df_with_bind_rows <- function(list_of_result_lists_input, method_name_str, template_for_cols, target_cols_config) {
    list_of_1row_dfs <- lapply(1:length(list_of_result_lists_input), function(j) {
      current_row_data_list <- list_of_result_lists_input[[j]]
      final_df_row_list <- vector("list", length(target_cols_config)) # Initialize with NAs or template
      names(final_df_row_list) <- target_cols_config
      
      # Set defaults from template for all target columns first
      for(col_name_target in target_cols_config) {
        if (col_name_target %in% names(template_for_cols)) {
          final_df_row_list[[col_name_target]] <- template_for_cols[[col_name_target]]
        } else if (col_name_target %in% c("method", "sim_id")) {
          # Handled next
        } else {
          final_df_row_list[[col_name_target]] <- NA # Default for any other unexpected columns
        }
      }
      
      final_df_row_list[["method"]] <- method_name_str
      final_df_row_list[["sim_id"]] <- as.integer(j)
      
      # Populate with actual data
      for(col_name_data in names(current_row_data_list)){ # Iterate over names in current_row_data_list
        if(col_name_data %in% target_cols_config && !is.null(current_row_data_list[[col_name_data]])){ # Check if it's a target column
          val_to_assign <- current_row_data_list[[col_name_data]]
          # Ensure atomic values for data frame cells, take first element if it's a vector unexpectedly
          if(length(val_to_assign) > 1 && !is.list(val_to_assign) && !is.data.frame(val_to_assign)) {
            val_to_assign <- val_to_assign[1]
          }
          final_df_row_list[[col_name_data]] <- val_to_assign
        }
      }
      
      # Special handling for status and error_message to ensure they are not NA
      if(is.na(final_df_row_list[["error_message"]]) || is.null(final_df_row_list[["error_message"]])) {
        final_df_row_list[["error_message"]] <- ""
      }
      if(is.na(final_df_row_list[["status"]]) || is.null(final_df_row_list[["status"]])) {
        final_df_row_list[["status"]] <- "Undefined"
      }
      
      # Convert any remaining list columns to character to avoid bind_rows issues
      # (Should ideally not happen if template and populate_results are robust)
      is_list_col <- sapply(final_df_row_list, is.list)
      if(any(is_list_col)){
        warning(paste("List column(s) being converted to character for sim_id", j, "method", method_name_str, ":", paste(names(final_df_row_list)[is_list_col], collapse=", ")))
        final_df_row_list[is_list_col] <- lapply(final_df_row_list[is_list_col], function(x) paste(unlist(x), collapse=";"))
      }
      
      return(as.data.frame(final_df_row_list, stringsAsFactors = FALSE))
    })
    
    list_of_1row_dfs_filtered <- list_of_1row_dfs[!sapply(list_of_1row_dfs, is.null)] # Should not be needed if all iterations return a df
    if (length(list_of_1row_dfs_filtered) == 0) { # If all iterations failed at a very early stage
      empty_df_structure <- template_for_cols
      empty_df_structure$method <- character(0)
      empty_df_structure$sim_id <- integer(0)
      # Ensure all target_cols_config are present
      for(col_conf in target_cols_config){
        if(!col_conf %in% names(empty_df_structure)) empty_df_structure[[col_conf]] <- NA
      }
      empty_df_final_structure <- empty_df_structure[intersect(names(empty_df_structure), target_cols_config)]
      return(as.data.frame(empty_df_final_structure))
    }
    
    # Robust column alignment before bind_rows
    all_col_names_present <- unique(unlist(lapply(list_of_1row_dfs_filtered, names)))
    target_cols_final_order <- intersect(target_cols_config, all_col_names_present) # Use defined order, only existing
    
    list_of_1row_dfs_aligned <- lapply(list_of_1row_dfs_filtered, function(df_row) {
      # Add missing columns with NA from template or pure NA
      missing_cols <- setdiff(target_cols_final_order, names(df_row))
      for (m_col in missing_cols) {
        df_row[[m_col]] <- if(m_col %in% names(template_for_cols)) template_for_cols[[m_col]] else NA
      }
      # Select and reorder to target_cols_final_order
      return(df_row[, target_cols_final_order, drop = FALSE])
    })
    
    return(suppressWarnings(dplyr::bind_rows(list_of_1row_dfs_aligned))) # Suppress type coercion warnings
  }
  
  stepwise_results_df <- convert_to_df_with_bind_rows(results_list_stepwise, "stepwise", results_template, target_column_names)
  grid_results_df     <- convert_to_df_with_bind_rows(results_list_grid, "grid_search", results_template, target_column_names)
  
  results_df <- dplyr::bind_rows(stepwise_results_df, grid_results_df)
  
  output_dir_name <- "ingarch_no_covariates_results_parallel" # Changed output directory name
  if (!dir.exists(output_dir_name)) {
    dir.create(output_dir_name, recursive = TRUE)
  }
  
  # Define file paths using the new directory
  results_csv_path <- file.path(output_dir_name, "ingarch_no_covariates_results.csv")
  sims_rds_path <- file.path(output_dir_name, "ingarch_no_covariates_simulations.rds")
  summary_stats_csv_path <- file.path(output_dir_name, "ingarch_no_covariates_summary_stats.csv")
  order_freq_csv_path <- file.path(output_dir_name, "ingarch_no_covariates_order_frequencies.csv")
  
  if (nrow(results_df) > 0) {
    # Ensure status and error_message are character, not factor, before writing
    if("status" %in% names(results_df) && is.factor(results_df$status)) results_df$status <- as.character(results_df$status)
    if("error_message" %in% names(results_df) && is.factor(results_df$error_message)) results_df$error_message <- as.character(results_df$error_message)
    
    write.csv(results_df, results_csv_path, row.names = FALSE, na = "") # Use na = "" for blank NAs
    cat("  Detailed results saved to:", results_csv_path, "\n")
  } else {
    cat("  No results to save to CSV.\n")
  }
  
  saveRDS(sims, sims_rds_path) # Save the list of simulated series
  cat("  Simulated datasets saved to:", sims_rds_path, "\n")
  
  # --- Summarize Results (logic largely unchanged) ---
  if (nrow(results_df) > 0) {
    # Ensure numeric columns are indeed numeric for summary stats
    numeric_cols_for_summary <- c("p_order", "q_order", "time", "n_models_tested", "aic", "bic", "intercept", "sigmasq")
    # Add beta/alpha columns that actually exist in results_df (up to task_max.p/q which defined the template)
    if(task_max.p > 0) numeric_cols_for_summary <- c(numeric_cols_for_summary, paste0("beta", 1:task_max.p))
    if(task_max.q > 0) numeric_cols_for_summary <- c(numeric_cols_for_summary, paste0("alpha", 1:task_max.q))
    
    # Filter to only include columns that are actually in results_df
    numeric_cols_for_summary <- intersect(numeric_cols_for_summary, names(results_df))
    
    for(col_s in numeric_cols_for_summary){
      if(col_s %in% names(results_df)){ # This check is now redundant due to intersect
        results_df[[col_s]] <- suppressWarnings(as.numeric(as.character(results_df[[col_s]])))
      }
    }
    
    results_df$status <- as.character(results_df$status) # Ensure character for grouping/filtering
    
    summary_stats <- results_df %>%
      filter(!is.na(method)) %>% # Ensure method is not NA
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
        mean_fitted_models = mean(n_models_tested[status == "Success"], na.rm = TRUE), # Ensure this is correctly populated
        median_fitted_models = median(n_models_tested[status == "Success"], na.rm = TRUE),
        sd_fitted_models = sd(n_models_tested[status == "Success"], na.rm = TRUE),
        failure_count = sum(status != "Success" & !(is.na(status) | status == ""), na.rm = TRUE), # Count non-success, non-empty/NA statuses
        .groups = 'drop'
      )
    
    order_freq <- results_df %>%
      filter(status == "Success", !is.na(p_order), !is.na(q_order)) %>%
      group_by(method, p_order, q_order) %>%
      summarize(count = n(), .groups = 'drop') %>%
      group_by(method) %>% # Group again to calculate frequency within each method
      mutate(freq = count / sum(count) * 100) %>%
      arrange(method, desc(freq))
    
    write.csv(summary_stats, summary_stats_csv_path, row.names = FALSE, na = "")
    write.csv(order_freq, order_freq_csv_path, row.names = FALSE, na = "")
    cat("  Summary statistics and order frequencies saved to CSV files in '", output_dir_name, "'.\n")
    
    cat("\n===== Summary Statistics (Parallel, No Covariates) =====\n"); print(summary_stats) # MODIFIED Title
    cat("\n===== Most Frequent (p_order,q_order) Orders (Top 5 per method - Parallel, No Covariates) =====\n") # MODIFIED Title
    if(nrow(order_freq) > 0) {
      top_orders <- order_freq %>% group_by(method) %>% slice_head(n = 5)
      print(as.data.frame(top_orders)) # Print as data frame for better console view
    } else {
      cat("No successful fits to determine frequent orders.\n")
    }
    cat("\n===== Status Summary (Parallel, No Covariates) =====\n") # MODIFIED Title
    print(as.data.frame(results_df %>% group_by(method, status) %>% summarise(count = n(), .groups = 'drop') %>% arrange(method, desc(count))))
    
  } else {
    cat("  No results available to summarize.\n")
  }
  
  return(results_df)
}

# Main function to call the parallel simulation study
main_parallel_no_cov <- function() {
  # Clean up any old flags if they exist (from previous script versions)
  if(exists(".coef_names_printed_flag_v4", envir = .GlobalEnv)){
    rm(".coef_names_printed_flag_v4", envir = .GlobalEnv)
  }
  
  overall_start_time <- Sys.time()
  results_data_frame <- tryCatch(run_simulation_study_no_covariates_parallel(),
                                 error = function(e) {
                                   cat("\nERROR during PARALLEL simulation study run (no covariates):\n")
                                   print(e)
                                   # Attempt to stop cluster if it exists from the error scope.
                                   # 'cl' might not be in this scope; better to handle in run_simulation_study...
                                   # if(exists("cl") && inherits(cl, "cluster")) try(stopCluster(cl), silent=TRUE)
                                   return(NULL)
                                 })
  overall_end_time <- Sys.time()
  overall_duration <- overall_end_time - overall_start_time
  
  if(!is.null(results_data_frame)){
    cat("\nParallel simulation study (no covariates) completed!\n")
  } else {
    cat("\nParallel simulation study (no covariates) encountered an error and did not complete successfully.\n")
  }
  cat("Total duration of the study:", format(overall_duration), "\n")
}

# Execute the main function
main_parallel_no_cov()