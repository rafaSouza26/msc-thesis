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

run_simulation_study_no_covariates_parallel <- function() {
  set.seed(12345)
  
  num_simulations <- 10
  sim_length <- 10
  progress_print_frequency <- max(1, floor(num_simulations / 10))
  
  cat("Starting INGARCH simulation study (without covariates) - PARALLEL EXECUTION.\n")
  cat("Configuration: num_simulations =", num_simulations, ", sim_length =", sim_length, "\n")
  
  cat("Step 1: Extracting model parameters from Excel file...\n")
  params <- extract_model_params("./data/modelosAveiro.xlsx", model_row = 1)
  cat("  True p_order=", params$p, ", True q_order=", params$q, ", True sigmasq=", params$sigmasq, "\n")
  
  ingarch_params <- list(
    intercept = params$intercept,
    past_obs = if(params$p > 0 && length(params$betas) >= params$p) params$betas[1:params$p] else NULL,
    past_mean = if(params$q > 0 && length(params$alphas) >= params$q) params$alphas[1:params$q] else NULL
  )
  ingarch_model <- list(
    past_obs = if(params$p > 0) 1:params$p else NULL,
    past_mean = if(params$q > 0) 1:params$q else NULL,
    external = FALSE
  )
  size_param <- if(!is.na(params$sigmasq) && params$sigmasq > 0) 1/params$sigmasq else {
    warning("sigmasq from Excel is NA, zero, or non-positive; cannot calculate 'size' for nbinom.")
    NA
  }
  if(is.na(size_param)) stop("size_param for nbinom is NA. Halting simulation.")
  cat("  Calculated 'size' parameter for nbinom simulation (1/sigmasq):", size_param, "\n")
  
  cat(paste("\nStep 2: Simulating", num_simulations, "INGARCH realizations of length", sim_length, "...\n"))
  sims <- vector("list", num_simulations)
  for(i in 1:num_simulations) {
    if(i == 1 || i == num_simulations || (i %% progress_print_frequency == 0 && num_simulations > 10)) {
      cat("  Generating simulation", i, "of", num_simulations, "...\n")
    }
    sim_result <- ingarch.sim(n = sim_length, param = ingarch_params, model = ingarch_model, link = "log", distr = "nbinom", size = size_param, n_start = 100)
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
  cat("Finished generating", actual_num_simulations, "valid simulations.\n")
  
  cat("\nStep 3: Running model selection in PARALLEL...\n")
  
  task_max.p <- 2
  task_max.q <- 2
  task_max.order_stepwise <- 5
  task_max.order_grid <- 4
  task_distribution <- "nbinom"
  task_link <- "log"
  task_ic <- "aic"
  
  define_results_template <- function(max_p_cols, max_q_cols) {
    na_betas_list <- list(); if(max_p_cols > 0) { na_betas_list <- as.list(rep(NA_real_, max_p_cols)); names(na_betas_list) <- paste0("beta", 1:max_p_cols) }
    na_alphas_list <- list(); if(max_q_cols > 0) { na_alphas_list <- as.list(rep(NA_real_, max_q_cols)); names(na_alphas_list) <- paste0("alpha", 1:max_q_cols) }
    template <- c(list(p_order=NA_integer_, q_order=NA_integer_, time=NA_real_, n_models_tested=NA_integer_, aic=NA_real_, bic=NA_real_, intercept=NA_real_, sigmasq=NA_real_), na_betas_list, na_alphas_list, list(status="Not Run", error_message=""))
    template[sapply(template, is.null)] <- NA
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
    populated_data <- template_data
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
        coeffs_vector <- tryCatch(stats::coef(fit_object), error = function(e) NULL)
      }
      
      current_p_order <- 0L; current_q_order <- 0L
      if (!is.null(fit_object$model)) {
        current_p_order <- if(is.null(fit_object$model$past_obs) || length(fit_object$model$past_obs)==0) 0L else as.integer(max(fit_object$model$past_obs))[1]
        current_q_order <- if(is.null(fit_object$model$past_mean) || length(fit_object$model$past_mean)==0) 0L else as.integer(max(fit_object$model$past_mean))[1]
      }
      populated_data$p_order <- current_p_order
      populated_data$q_order <- current_q_order
      populated_data$n_models_tested <- if(!is.null(fit_object$n_models_evaluated)) as.integer(fit_object$n_models_evaluated)[1] else NA_integer_
      populated_data$aic <- as.numeric(tryCatch(stats::AIC(fit_object), error = function(e) NA_real_))[1]
      populated_data$bic <- as.numeric(tryCatch(stats::BIC(fit_object), error = function(e) NA_real_))[1]
      
      if (!is.null(coeffs_vector)) {
        if ("(Intercept)" %in% names(coeffs_vector)) { populated_data$intercept <- as.numeric(coeffs_vector["(Intercept)"])[1]
        } else if ("intercept" %in% names(coeffs_vector)) { populated_data$intercept <- as.numeric(coeffs_vector["intercept"])[1] }
        
        if (max_p_cols_in_template > 0 && current_p_order > 0) {
          for (k_beta in 1:current_p_order) {
            if (k_beta > max_p_cols_in_template) break
            target_col_name_in_csv <- paste0("beta", k_beta)
            name_from_model_beta_underscore <- paste0("beta_", k_beta)
            name_from_model_past_obs <- paste0("past_obs", k_beta)
            if (name_from_model_beta_underscore %in% names(coeffs_vector)) {
              populated_data[[target_col_name_in_csv]] <- as.numeric(coeffs_vector[name_from_model_beta_underscore])[1]
            } else if (name_from_model_past_obs %in% names(coeffs_vector)) {
              populated_data[[target_col_name_in_csv]] <- as.numeric(coeffs_vector[name_from_model_past_obs])[1]
            }
          }
        }
        if (max_q_cols_in_template > 0 && current_q_order > 0) {
          for (k_alpha in 1:current_q_order) {
            if (k_alpha > max_q_cols_in_template) break
            target_col_name_in_csv <- paste0("alpha", k_alpha)
            name_from_model_alpha_underscore <- paste0("alpha_", k_alpha)
            name_from_model_past_mean <- paste0("past_mean", k_alpha)
            if (name_from_model_alpha_underscore %in% names(coeffs_vector)) {
              populated_data[[target_col_name_in_csv]] <- as.numeric(coeffs_vector[name_from_model_alpha_underscore])[1]
            } else if (name_from_model_past_mean %in% names(coeffs_vector)) {
              populated_data[[target_col_name_in_csv]] <- as.numeric(coeffs_vector[name_from_model_past_mean])[1]
            }
          }
        }
      }
      if (!is.null(fit_object$distr) && fit_object$distr == "nbinom") {
        if(!is.null(fit_object$sigmasq)) {
          populated_data$sigmasq <- as.numeric(fit_object$sigmasq)[1]
        } else if (!is.null(fit_object$distrcoefs) && "size" %in% names(fit_object$distrcoefs) && is.numeric(fit_object$distrcoefs$size) && fit_object$distrcoefs$size > 0) {
          populated_data$sigmasq <- 1 / as.numeric(fit_object$distrcoefs$size)[1]
        } else if (!is.null(model_summary_obj) && !is.null(model_summary_obj$dispersion) && "estimate" %in% names(model_summary_obj$dispersion)){
          alpha_disp <- as.numeric(model_summary_obj$dispersion["estimate"])[1]
          if(!is.na(alpha_disp)) populated_data$sigmasq <- alpha_disp
        }
      } else {
        populated_data$sigmasq <- NA_real_
      }
    } else {
      populated_data$status <- if(!is.null(fit_object$status_message)) fit_object$status_message else "Fit Error/Null"
      populated_data$error_message <- if(!is.null(fit_object$message)) as.character(fit_object$message)[1] else "Fit object error or NULL"
    }
    for(name_iter in names(template_data)){
      if(!name_iter %in% names(populated_data) || is.null(populated_data[[name_iter]])){
        populated_data[[name_iter]] <- template_data[[name_iter]]
      }
      if(name_iter %in% c("error_message", "status") && (is.na(populated_data[[name_iter]]) || is.null(populated_data[[name_iter]]))) {
        populated_data[[name_iter]] <- if(name_iter == "error_message") "" else "Undefined"
      }
    }
    return(populated_data)
  }
  
  num_cores_detected <- detectCores(logical = FALSE)
  cores_to_use <- max(1, num_cores_detected - 1)
  if(actual_num_simulations < cores_to_use) cores_to_use <- actual_num_simulations
  
  cl <- makeCluster(cores_to_use)
  registerDoParallel(cl)
  cat(paste("Registered parallel backend with", getDoParWorkers(), "cores for", actual_num_simulations, "simulations.\n"))
  
  custom_funcs_to_export <- c("auto.ingarch", "ingarch.sim", "newmodel",
                              "ingarch.string", "search.ingarch", "myingarch",
                              "populate_results_from_fit", "define_results_template")
  custom_funcs_to_export <- custom_funcs_to_export[sapply(custom_funcs_to_export, exists, envir = .GlobalEnv, mode = "function")]
  
  parallel_results_list <- foreach(
    i = 1:actual_num_simulations,
    .packages = c("tscount", "stats", "dplyr"),
    .export = c(custom_funcs_to_export,
                "sims", "task_max.p", "task_max.q", "task_max.order_stepwise", "task_max.order_grid",
                "task_distribution", "task_link", "task_ic"),
    .errorhandling = 'pass'
  ) %dopar% {
    worker_results_template <- define_results_template(task_max.p, task_max.q)
    sim_data_ts <- sims[[i]]
    
    current_stepwise_template <- worker_results_template
    stepwise_fit_object <- NULL
    stepwise_time_taken <- system.time({
      stepwise_fit_object <- tryCatch({
        auto.ingarch(y = sim_data_ts, max.p = task_max.p, max.q = task_max.q, max.order = task_max.order_stepwise,
                     distribution = task_distribution, link = task_link, ic = task_ic,
                     stepwise = TRUE, trace = FALSE, show_warnings = FALSE, parallel = FALSE)
      }, error = function(e) {
        return(list(error = TRUE, message = conditionMessage(e), status_message = "Stepwise Fit Error"))
      })
    })["elapsed"]
    current_stepwise_template$time <- as.numeric(stepwise_time_taken)
    iter_stepwise_result <- populate_results_from_fit(stepwise_fit_object, current_stepwise_template, task_max.p, task_max.q)
    iter_stepwise_result$time <- current_stepwise_template$time
    
    current_grid_template <- worker_results_template
    grid_fit_object <- NULL
    grid_time_taken <- system.time({
      grid_fit_object <- tryCatch({
        auto.ingarch(y = sim_data_ts, max.p = task_max.p, max.q = task_max.q, max.order = task_max.order_grid,
                     distribution = task_distribution, link = task_link, ic = task_ic,
                     stepwise = FALSE, trace = FALSE, show_warnings = FALSE, parallel = FALSE)
      }, error = function(e) {
        return(list(error = TRUE, message = conditionMessage(e), status_message = "Grid Fit Error"))
      })
    })["elapsed"]
    current_grid_template$time <- as.numeric(grid_time_taken)
    iter_grid_result <- populate_results_from_fit(grid_fit_object, current_grid_template, task_max.p, task_max.q)
    iter_grid_result$time <- current_grid_template$time
    if(iter_grid_result$status == "Success" && (is.na(iter_grid_result$n_models_tested) || iter_grid_result$n_models_tested == 0)){
      iter_grid_result$n_models_tested <- (task_max.p + 1) * (task_max.q + 1)
    }
    
    return(list(stepwise = iter_stepwise_result, grid = iter_grid_result))
  }
  
  stopCluster(cl)
  cat("Parallel processing finished.\n")
  
  results_list_stepwise <- vector("list", actual_num_simulations)
  results_list_grid <- vector("list", actual_num_simulations)
  
  for(i in 1:actual_num_simulations) {
    if (inherits(parallel_results_list[[i]], "error")) {
      cat("Error in parallel iteration", i, ":", conditionMessage(parallel_results_list[[i]]), "\n")
      error_template <- define_results_template(task_max.p, task_max.q)
      error_template$status <- "Parallel Iteration Error"
      error_template$error_message <- conditionMessage(parallel_results_list[[i]])
      results_list_stepwise[[i]] <- error_template
      results_list_grid[[i]] <- error_template
    } else if (is.list(parallel_results_list[[i]]) && !is.null(parallel_results_list[[i]]$stepwise) && !is.null(parallel_results_list[[i]]$grid)) {
      results_list_stepwise[[i]] <- parallel_results_list[[i]]$stepwise
      results_list_grid[[i]] <- parallel_results_list[[i]]$grid
    } else {
      cat("Unexpected result structure from parallel iteration", i, "\n")
      error_template <- define_results_template(task_max.p, task_max.q)
      error_template$status <- "Worker Result Invalid Structure"
      error_template$error_message <- "Worker did not return the expected list structure."
      results_list_stepwise[[i]] <- error_template
      results_list_grid[[i]] <- error_template
    }
  }
  
  cat("\nStep 4: Summarizing and saving results...\n")
  target_column_names <- c("method", "sim_id", names(results_template))
  
  convert_to_df_with_bind_rows <- function(list_of_result_lists_input, method_name_str, template_for_cols, target_cols_config) {
    list_of_1row_dfs <- lapply(1:length(list_of_result_lists_input), function(j) {
      current_row_data_list <- list_of_result_lists_input[[j]]
      final_df_row_list <- vector("list", length(target_cols_config))
      names(final_df_row_list) <- target_cols_config
      
      final_df_row_list[["method"]] <- method_name_str
      final_df_row_list[["sim_id"]] <- as.integer(j)
      
      for(col_name_template in names(template_for_cols)){
        if(col_name_template %in% names(current_row_data_list) && !is.null(current_row_data_list[[col_name_template]])){
          val_to_assign <- current_row_data_list[[col_name_template]]
          if(length(val_to_assign) > 1) val_to_assign <- val_to_assign[1]
          final_df_row_list[[col_name_template]] <- val_to_assign
        } else {
          final_df_row_list[[col_name_template]] <- template_for_cols[[col_name_template]]
        }
        if(col_name_template == "error_message" && (is.na(final_df_row_list[[col_name_template]]) || is.null(final_df_row_list[[col_name_template]]))) {
          final_df_row_list[[col_name_template]] <- ""
        }
        if(col_name_template == "status" && (is.na(final_df_row_list[[col_name_template]]) || is.null(final_df_row_list[[col_name_template]]))) {
          final_df_row_list[[col_name_template]] <- "Undefined"
        }
      }
      for(target_col in target_cols_config){
        if(!target_col %in% names(final_df_row_list) || is.null(final_df_row_list[[target_col]])){
          if(target_col %in% names(template_for_cols)) final_df_row_list[[target_col]] <- template_for_cols[[target_col]]
          else final_df_row_list[[target_col]] <- NA
        }
      }
      return(as.data.frame(final_df_row_list, stringsAsFactors = FALSE))
    })
    
    list_of_1row_dfs_filtered <- list_of_1row_dfs[!sapply(list_of_1row_dfs, is.null)]
    if (length(list_of_1row_dfs_filtered) > 0) {
      all_names <- unique(unlist(lapply(list_of_1row_dfs_filtered, names)))
      list_of_1row_dfs_aligned <- lapply(list_of_1row_dfs_filtered, function(df) {
        missing_cols <- setdiff(all_names, names(df))
        for (col in missing_cols) {
          if(col %in% names(template_for_cols)) df[[col]] <- template_for_cols[[col]]
          else df[[col]] <- NA
        }
        df[, all_names, drop = FALSE]
      })
      return(dplyr::bind_rows(list_of_1row_dfs_aligned))
    } else {
      empty_df_structure <- template_for_cols
      empty_df_structure$method <- character(0)
      empty_df_structure$sim_id <- integer(0)
      empty_df_final_structure <- empty_df_structure[target_cols_config]
      return(as.data.frame(empty_df_final_structure))
    }
  }
  
  stepwise_results_df <- convert_to_df_with_bind_rows(results_list_stepwise, "stepwise", results_template, target_column_names)
  grid_results_df     <- convert_to_df_with_bind_rows(results_list_grid, "grid_search", results_template, target_column_names)
  
  results_df <- dplyr::bind_rows(stepwise_results_df, grid_results_df)
  
  output_dir_name <- "ingarch_no_covariates_results_parallel"
  if (!dir.exists(output_dir_name)) {
    dir.create(output_dir_name, recursive = TRUE)
  }
  
  results_csv_path <- file.path(output_dir_name, "ingarch_no_covariates_results.csv")
  sims_rds_path <- file.path(output_dir_name, "ingarch_no_covariates_simulations.rds")
  summary_stats_csv_path <- file.path(output_dir_name, "ingarch_no_covariates_summary_stats.csv")
  order_freq_csv_path <- file.path(output_dir_name, "ingarch_no_covariates_order_frequencies.csv")
  
  if (nrow(results_df) > 0) {
    if("status" %in% names(results_df) && is.factor(results_df$status)) results_df$status <- as.character(results_df$status)
    if("error_message" %in% names(results_df) && is.factor(results_df$error_message)) results_df$error_message <- as.character(results_df$error_message)
    
    write.csv(results_df, results_csv_path, row.names = FALSE, na = "")
    cat("  Detailed results saved to:", results_csv_path, "\n")
  } else {
    cat("  No results to save to CSV.\n")
  }
  
  saveRDS(sims, sims_rds_path)
  cat("  Simulated datasets saved to:", sims_rds_path, "\n")
  
  if (nrow(results_df) > 0) {
    numeric_cols_for_summary <- c("p_order", "q_order", "time", "n_models_tested", "aic", "bic", "intercept", "sigmasq")
    if(task_max.p > 0) numeric_cols_for_summary <- c(numeric_cols_for_summary, paste0("beta", 1:task_max.p))
    if(task_max.q > 0) numeric_cols_for_summary <- c(numeric_cols_for_summary, paste0("alpha", 1:task_max.q))
    
    for(col_s in numeric_cols_for_summary){
      if(col_s %in% names(results_df)){
        results_df[[col_s]] <- suppressWarnings(as.numeric(as.character(results_df[[col_s]])))
      }
    }
    
    results_df$status <- as.character(results_df$status)
    
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
    
    cat("\n===== Summary Statistics (Parallel) =====\n"); print(summary_stats)
    cat("\n===== Most Frequent (p_order,q_order) Orders (Top 5 per method - Parallel) =====\n")
    if(nrow(order_freq) > 0) {
      top_orders <- order_freq %>% group_by(method) %>% slice_head(n = 5)
      print(top_orders)
    } else {
      cat("No successful fits to determine frequent orders.\n")
    }
    cat("\n===== Status Summary (Parallel) =====\n")
    print(results_df %>% group_by(method, status) %>% summarise(count = n(), .groups = 'drop') %>% arrange(method, desc(count)))
    
  } else {
    cat("  No results available to summarize.\n")
  }
  
  return(results_df)
}

main_parallel_no_cov <- function() {
  if(exists(".coef_names_printed_flag_v4", envir = .GlobalEnv)){
    rm(".coef_names_printed_flag_v4", envir = .GlobalEnv)
  }
  
  overall_start_time <- Sys.time()
  results_data_frame <- tryCatch(run_simulation_study_no_covariates_parallel(),
                                 error = function(e) {
                                   cat("\nERROR during PARALLEL simulation study run (no covariates):\n")
                                   print(e)
                                   if(exists("cl") && inherits(cl, "cluster")) try(stopCluster(cl), silent=TRUE)
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

main_parallel_no_cov()