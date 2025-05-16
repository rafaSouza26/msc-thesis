# INGARCH Model Simulation Study - Without Covariates
# Comparing model selection methods: auto.ingarch with stepwise vs grid search

# Load required packages
library(tscount)
library(dplyr)   
library(readxl)  
# library(ggplot2) # Only if plotting is added later

# Source custom functions
# Ensure these paths are correct relative to your script's location
source("./custom/auto.ingarch.R")
source("./custom/ingarch.sim.R")
source("./custom/newmodel.R")
source("./custom/ingarch.string.R")
source("./custom/search.ingarch.R")

# Function to extract model parameters from Excel file
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

# Main simulation study - Without Covariates
run_simulation_study_no_covariates <- function() {
  set.seed(12345) 
  
  num_simulations <- 1000 # Number of datasets to simulate
  sim_length <- 1000     # Length of each time series
  progress_print_frequency <- max(1, floor(num_simulations / 10)) # Print progress roughly 10 times + start/end
  
  cat("Starting INGARCH simulation study (without covariates).\n")
  cat("Configuration: num_simulations =", num_simulations, ", sim_length =", sim_length, "\n")
  
  cat("Step 1: Extracting model parameters from Excel file...\n")
  params <- extract_model_params("./data/modelosAveiro.xlsx", model_row = 1) # Ensure this path is correct
  # Minimal print of true parameters for verification
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
      stop(paste("Cannot extract numeric time series from ingarch.sim output in simulation", i))
    }
  }
  cat("Finished generating simulations.\n")
  
  cat("\nStep 3: Running model selection sequentially...\n")
  
  results_list_stepwise <- vector("list", num_simulations)
  results_list_grid <- vector("list", num_simulations)
  
  task_max.p <- 7 
  task_max.q <- 7 
  task_max.order_stepwise <- 5 
  task_max.order_grid <- 14    
  task_distribution <- "nbinom"
  task_link <- "log"
  
  define_results_template <- function(max_p_cols, max_q_cols) {
    na_betas_list <- list(); if(max_p_cols > 0) { na_betas_list <- as.list(rep(NA_real_, max_p_cols)); names(na_betas_list) <- paste0("beta", 1:max_p_cols) }
    na_alphas_list <- list(); if(max_q_cols > 0) { na_alphas_list <- as.list(rep(NA_real_, max_q_cols)); names(na_alphas_list) <- paste0("alpha", 1:max_q_cols) }
    template <- c(list(p_order=NA_integer_, q_order=NA_integer_, time=NA_real_, n_models_tested=NA_integer_, aic=NA_real_, bic=NA_real_, intercept=NA_real_, sigmasq=NA_real_), na_betas_list, na_alphas_list, list(error_message=""))
    template[sapply(template, is.null)] <- NA # Ensure NULLs become NAs
    # Explicitly type NAs
    for (name in c("p_order", "q_order", "n_models_tested")) template[[name]] <- NA_integer_
    for (name in c("time", "aic", "bic", "intercept", "sigmasq")) template[[name]] <- NA_real_
    if(max_p_cols > 0) for(k in 1:max_p_cols) template[[paste0("beta",k)]] <- NA_real_
    if(max_q_cols > 0) for(k in 1:max_q_cols) template[[paste0("alpha",k)]] <- NA_real_
    template$error_message <- ""
    return(template)
  }
  results_template <- define_results_template(task_max.p, task_max.q)
  target_column_names <- c("method", "sim_id", names(results_template))
  
  populate_results_from_fit <- function(fit_object, template_data, 
                                        max_p_cols_in_template, max_q_cols_in_template) { # Removed unused debug args
    populated_data <- template_data 
    if (!is.null(fit_object) && is.null(fit_object$error)) {
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
      if (!is.null(fit_object$distr) && fit_object$distr == "nbinom" && !is.null(fit_object$sigmasq)) {
        populated_data$sigmasq <- as.numeric(fit_object$sigmasq)[1]
      } else if (!is.null(coeffs_vector)) { 
        if ("sigmasq" %in% names(coeffs_vector)) { populated_data$sigmasq <- as.numeric(coeffs_vector["sigmasq"])[1] 
        } else if ("alpha" %in% names(coeffs_vector)) { populated_data$sigmasq <- as.numeric(coeffs_vector["alpha"])[1] 
        } else if ("dispersion" %in% names(coeffs_vector)) { populated_data$sigmasq <- as.numeric(coeffs_vector["dispersion"])[1] }
      }
    } else { 
      populated_data$error_message <- if(!is.null(fit_object$message)) as.character(fit_object$message)[1] else "Fit object error or NULL"
    }
    for(name_iter in names(results_template)){ 
      current_val <- populated_data[[name_iter]]
      if(length(current_val) != 1 || is.null(current_val)) { 
        original_template_val <- results_template[[name_iter]] 
        if(identical(class(original_template_val), "integer")) populated_data[[name_iter]] <- NA_integer_
        else if(identical(class(original_template_val), "numeric")) populated_data[[name_iter]] <- NA_real_
        else if(identical(class(original_template_val), "character")) {
          populated_data[[name_iter]] <- if(name_iter == "error_message") "" else NA_character_
        }
        else populated_data[[name_iter]] <- NA 
      }
      if(name_iter == "error_message" && is.na(populated_data[[name_iter]])) {
        populated_data[[name_iter]] <- "" 
      }
    }
    return(populated_data)
  }
  
  # Definition of convert_to_df_with_bind_rows (moved earlier)
  convert_to_df_with_bind_rows <- function(list_of_result_lists_input, method_name_str) {
    list_of_1row_dfs <- lapply(1:length(list_of_result_lists_input), function(j) {
      current_row_data_list <- list_of_result_lists_input[[j]] 
      
      final_df_row_list <- vector("list", length(target_column_names)) 
      names(final_df_row_list) <- target_column_names
      
      final_df_row_list[["method"]] <- method_name_str
      final_df_row_list[["sim_id"]] <- as.integer(j)
      
      for(col_name_template in names(results_template)){ 
        if(col_name_template %in% names(current_row_data_list)){
          final_df_row_list[[col_name_template]] <- current_row_data_list[[col_name_template]]
        } else {
          final_df_row_list[[col_name_template]] <- results_template[[col_name_template]] 
          if(col_name_template == "error_message" && (is.na(final_df_row_list[[col_name_template]]))) { # Ensure empty string for error if NA
            final_df_row_list[[col_name_template]] <- "" 
          }
        }
      }
      return(as.data.frame(final_df_row_list, stringsAsFactors = FALSE))
    })
    
    list_of_1row_dfs_filtered <- list_of_1row_dfs[!sapply(list_of_1row_dfs, is.null)] # Should not have NULLs if template is always used
    
    if (length(list_of_1row_dfs_filtered) > 0) {
      return(dplyr::bind_rows(list_of_1row_dfs_filtered))
    } else {
      empty_df <- data.frame(matrix(ncol = length(target_column_names), nrow = 0))
      tryCatch(colnames(empty_df) <- target_column_names, error = function(e){})
      return(empty_df)
    }
  } # End of convert_to_df_with_bind_rows definition
  
  # Sequential loop for simulations
  for(i in 1:num_simulations) {
    if(i == 1 || i == num_simulations || (i %% progress_print_frequency == 0 && num_simulations > 10)) {
      cat("  Processing simulation", i, "of", num_simulations, "...\n")
    }
    sim_data_ts <- sims[[i]]
    
    current_stepwise_template <- results_template 
    stepwise_time_taken <- system.time({
      stepwise_fit_object <- tryCatch({
        auto.ingarch(y = sim_data_ts, max.p = task_max.p, max.q = task_max.q, max.order = task_max.order_stepwise,
                     distribution = task_distribution, link = task_link, ic = "aic",
                     stepwise = TRUE, trace = FALSE, show_warnings = FALSE, parallel = FALSE)
      }, error = function(e) list(error = TRUE, message = conditionMessage(e)))
    })["elapsed"]
    current_stepwise_template$time <- as.numeric(stepwise_time_taken)
    results_list_stepwise[[i]] <- populate_results_from_fit(stepwise_fit_object, current_stepwise_template, 
                                                            task_max.p, task_max.q) # Removed debug args
    
    current_grid_template <- results_template 
    grid_time_taken <- system.time({
      grid_fit_object <- tryCatch({
        auto.ingarch(y = sim_data_ts, max.p = task_max.p, max.q = task_max.q, max.order = task_max.order_grid,
                     distribution = task_distribution, link = task_link, ic = "aicc",
                     stepwise = FALSE, trace = FALSE, show_warnings = FALSE, parallel = FALSE)
      }, error = function(e) list(error = TRUE, message = conditionMessage(e)))
    })["elapsed"]
    current_grid_template$time <- as.numeric(grid_time_taken)
    results_list_grid[[i]] <- populate_results_from_fit(grid_fit_object, current_grid_template, 
                                                        task_max.p, task_max.q) # Removed debug args
  }
  cat("Finished model selection processing.\n")
  
  cat("\nStep 4: Summarizing and saving results...\n")
  
  stepwise_results_df <- convert_to_df_with_bind_rows(results_list_stepwise, "stepwise")
  grid_results_df     <- convert_to_df_with_bind_rows(results_list_grid, "grid_search")
  
  results_df <- dplyr::bind_rows(stepwise_results_df, grid_results_df)
  
  if (nrow(results_df) > 0) {
    write.csv(results_df, "ingarch_no_covariates_results.csv", row.names = FALSE)
    cat("  Detailed results saved to: ingarch_no_covariates_results.csv\n")
  } else {
    cat("  No results to save to CSV.\n")
  }
  
  saveRDS(sims, "ingarch_no_covariates_simulations.rds")
  cat("  Simulated datasets saved to: ingarch_no_covariates_simulations.rds\n")
  
  if (nrow(results_df) > 0) {
    summary_stats <- results_df %>%
      group_by(method) %>%
      summarize(
        mean_p_order = mean(p_order, na.rm = TRUE), median_p_order = median(p_order, na.rm = TRUE), sd_p_order = sd(p_order, na.rm = TRUE),
        mean_q_order = mean(q_order, na.rm = TRUE), median_q_order = median(q_order, na.rm = TRUE), sd_q_order = sd(q_order, na.rm = TRUE),
        mean_time_secs = mean(time, na.rm = TRUE), median_time_secs = median(time, na.rm = TRUE),
        mean_fitted_models = mean(n_models_tested, na.rm = TRUE), median_fitted_models = median(n_models_tested, na.rm = TRUE),
        successful_fits = sum(is.na(error_message) | error_message == ""),
        failed_fits = sum(!is.na(error_message) & error_message != ""),
        .groups = 'drop'
      )
    
    order_freq <- results_df %>%
      filter(is.na(error_message) | error_message == "") %>% 
      filter(!is.na(p_order) & !is.na(q_order)) %>%
      group_by(method, p_order, q_order) %>%
      summarize(count = n(), .groups = 'drop') %>% 
      group_by(method) %>%
      mutate(freq = count / sum(count) * 100) %>%
      arrange(method, desc(freq))
    
    write.csv(summary_stats, "ingarch_no_covariates_summary_stats.csv", row.names = FALSE)
    write.csv(order_freq, "ingarch_no_covariates_order_frequencies.csv", row.names = FALSE)
    cat("  Summary statistics and order frequencies saved to CSV files.\n")
    
    cat("\n===== Summary Statistics =====\n"); print(summary_stats)
    cat("\n===== Most Frequent (p_order,q_order) Orders (Top 5 per method) =====\n")
    top_orders <- order_freq %>% group_by(method) %>% slice_head(n = 5)
    print(top_orders)
  } else {
    cat("  No results available to summarize.\n")
  }
  
  return(results_df)
}

# Main function to run the study
main <- function() {
  # Clear the global flag for diagnostic print before starting, if it was used
  if(exists(".coef_names_printed_flag_v4", envir = .GlobalEnv)){
    rm(".coef_names_printed_flag_v4", envir = .GlobalEnv)
  }
  
  overall_start_time <- Sys.time()
  results_data_frame <- tryCatch(run_simulation_study_no_covariates(),
                                 error = function(e) {
                                   cat("\nERROR during simulation study run:\n")
                                   print(e)
                                   return(NULL)
                                 })
  overall_end_time <- Sys.time()
  overall_duration <- overall_end_time - overall_start_time
  
  if(!is.null(results_data_frame)){
    cat("\nSimulation study completed!\n")
  } else {
    cat("\nSimulation study encountered an error and did not complete successfully.\n")
  }
  cat("Total duration of the study:", format(overall_duration), "\n")
}

# Execute main function
main()