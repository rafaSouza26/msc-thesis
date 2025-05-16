# INGARCH Model Simulation Study - Without Covariates
# Comparing model selection methods: auto.ingarch with stepwise vs grid search

# Load required packages
library(tscount)
library(ggplot2) 
library(dplyr)   
library(readxl)  

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
  
  cat("Extracting model parameters from Excel file...\n")
  params <- extract_model_params("./data/modelosAveiro.xlsx", model_row = 1) # Ensure this path is correct
  cat("Model parameters for simulation:\n")
  cat("True p =", params$p, "\n") 
  cat("True q =", params$q, "\n")
  cat("True sigmasq (alpha for nbinom) =", params$sigmasq, "\n")
  cat("True intercept =", params$intercept, "\n")
  if(params$p > 0 && length(params$betas) >= params$p) cat("True beta coefficients:", params$betas[1:params$p], "\n")
  if(params$q > 0 && length(params$alphas) >= params$q) cat("True alpha coefficients:", params$alphas[1:params$q], "\n")
  
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
  cat("Calculated 'size' parameter for nbinom simulation (1/sigmasq):", size_param, "\n")
  
  num_simulations <- 5 # Using 1 from your traceback, for faster testing. Increase for full study.
  sim_length <- 100    # Using 10 from your traceback, for faster testing. Increase for full study.
  cat(paste("Simulating", num_simulations, "INGARCH realizations of length", sim_length, "without covariates...\n"))
  sims <- vector("list", num_simulations)
  for(i in 1:num_simulations) {
    if(i %% 1 == 0 || num_simulations <= 10) cat("Generating simulation", i, "of", num_simulations, "\n")
    sim_result <- ingarch.sim(n = sim_length, param = ingarch_params, model = ingarch_model, link = "log", distr = "nbinom", size = size_param, n_start = 100)
    if (!is.null(sim_result) && is.list(sim_result) && "ts" %in% names(sim_result)) {
      sims[[i]] <- as.numeric(sim_result$ts)
    } else {
      stop(paste("Cannot extract numeric time series from ingarch.sim output in simulation", i, "- expected a list with a '$ts' component."))
    }
  }
  cat("Finished generating simulations.\n")
  
  cat("Running model selection sequentially for INGARCH without covariates...\n")
  
  results_list_stepwise <- vector("list", num_simulations)
  results_list_grid <- vector("list", num_simulations)
  
  task_max.p <- 7 
  task_max.q <- 7 
  task_max.order_stepwise <- 5 
  task_max.order_grid <- 5    
  task_distribution <- "nbinom"
  task_link <- "log"
  
  define_results_template <- function(max_p_cols, max_q_cols) {
    na_betas_list <- list(); if(max_p_cols > 0) { na_betas_list <- as.list(rep(NA_real_, max_p_cols)); names(na_betas_list) <- paste0("beta", 1:max_p_cols) }
    na_alphas_list <- list(); if(max_q_cols > 0) { na_alphas_list <- as.list(rep(NA_real_, max_q_cols)); names(na_alphas_list) <- paste0("alpha", 1:max_q_cols) }
    
    template <- c(
      list(p_order = NA_integer_, q_order = NA_integer_, 
           time = NA_real_, n_models_tested = NA_integer_, 
           aic = NA_real_, bic = NA_real_,
           intercept = NA_real_, sigmasq = NA_real_),
      na_betas_list, 
      na_alphas_list, 
      list(error_message = "") 
    )
    template$p_order <- NA_integer_ # Ensure types after c()
    template$q_order <- NA_integer_
    template$time <- NA_real_
    template$n_models_tested <- NA_integer_
    template$aic <- NA_real_
    template$bic <- NA_real_
    template$intercept <- NA_real_
    template$sigmasq <- NA_real_
    template$error_message <- "" 
    if(max_p_cols > 0) for(k_template in 1:max_p_cols) template[[paste0("beta",k_template)]] <- NA_real_
    if(max_q_cols > 0) for(k_template in 1:max_q_cols) template[[paste0("alpha",k_template)]] <- NA_real_
    return(template)
  }
  results_template <- define_results_template(task_max.p, task_max.q)
  target_column_names <- c("method", "sim_id", names(results_template))
  
  populate_results_from_fit <- function(fit_object, template_data, 
                                        max_p_for_columns, max_q_for_columns, # These are task_max.p/q for template size
                                        sim_idx_for_debug = NA, method_for_debug = NA) {
    populated_data <- template_data # Start with a fresh copy of the results_template
    
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
      
      # --- DIAGNOSTIC PRINT (prints once for the first successful fit encountered) ---
      if (!is.na(sim_idx_for_debug) && !exists(".coef_names_printed_flag_v4", envir = .GlobalEnv)) {
        if (!is.null(coeffs_vector) && length(coeffs_vector) > 0) {
          cat(paste0("\n[DIAGNOSTIC PRINT populate_results_from_fit] Method: ", method_for_debug, ", Sim_ID: ", sim_idx_for_debug, "\n"))
          cat("  Actual coefficient names found in fitted model object:\n")
          print(names(coeffs_vector))
          cat("  Corresponding values:\n")
          print(coeffs_vector)
          cat("  >>> IMPORTANT: Script is looking for names like '(Intercept)', 'beta_K', 'alpha_K', 'past_obsK', 'past_meanK'. Adjust if needed.\n")
        } else {
          cat(paste0("\n[DIAGNOSTIC PRINT populate_results_from_fit] Method: ", method_for_debug, ", Sim_ID: ", sim_idx_for_debug, " - No coefficients retrieved.\n"))
        }
        assign(".coef_names_printed_flag_v4", TRUE, envir = .GlobalEnv) 
      }
      # --- END DIAGNOSTIC PRINT ---
      
      current_p_order <- 0L
      current_q_order <- 0L
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
        # Intercept
        if ("(Intercept)" %in% names(coeffs_vector)) { populated_data$intercept <- as.numeric(coeffs_vector["(Intercept)"])[1]
        } else if ("intercept" %in% names(coeffs_vector)) { populated_data$intercept <- as.numeric(coeffs_vector["intercept"])[1] }
        
        # Betas (past_obs coefficients)
        if (max_p_for_columns > 0 && current_p_order > 0) {
          for (k_beta in 1:current_p_order) { # Iterate up to the actual fitted p_order
            if (k_beta > max_p_for_columns) break 
            target_col_name_in_csv <- paste0("beta", k_beta) # e.g. "beta1"
            
            # ** Coefficient names to search for in the model object **
            # Based on diagnostic print, "beta_K" is a primary candidate.
            name_from_model_beta_underscore <- paste0("beta_", k_beta) 
            name_from_model_past_obs <- paste0("past_obs", k_beta)
            
            if (name_from_model_beta_underscore %in% names(coeffs_vector)) {
              populated_data[[target_col_name_in_csv]] <- as.numeric(coeffs_vector[name_from_model_beta_underscore])[1]
            } else if (name_from_model_past_obs %in% names(coeffs_vector)) {
              populated_data[[target_col_name_in_csv]] <- as.numeric(coeffs_vector[name_from_model_past_obs])[1]
            } 
          }
        }
        # Alphas (past_mean coefficients)
        if (max_q_for_columns > 0 && current_q_order > 0) {
          for (k_alpha in 1:current_q_order) { 
            if (k_alpha > max_q_for_columns) break
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
      
      # Sigmasq (overdispersion for nbinom)
      if (!is.null(fit_object$distr) && fit_object$distr == "nbinom" && !is.null(fit_object$sigmasq)) {
        populated_data$sigmasq <- as.numeric(fit_object$sigmasq)[1]
      } else if (!is.null(coeffs_vector)) { 
        if ("sigmasq" %in% names(coeffs_vector)) { populated_data$sigmasq <- as.numeric(coeffs_vector["sigmasq"])[1] 
        } else if ("alpha" %in% names(coeffs_vector)) { # Common name for overdispersion
          populated_data$sigmasq <- as.numeric(coeffs_vector["alpha"])[1] 
        } else if ("dispersion" %in% names(coeffs_vector)) { # Another possible name
          populated_data$sigmasq <- as.numeric(coeffs_vector["dispersion"])[1]
        }
      }
    } else { 
      populated_data$error_message <- if(!is.null(fit_object$message)) as.character(fit_object$message)[1] else "Fit object error or NULL"
    }
    
    # Final check: ensure all elements of populated_data are scalar, otherwise NA of template type
    for(name_iter in names(results_template)){ 
      val_to_check <- populated_data[[name_iter]]
      if(is.null(val_to_check) || length(val_to_check) != 1) { 
        original_template_val <- results_template[[name_iter]] 
        if(identical(class(original_template_val), "integer")) populated_data[[name_iter]] <- NA_integer_
        else if(identical(class(original_template_val), "numeric")) populated_data[[name_iter]] <- NA_real_
        else if(identical(class(original_template_value), "character")) { # Typo fix: original_template_val
          populated_data[[name_iter]] <- if(name_iter == "error_message") "" else NA_character_
        }
        else populated_data[[name_iter]] <- NA 
      }
      if(name_iter == "error_message" && is.na(populated_data[[name_iter]])) {
        populated_data[[name_iter]] <- "" # Ensure error_message is "" if NA
      }
    }
    return(populated_data)
  }
  
  # Sequential loop for simulations
  for(i in 1:num_simulations) {
    if(i %% 1 == 0 || num_simulations <= 10 ) cat("Processing simulation", i, "of", num_simulations, "\n")
    sim_data_ts <- sims[[i]]
    
    # Get fresh templates for each fit result
    current_stepwise_template <- results_template 
    stepwise_time_taken <- system.time({
      stepwise_fit_object <- tryCatch({
        auto.ingarch(y = sim_data_ts, max.p = task_max.p, max.q = task_max.q, max.order = task_max.order_stepwise,
                     distribution = task_distribution, link = task_link, ic = "aic",
                     stepwise = TRUE, trace = FALSE, show_warnings = FALSE, parallel = FALSE)
      }, error = function(e) list(error = TRUE, message = conditionMessage(e)))
    })["elapsed"]
    current_stepwise_template$time <- as.numeric(stepwise_time_taken) # Add time before passing to populate
    results_list_stepwise[[i]] <- populate_results_from_fit(stepwise_fit_object, current_stepwise_template, 
                                                            task_max.p, task_max.q, 
                                                            sim_idx_for_debug = i, method_for_debug = "stepwise")
    
    current_grid_template <- results_template 
    grid_time_taken <- system.time({
      grid_fit_object <- tryCatch({
        auto.ingarch(y = sim_data_ts, max.p = task_max.p, max.q = task_max.q, max.order = task_max.order_grid,
                     distribution = task_distribution, link = task_link, ic = "aicc",
                     stepwise = FALSE, trace = FALSE, show_warnings = FALSE, parallel = FALSE)
      }, error = function(e) list(error = TRUE, message = conditionMessage(e)))
    })["elapsed"]
    current_grid_template$time <- as.numeric(grid_time_taken) # Add time before passing to populate
    results_list_grid[[i]] <- populate_results_from_fit(grid_fit_object, current_grid_template, 
                                                        task_max.p, task_max.q, 
                                                        sim_idx_for_debug = i, method_for_debug = "grid_search")
  }
  cat("Finished model selection processing.\n")
  
  cat("Summarizing results for INGARCH without covariates...\n")
  
  convert_to_df_with_bind_rows <- function(list_of_result_lists_input, method_name_str) {
    list_of_1row_dfs <- lapply(1:length(list_of_result_lists_input), function(j) {
      current_row_data_list <- list_of_result_lists_input[[j]] 
      
      final_df_row_list <- vector("list", length(target_column_names)) # target_column_names from outer scope
      names(final_df_row_list) <- target_column_names
      
      final_df_row_list[["method"]] <- method_name_str
      final_df_row_list[["sim_id"]] <- as.integer(j)
      
      for(col_name_template in names(results_template)){ # Iterate over defined data fields
        if(col_name_template %in% names(current_row_data_list)){
          final_df_row_list[[col_name_template]] <- current_row_data_list[[col_name_template]]
        } else {
          final_df_row_list[[col_name_template]] <- results_template[[col_name_template]] # Fallback to NA from template
          if(col_name_template == "error_message" && (is.na(final_df_row_list[[col_name_template]]))) {
            final_df_row_list[[col_name_template]] <- "Data structure error" 
          }
        }
      }
      return(as.data.frame(final_df_row_list, stringsAsFactors = FALSE))
    })
    
    # Filter out any NULL data frames that might have been produced (should be rare now)
    list_of_1row_dfs_filtered <- list_of_1row_dfs[!sapply(list_of_1row_dfs, is.null)]
    
    if (length(list_of_1row_dfs_filtered) > 0) {
      return(dplyr::bind_rows(list_of_1row_dfs_filtered))
    } else {
      empty_df <- data.frame(matrix(ncol = length(target_column_names), nrow = 0))
      tryCatch(colnames(empty_df) <- target_column_names, error = function(e){})
      return(empty_df)
    }
  }
  
  stepwise_results_df <- convert_to_df_with_bind_rows(results_list_stepwise, "stepwise")
  grid_results_df     <- convert_to_df_with_bind_rows(results_list_grid, "grid_search")
  
  results_df <- dplyr::bind_rows(stepwise_results_df, grid_results_df)
  
  cat("Saving detailed results to CSV file...\n")
  if (nrow(results_df) > 0) {
    write.csv(results_df, "ingarch_no_covariates_results.csv", row.names = FALSE)
  } else {
    cat("No results to save to CSV.\n")
  }
  
  saveRDS(sims, "ingarch_no_covariates_simulations.rds")
  
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
    
    cat("\n===== Summary Statistics =====\n"); print(summary_stats)
    cat("\n===== Most Frequent (p_order,q_order) Orders (Top 5 per method) =====\n")
    top_orders <- order_freq %>% group_by(method) %>% slice_head(n = 5)
    print(top_orders)
  } else {
    cat("No results available to summarize.\n")
  }
  
  cat("\nResults primarily saved to: ingarch_no_covariates_results.csv\n")
  cat("Other .rds and summary .csv files also saved in your working directory.\n")
  
  return(results_df)
}

# Main function to run the study
main <- function() {
  cat("Starting INGARCH simulation study (without covariates)...\n")
  overall_start_time <- Sys.time()
  if(exists(".coef_names_printed_flag_v4", envir = .GlobalEnv)){ # Ensure flag name matches
    rm(".coef_names_printed_flag_v4", envir = .GlobalEnv)
  }
  results_data_frame <- run_simulation_study_no_covariates()
  overall_end_time <- Sys.time()
  overall_duration <- overall_end_time - overall_start_time
  cat("Simulation study completed!\n")
  cat("Total duration of the study:", format(overall_duration), "\n")
}

# Execute main function
main()