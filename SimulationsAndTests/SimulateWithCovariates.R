# INGARCH Model Simulation Study - With Covariates
# Comparing model selection methods: auto.ingarch with stepwise vs grid search
# Version: Refactored for parameter structure consistency & text progress

# Load required packages
library(tscount)
library(dplyr)
library(readxl)
library(doParallel) # Still needed for parallel execution
library(foreach)    # Still needed for parallel execution
# Removed: library(progressr)

# Source custom functions
source("./custom/auto.ingarch.R")
source("./custom/ingarch.sim.R")
source("./custom/newmodel.R")
source("./custom/ingarch.string.R")
source("./custom/search.ingarch.R")
source("./custom/myingarch.R")

# Helper function to define a template for storing detailed results
define_results_template_with_cov <- function(max_p_cols, max_q_cols, xreg_col_names) {
  na_betas_list <- list()
  if(max_p_cols > 0) {
    na_betas_list <- as.list(rep(NA_real_, max_p_cols))
    names(na_betas_list) <- paste0("beta", 1:max_p_cols)
  }
  
  na_alphas_list <- list()
  if(max_q_cols > 0) {
    na_alphas_list <- as.list(rep(NA_real_, max_q_cols))
    names(na_alphas_list) <- paste0("alpha", 1:max_q_cols)
  }
  
  na_xreg_coefs_list <- list()
  if(length(xreg_col_names) > 0) {
    na_xreg_coefs_list <- as.list(rep(NA_real_, length(xreg_col_names)))
    valid_xreg_col_names <- make.names(xreg_col_names, unique = TRUE)
    names(na_xreg_coefs_list) <- valid_xreg_col_names
  }
  
  template <- c(
    list(p_order=NA_integer_, q_order=NA_integer_, time=NA_real_, 
         n_models_tested=NA_integer_, aic=NA_real_, bic=NA_real_, 
         intercept=NA_real_, sigmasq=NA_real_),
    na_betas_list,
    na_alphas_list,
    na_xreg_coefs_list,
    list(status="Not Run", error_message="") 
  )
  
  for (name in c("p_order", "q_order", "n_models_tested")) template[[name]] <- NA_integer_
  for (name in c("time", "aic", "bic", "intercept", "sigmasq")) template[[name]] <- NA_real_
  if(max_p_cols > 0) for(k in 1:max_p_cols) template[[paste0("beta",k)]] <- NA_real_
  if(max_q_cols > 0) for(k in 1:max_q_cols) template[[paste0("alpha",k)]] <- NA_real_
  if(length(xreg_col_names) > 0) {
    for(col_name in names(na_xreg_coefs_list)) template[[col_name]] <- NA_real_
  }
  template$status <- "Not Run" 
  template$error_message <- ""
  
  return(template)
}

# Helper function to populate the results template from a fitted model
populate_results_from_fit_with_cov <- function(fit_object, template_data, 
                                               max_p_cols_in_template, max_q_cols_in_template,
                                               xreg_col_names_in_template) {
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
      if ("(Intercept)" %in% names(coeffs_vector)) { 
        populated_data$intercept <- as.numeric(coeffs_vector["(Intercept)"])[1]
      } else if ("intercept" %in% names(coeffs_vector)) { 
        populated_data$intercept <- as.numeric(coeffs_vector["intercept"])[1] 
      }
      
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
      if (length(xreg_col_names_in_template) > 0 && !is.null(fit_object$xreg)) {
        original_xreg_names_from_model <- colnames(fit_object$xreg) 
        for (idx_xreg in 1:length(xreg_col_names_in_template)) {
          original_coef_name <- xreg_col_names_in_template[idx_xreg] 
          valid_template_col_name <- make.names(original_coef_name, unique = TRUE)
          if (original_coef_name %in% names(coeffs_vector)) {
            if(valid_template_col_name %in% names(populated_data)){ 
              populated_data[[valid_template_col_name]] <- as.numeric(coeffs_vector[original_coef_name])[1]
            }
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
    current_status <- if(!is.null(template_data$status) && template_data$status != "Not Run") template_data$status else "Fit Error/Null"
    current_error_message <- if(!is.null(fit_object$message)) as.character(fit_object$message)[1] else if(!is.null(template_data$error_message) && template_data$error_message != "") template_data$error_message else "Fit object error or NULL"
    
    populated_data <- template_data 
    populated_data$time <- template_data$time 
    populated_data$status <- current_status
    populated_data$error_message <- current_error_message
  }
  
  for(name_iter in names(results_template)) { 
    if(!name_iter %in% names(populated_data) || is.null(populated_data[[name_iter]]) || length(populated_data[[name_iter]]) == 0) {
      original_template_val <- results_template[[name_iter]] 
      if(identical(class(original_template_val), "integer")) populated_data[[name_iter]] <- NA_integer_
      else if(identical(class(original_template_val), "numeric")) populated_data[[name_iter]] <- NA_real_
      else if(identical(class(original_template_val), "character")) {
        populated_data[[name_iter]] <- if(name_iter %in% c("status", "error_message")) "" else NA_character_
      }
      else populated_data[[name_iter]] <- NA 
    }
    if(name_iter == "error_message" && (is.na(populated_data[[name_iter]]) || is.null(populated_data[[name_iter]]))) {
      populated_data[[name_iter]] <- "" 
    }
    if(name_iter == "status" && (is.na(populated_data[[name_iter]]) || is.null(populated_data[[name_iter]]))) {
      populated_data[[name_iter]] <- "Undefined"
    }
  }
  return(populated_data)
}

extract_model_params_with_covariates <- function(data_path, model_row = 2,
                                                 covariate_cols = c("Temp", "DTemp", "PM", "NO2")) {
  excel_data <- readxl::read_excel(data_path)
  cat("Excel file columns:", paste(colnames(excel_data), collapse=", "), "\n")
  cat("Excel file rows:", nrow(excel_data), "\n")
  if (model_row > nrow(excel_data)) {
    stop("model_row index out of bounds for the loaded Excel file.")
  }
  model_data <- excel_data[model_row, ]
  
  p <- as.numeric(model_data$p)
  q <- as.numeric(model_data$q)
  sigmasq <- as.numeric(model_data$sigmasq)
  intercept <- as.numeric(model_data$Intercept)
  
  beta_cols <- grep("^beta", colnames(excel_data), value = TRUE) 
  betas <- if(p > 0 && length(beta_cols) >= p) as.numeric(unlist(model_data[, beta_cols[1:p]])) else if (p == 0) NULL else {
    warning(paste("Found fewer beta columns (", length(beta_cols), ") than p=", p, ". Using available or NULL."))
    if(length(beta_cols) > 0) as.numeric(unlist(model_data[, beta_cols])) else NULL
  }
  
  alpha_cols <- grep("^alpha", colnames(excel_data), value = TRUE) 
  alphas <- if(q > 0 && length(alpha_cols) >= q) as.numeric(unlist(model_data[, alpha_cols[1:q]])) else if (q == 0) NULL else {
    warning(paste("Found fewer alpha columns (", length(alpha_cols), ") than q=", q, ". Using available or NULL."))
    if(length(alpha_cols) > 0) as.numeric(unlist(model_data[, alpha_cols])) else NULL
  }
  
  external_coefs <- numeric(length(covariate_cols))
  names(external_coefs) <- covariate_cols
  missing_cov_cols <- c()
  for(col_name in covariate_cols) {
    if(col_name %in% colnames(model_data)) {
      external_coefs[col_name] <- as.numeric(model_data[[col_name]])
    } else {
      warning(paste("Covariate coefficient column '", col_name, "' not found in Excel file.", sep=""))
      external_coefs[col_name] <- NA 
      missing_cov_cols <- c(missing_cov_cols, col_name)
    }
  }
  if(any(is.na(external_coefs))) {
    stop("Could not extract all required covariate coefficients: ", paste(missing_cov_cols, collapse=", "))
  }
  
  return(list(
    p = p, q = q, sigmasq = sigmasq, intercept = intercept,
    betas = betas, alphas = alphas, external_coefs = external_coefs
  ))
}

run_simulation_study_with_covariates_parallel <- function() {
  set.seed(54321)
  
  num_simulations <- 1000
  sim_length <- 1000     
  progress_print_frequency <- max(1, floor(num_simulations / 10))
  
  task_max.p <- 7 
  task_max.q <- 7 
  task_max.order_stepwise <- 5  
  task_max.order_grid <- 14     
  task_distribution <- "nbinom"
  task_link <- "log"
  task_ic <- "aic" 
  
  cat("Starting INGARCH simulation study (WITH covariates).\n")
  cat("Configuration: num_simulations =", num_simulations, ", sim_length =", sim_length, "\n")
  cat("  task_max.p =", task_max.p, ", task_max.q =", task_max.q, "\n")
  cat("  task_max.order_stepwise =", task_max.order_stepwise, ", task_max.order_grid =", task_max.order_grid, "\n")
  cat("  distribution =", task_distribution, ", link =", task_link, ", ic =", task_ic, "\n")
  
  cat("Loading covariate data...\n")
  covariate_data_path <- "./data/count_covariates_data.RData"
  if (!file.exists(covariate_data_path)) stop("Covariate data file not found: ", covariate_data_path)
  load(covariate_data_path) 
  if (!exists("d_1_data")) stop("Data frame 'd_1_data' not found in RData file.")
  
  cov_names_r <- c("mean.Temp.Lag5", "mean.DTemp.Lag0", "PM25.Lag0", "mean.NO2.Lag2")
  if (!all(cov_names_r %in% colnames(d_1_data))) stop("Missing required covariate columns in d_1_data.")
  
  if (nrow(d_1_data) < sim_length) stop(paste("Not enough covariate data rows. Need", sim_length, "have", nrow(d_1_data)))
  xreg_matrix <- as.matrix(d_1_data[1:sim_length, cov_names_r, drop = FALSE])
  if(anyNA(xreg_matrix)) warning("NA values found in xreg_matrix.", immediate. = TRUE)
  
  cat("Extracting true model parameters from Excel file...\n")
  excel_cov_names <- c("Temp", "DTemp", "PM", "NO2") 
  true_params <- extract_model_params_with_covariates("./data/modelosAveiro.xlsx", model_row = 2, covariate_cols = excel_cov_names)
  
  cat("True Model parameters from Excel:\n")
  cat("  p =", true_params$p, ", q =", true_params$q, ", sigmasq =", true_params$sigmasq, ", intercept =", true_params$intercept, "\n")
  cat("  Beta coefficients:", if(is.null(true_params$betas)) "None" else paste(round(true_params$betas,4), collapse=" "), "\n")
  cat("  Alpha coefficients:", if(is.null(true_params$alphas)) "None" else paste(round(true_params$alphas,4), collapse=" "), "\n")
  cat("  External coefficients (Excel names -> R names for xreg_matrix):\n")
  
  sim_external_coefs_for_ingarch.sim <- true_params$external_coefs 
  names(sim_external_coefs_for_ingarch.sim) <- cov_names_r 
  
  for(i in 1:length(cov_names_r)){
    cat("    ", excel_cov_names[i], " (simulated as ", cov_names_r[i], ") = ", round(true_params$external_coefs[excel_cov_names[i]], 4), "\n")
  }
  
  ingarch_params_for_sim <- list(
    intercept = true_params$intercept, 
    past_obs = true_params$betas,
    past_mean = true_params$alphas, 
    external = sim_external_coefs_for_ingarch.sim 
  )
  ingarch_model_for_sim <- list(
    past_obs = if(true_params$p > 0) 1:true_params$p else NULL,
    past_mean = if(true_params$q > 0) 1:true_params$q else NULL,
    external = TRUE 
  )
  if (is.na(true_params$sigmasq) || true_params$sigmasq <= 0) stop("Invalid sigmasq value.")
  size_param_for_sim <- 1/true_params$sigmasq
  cat("  Calculated 'size' parameter for nbinom simulation (1/sigmasq):", size_param_for_sim, "\n")
  
  cat(paste("Simulating", num_simulations, "INGARCH realizations of length", sim_length, "...\n"))
  sims <- vector("list", num_simulations)
  for(i in 1:num_simulations) {
    if(i == 1 || i == num_simulations || (i %% progress_print_frequency == 0 && num_simulations > 10)) { # Text progress
      cat("  Generating simulation", i, "of", num_simulations, "...\n")
    }
    sim_result <- tryCatch(
      ingarch.sim(n = sim_length, param = ingarch_params_for_sim, model = ingarch_model_for_sim,
                  xreg = xreg_matrix, link = task_link, distr = task_distribution,
                  size = size_param_for_sim, n_start = 100),
      error = function(e) { cat("Error sim", i, ":", conditionMessage(e), "\n"); NULL }
    )
    if(is.null(sim_result) || !("ts" %in% names(sim_result))) sims[[i]] <- NULL else sims[[i]] <- as.numeric(sim_result$ts)
  }
  cat("Finished generating simulations.\n")
  
  cat("Running PARALLEL model selection for INGARCH with covariates...\n")
  
  results_template <<- define_results_template_with_cov(task_max.p, task_max.q, cov_names_r) 
  target_column_names <<- c("method", "sim_id", names(results_template))
  
  num_cores_detected <- detectCores(logical = FALSE) # Physical cores
  cores_to_use <- max(1, num_cores_detected - 2) # Leave 2 cores free or use 1 if only 1-2 available
  if(num_simulations < cores_to_use) cores_to_use <- num_simulations # Don't use more cores than tasks
  
  cl <- makeCluster(cores_to_use)
  registerDoParallel(cl)
  cat(paste("Registered parallel backend with", getDoParWorkers(), "cores.\n"))
  
  # Removed progressr handlers
  
  custom_funcs_to_export <- c("auto.ingarch", "ingarch.sim", "newmodel",
                              "ingarch.string", "search.ingarch", "myingarch", 
                              "populate_results_from_fit_with_cov", "define_results_template_with_cov",
                              "results_template", "target_column_names") 
  custom_funcs_to_export <- custom_funcs_to_export[sapply(custom_funcs_to_export, exists, envir = .GlobalEnv, inherits = FALSE)]
  
  # Main parallel loop WITHOUT with_progress
  parallel_results_list <- foreach(
    i = 1:length(sims),
    .packages = c("tscount", "stats", "dplyr"), 
    .export = c("xreg_matrix", "task_max.p", "task_max.q", "task_max.order_stepwise", "task_max.order_grid", 
                "task_distribution", "task_link", "task_ic", "cov_names_r", custom_funcs_to_export),
    .errorhandling = 'pass' 
  ) %dopar% {
    
    # Optional: Print a message from worker - this might be messy
    # if (i %% (length(sims)/10) == 0 || i == 1 || i == length(sims)) { # Rough progress from workers
    #   cat(sprintf("Worker %d processing task %d\n", Sys.getpid(), i))
    # }
    
    sim_data <- sims[[i]]
    worker_results_template <- define_results_template_with_cov(task_max.p, task_max.q, cov_names_r)
    
    current_stepwise_template <- worker_results_template 
    current_grid_template     <- worker_results_template
    
    current_stepwise_template$time <- 0 
    current_grid_template$time   <- 0   
    
    if (is.null(sim_data)) {
      current_stepwise_template$status <- "Input Sim Failed"
      current_stepwise_template$error_message <- "Simulated data was NULL"
      current_grid_template$status <- "Input Sim Failed"
      current_grid_template$error_message <- "Simulated data was NULL"
      return(list(stepwise = current_stepwise_template, grid_search = current_grid_template))
    }
    if (length(sim_data) != nrow(xreg_matrix)) {
      current_stepwise_template$status <- "Input Length Mismatch"
      current_stepwise_template$error_message <- "Sim_data length != xreg_matrix rows"
      current_grid_template$status <- "Input Length Mismatch"
      current_grid_template$error_message <- "Sim_data length != xreg_matrix rows"
      return(list(stepwise = current_stepwise_template, grid_search = current_grid_template))
    }
    
    start_time_stepwise <- Sys.time()
    fit_stepwise <- tryCatch({
      model <- auto.ingarch(y = sim_data, xreg = xreg_matrix, 
                            max.p = task_max.p, max.q = task_max.q, max.order = task_max.order_stepwise,
                            distribution = task_distribution, link = task_link, ic = task_ic,
                            stepwise = TRUE, trace = FALSE, show_warnings = FALSE)
      model 
    }, error = function(e) {
      current_stepwise_template$status <- paste("Error") 
      current_stepwise_template$error_message <- conditionMessage(e)
      return(list(error = TRUE, message = conditionMessage(e))) 
    })
    end_time_stepwise <- Sys.time()
    current_stepwise_template$time <- as.numeric(difftime(end_time_stepwise, start_time_stepwise, units = "secs"))
    
    iter_stepwise_result <- populate_results_from_fit_with_cov(
      fit_stepwise, current_stepwise_template, task_max.p, task_max.q, cov_names_r)
    iter_stepwise_result$time <- current_stepwise_template$time 
    if (inherits(fit_stepwise, "error") || is.null(fit_stepwise) || !inherits(fit_stepwise, "tsglm")) {
      iter_stepwise_result$status <- current_stepwise_template$status 
      iter_stepwise_result$error_message <- current_stepwise_template$error_message
    }
    
    start_time_grid <- Sys.time()
    fit_grid <- tryCatch({
      model <- auto.ingarch(y = sim_data, xreg = xreg_matrix, 
                            max.p = task_max.p, max.q = task_max.q, max.order = task_max.order_grid,
                            distribution = task_distribution, link = task_link, ic = task_ic,
                            stepwise = FALSE, trace = FALSE, show_warnings = FALSE)
      model 
    }, error = function(e) {
      current_grid_template$status <- paste("Error") 
      current_grid_template$error_message <- conditionMessage(e)
      return(list(error = TRUE, message = conditionMessage(e))) 
    })
    end_time_grid <- Sys.time()
    current_grid_template$time <- as.numeric(difftime(end_time_grid, start_time_grid, units = "secs"))
    
    iter_grid_result <- populate_results_from_fit_with_cov(
      fit_grid, current_grid_template, task_max.p, task_max.q, cov_names_r)
    iter_grid_result$time <- current_grid_template$time 
    if (inherits(fit_grid, "error") || is.null(fit_grid) || !inherits(fit_grid, "tsglm")) {
      iter_grid_result$status <- current_grid_template$status 
      iter_grid_result$error_message <- current_grid_template$error_message
    }
    if(iter_grid_result$status == "Success" && (is.na(iter_grid_result$n_models_tested) || iter_grid_result$n_models_tested == 0)) { 
      iter_grid_result$n_models_tested <- (task_max.p + 1) * (task_max.q + 1)
    }
    
    # Removed p_progress call
    return(list(stepwise = iter_stepwise_result, grid_search = iter_grid_result))
  } 
  
  stopCluster(cl)
  cat("Parallel cluster stopped.\n")
  
  results <- list(stepwise = vector("list", length(sims)), grid_search = vector("list", length(sims)))
  for (i in 1:length(parallel_results_list)) {
    if (inherits(parallel_results_list[[i]], "error")) {
      error_msg <- conditionMessage(parallel_results_list[[i]])
      error_entry <- results_template 
      error_entry$status <- paste("Outer Error:", substr(error_msg,1,100))
      error_entry$error_message <- substr(error_msg,1,250)
      results$stepwise[[i]] <- error_entry
      results$grid_search[[i]] <- error_entry
    } else if (is.list(parallel_results_list[[i]]) && !is.null(parallel_results_list[[i]]$stepwise) && !is.null(parallel_results_list[[i]]$grid_search)) {
      results$stepwise[[i]] <- parallel_results_list[[i]]$stepwise
      results$grid_search[[i]] <- parallel_results_list[[i]]$grid_search
    } else {
      default_fail_entry <- results_template
      default_fail_entry$status <- "Unexpected Result Structure"
      default_fail_entry$error_message <- "Parallel result was not an error, but not the expected list."
      results$stepwise[[i]] <- default_fail_entry
      results$grid_search[[i]] <- default_fail_entry
    }
  }
  cat("Results restructured.\n")
  
  results_to_df <- function(list_of_result_lists_input, method_name_str) {
    if (length(list_of_result_lists_input) == 0 || all(sapply(list_of_result_lists_input, is.null))) {
      df_cols <- lapply(target_column_names, function(col_name) { 
        template_col_class <- class(results_template[[col_name]]) 
        if(length(template_col_class) > 1) template_col_class <- template_col_class[1]
        if (template_col_class == "integer") return(integer(0))
        if (template_col_class == "numeric") return(numeric(0))
        return(character(0)) 
      })
      names(df_cols) <- target_column_names
      return(as.data.frame(df_cols, stringsAsFactors = FALSE))
    }
    df_list_for_rbind <- lapply(1:length(list_of_result_lists_input), function(i) {
      single_sim_result_list <- list_of_result_lists_input[[i]]
      row_data <- results_template 
      if(!is.null(single_sim_result_list) && is.list(single_sim_result_list)){
        for (col_name in names(single_sim_result_list)) {
          if (col_name %in% names(row_data)) {
            if(length(single_sim_result_list[[col_name]]) == 1){
              row_data[[col_name]] <- single_sim_result_list[[col_name]]
            } else if (length(single_sim_result_list[[col_name]]) == 0 && (is.numeric(row_data[[col_name]]) || is.integer(row_data[[col_name]]))) {
              row_data[[col_name]] <- NA 
            } else if (length(single_sim_result_list[[col_name]]) == 0 && is.character(row_data[[col_name]])) {
              row_data[[col_name]] <- NA_character_
            } else {
              if (is.atomic(single_sim_result_list[[col_name]])) {
                row_data[[col_name]] <- single_sim_result_list[[col_name]][1] 
              } else { row_data[[col_name]] <- NA }
            }
          }
        }
      }
      final_row_list <- c(list(method = method_name_str, sim_id = i), row_data)
      final_df_row <- data.frame(matrix(ncol = length(target_column_names), nrow = 1))
      colnames(final_df_row) <- target_column_names
      for(col_name_target in target_column_names){
        if(col_name_target %in% names(final_row_list)){
          val_to_assign <- final_row_list[[col_name_target]]
          if(is.null(val_to_assign)){
            template_type <- class(results_template[[col_name_target]])
            if(length(template_type)>1) template_type <- template_type[1]
            if(template_type == "integer") val_to_assign <- NA_integer_
            else if(template_type == "numeric") val_to_assign <- NA_real_
            else val_to_assign <- NA_character_
          }
          final_df_row[1, col_name_target] <- val_to_assign
        } else {
          template_type <- class(results_template[[col_name_target]]) 
          if(length(template_type)>1) template_type <- template_type[1]
          if(template_type == "integer") final_df_row[1, col_name_target] <- NA_integer_
          else if(template_type == "numeric") final_df_row[1, col_name_target] <- NA_real_
          else final_df_row[1, col_name_target] <- NA_character_
        }
      }
      return(final_df_row)
    })
    if (length(df_list_for_rbind) == 0) { 
      df_cols <- lapply(target_column_names, function(col_name) {
        template_col_class <- class(results_template[[col_name]])
        if(length(template_col_class) > 1) template_col_class <- template_col_class[1]
        if (template_col_class == "integer") return(integer(0))
        if (template_col_class == "numeric") return(numeric(0))
        return(character(0)) 
      })
      names(df_cols) <- target_column_names
      return(as.data.frame(df_cols, stringsAsFactors = FALSE))
    }
    return(dplyr::bind_rows(df_list_for_rbind))
  }
  
  stepwise_results_df <- results_to_df(results$stepwise, "stepwise")
  grid_results_df     <- results_to_df(results$grid_search, "grid_search")
  results_df          <- dplyr::bind_rows(stepwise_results_df, grid_results_df)
  
  output_dir <- "./simulation_output_parallel" # Reverted to original for this specific output naming
  if (!dir.exists(output_dir)) { dir.create(output_dir, recursive=TRUE) } 
  sims_filename <- file.path(output_dir, "ingarch_with_covariates_simulations.rds")
  results_csv_filename <- file.path(output_dir, "ingarch_with_covariates_results_parallel.csv") # Original was _parallel.rds
  summary_filename <- file.path(output_dir, "ingarch_with_covariates_summary_parallel.csv")
  order_freq_filename <- file.path(output_dir, "ingarch_with_covariates_order_freq_parallel.csv")
  
  saveRDS(sims, sims_filename)
  write.csv(results_df, results_csv_filename, row.names = FALSE, na = "") 
  cat("Detailed results saved to:", results_csv_filename, "\n")
  
  results_df$status <- as.character(results_df$status) 
  
  numeric_cols_for_summary <- c("p_order", "q_order", "time", "n_models_tested", "aic", "bic", "intercept", "sigmasq", 
                                paste0("beta", 1:task_max.p), paste0("alpha", 1:task_max.q), 
                                make.names(cov_names_r, unique = TRUE))
  for(col_s in numeric_cols_for_summary){
    if(col_s %in% names(results_df)){
      results_df[[col_s]] <- as.numeric(as.character(results_df[[col_s]])) 
    }
  }
  
  summary_stats <- results_df %>%
    filter(!is.na(method)) %>% 
    group_by(method) %>%
    summarize(
      total_sims = n(),
      successful_fits = sum(status == "Success", na.rm = TRUE),
      mean_p = mean(p_order[status == "Success"], na.rm = TRUE),
      median_p = median(p_order[status == "Success"], na.rm = TRUE),
      sd_p = sd(p_order[status == "Success"], na.rm = TRUE),
      mean_q = mean(q_order[status == "Success"], na.rm = TRUE),
      median_q = median(q_order[status == "Success"], na.rm = TRUE),
      sd_q = sd(q_order[status == "Success"], na.rm = TRUE),
      mean_time_secs = mean(time, na.rm = TRUE),
      median_time_secs = median(time, na.rm = TRUE),
      sd_time_secs = sd(time, na.rm = TRUE),
      mean_models_tested = mean(n_models_tested[status == "Success"], na.rm = TRUE),
      median_models_tested = median(n_models_tested[status == "Success"], na.rm = TRUE),
      sd_models_tested = sd(n_models_tested[status == "Success"], na.rm = TRUE),
      failure_count = sum(status != "Success" & !(is.na(status) | status == ""), na.rm = TRUE),
      mean_intercept = mean(intercept[status == "Success"], na.rm = TRUE),
      mean_sigmasq = mean(sigmasq[status == "Success"], na.rm = TRUE),
      .groups = 'drop' 
    )
  
  order_freq <- results_df %>%
    filter(status == "Success") %>% 
    group_by(method, p_order, q_order) %>%
    summarize(count = n(), .groups = 'drop') %>%
    group_by(method) %>%
    mutate(freq = count / sum(count) * 100) %>%
    arrange(method, desc(freq))
  
  write.csv(summary_stats, summary_filename, row.names = FALSE, na = "")
  write.csv(order_freq, order_freq_filename, row.names = FALSE, na = "")
  
  cat("\n===== Summary Statistics (Parallel Run with Detailed Output) =====\n"); print(summary_stats)
  cat("\n===== Most Frequent (p,q) Orders (Parallel Run) =====\n")
  if(nrow(order_freq) > 0) print(order_freq %>% group_by(method) %>% slice_head(n = 5))
  else cat("No successful fits to determine frequent orders.\n")
  
  cat("\n===== Failure/Status Summary (Parallel Run) =====\n")
  print(results_df %>% group_by(method, status) %>% summarise(count = n(), .groups = 'drop') %>% arrange(method, desc(count)))
  
  cat(paste0("\nResults saved to directory: ", output_dir, "\n"))
  # Removed explicit "Please run warnings()" as it might not be relevant if no warnings occurred.
  # User can run it if they suspect issues.
  return(results_df)
}

main_parallel <- function() {
  if(exists(".coef_names_printed_flag_v4", envir = .GlobalEnv)){
    rm(".coef_names_printed_flag_v4", envir = .GlobalEnv)
  }
  
  overall_start_time <- Sys.time()
  results_data_frame <- tryCatch(run_simulation_study_with_covariates_parallel(),
                                 error = function(e) {
                                   cat("\nERROR during simulation study run:\n")
                                   print(e)
                                   # Clean up parallel backend on error
                                   if(exists("cl") && inherits(cl, "cluster")) try(stopCluster(cl), silent=TRUE)
                                   return(NULL)
                                 })
  overall_end_time <- Sys.time()
  overall_duration <- overall_end_time - overall_start_time
  
  if(!is.null(results_data_frame)){
    cat("\nParallel simulation study with detailed output completed!\n")
  } else {
    cat("\nParallel simulation study encountered an error and did not complete successfully.\n")
  }
  cat("Total duration of the study:", format(overall_duration), "\n")
}

main_parallel()