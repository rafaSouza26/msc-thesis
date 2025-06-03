# RerunStepwiseWithoutCovariates.R

# --- 0. Load required packages ---
library(tscount)
library(dplyr)
library(doParallel)
library(foreach)

# --- 1. Source custom functions ---
# Functions from ACTS package
source("./ACTS/auto.ingarch.R")
source("./ACTS/myingarch.R") # Called by auto.ingarch
source("./ACTS/newmodel.R") # Dependency for auto.ingarch
source("./ACTS/ingarch.string.R") # Dependency for auto.ingarch
source("./ACTS/search.ingarch.R") # Dependency for auto.ingarch

# Helper functions 

# Define a template for storing detailed results (no covariates)

define_results_template_no_cov <- function(max_p_cols, max_q_cols) {
  na_betas_list <- list(); if(max_p_cols > 0) { na_betas_list <- as.list(rep(NA_real_, max_p_cols)); names(na_betas_list) <- paste0("beta", 1:max_p_cols) }
  na_alphas_list <- list(); if(max_q_cols > 0) { na_alphas_list <- as.list(rep(NA_real_, max_q_cols)); names(na_alphas_list) <- paste0("alpha", 1:max_q_cols) }
  template <- c(list(p_order=NA_integer_, q_order=NA_integer_, time=NA_real_, n_models_tested=NA_integer_, aic=NA_real_, bic=NA_real_, intercept=NA_real_, sigmasq=NA_real_), na_betas_list, na_alphas_list, list(status="Not Run", error_message=""))
  
  for (name in c("p_order", "q_order", "n_models_tested")) template[[name]] <- NA_integer_
  for (name in c("time", "aic", "bic", "intercept", "sigmasq")) template[[name]] <- NA_real_
  if(max_p_cols > 0) for(k in 1:max_p_cols) template[[paste0("beta",k)]] <- NA_real_
  if(max_q_cols > 0) for(k in 1:max_q_cols) template[[paste0("alpha",k)]] <- NA_real_
  template$status <- "Not Run"
  template$error_message <- ""
  return(template)
}

# Populate the results template from a fitted model (no covariates)

populate_results_from_fit_no_cov <- function(fit_object, template_data,
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
      
      if (max_p_cols_in_template > 0 && current_p_order > 0 && !is.null(fit_object$model$past_obs)) {
        for (lag_val in fit_object$model$past_obs) { 
          if (lag_val > max_p_cols_in_template) next 
          template_beta_col <- paste0("beta", lag_val)
          model_coef_name_past_obs <- paste0("past_obs", lag_val) 
          model_coef_name_beta <- paste0("beta_", lag_val) 
          
          if (model_coef_name_past_obs %in% names(coeffs_vector)) {
            populated_data[[template_beta_col]] <- as.numeric(coeffs_vector[model_coef_name_past_obs])[1]
          } else if (model_coef_name_beta %in% names(coeffs_vector)) {
            populated_data[[template_beta_col]] <- as.numeric(coeffs_vector[model_coef_name_beta])[1]
          }
        }
      }
      if (max_q_cols_in_template > 0 && current_q_order > 0 && !is.null(fit_object$model$past_mean)) {
        for (lag_val in fit_object$model$past_mean) { 
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
    } 
    
    if (!is.null(fit_object$distr) && fit_object$distr == "nbinom") {
      if(!is.null(fit_object$sigmasq)) { 
        populated_data$sigmasq <- as.numeric(fit_object$sigmasq)[1]
      } else if (!is.null(model_summary_obj) && !is.null(model_summary_obj$dispersion) && "estimate" %in% names(model_summary_obj$dispersion)){
        alpha_disp <- as.numeric(model_summary_obj$dispersion["estimate"])[1] 
        if(!is.na(alpha_disp)) populated_data$sigmasq <- alpha_disp 
      } else if (!is.null(fit_object$distrcoefs)) { 
        if (is.numeric(fit_object$distrcoefs) && fit_object$distrcoefs > 0) { 
          populated_data$sigmasq <- 1 / as.numeric(fit_object$distrcoefs)[1] 
        } else if (is.list(fit_object$distrcoefs) && "size" %in% names(fit_object$distrcoefs) && is.numeric(fit_object$distrcoefs$size) && fit_object$distrcoefs$size > 0) {
          populated_data$sigmasq <- 1 / as.numeric(fit_object$distrcoefs$size)[1] 
        }
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


# --- 2. Configuration ---
set.seed(12345) # For reproducibility

# Set the number of simulations to process.
# Set to NULL or a non-positive number to process all available simulations.
num_simulations_to_process <- 50 # Example: Process only the first 100 simulations

# Paths
simulations_rds_path <- "./Results/Simulation/no_covariates/ingarch_no_covariates_simulations.rds" 
output_dir <- "./rerun_stepwise_without_covariates"

# auto.ingarch parameters 
task_max.p <- 7 
task_max.q <- 7 
task_max.order_stepwise <- 5 
task_distribution <- "nbinom"
task_link <- "log"
task_ic <- "aic"

cat("Starting INGARCH model re-run (stepwise only, NO covariates).\n")
if (!is.null(num_simulations_to_process) && is.numeric(num_simulations_to_process) && num_simulations_to_process > 0) {
  cat("User has requested to process up to", num_simulations_to_process, "simulations.\n")
} else {
  cat("User has requested to process all available simulations.\n")
}

# --- 3. Load Simulations ---
if (!file.exists(simulations_rds_path)) {
  stop("Simulations RDS file not found: ", simulations_rds_path, 
       "\nPlease ensure the path is correct and the RDS file with 'no covariates' simulations exists.")
}
cat("Loading simulations from:", simulations_rds_path, "\n")
sims_list_all <- readRDS(simulations_rds_path) 

if (!is.list(sims_list_all) || length(sims_list_all) == 0) {
  stop("Loaded simulations object is not a list or is empty.")
}
total_sims_available <- length(sims_list_all)
cat(total_sims_available, "total simulations available in the RDS file.\n")

actual_num_to_process <- total_sims_available
if (!is.null(num_simulations_to_process) && is.numeric(num_simulations_to_process) && num_simulations_to_process > 0) {
  if (num_simulations_to_process < total_sims_available) {
    actual_num_to_process <- num_simulations_to_process
    cat("Will process the first", actual_num_to_process, "simulations.\n")
  } else {
    cat("Requested number (", num_simulations_to_process, ") is >= total available. Processing all", total_sims_available, "simulations.\n")
  }
} else {
  cat("Processing all", total_sims_available, "available simulations.\n")
}

sims_list_to_process <- sims_list_all[1:actual_num_to_process]
sim_length <- length(sims_list_to_process[[1]]) # Determine sim_length from first loaded sim
cat("Determined simulation length from first simulation to be processed:", sim_length, "\n")

# --- 4. Setup for Parallel Processing ---
results_template <<- define_results_template_no_cov(task_max.p, task_max.q)
target_column_names <- c("sim_id", names(results_template)) 

num_cores_detected <- detectCores(logical = FALSE) 
cores_to_use <- max(1, num_cores_detected - 1) 
if(actual_num_to_process < cores_to_use) cores_to_use <- actual_num_to_process

cl <- makeCluster(cores_to_use)
registerDoParallel(cl)
cat(paste("Registered parallel backend with", getDoParWorkers(), "cores for processing", actual_num_to_process, "simulations.\n"))

custom_funcs_to_export <- c("auto.ingarch", "myingarch", "newmodel", "ingarch.string", "search.ingarch",
                            "populate_results_from_fit_no_cov", "define_results_template_no_cov")
custom_funcs_to_export <- custom_funcs_to_export[sapply(custom_funcs_to_export, exists, envir = .GlobalEnv, inherits = FALSE)]

vars_to_export <- c("task_max.p", "task_max.q", "task_max.order_stepwise",
                    "task_distribution", "task_link", "task_ic")

cat("Running PARALLEL model selection (stepwise only, NO covariates)...\n")

# --- 5. Parallel Loop for Model Fitting ---
parallel_results_list <- foreach(
  i = 1:actual_num_to_process,  
  .packages = c("tscount", "stats", "dplyr"), 
  .export = c(custom_funcs_to_export, vars_to_export), 
  .errorhandling = 'pass' 
) %dopar% {
  
  sim_data_ts <- ts(sims_list_to_process[[i]])
  worker_iter_template <- define_results_template_no_cov(task_max.p, task_max.q)
  
  worker_iter_template$time <- 0 
  
  if (is.null(sim_data_ts)) {
    worker_iter_template$status <- "Input Sim Failed"
    worker_iter_template$error_message <- "Simulated data was NULL"
    return(worker_iter_template)
  }
  
  start_time_stepwise <- Sys.time()
  fit_stepwise <- tryCatch({
    model <- auto.ingarch(y = sim_data_ts,
                          max.p = task_max.p, max.q = task_max.q, max.order = task_max.order_stepwise,
                          distribution = task_distribution, link = task_link, ic = task_ic,
                          stepwise = TRUE, trace = FALSE, show_warnings = FALSE, parallel = FALSE)
    model 
  }, error = function(e) {
    return(list(error = TRUE, message = conditionMessage(e), status_message = "Stepwise Fit Error")) 
  })
  end_time_stepwise <- Sys.time()
  worker_iter_template$time <- as.numeric(difftime(end_time_stepwise, start_time_stepwise, units = "secs"))
  
  iter_stepwise_result <- populate_results_from_fit_no_cov(
    fit_stepwise, worker_iter_template, task_max.p, task_max.q)
  
  iter_stepwise_result$time <- worker_iter_template$time 
  
  if (inherits(fit_stepwise, "error") || (is.list(fit_stepwise) && !is.null(fit_stepwise$error) && fit_stepwise$error)) {
    iter_stepwise_result$status <- fit_stepwise$status_message # Use status from tryCatch
    iter_stepwise_result$error_message <- if (is.list(fit_stepwise)) fit_stepwise$message else conditionMessage(fit_stepwise)
  }
  
  return(iter_stepwise_result)
} 

# --- 6. Stop Parallel Cluster ---
stopCluster(cl)
cat("Parallel cluster stopped.\n")

# --- 7. Process and Save Results ---
cat("Processing results...\n")

for(i in 1:length(parallel_results_list)){
  if(is.list(parallel_results_list[[i]])){
    parallel_results_list[[i]]$sim_id <- i
  } else { 
    temp_error_entry <- results_template 
    temp_error_entry$sim_id <- i
    temp_error_entry$status <- "Foreach Loop Error"
    temp_error_entry$error_message <- as.character(parallel_results_list[[i]])
    parallel_results_list[[i]] <- temp_error_entry
  }
}

results_df <- dplyr::bind_rows(parallel_results_list)

final_column_order <- target_column_names 
for (col_name in final_column_order) {
  if (!col_name %in% names(results_df)) {
    results_df[[col_name]] <- NA
  }
}
results_df <- results_df[, final_column_order, drop = FALSE]


if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

results_csv_filename <- file.path(output_dir, "rerun_stepwise_no_covariates_results.csv")
summary_filename <- file.path(output_dir, "rerun_stepwise_no_covariates_summary.csv")
order_freq_filename <- file.path(output_dir, "rerun_stepwise_no_covariates_order_freq.csv")

write.csv(results_df, results_csv_filename, row.names = FALSE, na = "") 
cat("Detailed results saved to:", results_csv_filename, "\n")

results_df$status <- as.character(results_df$status)
# Columns for numeric conversion (no covariate coefficient columns)
numeric_cols_for_summary <- c("p_order", "q_order", "time", "n_models_tested", "aic", "bic", "intercept", "sigmasq", 
                              paste0("beta", 1:task_max.p), paste0("alpha", 1:task_max.q))
numeric_cols_for_summary <- intersect(numeric_cols_for_summary, names(results_df)) # Ensure columns exist

for(col_s in numeric_cols_for_summary){
  if(col_s %in% names(results_df)){
    results_df[[col_s]] <- as.numeric(as.character(results_df[[col_s]])) 
  }
}

summary_stats <- results_df %>%
  summarize(
    total_sims_processed = n(),
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
    mean_sigmasq = mean(sigmasq[status == "Success"], na.rm = TRUE), # This is 1/size or dispersion alpha
    .groups = 'drop' 
  )

order_freq <- results_df %>%
  filter(status == "Success") %>% 
  group_by(p_order, q_order) %>%
  summarize(count = n(), .groups = 'drop') %>%
  mutate(freq = count / sum(count) * 100) %>%
  arrange(desc(freq))

write.csv(summary_stats, summary_filename, row.names = FALSE, na = "")
cat("Summary statistics saved to:", summary_filename, "\n")
write.csv(order_freq, order_freq_filename, row.names = FALSE, na = "")
cat("Order frequencies saved to:", order_freq_filename, "\n")

# --- 8. Print Summary to Console ---
cat("\n===== Summary Statistics (Stepwise Re-run, NO Covariates) =====\n"); print(summary_stats)
cat("\n===== Most Frequent (p,q) Orders (Stepwise Re-run, NO Covariates) =====\n")
if(nrow(order_freq) > 0) {
  print(order_freq %>% slice_head(n = 10))
} else {
  cat("No successful fits to determine frequent orders.\n")
}

cat("\n===== Failure/Status Summary (Stepwise Re-run, NO Covariates) =====\n")
print(results_df %>% group_by(status) %>% summarise(count = n(), .groups = 'drop') %>% arrange(desc(count)))

cat(paste0("\nScript completed. Results saved to directory: ", output_dir, "\n"))