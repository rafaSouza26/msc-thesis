# --------------------------------------------------------------------------
# R SCRIPT: Multiple auto.ingarch runs PER INPUT TIME SERIES by varying
# max.p/q settings ONLY. Uses fixed AIC, Nbinom dist, Log link, and
# default start.p/q for auto.ingarch.
# Simplified numeric Run IDs. Results saved for Python.
# --------------------------------------------------------------------------
cat("Starting R script: Multiple runs/sim, fixed AIC/Nbinom/Log, varying Max P/Q...\n")
script_start_time <- Sys.time()

# STEP 0: Load custom functions and 'jsonlite'
# --------------------------------------------------------------------------
cat("Step 0: Loading custom functions and 'jsonlite' package...\n")
tryCatch({
  source("./custom/newmodel.R")
  source("./custom/ingarch.string.R")
  source("./custom/search.ingarch.R")
  cat("Custom functions loaded successfully.\n")
}, error = function(e) {
  stop("Error loading custom functions: ", e$message)
})

if (!require(jsonlite)) {
  install.packages("jsonlite", repos = "http://cran.us.r-project.org")
  library(jsonlite)
} else {
  library(jsonlite)
}
cat("jsonlite package is ready.\n")

# STEP 1: Load all simulated time series
# --------------------------------------------------------------------------
cat("Step 1: Loading all simulated time series...\n")
simulations_filepath <- "C:/Users/Rafael/Desktop/msc-thesis/results/simulatedDataResults/ingarch_no_covariates_simulations.rds" # ADJUST THIS PATH

if (!file.exists(simulations_filepath)) {
  stop("Simulations file not found: ", simulations_filepath)
}
simulations_list_all <- readRDS(simulations_filepath)
if (!is.list(simulations_list_all) || length(simulations_list_all) == 0) {
  stop("Loaded data is not a list or is empty.")
}
cat("Loaded", length(simulations_list_all), "simulated time series.\n")

# STEP 2: Define parameter combinations for multiple runs
# --------------------------------------------------------------------------
cat("Step 2: Defining parameter combinations for multiple runs...\n")

# --- USER CONFIGURATION ---
num_original_simulations_to_process <- 1 # How many of your 1000 simulations to use as input
# --------------------------

if (num_original_simulations_to_process > length(simulations_list_all)) {
  num_original_simulations_to_process <- length(simulations_list_all)
  cat("Warning: num_original_simulations_to_process adjusted to available:", num_original_simulations_to_process, "\n")
}
cat("Will process the first", num_original_simulations_to_process, "simulations from the loaded list.\n")

# Fixed parameters for all runs
fixed_ic         <- "aic"
fixed_distribution <- "nbinom"
fixed_link         <- "log"
cat("Using fixed IC:", fixed_ic, ", Distribution:", fixed_distribution,", Link:", fixed_link, "\n")

# Max P/Q settings to vary across runs
max_pq_settings_to_try <- list(
  c(max_p=5, max_q=5),
  c(max_p=7, max_q=7),
  c(max_p=7, max_q=4),
  c(max_p=4, max_q=7)
)

total_runs_planned_per_simulation <- length(max_pq_settings_to_try)
total_overall_runs_planned <- num_original_simulations_to_process * total_runs_planned_per_simulation

cat("Runs planned per input simulation:", total_runs_planned_per_simulation, "(varying max.p/q settings)\n")
cat("Total overall auto.ingarch runs planned:", total_overall_runs_planned, "\n")
if (total_overall_runs_planned == 0) {
  stop("No runs planned. Check max_pq_settings_to_try.")
}

# Base parameters for auto.ingarch call (some will be overridden in loops)
auto_ingarch_call_params <- list(
  # 'y', 'max.p', 'max.q' will be set in loops
  ic = fixed_ic,
  distribution = fixed_distribution,
  link = fixed_link,
  stepwise = TRUE,
  trace = FALSE, # Set TRUE to debug individual runs
  show_warnings = FALSE
  # start.p and start.q are NOT specified, so auto.ingarch uses its defaults (2,2)
)

# STEP 3: Loop through selected simulations and Max P/Q settings
# --------------------------------------------------------------------------
cat("Step 3: Starting auto.ingarch runs...\n")
all_run_results <- list() 
overall_run_counter <- 0 

for (sim_idx in 1:num_original_simulations_to_process) {
  current_ts_data <- as.ts(simulations_list_all[[sim_idx]])
  cat("\n--- Processing Original Simulation Index:", sim_idx, "---\n")
  
  for (mpq_setting in max_pq_settings_to_try) {
    overall_run_counter <- overall_run_counter + 1
    current_run_id_numeric_string <- as.character(overall_run_counter)
    
    current_max_p <- mpq_setting[1]
    current_max_q <- mpq_setting[2]
    
    descriptive_params <- paste0("SimIdx=", sim_idx, ", MaxP=", current_max_p, ", MaxQ=", current_max_q, 
                                 " (AIC, Nbinom, Log, Default StartP/Q)")
    cat("  Run ID:", current_run_id_numeric_string, "(", overall_run_counter, "/", total_overall_runs_planned, ") Params: ", descriptive_params, "\n", sep="")
    
    run_params <- auto_ingarch_call_params
    run_params$y <- current_ts_data
    run_params$max.p <- current_max_p
    run_params$max.q <- current_max_q
    
    ingarch_run_result_obj <- NULL; error_message_for_run <- NULL
    stepwise_path_matrix_run <- NULL; final_p_run <- NA; final_q_run <- NA
    final_ic_value_run <- NA; final_coeffs_run <- NULL
    actual_start_p_used <- 2; actual_start_q_used <- 2 # Documenting default
    
    tryCatch({
      ingarch_run_result_obj <- do.call(auto.ingarch, run_params)
      if (!is.null(ingarch_run_result_obj)) {
        stepwise_path_matrix_run <- ingarch_run_result_obj$results
        if (!is.null(ingarch_run_result_obj$model)) {
          final_p_run <- if (!is.null(ingarch_run_result_obj$model$past_obs)) length(ingarch_run_result_obj$model$past_obs) else 0
          final_q_run <- if (!is.null(ingarch_run_result_obj$model$past_mean)) length(ingarch_run_result_obj$model$past_mean) else 0
        }
        if (!is.null(stepwise_path_matrix_run) && is.matrix(stepwise_path_matrix_run) && nrow(stepwise_path_matrix_run) > 0 && !is.na(final_p_run) && !is.na(final_q_run)) {
          df_temp <- as.data.frame(stepwise_path_matrix_run)
          matching_rows <- df_temp[df_temp$p == final_p_run & df_temp$q == final_q_run, ]
          if (nrow(matching_rows) > 0) final_ic_value_run <- min(matching_rows$ic, na.rm = TRUE)
        }
        final_coeffs_run <- ingarch_run_result_obj$coefficients
        cat("    Run ID:", current_run_id_numeric_string, "completed. Final (p=", final_p_run, ",q=", final_q_run, "), AIC=", final_ic_value_run, "\n")
      }
    }, error = function(e_run) {
      cat("    ERROR in Run ID:", current_run_id_numeric_string, ":", e_run$message, "\n")
      error_message_for_run <<- e_run$message
    })
    
    all_run_results[[current_run_id_numeric_string]] <- list(
      run_id = current_run_id_numeric_string,
      original_simulation_index = sim_idx,
      parameters_used = list(
        ic = fixed_ic, distribution = fixed_distribution, link = fixed_link,
        max_p_setting = current_max_p, max_q_setting = current_max_q,
        start_p_function_default = actual_start_p_used, 
        start_q_function_default = actual_start_q_used
      ),
      final_model_orders = list(p = final_p_run, q = final_q_run),
      final_model_ic_value = final_ic_value_run,
      final_model_coefficients = final_coeffs_run,
      stepwise_path_evaluated = stepwise_path_matrix_run,
      error_message = error_message_for_run
    )
  } # end max_pq_settings loop
} # end original simulations loop
cat("---\nAll planned auto.ingarch runs completed.\n")

# STEP 4: Format and Save all collected paths and summaries
# (This part is largely the same, just ensure variable names and filenames match)
# --------------------------------------------------------------------------
cat("Step 4: Formatting and saving all collected results...\n")
all_paths_list_for_df <- list()
all_summaries_list_for_json <- list()

for (run_id_key in names(all_run_results)) { # "1", "2", "3", ...
  run_data <- all_run_results[[run_id_key]]
  
  summary_entry <- run_data 
  summary_entry$stepwise_path_evaluated <- NULL 
  all_summaries_list_for_json[[length(all_summaries_list_for_json) + 1]] <- summary_entry
  
  if (!is.null(run_data$stepwise_path_evaluated) && is.matrix(run_data$stepwise_path_evaluated) && nrow(run_data$stepwise_path_evaluated) > 0) {
    path_matrix_run <- run_data$stepwise_path_evaluated
    df_run <- as.data.frame(path_matrix_run)
    if (nrow(df_run) > 0) {
      df_run$step <- 1:nrow(df_run)
      df_run$run_id <- run_data$run_id 
      df_run$original_simulation_index <- run_data$original_simulation_index
      df_run$ic_used <- run_data$parameters_used$ic # Will always be "aic"
      df_run$dist_used <- run_data$parameters_used$distribution # Will always be "nbinom"
      df_run$link_used <- run_data$parameters_used$link # Will always be "log"
      df_run$max_p_run_setting <- run_data$parameters_used$max_p_setting
      df_run$max_q_run_setting <- run_data$parameters_used$max_q_setting
      all_paths_list_for_df[[length(all_paths_list_for_df) + 1]] <- df_run
    }
  }
}

if (length(all_paths_list_for_df) > 0) {
  combined_paths_df <- do.call(rbind, all_paths_list_for_df)
  if (!is.null(combined_paths_df)) rownames(combined_paths_df) <- NULL
} else {
  combined_paths_df <- data.frame()
  cat("Warning: No valid stepwise paths were collected for CSV.\n")
}

output_csv_all_paths <- "all_paths_fixed_aic_nbinom_log_varied_maxpq.csv" # Descriptive filename
if (nrow(combined_paths_df) > 0) {
  tryCatch({
    write.csv(combined_paths_df, output_csv_all_paths, row.names = FALSE)
    cat("\nAll collected stepwise paths saved to CSV:", output_csv_all_paths, "\n")
  }, error = function(e_csv) {
    cat("\nError saving CSV file:", e_csv$message, "\n")
  })
} else {
  cat("\nNo data to save in combined paths CSV.\n")
}

output_json_all_summaries <- "all_summaries_fixed_aic_nbinom_log_varied_maxpq.json" # Descriptive filename
cleaned_summaries_list_for_json <- lapply(all_summaries_list_for_json, function(item) {
  # (Same cleaning logic as before)
  if (!is.null(item$final_model_coefficients) && !is.list(item$final_model_coefficients)) {
    if(is.numeric(item$final_model_coefficients) && !is.null(names(item$final_model_coefficients))) { item$final_model_coefficients <- as.list(item$final_model_coefficients)
    } else if (!is.atomic(item$final_model_coefficients)) { item$final_model_coefficients <- tryCatch(as.list(item$final_model_coefficients), error = function(e) "conversion_failed") }
  }
  if (!is.null(item$final_model_orders$p) && length(item$final_model_orders$p) > 1) item$final_model_orders$p <- item$final_model_orders$p[1]
  if (!is.null(item$final_model_orders$q) && length(item$final_model_orders$q) > 1) item$final_model_orders$q <- item$final_model_orders$q[1]
  return(item)
})

if (length(cleaned_summaries_list_for_json) > 0) {
  tryCatch({
    write_json(cleaned_summaries_list_for_json, output_json_all_summaries, auto_unbox = TRUE, pretty = TRUE, na = "string")
    cat("All model summaries saved to JSON:", output_json_all_summaries, "\n")
  }, error = function(e_json) {
    cat("\nError saving JSON file:", e_json$message, "\n")
  })
} else {
  cat("\nNo summaries collected to save in JSON.\n")
}

script_end_time <- Sys.time()
total_time_taken <- script_end_time - script_start_time
cat("\nScript finished. Total time taken:", format(total_time_taken), "\n")
cat("Please inspect '", output_csv_all_paths, "' and '", output_json_all_summaries, "'.\n", sep="")
cat("The JSON file contains results for runs '1', '2', '3', '4' (etc.) for each simulation processed.\n")
cat("Each run used AIC, nbinom dist, log link, default start p/q, but varied max p/q search limits.\n")
cat("Check 'final_model_orders' in the JSON to find run IDs for your target (2,6), 'close', and 'off' models.\n")
cat("Then use these simple numeric string run_ids in the Python plotting script.\n")

# END OF R SCRIPT