# --------------------------------------------------------------------------
# Script to apply auto.ingarch to a single time series multiple times
# with varying start.p and start.q to find diverse stepwise paths.
# Results are saved for Python.
# --------------------------------------------------------------------------

# STEP 0: Load custom functions and 'jsonlite'
# ... (same as your provided auto.ingarch, ensure paths are correct) ...
tryCatch({
  source("./custom/auto.ingarch.R")
  source("./custom/newmodel.R")
  source("./custom/ingarch.string.R")
  source("./custom/search.ingarch.R")
  cat("Custom functions loaded successfully.\n")
}, error = function(e) {
  stop("Error loading custom functions: ", e$message)
})

if (!require(jsonlite)) {
  install.packages("jsonlite")
  library(jsonlite)
} else {
  library(jsonlite)
}

# STEP 1: Load simulated data (first time series only)
# --------------------------------------------------------------------------
simulations_filepath <- "C:/Users/Rafael/Desktop/msc-thesis/results/simulatedDataResults/ingarch_no_covariates_simulations.rds" # ADJUST THIS PATH

if (!file.exists(simulations_filepath)) {
  stop("Simulations file not found: ", simulations_filepath)
}
simulations_list <- readRDS(simulations_filepath)
if (!is.list(simulations_list) || length(simulations_list) == 0) {
  stop("Loaded data is not a list or is empty.")
}
selected_ts_data <- as.ts(simulations_list[[1]])
cat("Using the first time series (", length(selected_ts_data), " obs) from: ", basename(simulations_filepath), "\n")


# STEP 2: Define starting parameter combinations and run auto.ingarch
# --------------------------------------------------------------------------
# Define different start.p and start.q combinations to try
# Adjust these based on how wide you want to search.
# Keep start.p <= max.p and start.q <= max.q in auto_ingarch_base_params
start_pq_combinations <- list(
  c(p=0, q=0),
  c(p=1, q=0),
  c(p=0, q=1),
  c(p=1, q=1),
  c(p=2, q=1),
  c(p=1, q=2),
  c(p=2, q=2),
  c(p=3, q=0), # Example: Higher start p
  c(p=0, q=3), # Example: Higher start q
  c(p=3, q=3)
  # Add more combinations if needed, up to max.p/max.q
)

# Base parameters for auto.ingarch (start.p and start.q will be overridden)
auto_ingarch_base_params <- list(
  max.p = 7,
  max.q = 7,
  distribution = "poisson", # or "nbinom"
  link = "log",           # or "identity"
  ic = "aic",             # or "aicc", "bic", "qic"
  stepwise = TRUE,
  trace = FALSE,         # Set to TRUE to see detailed output for each run
  show_warnings = FALSE
)

all_run_results <- list() # To store results from all runs

cat("\nStarting multiple auto.ingarch runs with varying start.p and start.q...\n")

for (i in 1:length(start_pq_combinations)) {
  current_start_p <- start_pq_combinations[[i]]['p']
  current_start_q <- start_pq_combinations[[i]]['q']
  
  run_id <- paste0("start_p", current_start_p, "_q", current_start_q)
  cat("---\nProcessing Run ID:", run_id, "---\n")
  
  # Ensure start.p and start.q are not greater than max.p/max.q
  if (current_start_p > auto_ingarch_base_params$max.p) {
    cat("Skipping run ", run_id, ": start.p (", current_start_p, ") > max.p (", auto_ingarch_base_params$max.p, ")\n")
    next
  }
  if (current_start_q > auto_ingarch_base_params$max.q) {
    cat("Skipping run ", run_id, ": start.q (", current_start_q, ") > max.q (", auto_ingarch_base_params$max.q, ")\n")
    next
  }
  
  current_params <- auto_ingarch_base_params
  current_params$start.p <- current_start_p
  current_params$start.q <- current_start_q
  
  tryCatch({
    ingarch_run_result <- do.call(auto.ingarch, c(list(y = selected_ts_data), current_params))
    
    # Extract information
    stepwise_path_matrix_run <- ingarch_run_result$results
    
    final_p_run <- length(ingarch_run_result$model$past_obs)
    final_q_run <- length(ingarch_run_result$model$past_mean)
    if(is.null(ingarch_run_result$model$past_obs)) final_p_run <- 0
    if(is.null(ingarch_run_result$model$past_mean)) final_q_run <- 0
    
    final_ic_run <- NA
    if (!is.null(stepwise_path_matrix_run) && nrow(stepwise_path_matrix_run) > 0) {
      final_model_row_run <- stepwise_path_matrix_run[
        stepwise_path_matrix_run[, "p"] == final_p_run & stepwise_path_matrix_run[, "q"] == final_q_run, , drop = FALSE
      ]
      if (nrow(final_model_row_run) > 0) {
        final_ic_run <- final_model_row_run[nrow(final_model_row_run), "ic"]
      }
    }
    
    all_run_results[[run_id]] <- list(
      run_id = run_id,
      start_params = list(p = current_start_p, q = current_start_q),
      final_model_orders = list(p = final_p_run, q = final_q_run),
      final_model_ic_type = current_params$ic,
      final_model_ic_value = final_ic_run,
      final_model_coefficients = ingarch_run_result$coefficients,
      stepwise_path = stepwise_path_matrix_run,
      full_ingarch_object = ingarch_run_result # Optionally save the full object if memory allows
    )
    cat("Run ID:", run_id, "completed. Final model: (p=", final_p_run, ", q=", final_q_run, "), IC=", final_ic_run, "\n")
    
  }, error = function(e_run) {
    cat("Error in Run ID:", run_id, ":", e_run$message, "\n")
    all_run_results[[run_id]] <- list(
      run_id = run_id,
      start_params = list(p = current_start_p, q = current_start_q),
      error_message = e_run$message,
      stepwise_path = NULL
    )
  })
}
cat("---\nAll runs attempted.\n")

# STEP 3: Format and Save all collected paths and summaries
# --------------------------------------------------------------------------
all_paths_list_for_df <- list()
all_summaries_list <- list()

for (run_id_name in names(all_run_results)) {
  run_data <- all_run_results[[run_id_name]]
  
  # Add to summaries list
  summary_entry <- run_data # Copy most fields
  summary_entry$stepwise_path <- NULL # Remove bulky path from summary list
  summary_entry$full_ingarch_object <- NULL # Remove full object
  all_summaries_list[[run_id_name]] <- summary_entry
  
  # Format path for combined CSV
  if (!is.null(run_data$stepwise_path) && nrow(run_data$stepwise_path) > 0) {
    path_matrix_run <- run_data$stepwise_path
    df_run <- as.data.frame(path_matrix_run)
    df_run$step <- 1:nrow(df_run)
    df_run$run_id <- run_id_name
    df_run$start_p <- run_data$start_params$p
    df_run$start_q <- run_data$start_params$q
    all_paths_list_for_df[[run_id_name]] <- df_run
  }
}

# Combine all paths into one data frame
combined_paths_df <- do.call(rbind, all_paths_list_for_df)
if (!is.null(combined_paths_df)) rownames(combined_paths_df) <- NULL


# Save combined paths CSV
output_csv_all_paths <- "all_stepwise_paths_varied_starts.csv"
if (!is.null(combined_paths_df) && nrow(combined_paths_df) > 0) {
  tryCatch({
    write.csv(combined_paths_df, output_csv_all_paths, row.names = FALSE)
    cat("\nAll collected stepwise paths saved to CSV:", output_csv_all_paths, "\n")
  }, error = function(e) {
    cat("\nError saving combined paths CSV file:", e$message, "\n")
  })
} else {
  cat("\nNo valid paths collected to save in combined CSV.\n")
}

# Save all summaries JSON
output_json_all_summaries <- "all_model_summaries_varied_starts.json"
if (length(all_summaries_list) > 0) {
  tryCatch({
    write_json(all_summaries_list, output_json_all_summaries, auto_unbox = TRUE, pretty = TRUE)
    cat("All model summaries saved to JSON:", output_json_all_summaries, "\n")
  }, error = function(e) {
    cat("\nError saving all summaries JSON file:", e$message, "\n")
  })
} else {
  cat("\nNo summaries collected to save in JSON.\n")
}

cat("\nScript finished. Multiple runs completed and results saved.\n")
cat("You should now inspect '", output_csv_all_paths, "' and '", output_json_all_summaries, "'.\n", sep="")
cat("From these files, select three 'run_id' entries whose paths and final models represent your desired 'correct', 'close', and 'wrong' scenarios.\n")
cat("You can then filter the CSV by those 'run_id's for plotting or further analysis in Python.\n")


# --------------------------------------------------------------------------
# STEP 4: Next Steps - User Selection and Plotting
# --------------------------------------------------------------------------
# After running this script, you will need to:
# 1. Open 'all_stepwise_paths_varied_starts.csv' and 'all_model_summaries_varied_starts.json'.
# 2. Examine the `final_model_orders` in the JSON (or the last p,q for each run_id in the CSV)
#    to identify three `run_id`s that give you a "correct", "close", and "wrong" outcome
#    based on your own criteria for the `selected_ts_data`.
#
# Example: Let's say you identify these run_ids (replace with your actual findings):
#   selected_run_id_correct <- "start_p2_q2" # Assume this led to your 'correct' (e.g. 2,6)
#   selected_run_id_close   <- "start_p1_q1" # Assume this led to your 'close' (e.g. 1,5)
#   selected_run_id_wrong   <- "start_p0_q3" # Assume this led to your 'wrong' (e.g. 5,5)
#
# Then, you can filter `combined_paths_df` for these specific runs for plotting:
#
# if (!is.null(combined_paths_df) && nrow(combined_paths_df) > 0 &&
#     exists("selected_run_id_correct") && exists("selected_run_id_close") && exists("selected_run_id_wrong")) {
#
#   ids_for_plotting <- c(selected_run_id_correct, selected_run_id_close, selected_run_id_wrong)
#   df_for_plotting <- combined_paths_df[combined_paths_df$run_id %in% ids_for_plotting, ]
#
#   # Add a scenario label for clearer plotting
#   df_for_plotting$scenario_type <- factor(df_for_plotting$run_id,
#                                          levels = ids_for_plotting,
#                                          labels = c("Path to 'Correct' Model",
#                                                     "Path to 'Close' Model",
#                                                     "Path to 'Wrong' Model"))
#
#   if (!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)} else {library(ggplot2)}
#   if (!require(ggrepel)) {install.packages("ggrepel"); library(ggrepel)} else {library(ggrepel)}
#
#   max_p_plot <- max(0, df_for_plotting$p, na.rm = TRUE)
#   max_q_plot <- max(0, df_for_plotting$q, na.rm = TRUE)
#
#   plot_selected_paths <- ggplot(df_for_plotting, aes(x = p, y = q)) +
#     geom_path(aes(group = run_id, color = scenario_type), # Color by scenario_type
#               arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
#               linewidth = 0.7) +
#     geom_point(aes(color = scenario_type), size = 2.5, alpha = 0.8) +
#     geom_text_repel(aes(label = paste0("(", p, ",", q, ")")),
#                     size = 2.5,
#                     segment.alpha = 0.4,
#                     box.padding = unit(0.4, "lines"),
#                     max.overlaps = Inf) + # Allow more overlaps if needed
#     facet_wrap(~ scenario_type, scales = "fixed", ncol = 1) +
#     labs(title = "Selected Stepwise Paths from Varied Starts",
#          subtitle = paste("Using first series from:", basename(simulations_filepath)),
#          x = "Order p", y = "Order q", color = "Path Type") +
#     scale_x_continuous(breaks = 0:max_p_plot, minor_breaks = NULL) +
#     scale_y_continuous(breaks = 0:max_q_plot, minor_breaks = NULL) +
#     theme_minimal(base_size = 10) +
#     theme(legend.position = "none", # Legend is redundant due to facet titles
#           strip.text = element_text(face = "bold"))
#
#   print(plot_selected_paths)
#   # ggsave("plot_selected_3paths.png", plot_selected_paths, width = 7, height = 9, dpi = 300)
# } else {
#   cat("\nCannot create plot: 'combined_paths_df' is empty or selected run_ids are not defined.\n")
#   cat("Please run the script, inspect results, and then define selected_run_id_correct/close/wrong to plot.\n")
# }
# --------------------------------------------------------------------------