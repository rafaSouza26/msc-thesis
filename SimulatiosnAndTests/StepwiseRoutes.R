# Main script: run_simulations_no_ts_validation.R

# --- User Instructions ---
# 1. CRITICAL: Set your R working directory to the PARENT directory of your 'custom' folder
#    BEFORE running this script. (See detailed example in instructions above).
#
# 2. CRITICAL: Ensure that within your custom .R files (e.g., auto.ingarch.R),
#    any `source()` calls to other files ALSO IN THE 'custom' FOLDER use the
#    relative path starting with "./custom/".
#    Example: inside 'custom/auto.ingarch.R', use source("./custom/newmodel.R")

cat("Starting R script for stepwise path extraction (no explicit TS validation in main script)...\n")
script_start_time <- Sys.time()

# --- Load necessary packages ---
cat("Loading 'tscount' package...\n")
if (!require(tscount)) {
  install.packages("tscount", repos = "http://cran.us.r-project.org")
  library(tscount)
}
cat("'tscount' package is ready.\n")

# --- Define paths and Sourcing ---
# All paths for sourcing custom functions will be relative to the current working directory.
# The user MUST set the working directory to be the parent of the 'custom' folder.

current_wd <- getwd()
cat("Current working directory (should be parent of 'custom'):", current_wd, "\n")
cat("The 'custom' folder is expected to be at:", file.path(current_wd, "custom"), "\n")

# Files to source. Paths are relative to the current working directory.
# Order matters for dependencies.
custom_files_to_source <- c(
  "./custom/myingarch.R",      # Dependency for search.ingarch.R
  "./custom/newmodel.R",       # Dependency for auto.ingarch.R
  "./custom/ingarch.string.R", # Dependency for auto.ingarch.R
  "./custom/search.ingarch.R", # Dependency for auto.ingarch.R, sources myingarch.R internally
  "./custom/auto.ingarch.R"    # Main function we'll be calling, sources its dependencies
)

cat("Sourcing custom functions using paths relative to current WD:\n")
for (path_to_script in custom_files_to_source) {
  if (!file.exists(path_to_script)) {
    stop(paste("File '", path_to_script, "' not found. ",
               "Please ensure your R working directory is set to the PARENT of the 'custom' folder, ",
               "and the file '", basename(path_to_script), "' exists inside 'custom'. Current WD: ", current_wd, sep=""))
  }
  cat("  Sourcing:", path_to_script, "\n")
  tryCatch({
    source(path_to_script) # e.g., source("./custom/myingarch.R")
  }, error = function(e) {
    stop(paste("Error sourcing '", path_to_script, "': ", e$message, sep=""))
  })
}

# Check if the main function is loaded
if (!exists("auto.ingarch") || !is.function(auto.ingarch)) {
  stop("The 'auto.ingarch' function was not loaded correctly. Check sourcing steps, file paths, and internal source calls within your custom scripts.")
}
cat("Custom functions sourced successfully.\n")


# --- User-defined parameters for the analysis ---
# Path to your .RDS file containing the time series simulations list
# This path should also be relative to your current working directory, or an absolute path.
rds_file_path <- "./results/simulatedDataResults/ingarch_no_covariates_simulations.rds" # Adjust if needed

# Number of simulations to process from the list
num_simulations_to_process <- 100 # For testing; set to length(simulations_list) for all.

# Path for the output CSV file (will be saved in your current working directory)
output_csv_path <- "auto_ingarch_stepwise_paths.csv"

# --- Load the simulations data from .rds file ---
cat("Loading simulations data from .rds file:", rds_file_path, "\n")
if (!file.exists(rds_file_path)) {
  stop(paste("Error: Cannot find .rds file at '", rds_file_path, "'. This path is relative to '", current_wd, "'. Please check.", sep=""))
}
simulations_list <- readRDS(rds_file_path)

if (!is.list(simulations_list)) {
  stop(paste("Error: The object loaded from '", rds_file_path, "' is not a list. It is a '", class(simulations_list)[1], "'. Expected a list of time series.", sep=""))
}
cat("Successfully loaded", length(simulations_list), "simulations/elements from '", rds_file_path, "'.\n")


# --- Initialize a list to store the results ---
results_list <- list()

# --- Loop through the specified number of simulations ---
if (num_simulations_to_process > length(simulations_list)) {
  warning(paste("Requested 'num_simulations_to_process' (", num_simulations_to_process,
                ") is greater than available elements (", length(simulations_list),
                "). Processing all ", length(simulations_list), " available elements.", sep=""))
  num_simulations_to_process <- length(simulations_list)
}
if (num_simulations_to_process == 0 && length(simulations_list) > 0) {
  warning("num_simulations_to_process is 0. No elements will be processed.")
}

cat("\nStarting processing of up to", num_simulations_to_process, "elements from the list...\n")

for (i in 1:num_simulations_to_process) {
  current_element_id <- as.character(i)
  sim_id_label <- paste0("run_", current_element_id)
  
  cat("-------------------------------------------------\n")
  cat("Processing element #", current_element_id, " (Output ID: ", sim_id_label, ")...\n", sep="")
  
  current_ts <- simulations_list[[i]]
  
  # Removed explicit time series validation here, as requested.
  # auto.ingarch.R is expected to handle validation.
  
  path_taken_for_current_ts <- NA_character_
  
  # Call auto.ingarch with specified parameters
  model_output_object <- tryCatch({
    cat("  Applying auto.ingarch (stepwise=TRUE, max.p=7, max.q=7, link='log', distr='nbinom')...\n")
    # Ensure auto.ingarch is called from the global environment where it was sourced
    # No need to change WD here if sourcing setup is correct
    auto.ingarch(
      y = current_ts,
      stepwise = TRUE,
      max.p = 7,
      max.q = 7,
      link = "log",
      distribution = "nbinom",
      xreg = NULL,
      trace = FALSE,          # As per your original request
      show_warnings = FALSE,  # As per your original request
      ic = "aicc"             # Default from auto.ingarch, can be changed
    )
  }, error = function(e) {
    # This catches errors from the auto.ingarch call itself (e.g., if auto.ingarch stops due to invalid data)
    warning(paste("  Error during auto.ingarch call for element #", current_element_id, ": ", e$message, sep=""), call. = FALSE)
    return(NULL) # Return NULL to indicate an error occurred during the main call
  })
  
  # --- Extract the path information from model_output_object$results ---
  if (!is.null(model_output_object)) {
    # Check if auto.ingarch returned an error structure (e.g., from myingarch.R)
    if (!is.null(model_output_object$ic) && is.infinite(model_output_object$ic) && !is.null(model_output_object$error)) {
      cat("  auto.ingarch indicated an internal error for element #", current_element_id, ": ", model_output_object$error, "\n")
      path_taken_for_current_ts <- paste("Error in auto.ingarch (reported by function):", model_output_object$error)
    } else if (!is.null(model_output_object$results) && is.matrix(model_output_object$results) && ncol(model_output_object$results) >= 3) {
      cat("  auto.ingarch call completed. Attempting to extract path from model_output_object$results...\n")
      path_steps <- character(nrow(model_output_object$results))
      # Ensure column names "p", "q", "ic" exist
      if (all(c("p", "q", "ic") %in% colnames(model_output_object$results))) {
        for (step_idx in 1:nrow(model_output_object$results)) {
          p_val <- model_output_object$results[step_idx, "p"]
          q_val <- model_output_object$results[step_idx, "q"]
          ic_val <- model_output_object$results[step_idx, "ic"]
          ic_display <- if (is.na(ic_val) || is.infinite(ic_val)) as.character(ic_val) else sprintf("%.2f", ic_val)
          path_steps[step_idx] <- sprintf("(%s,%s) IC:%s",
                                          ifelse(is.na(p_val), "NA", p_val),
                                          ifelse(is.na(q_val), "NA", q_val),
                                          ic_display)
        }
        path_taken_for_current_ts <- paste(path_steps, collapse = " -> ")
        cat("    Path extracted successfully.\n")
      } else {
        warning(paste("  Could not extract path for element #", current_element_id,
                      ". 'model_output_object$results' matrix is missing expected columns (p, q, ic).", sep=""), call. = FALSE)
        path_taken_for_current_ts <- "Path extraction failed: model_output_object$results missing p,q,ic columns"
      }
    } else {
      # This case handles if model_output_object is not NULL, but doesn't fit the error pattern
      # and also doesn't have the expected $results structure.
      warning(paste("  Could not extract path for element #", current_element_id,
                    ". 'model_output_object$results' was NULL, not a matrix, or had too few columns, or an unrecognized structure.", sep=""), call. = FALSE)
      cat("    Structure of model_output_object for element #", current_element_id, ":\n")
      utils::str(model_output_object, max.level = 2) # Print structure to help debug
      path_taken_for_current_ts <- "Path extraction failed: model_output_object$results has unexpected structure or object is not a standard model output"
    }
  } else { # This 'else' corresponds to the tryCatch for auto.ingarch call returning NULL
    path_taken_for_current_ts <- "Error: auto.ingarch main call failed (tryCatch returned NULL)"
    # The specific error message from auto.ingarch (if it stopped) would have been printed by the warning in tryCatch
    cat("  auto.ingarch main call failed, or the function stopped and returned NULL via tryCatch.\n")
  }
  # --- End of path extraction ---
  
  results_list[[sim_id_label]] <- list(simulation_id = sim_id_label, path_taken = path_taken_for_current_ts)
  cat("  Finished processing element #", current_element_id, ". Path/Status recorded.\n", sep="")
}
cat("-------------------------------------------------\n")
cat("All specified elements processed.\n")

# --- Convert the list of results to a data frame ---
cat("\nConverting results to data frame...\n")
if (length(results_list) > 0) {
  paths_df <- data.frame(
    Simulation_ID = sapply(results_list, function(x) if(!is.null(x$simulation_id)) x$simulation_id else NA_character_),
    Path_Taken = sapply(results_list, function(x) if(!is.null(x$path_taken)) x$path_taken else NA_character_)
  )
  rownames(paths_df) <- NULL
} else {
  paths_df <- data.frame(Simulation_ID = character(0), Path_Taken = character(0))
  cat("No results to convert to data frame.\n")
}

# --- Save the data frame to a CSV file ---
cat("Saving results to CSV file:", output_csv_path, "\n")
tryCatch({
  write.csv(paths_df, file = output_csv_path, row.names = FALSE, quote = TRUE)
  cat("Successfully saved paths to:", output_csv_path, "\n")
}, error = function(e) {
  cat("Error saving CSV file '", output_csv_path, "': ", e$message, "\n", sep="")
  cat("The data frame to be saved looks like this (first few rows):\n")
  print(utils::head(paths_df))
})

script_end_time <- Sys.time()
total_time_taken <- script_end_time - script_start_time
cat("\nScript finished. Total time taken:", format(total_time_taken), "\n")
cat("Please inspect '", output_csv_path, "'.\n", sep="")

