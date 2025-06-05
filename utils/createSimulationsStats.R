# --- Configuration ---
# Define your input RDS files and their associated names.
# Each element in the list is another list containing 'path' and 'name'.
# IMPORTANT: For Windows paths in R strings:
#   - Use forward slashes: "C:/Users/YourUser/Documents/Sim1.rds" (Recommended)
#   - OR use escaped backslashes: "C:\\Users\\YourUser\\Documents\\Sim1.rds"

files_to_process <- list(
  # Example 1:
  list(
    path = "./Results/Simulation/no_covariates/Model1/ingarch_no_covariates_simulations.rds",
    name = "M1"
  ),
  # Example 2:
  list(
    path = "./Results/Simulation/no_covariates/Model2/ingarch_no_covariates_simulations.rds",
    name = "M2"
  )
)

# Define the full path (including filename) for your output CSV file.
# Use forward slashes or escaped backslashes.
output_csv_path <- "./Utils/simulationStats.csv"

# --- End Configuration ---


# --- Initialization ---
# This list will store a data frame of results for each processed file
results_collector_list <- list()

# Define expected structure (for warnings/checks)
EXPECTED_SIM_COUNT <- 1000
EXPECTED_SIM_LENGTH <- 1000

# --- Processing Loop ---
if (length(files_to_process) > 0) {
  for (i in 1:length(files_to_process)) {
    file_info <- files_to_process[[i]]
    rds_file_path <- file_info$path
    simulation_set_name <- file_info$name
    
    cat("Processing:", simulation_set_name, "from file:", rds_file_path, "\n")
    
    # Initialize a data frame row for the current file's statistics (with NAs)
    # This ensures a row is added even if processing fails for this file.
    current_file_stats_df <- data.frame(
      `Simulation Name` = simulation_set_name, # Backticks for space in name
      mean_of_counts = NA_real_,
      median_of_counts = NA_real_,
      sd_of_counts = NA_real_,
      mean_under_10 = NA_real_,
      sd_under_10 = NA_real_,
      stringsAsFactors = FALSE,
      check.names = FALSE # Allows non-standard column names (like with spaces)
    )
    
    tryCatch({
      # Read the RDS file
      simulations_list <- readRDS(rds_file_path)
      
      # --- Basic Validation ---
      if (!is.list(simulations_list)) {
        stop("Data read from RDS is not a list.")
      }
      if (length(simulations_list) == 0) {
        stop("The list read from RDS is empty.")
      }
      if (length(simulations_list) != EXPECTED_SIM_COUNT) {
        warning(paste("Expected", EXPECTED_SIM_COUNT, "simulations, but found", length(simulations_list), "in file:", rds_file_path))
      }
      
      # --- Calculate statistics for each simulation within the list ---
      all_simulation_means <- sapply(simulations_list, function(single_simulation) {
        if (!is.numeric(single_simulation)) {
          warning(paste("A simulation in", rds_file_path, "is not numeric. Skipping its mean calculation."))
          return(NA_real_)
        }
        if (length(single_simulation) != EXPECTED_SIM_LENGTH) {
          warning(paste("A simulation in", rds_file_path, "does not have expected length", EXPECTED_SIM_LENGTH, ". Actual:", length(single_simulation)))
        }
        return(mean(single_simulation, na.rm = TRUE))
      })
      
      all_counts_under_10 <- sapply(simulations_list, function(single_simulation) {
        if (!is.numeric(single_simulation)) {
          warning(paste("A simulation in", rds_file_path, "is not numeric. Skipping its 'counts under 10' calculation."))
          return(NA_integer_)
        }
        return(sum(single_simulation < 10, na.rm = TRUE))
      })
      
      # Remove NAs that might have resulted from non-numeric elements before aggregate stats
      all_simulation_means <- na.omit(all_simulation_means)
      all_counts_under_10 <- na.omit(all_counts_under_10)
      
      if(length(all_simulation_means) == 0 || length(all_counts_under_10) == 0){
        stop("No valid simulation data found after initial processing to calculate aggregate statistics.")
      }
      
      # --- Calculate overall statistics for the current file ---
      current_file_stats_df$mean_of_counts <- mean(all_simulation_means, na.rm = TRUE)
      current_file_stats_df$median_of_counts <- median(all_simulation_means, na.rm = TRUE)
      current_file_stats_df$sd_of_counts <- sd(all_simulation_means, na.rm = TRUE)
      
      current_file_stats_df$mean_under_10 <- mean(all_counts_under_10, na.rm = TRUE)
      current_file_stats_df$sd_under_10 <- sd(all_counts_under_10, na.rm = TRUE)
      
      cat("Successfully processed:", simulation_set_name, "\n")
      
    }, error = function(e) {
      # Error occurred while processing this file
      cat("Error processing file:", rds_file_path, "\n")
      cat("Error message:", e$message, "\n")
      # The current_file_stats_df will retain its NA values for this file
    })
    
    # Add the statistics (or NAs if error) for the current file to our collector list
    results_collector_list[[i]] <- current_file_stats_df
  }
} else {
  cat("The 'files_to_process' list is empty. No files to process.\n")
}

# --- Combine results and write to CSV ---
if (length(results_collector_list) > 0) {
  # Combine all the individual data frame rows into one data frame
  final_results_df <- do.call(rbind, results_collector_list)
  
  # Ensure the output directory exists
  output_dir <- dirname(output_csv_path)
  if (!dir.exists(output_dir)) {
    cat("Output directory '", output_dir, "' does not exist. Attempting to create it.\n", sep="")
    tryCatch({
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE) # showWarnings = FALSE to suppress warning if dir already exists
      cat("Output directory created successfully or already existed.\n")
    }, error = function(e_dir){
      cat("Failed to create output directory '", output_dir, "': ", e_dir$message, "\nProceeding to write in current directory if path was relative, or fail if absolute.\n", sep="")
    })
  }
  
  # Write the final data frame to a CSV file
  tryCatch({
    write.csv(final_results_df, file = output_csv_path, row.names = FALSE, na = "NA")
    cat("Results successfully saved to:", output_csv_path, "\n")
  }, error = function(e_write) {
    cat("Error writing CSV file to:", output_csv_path, "\n")
    cat("Error message:", e_write$message, "\n")
  })
  
} else {
  cat("No results were generated to save to CSV.\n")
}

cat("Script finished.\n")