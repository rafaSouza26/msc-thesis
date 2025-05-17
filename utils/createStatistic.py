import pandas as pd
import pyreadr
import numpy as np
import os

def analyze_rds_file(input_file, output_file=None):
    # Read the RDS file
    print(f"Reading file: {input_file}")
    result = pyreadr.read_r(input_file)
    
    # RDS files are dictionaries with one key, get the dataframe
    df = result[None]  # None is the key for the default dataframe
    
    # Group by method and calculate statistics
    stats = df.groupby('method').agg({
        'p': ['mean', 'median', 'std'],
        'q': ['mean', 'median', 'std'],
        'time': ['mean', 'median', 'std'],
        'n_models': ['mean', 'median', 'std']
    }).reset_index()
    
    # Flatten the multi-level column names
    stats.columns = [
        'method' if col[0] == 'method' else f"{col[0]}_{col[1]}"
        for col in stats.columns
    ]
    
    # Rename columns to match requested names
    stats = stats.rename(columns={
        'p_mean': 'mean_p',
        'p_median': 'median_p',
        'p_std': 'sd_p',
        'q_mean': 'mean_q',
        'q_median': 'median_q',
        'q_std': 'sd_q',
        'time_mean': 'mean_time',
        'time_median': 'median_time',
        'time_std': 'sd_time',
        'n_models_mean': 'mean_models_tested',
        'n_models_median': 'median_models',
        'n_models_std': 'sd_models_tested'
    })
    
    # If output file is not specified, create one based on the input filename
    if output_file is None:
        input_dir = os.path.dirname(input_file)
        input_filename = os.path.splitext(os.path.basename(input_file))[0]
        output_file = os.path.join(input_dir, f"{input_filename}_stats.csv")
    
    # Save to CSV
    print(f"Saving results to: {output_file}")
    stats.to_csv(output_file, index=False)
    print("Analysis complete!")
    
    return stats

if __name__ == "__main__":
    # ====================================================================
    # EDIT THESE PATHS TO MATCH YOUR WINDOWS FOLDER AND FILE STRUCTURE
    # ====================================================================
    
    # Input RDS file path - EDIT THIS LINE with your Windows path
    input_file = r"C:\Users\rafas\Desktop\msc-thesis\old\results\simulatedDataResults\ingarch_no_covariates_results.rds"
    
    # Output CSV file path - EDIT THIS LINE with your Windows path
    # If you want the output in the same folder as the input with an automatic name,
    # you can set this to None
    output_file = r"C:\Users\rafas\Desktop\msc-thesis\utils\no_covariates_stats.csv"
    
    # ====================================================================
    # END OF PATH CONFIGURATION
    # ====================================================================
    
    # Process the file
    result_stats = analyze_rds_file(input_file, output_file)
    
    # Display the first few rows of the results
    print("\nResult Preview:")
    print(result_stats.head())