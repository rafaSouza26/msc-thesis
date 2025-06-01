import pandas as pd
import numpy as np
import os

def analyze_csv_file(input_file, output_file=None):
    """
    Reads a CSV file, calculates statistics for specified columns grouped by 'method',
    and saves the results to a new CSV file.
    """
    # Read the CSV file
    print(f"Reading CSV file: {input_file}")
    try:
        df = pd.read_csv(input_file)
    except FileNotFoundError:
        print(f"Error: The file '{input_file}' was not found.")
        return None
    except Exception as e:
        print(f"Error reading CSV file '{input_file}': {e}")
        return None
    
    if df.empty:
        print(f"Warning: The CSV file '{input_file}' is empty or could not be read into a DataFrame.")
        return None
        
    # Check for required columns for aggregation
    required_cols_for_agg = ['method', 'p_order', 'q_order', 'time', 'n_models_tested']
    missing_cols = [col for col in required_cols_for_agg if col not in df.columns]
    if missing_cols:
        print(f"Error: The CSV file '{input_file}' is missing the following required columns for analysis: {', '.join(missing_cols)}")
        print(f"Available columns: {df.columns.tolist()}")
        return None

    # Group by method and calculate statistics
    print("Calculating statistics...")
    try:
        stats = df.groupby('method').agg({
            'p_order': ['mean', 'median', 'std'],
            'q_order': ['mean', 'median', 'std'],
            'time': ['mean', 'median', 'std'],
            'n_models_tested': ['mean', 'median', 'std']
        }).reset_index()
    except Exception as e:
        print(f"Error during statistics calculation: {e}")
        return None
    
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
        if not input_dir: # If input_file is just a filename without a path
            input_dir = "." # Save in current directory
        input_filename = os.path.splitext(os.path.basename(input_file))[0]
        output_file = os.path.join(input_dir, f"{input_filename}_stats.csv")
    
    # Ensure output directory exists
    output_dir_path = os.path.dirname(output_file)
    if output_dir_path and not os.path.exists(output_dir_path):
        try:
            os.makedirs(output_dir_path)
            print(f"Created output directory: {output_dir_path}")
        except Exception as e:
            print(f"Error creating output directory '{output_dir_path}': {e}")
            # Optionally, default to saving in the current directory or handle error
            output_file = os.path.basename(output_file) # Save in CWD as fallback
            print(f"Attempting to save in current directory as: {output_file}")


    # Save to CSV
    print(f"Saving results to: {output_file}")
    try:
        stats.to_csv(output_file, index=False)
        print("Analysis complete! 🎉")
    except Exception as e:
        print(f"Error saving results to '{output_file}': {e}")
        return None
        
    return stats

if __name__ == "__main__":
    # ====================================================================
    # EDIT THESE PATHS TO MATCH YOUR WINDOWS FOLDER AND FILE STRUCTURE
    # ====================================================================
    
    # Input CSV file path - EDIT THIS LINE with your Windows path
    # Ensure this file has columns like 'method', 'p', 'q', 'time', 'n_models'
    input_file = r"C:\Users\Rafael\Desktop\msc-thesis\Results\Simulation\with_covariates\ingarch_with_covariates_results.csv"
    
    # Output CSV file path - EDIT THIS LINE with your Windows path
    # If you want the output in the same folder as the input with an automatic name,
    # you can set this to None
    output_file = r"C:\Users\Rafael\Desktop\msc-thesis\utils\ingarch_with_covariates_stats.csv"
    # output_file = None # Uncomment this line to auto-generate output filename
    
    # ====================================================================
    # END OF PATH CONFIGURATION
    # ====================================================================
    
    # Process the file
    result_stats = analyze_csv_file(input_file, output_file)
    
    # Display the first few rows of the results if successful
    if result_stats is not None:
        print("\nResult Preview:")
        print(result_stats.head())
    else:
        print("\nAnalysis did not complete successfully.")