# Install pyreadr: pip install pyreadr
# Install pandas: pip install pandas
# Install numpy: pip install numpy

import pyreadr
import pandas as pd
import numpy as np
import os

def analyze_ingarch_results(rds_file_path, target_orders):
    """
    Reads an RDS file containing INGARCH model results using pyreadr,
    filters by method, calculates statistics for specified (p,q) orders,
    and their selection frequency (as a percentage).

    Args:
        rds_file_path (str): The path to the .rds file.
        target_orders (list of tuples): A list of (p,q) tuples to analyze.
                                         Example: [(2,6), (1,4)]
    Returns:
        pandas.DataFrame: A DataFrame with the results, including (p,q) order,
                          AVG_n_models, STD_n_models, Count, and Selection_Frequency (%).
                          Returns an empty DataFrame if the RDS file cannot be read
                          or if the expected data structure is not found.
    """
    try:
        print(f"Reading RDS file: {rds_file_path} using pyreadr...")
        r_data_dict = pyreadr.read_r(rds_file_path)

        df = None
        if not r_data_dict:
            print("Error: pyreadr returned an empty dictionary. No data found in RDS file.")
            return pd.DataFrame()

        if len(r_data_dict) == 1:
            df = list(r_data_dict.values())[0]
        else:
            found_df_by_name = False
            for key, value in r_data_dict.items():
                if isinstance(value, pd.DataFrame):
                    if "ingarch" in key.lower() or "results" in key.lower():
                        df = value
                        print(f"Found DataFrame named '{key}' in RDS file.")
                        found_df_by_name = True
                        break
            if not found_df_by_name:
                for key, value in r_data_dict.items():
                    if isinstance(value, pd.DataFrame):
                        df = value
                        print(f"Warning: Multiple objects found. Using the first DataFrame encountered, named '{key}'.")
                        break
        
        if df is None:
            print("Error: Could not find a pandas DataFrame in the RDS file contents.")
            print("Contents found:", list(r_data_dict.keys()))
            return pd.DataFrame()

        if not isinstance(df, pd.DataFrame):
            print(f"Error: Extracted object is not a pandas DataFrame. Type: {type(df)}")
            return pd.DataFrame()

        print("Successfully loaded RDS file into pandas DataFrame. Columns:", df.columns.tolist())

    except Exception as e:
        print(f"Error reading RDS file with pyreadr: {e}")
        print("Please ensure the RDS file path is correct and the file is not corrupted.")
        return pd.DataFrame()

    required_columns = ['method', 'p', 'q', 'n_models']
    if not all(col in df.columns for col in required_columns):
        print(f"Error: DataFrame is missing one or more required columns. Found: {df.columns.tolist()}")
        print(f"Required columns are: {required_columns}")
        return pd.DataFrame()

    df_stepwise = df[df['method'] == 'stepwise'].copy()

    if df_stepwise.empty:
        print("No data found for method 'stepwise'.")
        return pd.DataFrame(columns=['(p,q)', 'AVG_n_models', 'STD_n_models', 'Count', 'Selection_Frequency'])

    # Total number of simulations using the stepwise method
    total_stepwise_simulations = df_stepwise.shape[0]
    print(f"Total simulations using 'stepwise' method: {total_stepwise_simulations}")

    # Convert p, q, and n_models to numeric, coercing errors to NaN
    df_stepwise.loc[:, 'p'] = pd.to_numeric(df_stepwise['p'], errors='coerce')
    df_stepwise.loc[:, 'q'] = pd.to_numeric(df_stepwise['q'], errors='coerce')
    df_stepwise.loc[:, 'n_models'] = pd.to_numeric(df_stepwise['n_models'], errors='coerce')

    # Drop rows where p, q, or n_models could not be converted to numeric
    df_stepwise.dropna(subset=['p', 'q', 'n_models'], inplace=True)
    
    # Check if empty again after dropna and before astype
    if df_stepwise.empty: 
        print("No valid data found for method 'stepwise' after coercing p, q, n_models to numeric and dropping NaNs.")
        return pd.DataFrame(columns=['(p,q)', 'AVG_n_models', 'STD_n_models', 'Count', 'Selection_Frequency'])

    # Ensure p and q are integers for matching. Explicitly cast to np.int32 to address the FutureWarning.
    try:
        df_stepwise.loc[:, 'p'] = df_stepwise['p'].astype(np.int32)
        df_stepwise.loc[:, 'q'] = df_stepwise['q'].astype(np.int32)
    except Exception as e:
        print(f"Error converting p or q to np.int32: {e}")
        return pd.DataFrame(columns=['(p,q)', 'AVG_n_models', 'STD_n_models', 'Count', 'Selection_Frequency'])


    results = []
    for p_val, q_val in target_orders:
        subset = df_stepwise[(df_stepwise['p'] == p_val) & (df_stepwise['q'] == q_val)]

        count = subset.shape[0]
        selection_frequency_pct = 0.0 # Default to 0

        if not subset.empty:
            avg_n_models = subset['n_models'].mean()
            std_n_models = subset['n_models'].std()
            std_n_models = std_n_models if not pd.isna(std_n_models) else 0.0
        else:
            avg_n_models = np.nan
            std_n_models = np.nan
            # count is already 0 if subset is empty

        if total_stepwise_simulations > 0:
            selection_frequency_pct = (count / total_stepwise_simulations) * 100 # Calculate as percentage
        else: 
            selection_frequency_pct = np.nan if count > 0 else 0.0


        results.append({
            '(p,q)': f"({p_val},{q_val})",
            'AVG_n_models': avg_n_models,
            'STD_n_models': std_n_models,
            'Count': count,
            'Selection_Frequency': selection_frequency_pct # Store percentage
        })

    results_df = pd.DataFrame(results)
    return results_df

if __name__ == "__main__":
    # --- Configuration ---
    rds_file_path = r'C:\Users\Rafael\Desktop\msc-thesis\old\results\simulatedDataResults\ingarch_no_covariates_results.rds' # EXAMPLE PATH
    if rds_file_path == 'path/to/your/ingarch_no_covariates_results.rds' or not os.path.exists(rds_file_path):
        print(f"WARNING: The RDS file path '{rds_file_path}' is a placeholder or does not exist.")
        print("Please update the 'rds_file_path' variable in the script with the correct path to your file.")
        if not os.path.exists(rds_file_path):
            print("Script will likely fail or produce no results as RDS file was not found and no dummy data loading is implemented here.")
            pass


    target_pq_orders = [(2,6), (1,4), (1,5), (1,7), (2,5), (1,1)]
    # --- End Configuration ---

    analysis_results_df = analyze_ingarch_results(rds_file_path, target_pq_orders)

    if not analysis_results_df.empty:
        print("\n--- Analysis Results ---")
        print("Target (p,q) orders, 'n_models' statistics, and selection frequency:")

        # Adjusted header for Sel_Freq (%)
        header = f"{'(p,q)':<10} | {'AVG_n_models':<15} | {'STD_n_models':<15} | {'Count':<7} | {'Sel_Freq (%)':<12}"
        print(header)
        print("-" * len(header))

        for index, row in analysis_results_df.iterrows():
            avg_str = f"{row['AVG_n_models']:.2f}" if not pd.isna(row['AVG_n_models']) else "N/A"
            std_str = f"{row['STD_n_models']:.2f}" if not pd.isna(row['STD_n_models']) else "N/A"
            count_str = f"{int(row['Count'])}"
            # Format frequency as percentage with 2 decimal places
            freq_str = f"{row['Selection_Frequency']:.2f}%" if not pd.isna(row['Selection_Frequency']) else "N/A"
            print(f"{row['(p,q)']:<10} | {avg_str:<15} | {std_str:<15} | {count_str:<7} | {freq_str:<12}")
        
        print("\nNote: 'N/A' means no data was found for that specific (p,q) order, or only one data point was found for STD calculation.")
        print("Sel_Freq (%) is the selection frequency as a percentage relative to all 'stepwise' method simulations.")
    else:
        print("No results to display. This could be due to an issue with the RDS file, the path, or no matching data.")

