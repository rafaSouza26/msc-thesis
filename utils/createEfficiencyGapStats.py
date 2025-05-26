# Install pandas: pip install pandas
# Install numpy: pip install numpy

import pandas as pd
import numpy as np
import os

def analyze_ingarch_csv_results(csv_file_path, target_orders):
    """
    Reads a CSV file containing INGARCH model results,
    filters by method 'stepwise', calculates statistics for specified (p_order,q_order) orders,
    and their selection frequency (as a percentage).
    Assumes input columns: 'method', 'p_order', 'q_order', 'n_models_tested'.

    Args:
        csv_file_path (str): The path to the .csv file.
        target_orders (list of tuples): A list of (p,q) tuples to analyze.
                                      Example: [(2,6), (1,4)]
    Returns:
        pandas.DataFrame: A DataFrame with the results, including (p,q) order,
                          AVG_n_models, STD_n_models, Count, and Selection_Frequency (%).
                          Returns an empty DataFrame if the CSV file cannot be read,
                          is empty, or if the expected data structure is not found.
    """
    try:
        print(f"Reading CSV file: {csv_file_path} using pandas...")
        df = pd.read_csv(csv_file_path)
        print("Successfully loaded CSV file into pandas DataFrame. Columns:", df.columns.tolist())

    except FileNotFoundError:
        print(f"Error: The file '{csv_file_path}' was not found.")
        print("Please ensure the CSV file path is correct.")
        return pd.DataFrame()
    except pd.errors.EmptyDataError:
        print(f"Error: The CSV file '{csv_file_path}' is empty.")
        return pd.DataFrame()
    except Exception as e:
        print(f"Error reading CSV file with pandas: {e}")
        print("Please ensure the CSV file path is correct and the file is not corrupted and is a valid CSV.")
        return pd.DataFrame()

    if df.empty:
        print("Warning: The loaded DataFrame is empty.")
        return pd.DataFrame()

    # Updated required columns
    required_columns = ['method', 'p_order', 'q_order', 'n_models_tested']
    if not all(col in df.columns for col in required_columns):
        print(f"Error: DataFrame is missing one or more required columns. Found: {df.columns.tolist()}")
        print(f"Required columns are: {required_columns}")
        return pd.DataFrame()

    # Standardize 'method' column: convert to lowercase and strip whitespace
    if 'method' in df.columns:
        df['method'] = df['method'].astype(str).str.lower().str.strip()
    else: 
        print("Error: 'method' column is missing, which is critical for filtering.")
        return pd.DataFrame()


    df_stepwise = df[df['method'] == 'stepwise'].copy() 

    if df_stepwise.empty:
        print("No data found for method 'stepwise'.")
        return pd.DataFrame(columns=['(p,q)', 'AVG_n_models', 'STD_n_models', 'Count', 'Selection_Frequency'])
    
    # Convert p_order, q_order, and n_models_tested to numeric, coercing errors to NaN
    df_stepwise.loc[:, 'p_order'] = pd.to_numeric(df_stepwise['p_order'], errors='coerce')
    df_stepwise.loc[:, 'q_order'] = pd.to_numeric(df_stepwise['q_order'], errors='coerce')
    df_stepwise.loc[:, 'n_models_tested'] = pd.to_numeric(df_stepwise['n_models_tested'], errors='coerce')

    # Drop rows where p_order, q_order, or n_models_tested could not be converted to numeric (became NaN)
    df_stepwise.dropna(subset=['p_order', 'q_order', 'n_models_tested'], inplace=True)
    
    total_stepwise_simulations_valid_data = df_stepwise.shape[0]
    print(f"Total valid simulations for 'stepwise' method after cleaning: {total_stepwise_simulations_valid_data}")


    if df_stepwise.empty: 
        print("No valid data found for method 'stepwise' after coercing p_order, q_order, n_models_tested to numeric and dropping NaNs.")
        return pd.DataFrame(columns=['(p,q)', 'AVG_n_models', 'STD_n_models', 'Count', 'Selection_Frequency'])

    # Ensure p_order and q_order are integers for matching.
    try:
        df_stepwise.loc[:, 'p_order'] = df_stepwise['p_order'].astype(np.int32)
        df_stepwise.loc[:, 'q_order'] = df_stepwise['q_order'].astype(np.int32)
    except Exception as e:
        print(f"Error converting p_order or q_order to np.int32: {e}")
        return pd.DataFrame(columns=['(p,q)', 'AVG_n_models', 'STD_n_models', 'Count', 'Selection_Frequency'])


    results = []
    for p_val, q_val in target_orders:
        # Use updated column names for filtering
        subset = df_stepwise[(df_stepwise['p_order'] == p_val) & (df_stepwise['q_order'] == q_val)]

        count = subset.shape[0]
        selection_frequency_pct = 0.0 

        if not subset.empty:
            # Use updated column name for statistics
            avg_n_models = subset['n_models_tested'].mean()
            std_n_models = subset['n_models_tested'].std()
            std_n_models = std_n_models if not pd.isna(std_n_models) else 0.0 
        else:
            avg_n_models = np.nan 
            std_n_models = np.nan 
            
        if total_stepwise_simulations_valid_data > 0:
            selection_frequency_pct = (count / total_stepwise_simulations_valid_data) * 100
        elif count > 0 : 
             selection_frequency_pct = np.nan 
        else: 
            selection_frequency_pct = 0.0


        results.append({
            '(p,q)': f"({p_val},{q_val})",
            'AVG_n_models': avg_n_models,
            'STD_n_models': std_n_models,
            'Count': count,
            'Selection_Frequency': selection_frequency_pct
        })

    results_df = pd.DataFrame(results)
    return results_df

if __name__ == "__main__":
    # --- Configuration ---
    # EXAMPLE: Replace with the actual path to your CSV file
    csv_file_path = r'C:\Users\Rafael\Desktop\msc-thesis\SimulationResults\ingarch_with_covariates_results.csv' 
    
    if csv_file_path == r'C:\Users\YourUser\Desktop\path_to_your_data\ingarch_results_new_columns.csv' or \
       not os.path.exists(csv_file_path):
        print(f"WARNING: The CSV file path '{csv_file_path}' is a placeholder or does not exist.")
        print("Please update the 'csv_file_path' variable in the script with the correct path to your file.")
        if not os.path.exists(csv_file_path):
            print("Attempting to run with dummy data as file was not found.")
            # Updated dummy data with new column names
            dummy_data = {
                'method': ['stepwise', 'stepwise', 'grid_search', 'stepwise', 'stepwise', 'stepwise', 'stepwise'],
                'p_order': [2, 1, 2, 2, 1, 2, 3],         # Renamed
                'q_order': [6, 4, 6, 6, 5, 5, 1],         # Renamed
                'n_models_tested': [10, 12, 100, 11, 15, 13, 9], # Renamed
                'other_col': [1,2,3,4,5,6,7]
            }
            dummy_csv_path = "dummy_ingarch_results_new_columns.csv"
            pd.DataFrame(dummy_data).to_csv(dummy_csv_path, index=False)
            csv_file_path = dummy_csv_path 
            print(f"Created and using dummy CSV file: {dummy_csv_path}")
            
    target_pq_orders = [(2,6), (1,4), (1,5), (1,7), (2,5), (1,1), (3,3)] 
    # --- End Configuration ---

    analysis_results_df = analyze_ingarch_csv_results(csv_file_path, target_pq_orders)

    if not analysis_results_df.empty:
        print("\n--- Analysis Results ---")
        print("Target (p,q) orders, 'n_models_tested' statistics, and selection frequency (for 'stepwise' method):") # Updated description

        header = f"{'(p,q)':<10} | {'AVG_n_models':<15} | {'STD_n_models':<15} | {'Count':<7} | {'Sel_Freq (%)':<12}"
        print(header)
        print("-" * len(header))

        for index, row in analysis_results_df.iterrows():
            avg_str = f"{row['AVG_n_models']:.2f}" if not pd.isna(row['AVG_n_models']) else "N/A"
            std_str = f"{row['STD_n_models']:.2f}" if not pd.isna(row['STD_n_models']) else "N/A"
            count_str = f"{int(row['Count'])}"
            freq_str = f"{row['Selection_Frequency']:.2f}%" if not pd.isna(row['Selection_Frequency']) else "N/A"
            print(f"{row['(p,q)']:<10} | {avg_str:<15} | {std_str:<15} | {count_str:<7} | {freq_str:<12}")
        
        print("\nNote: 'N/A' means no data was found for that specific (p,q) order, or only one data point was found for STD calculation.")
        print("Sel_Freq (%) is the selection frequency as a percentage relative to all valid 'stepwise' method simulations in the input file.")
    else:
        print("No results to display. This could be due to an issue with the CSV file, the path, or no matching 'stepwise' data.")

    if 'dummy_csv_path' in locals() and os.path.exists(dummy_csv_path):
        os.remove(dummy_csv_path)
        print(f"Cleaned up dummy CSV file: {dummy_csv_path}")