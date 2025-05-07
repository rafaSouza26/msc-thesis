# %% Import Libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pyreadr # For reading .rds files

# %% Plotting Function (Difference Boxplots)
def create_difference_boxplots(filepath, 
                               case_label, 
                               real_p, 
                               real_q, 
                               output_dir="ingarch_order_differences_boxplots", 
                               method_rename_map=None, 
                               custom_palette=None):
    """
    Creates and saves box plots of the difference between estimated and real model orders (p, q)
    from an .rds file, allowing custom method display names, palettes, and adjusted box spacing.

    Args:
        filepath (str): The path to the .rds file.
        case_label (str): Label for the dataset (e.g., "No Covariates"). Used in title/filename.
        real_p (int): The real order p.
        real_q (int): The real order q.
        output_dir (str): The directory where plots will be saved.
        method_rename_map (dict, optional): Maps original method names to display names.
                                            e.g., {"stepwise": "Stepwise Strategy A"}.
                                            Defaults to None (use original names).
        custom_palette (dict, optional): Maps display names (after renaming) to colors.
                                         e.g., {"Stepwise Strategy A": "blue"}.
                                         Defaults to None.
    """
    # Create the output directory if it doesn't exist
    try:
        os.makedirs(output_dir, exist_ok=True)
        print(f"Output directory for difference boxplots set to: {output_dir}")
    except OSError as e:
        print(f"Error creating directory {output_dir}: {e}")
        return

    # Check if the file exists
    if not os.path.exists(filepath):
        print(f"Error: File not found at '{filepath}'. Please check the path.")
        return

    try:
        # --- Read .rds file ---
        print(f"Attempting to read RDS file: {filepath}")
        rds_data = pyreadr.read_r(filepath)
        if not rds_data:
            print(f"Error: RDS file '{filepath}' is empty or pyreadr could not extract any R objects.")
            return
        
        df_key = list(rds_data.keys())[0] 
        df = rds_data[df_key]
        print(f"Successfully read object '{df_key}' from RDS file: {os.path.basename(filepath)} with {df.shape[0]} rows and {df.shape[1]} columns.")

        # --- Column Name Handling (Focus on method, p, q) ---
        expected_cols_map = {
            'method': ['method', 'Method'],
            'p': ['p', 'P'],
            'q': ['q', 'Q']
        }
        rename_operation_map = {}
        current_columns_in_df = df.columns.tolist()
        
        for standard_col_name, possible_file_names in expected_cols_map.items():
            found_column_in_file = False
            for file_col_name_variant in possible_file_names:
                matching_df_cols = [col_in_df for col_in_df in current_columns_in_df if col_in_df.lower() == file_col_name_variant.lower()]
                if matching_df_cols:
                    actual_col_name_in_df = matching_df_cols[0]
                    if actual_col_name_in_df != standard_col_name:
                        if actual_col_name_in_df not in rename_operation_map and standard_col_name not in current_columns_in_df :
                             rename_operation_map[actual_col_name_in_df] = standard_col_name
                    found_column_in_file = True
                    break
            if not found_column_in_file:
                print(f"Error: Required column for '{standard_col_name}' (e.g., {possible_file_names}) not found in RDS object '{df_key}'.")
                print(f"   Available columns in RDS object: {current_columns_in_df}")
                return
        
        df.rename(columns=rename_operation_map, inplace=True)
        print(f"Columns after attempting initial standardization: {df.columns.tolist()}")

        required_standard_cols = ['method', 'p', 'q']
        for req_col in required_standard_cols:
            if req_col not in df.columns:
                print(f"Error: Standardized column '{req_col}' is missing after renaming attempts.")
                print(f"   Available columns: {df.columns.tolist()}")
                return

        # Standardize 'method' column content (e.g., to lowercase)
        df['method'] = df['method'].astype(str).str.lower().str.strip()

        # --- Apply custom method renaming ---
        if method_rename_map:
            # Create 'display_method' column, using mapped name or original if not in map
            df['display_method'] = df['method'].map(method_rename_map).fillna(df['method'])
            print(f"Applied method rename map. Unique display methods: {df['display_method'].unique()}")
        else:
            # If no map, display_method is the same as the standardized method
            df['display_method'] = df['method']
            print(f"No method rename map provided. Using original standardized methods: {df['display_method'].unique()}")


        # --- Data Cleaning (p, q) ---
        for col_to_clean in ['p', 'q']:
            df[col_to_clean] = pd.to_numeric(df[col_to_clean], errors='coerce')
        df.dropna(subset=['display_method', 'p', 'q'], inplace=True) # Use display_method here for consistency
        if df.empty:
            print(f"Warning: No valid data remaining after cleaning for file '{filepath}'. Skipping difference boxplot.")
            return
        df[['p', 'q']] = df[['p', 'q']].astype(int)

        individual_estimations_df = df.copy()

        individual_estimations_df['p_diff'] = individual_estimations_df['p'] - real_p
        individual_estimations_df['q_diff'] = individual_estimations_df['q'] - real_q

        # Melt the DataFrame for plotting, keeping 'display_method' as the identifier for hue
        plot_df = individual_estimations_df.melt(id_vars=['display_method'], 
                                                 value_vars=['p_diff', 'q_diff'],
                                                 var_name='difference_source',
                                                 value_name='difference_value')
        plot_df['Type'] = plot_df['difference_source'].apply(lambda x: 'P' if x == 'p_diff' else 'Q')

        # --- Determine Palette for Plotting ---
        active_palette = None
        if custom_palette:
            active_palette = custom_palette # User-provided palette, keyed by display_method values
            print("Using user-provided custom_palette.")
        elif not method_rename_map: # No renaming by user, and no specific custom_palette
            # Use a default internal palette keyed by original standardized method names
            default_internal_palette = {
                "grid_search": "mediumseagreen",
                "stepwise": "firebrick"
            }
            active_palette = default_internal_palette
            print("Using default internal palette for original method names.")
        # If method_rename_map IS used, but custom_palette is NOT provided,
        # active_palette remains None, and seaborn will assign default colors.
        else:
            print("Methods renamed, but no custom_palette provided. Seaborn will assign default colors.")

        # Filter the active_palette to include only display methods present in the current data
        if active_palette:
            present_display_methods = plot_df['display_method'].unique()
            active_palette = {k: v for k, v in active_palette.items() if k in present_display_methods}
            if not active_palette: # If filtering resulted in an empty palette
                print("Warning: Active palette became empty after filtering for present methods. Seaborn will assign default colors.")
                active_palette = None 
        
        # --- Create Box Plot ---
        plt.figure(figsize=(8, 7))
        order_x_axis = ['P', 'Q']

        # Adjust box width here to create more space between dodged boxes
        # Changed width from 0.6 to 0.5 for more spacing
        sns.boxplot(x='Type', y='difference_value', hue='display_method', data=plot_df,
                    palette=active_palette, order=order_x_axis, width=0.5, gap=0.2) 
        
        plt.title(f'Order Estimation Differences from True (P={real_p}, Q={real_q})\n{case_label}', fontsize=14)
        plt.xlabel('Order Parameter', fontsize=12)
        plt.ylabel('Estimated Order - True Order', fontsize=12)
        plt.axhline(0, color='grey', linestyle='--', linewidth=0.8)
        plt.legend(title='Method', loc='upper right') # Legend now uses display_method names
        plt.grid(axis='y', linestyle=':', alpha=0.7)
        plt.tight_layout()

        filename = f"order_differences_{case_label.replace(' ', '_').lower()}.rds_based.png"
        save_path = os.path.join(output_dir, filename)
        plt.savefig(save_path, dpi=300)
        plt.close()
        print(f"Difference boxplot saved successfully to: {save_path}")

    except Exception as e:
        print(f"An unexpected error occurred while processing {filepath} for difference boxplots: {e}")
        import traceback
        print(traceback.format_exc())

# %% --- Example Usage ---
# IMPORTANT: Remember to install pyreadr: pip install pyreadr

# --- Define your .rds file paths here ---
path_to_rds_no_covariates = r'C:\Users\Rafael\Desktop\msc-thesis\results\simulatedDataResults\ingarch_no_covariates_results.rds' 
path_to_rds_with_covariates = r'C:\Users\Rafael\Desktop\msc-thesis\results\simulatedDataResults\ingarch_with_covariates_results.rds' 

# --- Define Output Directory for Difference Boxplots ---
diff_plot_output_directory = "ingarch_order_differences_boxplots" 

# --- Define Real Model Orders (True P and Q) ---
REAL_P_ORDER = 2
REAL_Q_ORDER = 6

# --- Define Custom Method Names and Palettes (Optional) ---
# Example: map original names to new display names
user_defined_method_names = {
    "stepwise": "Stepwise",
    "grid_search": "Grid Search"
    # Add other original method names from your files if needed
}

# Example: define colors for these *new* display names
user_defined_palette = {
    "Stepwise": "mediumseagreen", 
    "Grid Search": "#377E7F" # A teal-like color
}


# --- Run the difference boxplot generation ---
print("--- Generating Difference Boxplots from .rds files ---")

# Process the "No Covariates" results file
print(f"\nProcessing 'No Covariates' data from: {path_to_rds_no_covariates}")
if os.path.exists(path_to_rds_no_covariates):
    create_difference_boxplots(path_to_rds_no_covariates, 
                               "No Covariates Data",
                               REAL_P_ORDER, REAL_Q_ORDER,
                               output_dir=diff_plot_output_directory,
                               method_rename_map=user_defined_method_names, # Pass the rename map
                               custom_palette=user_defined_palette)       # Pass the custom palette
else:
    print(f"File not found: {path_to_rds_no_covariates}")
    print("Please update 'path_to_rds_no_covariates' with the correct file path.")

# Process the "With Covariates" results file
print(f"\nProcessing 'With Covariates' data from: {path_to_rds_with_covariates}")
if os.path.exists(path_to_rds_with_covariates):
    create_difference_boxplots(path_to_rds_with_covariates, 
                               "With Covariates Data",
                               REAL_P_ORDER, REAL_Q_ORDER,
                               output_dir=diff_plot_output_directory,
                               method_rename_map=user_defined_method_names, # Pass the rename map
                               custom_palette=user_defined_palette)       # Pass the custom palette
else:
    print(f"File not found: {path_to_rds_with_covariates}")
    print("Please update 'path_to_rds_with_covariates' with the correct file path.")

# Example of plotting without custom names/palette (will use defaults or seaborn's choices)
# print(f"\nProcessing 'No Covariates' data (default names/colors) from: {path_to_rds_no_covariates}")
# if os.path.exists(path_to_rds_no_covariates):
#     create_difference_boxplots(path_to_rds_no_covariates, 
#                                "No Covariates Data (Default Names)",
#                                REAL_P_ORDER, REAL_Q_ORDER,
#                                output_dir=diff_plot_output_directory)
# else:
#     print(f"File not found: {path_to_rds_no_covariates}")


print("\nDifference boxplot generation complete.")