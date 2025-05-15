# %% Import Libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pyreadr # For reading .rds files

# %% Plotting Function (Execution Time Boxplots)
def create_time_boxplots(filepath, 
                         case_label, 
                         output_dir="ingarch_plots_time", 
                         method_rename_map=None, 
                         custom_palette=None):
    """
    Creates and saves box plots of execution times for different methods
    from an .rds file.

    Args:
        filepath (str): The path to the .rds file.
        case_label (str): Label for the dataset (e.g., "No Covariates"). Used in title/filename.
        output_dir (str): The directory where plots will be saved.
        method_rename_map (dict, optional): Maps original method names to display names.
                                            e.g., {"stepwise": "Stepwise Search"}.
                                            Defaults to None (use original names).
        custom_palette (dict, optional): Maps display names (after renaming) to colors.
                                         e.g., {"Stepwise Search": "blue"}.
                                         Defaults to None.
    """
    # Create the output directory if it doesn't exist
    try:
        os.makedirs(output_dir, exist_ok=True)
        print(f"Output directory for time boxplots set to: {output_dir}")
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

        # --- Column Name Handling (Focus on method and time) ---
        expected_cols_map = {
            'method': ['method', 'Method'], # Standard name: list of possible names in file
            'time': ['time', 'Time', 'execution_time', 'ExecutionTime'] # Standard name for time
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

        required_standard_cols = ['method', 'time']
        for req_col in required_standard_cols:
            if req_col not in df.columns:
                print(f"Error: Standardized column '{req_col}' is missing after renaming attempts.")
                print(f"   Available columns: {df.columns.tolist()}")
                return

        # Standardize 'method' column content (e.g., to lowercase)
        df['method'] = df['method'].astype(str).str.lower().str.strip()

        # --- Apply custom method renaming ---
        if method_rename_map:
            df['display_method'] = df['method'].map(method_rename_map).fillna(df['method'])
            print(f"Applied method rename map. Unique display methods: {df['display_method'].unique()}")
        else:
            df['display_method'] = df['method']
            print(f"No method rename map provided. Using original standardized methods: {df['display_method'].unique()}")

        # --- Data Cleaning (time) ---
        df['time'] = pd.to_numeric(df['time'], errors='coerce')
        df.dropna(subset=['display_method', 'time'], inplace=True)
        if df.empty:
            print(f"Warning: No valid data remaining after cleaning for file '{filepath}'. Skipping time boxplot.")
            return
        
        # Data for plotting is in 'df' with 'display_method' and 'time' columns

        # --- Determine Palette for Plotting ---
        active_palette = None
        if custom_palette:
            active_palette = custom_palette
            print("Using user-provided custom_palette.")
        elif not method_rename_map: 
            default_internal_palette = {
                "grid_search": "teal", # Adjusted default color
                "stepwise": "mediumseagreen"  # Adjusted default color
            }
            active_palette = default_internal_palette
            print("Using default internal palette for original method names.")
        else:
            print("Methods renamed, but no custom_palette provided. Seaborn will assign default colors.")

        if active_palette:
            present_display_methods = df['display_method'].unique()
            active_palette = {k: v for k, v in active_palette.items() if k in present_display_methods}
            if not active_palette:
                print("Warning: Active palette became empty after filtering. Seaborn will assign default colors.")
                active_palette = None 
        
        # --- Create Box Plot ---
        plt.figure(figsize=(8, 7)) # Adjust figure size as needed
        
        # Order for x-axis categories (methods)
        # If you want a specific order, define it here, e.g., based on user_defined_method_names.values()
        # For now, it will use the order of appearance or sorted order by default from Seaborn.
        # If you have user_defined_method_names, you might want to control the order:
        plot_order = None
        if method_rename_map:
             # Attempt to order based on the values in the rename map,
             # filtering by those actually present in the data.
             ordered_display_names = [name for name in method_rename_map.values() if name in df['display_method'].unique()]
             if len(ordered_display_names) == len(df['display_method'].unique()): # Ensure all present methods are in the map
                 plot_order = ordered_display_names

        # Create the box plot
        ax = sns.boxplot(x='display_method', y='time', data=df,
                    palette=active_palette, order=plot_order, width=0.5) # Using display_method for x-axis
        
        # Calculate means for each method and add red dots for them
        if plot_order:
            # If a specific order was defined and used
            for i, method in enumerate(plot_order):
                mean_time = df[df['display_method'] == method]['time'].mean()
                ax.plot(i, mean_time, 'ro', ms=8)  # Add red dot at mean
                
                # Add text label at bottom left of each box plot
                # Get the y-min for this group's data
                min_y = df[df['display_method'] == method]['time'].min()
                # Position text at the bottom left of the box
                ax.text(i - 0.4, min_y - (df['time'].max() - df['time'].min()) * 0.08,
                        f'Mean: {mean_time:.2f}', 
                        ha='left', va='top', fontsize=9, color='red')
        else:
            # If using default order
            for i, method in enumerate(df['display_method'].unique()):
                mean_time = df[df['display_method'] == method]['time'].mean()
                ax.plot(i, mean_time, 'ro', ms=8)  # Add red dot at mean
                
                # Add text label at bottom left of each box plot
                # Get the y-min for this group's data
                min_y = df[df['display_method'] == method]['time'].min()
                # Position text at the bottom left of the box
                ax.text(i - 0.4, min_y - (df['time'].max() - df['time'].min()) * 0.08,
                        f'Mean: {mean_time:.2f}', 
                        ha='left', va='top', fontsize=9, color='red')
        
        plt.title(f'Comparison of Execution Times\n{case_label}', fontsize=14)
        plt.xlabel('Search Method', fontsize=12)
        plt.ylabel('Time (milliseconds)', fontsize=12) # As per example image
        plt.grid(axis='y', linestyle='--', alpha=0.7) # Dashed grid lines
        plt.tight_layout()

        filename = f"execution_time_comparison_{case_label.replace(' ', '_').lower()}.png"
        save_path = os.path.join(output_dir, filename)
        plt.savefig(save_path, dpi=300)
        plt.close()
        print(f"Execution time boxplot saved successfully to: {save_path}")

    except Exception as e:
        print(f"An unexpected error occurred while processing {filepath} for time boxplots: {e}")
        import traceback
        print(traceback.format_exc())

# %% --- Example Usage ---
# IMPORTANT: Remember to install pyreadr: pip install pyreadr

# --- Define your .rds file paths here ---
path_to_rds_no_covariates = r'C:\Users\Rafael\Desktop\msc-thesis\results\simulatedDataResults\ingarch_no_covariates_results.rds' 
path_to_rds_with_covariates = r'C:\Users\Rafael\Desktop\msc-thesis\results\simulatedDataResults\ingarch_with_covariates_results.rds' 

# --- Define Output Directory for Time Boxplots ---
time_plot_output_directory = "ingarch_execution_times_boxplots" 

# --- Define Custom Method Names and Palettes (Optional) ---
# These will be used for the x-axis labels and box colors
user_defined_method_names = {
    "stepwise": "Stepwise",  # Changed to match example image
    "grid_search": "Grid"  # Changed to match example image
}

# Example: define colors for these *new* display names
# Using colors closer to the example image (teal/greenish)
user_defined_palette = {
    "Stepwise": "mediumseagreen", 
    "Grid": "#377E7F" # A teal-like color
}


# --- Run the execution time boxplot generation ---
print("--- Generating Execution Time Boxplots from .rds files ---")

# Process the "No Covariates" results file
print(f"\nProcessing 'No Covariates' data from: {path_to_rds_no_covariates}")
if os.path.exists(path_to_rds_no_covariates):
    create_time_boxplots(path_to_rds_no_covariates, 
                         "No Covariates Data",
                         output_dir=time_plot_output_directory,
                         method_rename_map=user_defined_method_names,
                         custom_palette=user_defined_palette)
else:
    print(f"File not found: {path_to_rds_no_covariates}")
    print("Please update 'path_to_rds_no_covariates' with the correct file path.")

# Process the "With Covariates" results file
print(f"\nProcessing 'With Covariates' data from: {path_to_rds_with_covariates}")
if os.path.exists(path_to_rds_with_covariates):
    create_time_boxplots(path_to_rds_with_covariates, 
                         "With Covariates Data",
                         output_dir=time_plot_output_directory,
                         method_rename_map=user_defined_method_names,
                         custom_palette=user_defined_palette)
else:
    print(f"File not found: {path_to_rds_with_covariates}")
    print("Please update 'path_to_rds_with_covariates' with the correct file path.")

print("\nExecution time boxplot generation complete.")