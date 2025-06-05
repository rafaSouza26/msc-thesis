# %% Import Libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
# import traceback # traceback is already imported in your provided script's try-except block if needed there

# %% Plotting Function (Execution Time Boxplots)
def create_time_boxplots(filepath,
                         case_label,
                         output_dir="ingarch_plots_time",
                         method_rename_map=None,
                         custom_palette=None,
                         show_mean_value=False): # Flag to show mean's numerical value
    """
    Creates and saves box plots of execution times for different methods
    from a .csv file.

    Args:
        filepath (str): The path to the .csv file.
        case_label (str): Label for the dataset (e.g., "No Covariates"). Used in title/filename.
        output_dir (str): The directory where plots will be saved.
        method_rename_map (dict, optional): Maps original method names to display names.
                                         e.g., {"stepwise": "Stepwise Search"}.
                                         Defaults to None (use original names).
        custom_palette (dict, optional): Maps display names (after renaming) to colors.
                                        e.g., {"Stepwise Search": "blue"}.
                                        Defaults to None.
        show_mean_value (bool, optional): If True, displays the numerical mean value as red text
                                          styled like "Mean: value" next to each mean marker.
                                          Defaults to False.
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
        # --- Read .csv file ---
        print(f"Attempting to read CSV file: {filepath}")
        df = pd.read_csv(filepath)
        if df.empty:
            print(f"Error: CSV file '{filepath}' is empty or pandas could not read it properly.")
            return
        print(f"Successfully read CSV file: {os.path.basename(filepath)} with {df.shape[0]} rows and {df.shape[1]} columns.")

        # --- Column Name Handling (Focus on method and time) ---
        expected_cols_map = {
            'method': ['method', 'Method'],
            'time': ['time', 'Time', 'execution_time', 'ExecutionTime']
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
                print(f"Error: Required column for '{standard_col_name}' (e.g., {possible_file_names}) not found in CSV file '{filepath}'.")
                print(f"   Available columns in CSV file: {current_columns_in_df}")
                return

        df.rename(columns=rename_operation_map, inplace=True)
        print(f"Columns after attempting initial standardization: {df.columns.tolist()}")

        required_standard_cols = ['method', 'time']
        for req_col in required_standard_cols:
            if req_col not in df.columns:
                print(f"Error: Standardized column '{req_col}' is missing after renaming attempts.")
                print(f"   Available columns: {df.columns.tolist()}")
                return

        df['method'] = df['method'].astype(str).str.lower().str.strip()

        if method_rename_map:
            df['display_method'] = df['method'].map(method_rename_map).fillna(df['method'])
            print(f"Applied method rename map. Unique display methods: {df['display_method'].unique()}")
        else:
            df['display_method'] = df['method']
            print(f"No method rename map provided. Using original standardized methods: {df['display_method'].unique()}")

        df['time'] = pd.to_numeric(df['time'], errors='coerce')
        df.dropna(subset=['display_method', 'time'], inplace=True)
        if df.empty:
            print(f"Warning: No valid data remaining after cleaning for file '{filepath}'. Skipping time boxplot.")
            return

        active_palette = None
        if custom_palette:
            active_palette = custom_palette
            print("Using user-provided custom_palette.")
        elif not method_rename_map:
            default_internal_palette = {"grid_search": "teal", "stepwise": "mediumseagreen"}
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

        plt.figure(figsize=(8, 7))
        plot_order = None
        if method_rename_map:
             ordered_display_names = [name for name in method_rename_map.values() if name in df['display_method'].unique()]
             if len(ordered_display_names) == len(df['display_method'].unique()): # Ensure all present methods are in the map
                 plot_order = ordered_display_names

        ax = sns.boxplot(x='display_method', y='time', data=df,
                         palette=active_palette, order=plot_order, width=0.5)

        plot_methods_in_order = [tick.get_text() for tick in ax.get_xticklabels()]
        legend_handles = []

        mean_points_x = range(len(plot_methods_in_order))
        mean_points_y = []
        for method_name in plot_methods_in_order:
            mean_val = df[df['display_method'] == method_name]['time'].mean()
            mean_points_y.append(mean_val)

        if any(pd.notnull(y) for y in mean_points_y): # Check if there are any valid mean values to plot
            # Plot the mean markers (red dots)
            mean_handle, = ax.plot(mean_points_x, mean_points_y,
                                   linestyle='', marker='o', color='red', ms=8, 
                                   label='Mean' if not show_mean_value else None) # Only label for legend if text isn't shown
            
            if not show_mean_value: # Add to legend handles only if text annotation is not shown
                 legend_handles.append(mean_handle)

            if show_mean_value:
                for i, mean_val in enumerate(mean_points_y):
                    if pd.notnull(mean_val):
                        # MODIFIED PART for red text "Mean: value"
                        ax.annotate(f'Mean: {mean_val:.2f}',       # Text format
                                    xy=(mean_points_x[i], mean_val),# Point to annotate
                                    xytext=(8, -6),                  # Offset text (8 points right, 0 points vertical)
                                    textcoords='offset points',     # Specify offset in points
                                    color='red',                    # Text color red
                                    ha='left',                      # Horizontal alignment
                                    va='center',                    # Vertical alignment
                                    fontsize=9)                     # Font size
                        # END OF MODIFIED PART
        
        if legend_handles: # Only show legend if there are specific handles (e.g. Mean dot when not annotated by text)
            ax.legend(handles=legend_handles, loc='best')


        plt.title(f'Comparison of Execution Times\n{case_label}', fontsize=14)
        plt.xlabel('Search Method', fontsize=12)
        plt.ylabel('Time (milliseconds)', fontsize=12)
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.tight_layout()

        filename = f"execution_time_comparison_{case_label.replace(' ', '_').lower()}.png"
        save_path = os.path.join(output_dir, filename)
        plt.savefig(save_path, dpi=300)
        plt.close()
        print(f"Execution time boxplot saved successfully to: {save_path}")

    except Exception as e:
        print(f"An unexpected error occurred while processing {filepath} for time boxplots: {e}")
        import traceback # Ensure traceback is available here
        print(traceback.format_exc())

# %% --- Example Usage ---
# IMPORTANT: Ensure your CSV files have columns like 'method' (or 'Method') and 'time' (or 'Time', 'execution_time').

# --- Define your .csv file paths here ---
# Replace these with the actual paths to your CSV files
path_to_csv_no_covariates_m1 = r'C:\Users\Rafael\Desktop\msc-thesis\Results\Simulation\no_covariates\Model1\ingarch_no_covariates_results.csv'
path_to_csv_no_covariates_m2 = r'C:\Users\Rafael\Desktop\msc-thesis\Results\Simulation\no_covariates\Model2\ingarch_no_covariates_results.csv'
path_to_csv_with_covariates = r'C:\Users\Rafael\Desktop\msc-thesis\Results\Simulation\with_covariates\Combined_cpu3\ingarch_with_covariates_results.csv'

# --- Define Output Directory for Time Boxplots ---
time_plot_output_directory = "ingarch_execution_times_boxplots"

# --- Define Custom Method Names and Palettes (Optional) ---
user_defined_method_names = {
    "stepwise": "Stepwise",
    "grid_search": "Grid"
}

# Palette to match the style of the image provided earlier
user_defined_palette = {
    "Stepwise": "mediumseagreen", 
    "Grid": "#377E7F" # Darker teal/blue-green
}

# --- Run the execution time boxplot generation ---
print("--- Generating Execution Time Boxplots from .csv files ---")

# Process the "No Covariates" m1 results file, showing mean values
print(f"\nProcessing 'No Covariates' data from: {path_to_csv_no_covariates_m1}")
if os.path.exists(path_to_csv_no_covariates_m1):
    create_time_boxplots(path_to_csv_no_covariates_m1,
                         "No Covariates Data M1",
                         output_dir=time_plot_output_directory,
                         method_rename_map=user_defined_method_names,
                         custom_palette=user_defined_palette,
                         show_mean_value=False) # Show styled mean numerical value for this plot
else:
    print(f"File not found: {path_to_csv_no_covariates_m1}")
    print("Please update 'path_to_csv_no_covariates' with the correct file path.")

# Process the "No Covariates" m2 results file, showing mean values
print(f"\nProcessing 'No Covariates' data from: {path_to_csv_no_covariates_m2}")
if os.path.exists(path_to_csv_no_covariates_m2):
    create_time_boxplots(path_to_csv_no_covariates_m2,
                         "No Covariates Data M2",
                         output_dir=time_plot_output_directory,
                         method_rename_map=user_defined_method_names,
                         custom_palette=user_defined_palette,
                         show_mean_value=False) # Show styled mean numerical value for this plot
else:
    print(f"File not found: {path_to_csv_no_covariates_m2}")
    print("Please update 'path_to_csv_no_covariates' with the correct file path.")

# Process the "With Covariates" results file, also showing mean values
print(f"\nProcessing 'With Covariates' data from: {path_to_csv_with_covariates}")
if os.path.exists(path_to_csv_with_covariates):
    create_time_boxplots(path_to_csv_with_covariates,
                         "With Covariates Data",
                         output_dir=time_plot_output_directory,
                         method_rename_map=user_defined_method_names,
                         custom_palette=user_defined_palette,
                         show_mean_value=False) # Show styled mean numerical value for this plot
else:
    print(f"File not found: {path_to_csv_with_covariates}")
    print("Please update 'path_to_csv_with_covariates' with the correct file path.")


print("\nExecution time boxplot generation complete.")