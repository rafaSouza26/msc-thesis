# %% Import Libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pyreadr # For reading .rds files
import traceback

# %% Plotting Function (AIC Difference Boxplots)
def create_aic_difference_boxplots(filepath,
                                   case_label,
                                   true_p_order,
                                   true_q_order,
                                   output_dir="ingarch_order_differences_boxplots", # Reverted to original default
                                   method_rename_map=None,
                                   custom_palette=None):
    """
    Creates and saves box plots of the AIC difference between models with true orders (p,q)
    and other models, from an .rds file.
    The reference AIC for true orders is calculated as the mean AIC for models
    with (true_p_order, true_q_order) within each estimation method.
    Differences are AIC(true_order_model_mean) - AIC(other_model).

    Args:
        filepath (str): The path to the .rds file.
        case_label (str): Label for the dataset (e.g., "No Covariates"). Used in title/filename.
        true_p_order (int): The true order p, used to identify reference models.
        true_q_order (int): The true order q, used to identify reference models.
        output_dir (str): The directory where plots will be saved.
        method_rename_map (dict, optional): Maps original method names to display names.
        custom_palette (dict, optional): Maps display names (after renaming) to colors.
    """
    try:
        os.makedirs(output_dir, exist_ok=True)
        print(f"Output directory for AIC difference boxplots set to: {output_dir}")
    except OSError as e:
        print(f"Error creating directory {output_dir}: {e}")
        return

    if not os.path.exists(filepath):
        print(f"Error: File not found at '{filepath}'. Please check the path.")
        return

    try:
        print(f"Attempting to read RDS file: {filepath}")
        try:
            rds_data = pyreadr.read_r(filepath)
        except Exception as read_error:
            print(f"Error reading RDS file '{filepath}' with pyreadr: {read_error}")
            print("Ensure 'pyreadr' is installed (pip install pyreadr) and the file path is correct.")
            print("If R is not installed or not found by pyreadr, you might need to install R and add it to your system's PATH.")
            return

        if not rds_data:
            print(f"Error: RDS file '{filepath}' is empty or pyreadr could not extract any R objects.")
            return

        df_key = list(rds_data.keys())[0]
        df = rds_data[df_key]
        print(f"Successfully read object '{df_key}' from RDS file: {os.path.basename(filepath)} with {df.shape[0]} rows and {df.shape[1]} columns.")

        # --- Column Name Handling (Focus on method, p, q, aic) ---
        expected_cols_map = {
            'method': ['method', 'Method'],
            'p': ['p', 'P'],
            'q': ['q', 'Q'],
            'aic': ['aic', 'AIC', 'aic_value', 'aicvalue'] # Added AIC variants
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
                print(f"  Available columns in RDS object: {current_columns_in_df}")
                return

        df.rename(columns=rename_operation_map, inplace=True)
        print(f"Columns after attempting initial standardization: {df.columns.tolist()}")

        required_standard_cols = ['method', 'p', 'q', 'aic']
        for req_col in required_standard_cols:
            if req_col not in df.columns:
                print(f"Error: Standardized column '{req_col}' is missing after renaming attempts.")
                print(f"  Available columns: {df.columns.tolist()}")
                return

        df['method'] = df['method'].astype(str).str.lower().str.strip()

        if method_rename_map:
            df['display_method'] = df['method'].map(method_rename_map).fillna(df['method'])
            print(f"Applied method rename map. Unique display methods: {df['display_method'].unique()}")
        else:
            df['display_method'] = df['method']
            print(f"No method rename map provided. Using original standardized methods: {df['display_method'].unique()}")

        for col_to_clean in ['p', 'q', 'aic']:
            df[col_to_clean] = pd.to_numeric(df[col_to_clean], errors='coerce')

        df.dropna(subset=['display_method', 'p', 'q', 'aic'], inplace=True)
        if df.empty:
            print(f"Warning: No valid data remaining after cleaning for file '{filepath}'. Skipping AIC difference boxplot.")
            return

        df[['p', 'q']] = df[['p', 'q']].astype(int)

        plot_data_rows = []
        unique_display_methods_in_data = sorted(df['display_method'].unique())

        for method_name in unique_display_methods_in_data:
            method_df = df[df['display_method'] == method_name].copy()
            true_order_models_for_method = method_df[
                (method_df['p'] == true_p_order) & (method_df['q'] == true_q_order)
            ]

            if true_order_models_for_method.empty or true_order_models_for_method['aic'].isnull().all():
                print(f"Warning: No valid reference models with orders (P={true_p_order}, Q={true_q_order}) found for method '{method_name}'. This method will be skipped for AIC diff.")
                continue

            mean_reference_aic = true_order_models_for_method['aic'].mean()
            if pd.isna(mean_reference_aic):
                print(f"Warning: Mean reference AIC is NaN for method '{method_name}' with orders (P={true_p_order}, Q={true_q_order}). This method will be skipped.")
                continue

            print(f"Method '{method_name}': Mean AIC for (P={true_p_order},Q={true_q_order}) models = {mean_reference_aic:.2f}")

            other_models_for_method = method_df[
                ~((method_df['p'] == true_p_order) & (method_df['q'] == true_q_order))
            ]

            if other_models_for_method.empty:
                print(f"Info: No 'other' models (P,Q different from ({true_p_order},{true_q_order})) to compare against for method '{method_name}'.")
                # If you want to plot methods that ONLY have reference models, you might add placeholder differences (e.g., zero)
                # or ensure the plotting logic can handle cases where a category in 'plot_df' might be empty or have specific values.
                # For now, we stick to only adding rows if 'other_models_for_method' exist to calculate differences.
                continue


            for _, row in other_models_for_method.iterrows():
                other_model_aic = row['aic']
                if not pd.isna(other_model_aic):
                    aic_diff = mean_reference_aic - other_model_aic
                    plot_data_rows.append({
                        'display_method': method_name,
                        'aic_difference': aic_diff
                    })

        if not plot_data_rows:
            print(f"Warning: No data available for plotting AIC differences for file '{filepath}'.")
            # This could happen if all methods were skipped, or no 'other' models were found.
            return

        plot_df = pd.DataFrame(plot_data_rows)
        # Ensure plot_order uses only methods that are actually in plot_df after processing
        plot_order = sorted(plot_df['display_method'].unique())

        if not plot_order: # Double check if plot_df ended up empty or without valid methods
            print(f"Warning: No methods have data for AIC difference plotting after processing '{filepath}'. Skipping plot generation.")
            return


        active_palette = None
        if custom_palette:
            active_palette = {k: v for k, v in custom_palette.items() if k in plot_order}
            print("Using user-provided custom_palette, filtered for methods present in plot.")
            if len(active_palette) < len(plot_order):
                print("  Not all methods in the plot have a color in custom_palette. Seaborn will use default colors for missing ones.")
        elif not method_rename_map:
            default_internal_palette = {
                "grid_search": "mediumseagreen",
                "stepwise": "firebrick"
            }
            active_palette = {k: v for k, v in default_internal_palette.items() if k in plot_order}
            print("Using default internal palette for original method names, filtered for methods present in plot.")
            if len(active_palette) < len(plot_order) and plot_order: # Check plot_order is not empty
                print("  Not all methods in the plot have a color in default_internal_palette. Seaborn will use default colors for missing ones.")
        else:
             print("Methods renamed, but no custom_palette provided. Seaborn will assign default colors for all methods.")
             active_palette = None

        if active_palette and not active_palette :
            print("Warning: Active palette became empty after filtering. Seaborn will assign default colors.")
            active_palette = None


        fig_width = max(8, len(plot_order) * 1.2 if len(plot_order) > 1 else 8) # Ensure min width even for one box
        plt.figure(figsize=(fig_width, 7))

        ax = sns.boxplot(x='display_method', y='aic_difference', hue='display_method',
                         data=plot_df, palette=active_palette, order=plot_order,
                         width=0.6, dodge=False, legend=False)


        # --- Start of modifications for average dot and text ---
        for i, method_name_in_plot in enumerate(plot_order): # Iterate using the actual order of boxes
            # Filter data for the current box being plotted
            method_data = plot_df[plot_df['display_method'] == method_name_in_plot]['aic_difference']
            if not method_data.empty:
                mean_val = method_data.mean()
                # Plot the mean as a red dot
                # The x-coordinate for plot() within a boxplot context is the categorical index (0, 1, 2...)
                ax.plot(i, mean_val, 'ro', markersize=8, zorder=5) # 'ro' is red circle. zorder to be on top.
                # Add text for the mean value
                # Adjust text position slightly for better visibility (e.g., offset from the dot)
                ax.text(i + 0.05, mean_val, f'{mean_val:.2f}', color='red',
                        ha='left', va='center', fontsize=9, fontweight='bold')
        # --- End of modifications ---

        plt.title(f'AIC Difference relative to Model (P={true_p_order}, Q={true_q_order})\n{case_label}', fontsize=14)
        plt.xlabel('Method', fontsize=12)
        plt.ylabel(f'AIC(P={true_p_order},Q={true_q_order}) - AIC(Other Model)', fontsize=12)
        plt.axhline(0, color='grey', linestyle='--', linewidth=0.8)

        # Legend handling: Since hue=x, a default legend might be redundant or show multiple entries for colors.
        # We can create a custom legend if needed, or rely on x-axis labels if colors are distinct per box.
        # If custom_palette is provided and maps to display_method, it should color boxes accordingly.
        # If legend=False was set in sns.boxplot, this part might be optional or need specific handles.
        if active_palette : # Only add legend if a palette was used that might need explaining
            handles = [plt.Rectangle((0,0),1,1, color=active_palette.get(label, 'gray')) for label in plot_order if label in active_palette]
            labels = [label for label in plot_order if label in active_palette]
            if handles: # ensure there are handles to create a legend
                 ax.legend(handles, labels, title='Method', loc='upper right')
        elif ax.get_legend() is not None: # Fallback if sns created one
            ax.legend(title='Method', loc='upper right')


        plt.grid(axis='y', linestyle=':', alpha=0.7)
        if len(plot_order) > 4:
            plt.xticks(rotation=45, ha='right')
        else:
            plt.xticks(rotation=0)
        plt.tight_layout()

        filename = f"aic_differences_{case_label.replace(' ', '_').lower()}_p{true_p_order}q{true_q_order}.rds_based.png"
        save_path = os.path.join(output_dir, filename)
        plt.savefig(save_path, dpi=300)
        plt.close()
        print(f"AIC difference boxplot saved successfully to: {save_path}")

    except Exception as e:
        print(f"An unexpected error occurred while processing {filepath} for AIC difference boxplots: {e}")
        print(traceback.format_exc())

# %% --- Example Usage ---
# IMPORTANT: Remember to install pyreadr: pip install pyreadr pandas matplotlib seaborn

# --- Define your .rds file paths here ---
# Please ensure these paths are correct on your system.
# For demonstration, I'll use placeholder names. Replace with your actual paths.
path_to_rds_no_covariates = r'C:\Users\Rafael\Desktop\msc-thesis\results\simulatedDataResults\ingarch_no_covariates_results.rds'
path_to_rds_with_covariates = r'C:\Users\Rafael\Desktop\msc-thesis\results\simulatedDataResults\ingarch_with_covariates_results.rds'


# --- Define Output Directory for AIC Difference Boxplots ---
plot_output_directory = "ingarch_order_differences_boxplots_with_avg"

# --- Define True Model Orders (P and Q for reference AIC calculation) ---
TRUE_P_FOR_REF_AIC = 2
TRUE_Q_FOR_REF_AIC = 6

# --- Define Custom Method Names and Palettes (Optional) ---
user_defined_method_names = {
    "stepwise": "Stepwise",
    "grid_search": "Grid Search"
}

user_defined_palette = {
    "Stepwise": "mediumseagreen",
    "Grid Search": "#377E7F"
}


print("--- Generating AIC Difference Boxplots from .rds files ---")

# Check if files exist before attempting to process
# (This is good practice but the function itself also checks)

print(f"\nProcessing 'No Covariates' data from: {path_to_rds_no_covariates}")
if os.path.exists(path_to_rds_no_covariates):
    create_aic_difference_boxplots(filepath=path_to_rds_no_covariates,
                                   case_label="No Covariates Data",
                                   true_p_order=TRUE_P_FOR_REF_AIC,
                                   true_q_order=TRUE_Q_FOR_REF_AIC,
                                   output_dir=plot_output_directory,
                                   method_rename_map=user_defined_method_names,
                                   custom_palette=user_defined_palette)
else:
    print(f"File not found: {path_to_rds_no_covariates}")
    print(f"Please update 'path_to_rds_no_covariates' with the correct file path.")

print(f"\nProcessing 'With Covariates' data from: {path_to_rds_with_covariates}")
if os.path.exists(path_to_rds_with_covariates):
    create_aic_difference_boxplots(filepath=path_to_rds_with_covariates,
                                   case_label="With Covariates Data",
                                   true_p_order=TRUE_P_FOR_REF_AIC,
                                   true_q_order=TRUE_Q_FOR_REF_AIC,
                                   output_dir=plot_output_directory,
                                   method_rename_map=user_defined_method_names,
                                   custom_palette=user_defined_palette)
else:
    print(f"File not found: {path_to_rds_with_covariates}")
    print(f"Please update 'path_to_rds_with_covariates' with the correct file path.")

print("\nAIC difference boxplot generation complete.")
print(f"Check the '{plot_output_directory}' directory for the plots if processing was successful.")