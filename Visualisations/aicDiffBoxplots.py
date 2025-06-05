# %% Import Libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import traceback
from matplotlib.lines import Line2D # For custom legend

# %% Plotting Function (AIC Difference Boxplots)
def create_aic_difference_boxplots(filepath,
                                   case_label, # This will be the specific title part
                                   true_p_order, # Will be set to 2 for reference
                                   true_q_order, # Will be set to 6 for reference
                                   output_dir="ingarch_order_differences_boxplots",
                                   method_rename_map=None,
                                   custom_palette=None,
                                   show_mean_value=False):
    """
    Creates and saves box plots of the AIC difference between models with specified reference orders (p,q)
    and other models, from a .csv file.
    The reference AIC is taken from the first model instance found with (true_p_order, true_q_order)
    within each estimation method.
    Differences are AIC(reference_order_model_instance) - AIC(other_model).

    Args:
        filepath (str): The path to the .csv file.
        case_label (str): Custom label for the dataset (e.g., "Specific Analysis XYZ"). Used in title/filename.
        true_p_order (int): The reference order p (e.g., 2).
        true_q_order (int): The reference order q (e.g., 6).
        output_dir (str): The directory where plots will be saved.
        method_rename_map (dict, optional): Maps original method names to display names.
        custom_palette (dict, optional): Maps display names (after renaming) to colors for boxplots.
        show_mean_value (bool, optional): If True, displays the numerical mean value as text
                                         above each mean marker. Defaults to False.
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

    df = None
    file_extension = os.path.splitext(filepath)[1].lower()

    if file_extension != '.csv':
        print(f"Error: Unsupported file type '{file_extension}'. This script only supports .csv files.")
        return

    print(f"Attempting to read CSV file: {filepath}")
    try:
        df = pd.read_csv(filepath)
        print(f"Successfully read CSV file: {os.path.basename(filepath)} with {df.shape[0]} rows and {df.shape[1]} columns.")
    except Exception as read_error:
        print(f"Error reading CSV file '{filepath}': {read_error}")
        return

    if df is None or df.empty:
        print(f"Error: No data loaded from file '{filepath}'. The DataFrame is empty.")
        return

    try:
        # --- Column Name Handling (Focus on method, p, q, aic) ---
        standard_names = {'method': 'method', 'p': 'p', 'q': 'q', 'aic': 'aic'}
        possible_col_names_map = {
            'method': standard_names['method'], 'Method': standard_names['method'],
            'p': standard_names['p'], 'P': standard_names['p'], 'p_order': standard_names['p'],
            'q': standard_names['q'], 'Q': standard_names['q'], 'q_order': standard_names['q'],
            'aic': standard_names['aic'], 'AIC': standard_names['aic'],
            'aic_value': standard_names['aic'], 'aicvalue': standard_names['aic']
        }

        rename_map = {}
        current_columns_lower = {col.lower(): col for col in df.columns}

        for standard_name_key, standard_name_val in standard_names.items():
            found = False
            if standard_name_val.lower() in current_columns_lower:
                original_case_col_name = current_columns_lower[standard_name_val.lower()]
                if original_case_col_name != standard_name_val:
                    rename_map[original_case_col_name] = standard_name_val
                found = True
            if not found:
                for possible_name_variant_lower, mapped_standard_name in possible_col_names_map.items():
                    if mapped_standard_name == standard_name_val:
                        if possible_name_variant_lower in current_columns_lower:
                            original_case_col_name = current_columns_lower[possible_name_variant_lower]
                            if original_case_col_name != standard_name_val:
                                rename_map[original_case_col_name] = standard_name_val
                            found = True
                            break
            if not found:
                possible_options_for_error = [k for k,v in possible_col_names_map.items() if v == standard_name_val]
                if standard_name_val not in possible_options_for_error:
                    possible_options_for_error.insert(0, standard_name_val)
                print(f"Error: Required column for '{standard_name_val}' (e.g., {possible_options_for_error}) not found in file '{os.path.basename(filepath)}'.")
                print(f"  Available columns: {df.columns.tolist()}")
                return

        df.rename(columns=rename_map, inplace=True)
        print(f"Columns after attempting initial standardization: {df.columns.tolist()}")

        required_standard_cols = list(standard_names.values())
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
        unique_display_methods_in_data_df = sorted(df['display_method'].unique())

        for method_name in unique_display_methods_in_data_df:
            method_df = df[df['display_method'] == method_name].copy()
            
            reference_order_models_for_method = method_df[
                (method_df['p'] == true_p_order) & (method_df['q'] == true_q_order)
            ]

            if reference_order_models_for_method.empty:
                print(f"Info: No reference models with orders (P={true_p_order}, Q={true_q_order}) found for method '{method_name}'. This method will be skipped for AIC diff.")
                continue

            first_model_aic = reference_order_models_for_method['aic'].iloc[0]

            if pd.isna(first_model_aic):
                print(f"Warning: The first reference model (P={true_p_order}, Q={true_q_order}) found for method '{method_name}' has a NaN AIC value. This method will be skipped.")
                continue
            
            reference_aic_value = first_model_aic 

            print_msg_prefix = f"Method '{method_name}': Using reference AIC (P={true_p_order},Q={true_q_order})"
            
            if len(reference_order_models_for_method) > 1:
                unique_aics_in_ref_models = reference_order_models_for_method['aic'].unique()
                valid_unique_aics = [aic for aic in unique_aics_in_ref_models if pd.notna(aic)]

                if len(valid_unique_aics) > 1:
                    all_aics_list = reference_order_models_for_method['aic'].tolist()
                    print(f"{print_msg_prefix} = {reference_aic_value:.2f} (from the first instance).")
                    print(f"  Warning: Multiple reference models (P={true_p_order}, Q={true_q_order}) found for method '{method_name}' with *differing* valid AICs: {valid_unique_aics}. All AICs found: {all_aics_list}.")
                elif not valid_unique_aics or valid_unique_aics[0] != reference_aic_value:
                    all_aics_list = reference_order_models_for_method['aic'].tolist()
                    print(f"{print_msg_prefix} = {reference_aic_value:.2f} (from the first valid instance).")
                    print(f"  Note: Additional reference model instances for (P={true_p_order}, Q={true_q_order}) for method '{method_name}' may have non-matching or NaN AICs. All AICs found: {all_aics_list}.")
                else: 
                    print(f"{print_msg_prefix} = {reference_aic_value:.2f} (multiple instances found, all with this AIC).")
            else:
                print(f"{print_msg_prefix} = {reference_aic_value:.2f} (single instance found).")

            other_models_for_method = method_df[
                ~((method_df['p'] == true_p_order) & (method_df['q'] == true_q_order))
            ]

            if other_models_for_method.empty:
                print(f"Info: No 'other' models (P,Q different from ({true_p_order},{true_q_order})) to compare against for method '{method_name}'.")
                continue

            for _, row in other_models_for_method.iterrows():
                other_model_aic = row['aic']
                if not pd.isna(other_model_aic):
                    aic_diff = reference_aic_value - other_model_aic 
                    plot_data_rows.append({
                        'display_method': method_name,
                        'aic_difference': aic_diff
                    })

        if not plot_data_rows:
            print(f"Warning: No data available for plotting AIC differences for file '{filepath}' (perhaps only reference models existed or no other models had valid AICs).")
            return

        plot_df = pd.DataFrame(plot_data_rows)
        
        if plot_df.empty or 'display_method' not in plot_df.columns or plot_df['display_method'].nunique() == 0:
            print(f"Warning: No data available in plot_df for plotting AIC differences for file '{filepath}'.")
            return

        all_methods_in_plot_data = sorted(plot_df['display_method'].unique())

        # Define the desired primary order for display methods
        preferred_order_list = ["Stepwise", "Grid"] 
        
        plot_order = []
        # Add preferred methods first, in the specified order, if they exist in the data
        for method_name in preferred_order_list:
            if method_name in all_methods_in_plot_data:
                plot_order.append(method_name)
        
        # Add any other methods from the data that are not in the preferred list,
        # maintaining their original sorted order relative to each other.
        for method_name in all_methods_in_plot_data:
            if method_name not in plot_order: # Appends remaining sorted methods
                plot_order.append(method_name)

        if not plot_order:
            print(f"Warning: No methods determined for plot order for '{filepath}'. Skipping plot generation.")
            return

        active_palette = None
        if custom_palette:
            active_palette = {k: v for k, v in custom_palette.items() if k in plot_order}
            print("Using user-provided custom_palette, filtered for methods present in plot.")
            if len(active_palette) < len(plot_order):
                print("  Not all methods in the plot have a color in custom_palette. Seaborn will use default colors for missing ones.")
        elif not method_rename_map: 
            default_internal_palette = {
                "grid_search": "mediumseagreen", # Raw method name
                "stepwise": "firebrick"          # Raw method name
            }
            # If no rename map, display_method IS the raw method name.
            # plot_order contains display_method names.
            active_palette = {k: v for k, v in default_internal_palette.items() if k in plot_order}
            print("Using default internal palette for original method names, filtered for methods present in plot.")
            if active_palette and len(active_palette) < len(plot_order) and plot_order:
                print("  Not all methods in the plot have a color in default_internal_palette. Seaborn will use default colors for missing ones.")
        else:
            # This case implies method_rename_map IS provided, but custom_palette is NOT.
            print("Seaborn will assign default colors for all methods (custom_palette not provided but methods were renamed).")
            active_palette = None # Explicitly set to None

        if active_palette and not active_palette : # Check if it became an empty dict
             print("Warning: Active palette became empty after filtering. Seaborn will assign default colors.")
             active_palette = None

        fig_width = max(8, len(plot_order) * 1.2 if len(plot_order) > 1 else 8)
        plt.figure(figsize=(fig_width, 7))

        ax = sns.boxplot(x='display_method', y='aic_difference', hue='display_method',
                         data=plot_df, palette=active_palette, order=plot_order,
                         width=0.6, dodge=False, legend=False) 

        for i, method_name_in_plot in enumerate(plot_order):
            method_data = plot_df[plot_df['display_method'] == method_name_in_plot]['aic_difference']
            if not method_data.empty:
                mean_val = method_data.mean()
                ax.plot(i, mean_val, 'ro', markersize=8, zorder=5) 
                
                if show_mean_value:
                    if pd.notnull(mean_val):
                        ax.annotate(f'Mean: {mean_val:.2f}',
                                    xy=(i, mean_val),
                                    xytext=(8, -9), 
                                    textcoords='offset points',
                                    color='red',
                                    ha='left',
                                    va='center',
                                    fontsize=9)

        plt.title(f'AIC Difference relative to Model (P={true_p_order}, Q={true_q_order})\n{case_label}', fontsize=14)
        plt.xlabel('Method', fontsize=12) # X-axis title remains generic
        plt.ylabel(f'AIC(Ref. P={true_p_order},Q={true_q_order}) - AIC(Other Model)', fontsize=12)
        plt.axhline(0, color='grey', linestyle='--', linewidth=0.8)

        mean_dot_legend_entry = Line2D([0], [0], marker='o', color='w', 
                                       markerfacecolor='red', markersize=8,
                                       label='Mean')
        ax.legend(handles=[mean_dot_legend_entry], loc='best', fontsize=10)

        plt.grid(axis='y', linestyle=':', alpha=0.7)
        if len(plot_order) > 4:
            plt.xticks(rotation=45, ha='right')
        else:
            plt.xticks(rotation=0) # Tick labels will be "Stepwise", "Grid" etc. from plot_order
        plt.tight_layout()

        file_basename = os.path.basename(filepath)
        filename_prefix = os.path.splitext(file_basename)[0]
        safe_case_label = "".join(c if c.isalnum() or c in (' ', '_') else '_' for c in case_label).replace(' ', '_')
        filename = f"aic_differences_{filename_prefix}_{safe_case_label.lower()}_ref_p{true_p_order}q{true_q_order}.png"
        save_path = os.path.join(output_dir, filename)
        plt.savefig(save_path, dpi=300)
        plt.close()
        print(f"AIC difference boxplot saved successfully to: {save_path}")

    except Exception as e:
        print(f"An unexpected error occurred while processing {filepath} for AIC difference boxplots: {e}")
        print(traceback.format_exc())

# %% --- Example Usage ---
if __name__ == '__main__': # Ensures this runs only when script is executed directly
    # Define your file paths and corresponding specific titles
    files_and_titles = [
        {
            "path": r'C:\Users\Rafael\Desktop\msc-thesis\Results\Simulation\no_covariates\Model1\ingarch_no_covariates_results.csv',
            "title": "No Covariates M1"
        },
        {
            "path": r'C:\Users\Rafael\Desktop\msc-thesis\Results\Simulation\no_covariates\Model2\ingarch_no_covariates_results.csv',
            "title": "No Covariates M2"
        },
        {
            "path": r'C:\Users\Rafael\Desktop\msc-thesis\Results\Simulation\with_covariates\Combined\ingarch_with_covariates_results_combined.csv',
            "title": "With Covariates"
        }
    ]

    # Set reference P and Q orders to (2,6)
    REFERENCE_P_ORDER = 2
    REFERENCE_Q_ORDER = 6

    plot_output_directory = "ingarch_aic_difference_boxplots"

    # Optional: Define how method names from CSV should be displayed and their colors
    # Updated to "Stepwise" and "Grid"
    user_defined_method_names = {
        "stepwise": "Stepwise",  # Changed from "Stepwise Search"
        "grid_search": "Grid"    # Changed from "Grid Search"
    }

    # Updated palette keys to match the new display names
    user_defined_palette = {
        "Stepwise": "mediumseagreen", # Key changed
        "Grid": "#377E7F"           # Key changed (dark teal)
    }

    print("--- Generating AIC Difference Boxplots (Reference P=2, Q=6) from CSV data files ---")

    for item in files_and_titles:
        file_path = item["path"]
        current_case_label = item["title"] 

        print(f"\nProcessing file: {file_path} with custom title: '{current_case_label}'")

        # Before calling, check if the file exists to avoid errors for placeholder paths
        if not os.path.exists(file_path):
            print(f"File not found: {file_path}. Please update the path in 'files_and_titles'. Skipping this file.")
            continue
            
        create_aic_difference_boxplots(filepath=file_path,
                                         case_label=current_case_label,
                                         true_p_order=REFERENCE_P_ORDER, 
                                         true_q_order=REFERENCE_Q_ORDER, 
                                         output_dir=plot_output_directory,
                                         method_rename_map=user_defined_method_names,
                                         custom_palette=user_defined_palette,
                                         show_mean_value=False)

    print("\n--- AIC difference boxplot generation process complete. ---")
    print(f"Check the '{plot_output_directory}' directory for the plots if processing was successful.")