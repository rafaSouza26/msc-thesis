# %% Import Libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import traceback

# %% Plotting Function (AIC Difference Boxplots)
def create_aic_difference_boxplots(filepath,
                                     case_label,
                                     true_p_order,
                                     true_q_order,
                                     output_dir="ingarch_order_differences_boxplots",
                                     method_rename_map=None,
                                     custom_palette=None,
                                     show_mean_value=False): # New flag
    """
    Creates and saves box plots of the AIC difference between models with true orders (p,q)
    and other models, from a .csv file.
    The reference AIC for true orders is calculated as the mean AIC for models
    with (true_p_order, true_q_order) within each estimation method.
    Differences are AIC(true_order_model_mean) - AIC(other_model).

    Args:
        filepath (str): The path to the .csv file.
        case_label (str): Label for the dataset (e.g., "No Covariates"). Used in title/filename.
        true_p_order (int): The true order p, used to identify reference models.
        true_q_order (int): The true order q, used to identify reference models.
        output_dir (str): The directory where plots will be saved.
        method_rename_map (dict, optional): Maps original method names to display names.
        custom_palette (dict, optional): Maps display names (after renaming) to colors.
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

        for standard_name_key, standard_name_val in standard_names.items(): # Iterate through 'method', 'p', 'q', 'aic'
            found = False
            # 1. Check if the standard name itself (e.g., 'p') is already a column (case-insensitive check for its presence)
            if standard_name_val.lower() in current_columns_lower:
                original_case_col_name = current_columns_lower[standard_name_val.lower()]
                if original_case_col_name != standard_name_val: # If 'P' exists and standard is 'p', rename 'P' to 'p'
                    rename_map[original_case_col_name] = standard_name_val
                found = True
            # 2. If not found as is, check variants from possible_col_names_map
            if not found:
                for possible_name_variant_lower, mapped_standard_name in possible_col_names_map.items():
                    if mapped_standard_name == standard_name_val: # Ensure we are checking variants for the current standard_name_val
                        if possible_name_variant_lower in current_columns_lower:
                            original_case_col_name = current_columns_lower[possible_name_variant_lower]
                            # Only add to rename_map if it's not already the standard name
                            if original_case_col_name != standard_name_val:
                                rename_map[original_case_col_name] = standard_name_val
                            found = True
                            break # Found a variant for this standard_name_val
            
            if not found: # After checking standard name and its variants
                # Construct a list of possible names for the error message for this specific standard_name_val
                possible_options_for_error = [k for k,v in possible_col_names_map.items() if v == standard_name_val]
                # Add the standard name itself as a primary option
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
            print(f"Warning: No data available for plotting AIC differences for file '{filepath}' (perhaps only true models existed or no other models had valid AICs).")
            return

        plot_df = pd.DataFrame(plot_data_rows)
        plot_order = sorted(plot_df['display_method'].unique()) 

        if not plot_order: 
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
            if active_palette and len(active_palette) < len(plot_order) and plot_order: 
                print("  Not all methods in the plot have a color in default_internal_palette. Seaborn will use default colors for missing ones.")
        else: 
            print("Seaborn will assign default colors for all methods (custom_palette not provided or methods were renamed without one).")
            active_palette = None 

        if active_palette and not active_palette : 
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
                ax.plot(i, mean_val, 'ro', markersize=8, zorder=5) # Red dot for mean
                
                if show_mean_value: # Display mean value if flag is True
                    if pd.notnull(mean_val):
                        ax.annotate(f'Mean: {mean_val:.2f}', # Text format "Mean: value"
                                    xy=(i, mean_val),               # Point to annotate (center of the red dot)
                                    xytext=(8, -7),                  # Offset text 8 points to the right, 0 points vertical
                                    textcoords='offset points',     # Specify offset in points
                                    color='red',                    # Text color red
                                    ha='left',                      # Horizontal alignment of text
                                    va='center',                    # Vertical alignment of text
                                    fontsize=9)                     # Font size (adjust as needed)


        plt.title(f'AIC Difference relative to Model (P={true_p_order}, Q={true_q_order})\n{case_label}', fontsize=14)
        plt.xlabel('Method', fontsize=12)
        plt.ylabel(f'AIC(P={true_p_order},Q={true_q_order}) - AIC(Other Model)', fontsize=12)
        plt.axhline(0, color='grey', linestyle='--', linewidth=0.8) 

        if active_palette :
            handles = [plt.Rectangle((0,0),1,1, color=active_palette.get(label, 'gray')) for label in plot_order if label in active_palette]
            labels = [label for label in plot_order if label in active_palette]
            if handles: 
                ax.legend(handles, labels, title='Method', loc='best')

        plt.grid(axis='y', linestyle=':', alpha=0.7)
        if len(plot_order) > 4: 
            plt.xticks(rotation=45, ha='right')
        else:
            plt.xticks(rotation=0)
        plt.tight_layout() 

        file_basename = os.path.basename(filepath)
        filename_prefix = os.path.splitext(file_basename)[0]
        safe_case_label = "".join(c if c.isalnum() or c in (' ', '_') else '_' for c in case_label).replace(' ', '_')
        filename = f"aic_differences_{filename_prefix}_{safe_case_label.lower()}_p{true_p_order}q{true_q_order}.png"
        save_path = os.path.join(output_dir, filename)
        plt.savefig(save_path, dpi=300)
        plt.close() 
        print(f"AIC difference boxplot saved successfully to: {save_path}")

    except Exception as e:
        print(f"An unexpected error occurred while processing {filepath} for AIC difference boxplots: {e}")
        print(traceback.format_exc())

# %% --- Example Usage ---
# IMPORTANT: Pandas, Matplotlib, Seaborn are required.
file_paths_to_process = [
    r'C:\Users\Rafael\Desktop\msc-thesis\SimulationResults\ingarch_no_covariates_results.csv',
    r'C:\Users\Rafael\Desktop\msc-thesis\SimulationResults\ingarch_with_covariates_results.csv'
]

case_labels_for_files = [
    "No Covariates Data",
    "With Covariates Data"
]

if len(file_paths_to_process) != len(case_labels_for_files):
    print("Error: The number of file paths and case labels must match.")
    print(f"  Number of file paths: {len(file_paths_to_process)}")
    print(f"  Number of case labels: {len(case_labels_for_files)}")
else:
    plot_output_directory = "ingarch_aic_difference_boxplots"

    TRUE_P_FOR_REF_AIC = 2
    TRUE_Q_FOR_REF_AIC = 6

    user_defined_method_names = {
        "stepwise": "Stepwise",
        "grid_search": "Grid Search"
    }
    user_defined_palette = {
        "Stepwise": "mediumseagreen", # Example color
        "Grid Search": "#377E7F"   # Example color (dark teal)
    }

    print("--- Generating AIC Difference Boxplots from CSV data files ---")

    for i, file_path in enumerate(file_paths_to_process):
        current_case_label = case_labels_for_files[i]
        print(f"\nProcessing file: {file_path} with case label: '{current_case_label}'")
        
        if os.path.exists(file_path):
            create_aic_difference_boxplots(filepath=file_path,
                                           case_label=current_case_label,
                                           true_p_order=TRUE_P_FOR_REF_AIC,
                                           true_q_order=TRUE_Q_FOR_REF_AIC,
                                           output_dir=plot_output_directory,
                                           method_rename_map=user_defined_method_names,
                                           custom_palette=user_defined_palette,
                                           show_mean_value=True) # <<< Set to True to show styled mean values
        else:
            print(f"File not found: {file_path}")
            print(f"  Please update 'file_paths_to_process' with the correct CSV file path(s).")
    
    # Example of calling without showing mean value (default behavior)
    if len(file_paths_to_process) > 0 and os.path.exists(file_paths_to_process[0]):
        print(f"\nProcessing file (again for demo, mean value NOT shown): {file_paths_to_process[0]} with case label: '{case_labels_for_files[0]} (No Mean Value)'")
        create_aic_difference_boxplots(filepath=file_paths_to_process[0],
                                       case_label=f"{case_labels_for_files[0]} (No Mean Value)", # Modified label for demo
                                       true_p_order=TRUE_P_FOR_REF_AIC,
                                       true_q_order=TRUE_Q_FOR_REF_AIC,
                                       output_dir=plot_output_directory,
                                       method_rename_map=user_defined_method_names,
                                       custom_palette=user_defined_palette,
                                       show_mean_value=False) # <<< Explicitly False or omit

    print("\n--- AIC difference boxplot generation process complete. ---")
    print(f"Check the '{plot_output_directory}' directory for the plots if processing was successful.")