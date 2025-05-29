# %% Import Libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors # Import for colormap manipulation
import seaborn as sns
import os # Make sure os is imported

# %% Plotting Function (MODIFIED FOR SPECIAL ZERO COLORING)
def create_heatmap(data_matrix, row_sums, col_sums, total_for_percentage, calculated_total, title, save_path, zero_marker_color): # Added zero_marker_color
    """
    Creates and saves a heatmap visualization with specified parameters.
    Cells with an original count of 0 are colored with zero_marker_color.
    (Counts just above grid, large gap below title, saves plot to file)

    Args:
        data_matrix (pd.DataFrame): The 8x8 matrix of counts (p x q).
        row_sums (pd.Series): Marginal totals for rows (q orders).
        col_sums (pd.Series): Marginal totals for columns (p orders).
        total_for_percentage (int): The total number of models (e.g., 1000) to use for percentage calculation.
        calculated_total (int): The actual sum of counts from the data matrix.
        title (str): The title for the heatmap plot.
        save_path (str): The full path (directory + filename) where the plot should be saved.
        zero_marker_color (str): The color to use for cells with an original count of 0 (e.g., 'white', 'black').
    """
    # --- Figure Setup ---
    fig, ax = plt.subplots(figsize=(10, 10.5))

    # --- Identify original zero-count cells ---
    # These are cells that had 0 in the input data_matrix (raw counts)
    # or were filled with 0 during reindexing because they were missing.
    original_zero_mask = (data_matrix == 0)

    # --- Percentage Calculation (for annotations) ---
    # This matrix will be used for the text annotations in the cells.
    if total_for_percentage == 0:
        percentage_matrix_annotations = data_matrix * 0 # Should be all zeros
        print(f"Warning: Total sum for percentage calculation is zero for '{title}'. Percentages cannot be calculated meaningfully.")
    else:
        percentage_matrix_annotations = (data_matrix.astype(float) / total_for_percentage) * 100

    # --- Prepare data for heatmap coloring ---
    # Create a copy of the percentage matrix for coloring.
    # Values that were originally 0 count (and thus 0%) will be set to NaN.
    # These NaN values will then be colored by cmap.set_bad().
    percentage_matrix_for_coloring = percentage_matrix_annotations.copy()
    percentage_matrix_for_coloring[original_zero_mask] = np.nan

    # --- Create Custom Colormap ---
    # Get a copy of the 'viridis_r' colormap to modify
    custom_cmap = plt.cm.get_cmap('viridis_r').copy()
    # Set the color for NaN values (which correspond to our original zero-count cells)
    custom_cmap.set_bad(color=zero_marker_color)

    # --- Create Heatmap ---
    cbar_kws_params = {
        'orientation': 'horizontal',
        'location': 'bottom',
        'shrink': 0.7,
        'aspect': 30,
        'pad': 0.15
    }
    heatmap = sns.heatmap(percentage_matrix_for_coloring, # Data with NaNs for special coloring
                          ax=ax,
                          annot=percentage_matrix_annotations, # Annotate with original percentages (e.g., "0.0" for zeros)
                          fmt=".1f",
                          cmap=custom_cmap, # Use the modified colormap
                          linewidths=.5,
                          linecolor='lightgray',
                          cbar=True,
                          cbar_kws=cbar_kws_params,
                          vmin=0, # Percentages typically range from 0
                          vmax=100 if total_for_percentage > 0 else 1 # Ensure vmax is valid even if all are 0
                         )
    
    # If total_for_percentage is 0, all percentages are 0.
    # The colorbar might need adjustment if vmax was set to 1 and all values are 0 or NaN.
    # However, vmin=0, vmax=100 should handle 0% correctly by mapping it to the lowest color of viridis_r.
    # NaNs will get 'zero_marker_color'.

    # --- Adjust Colorbar Label Spacing ---
    cbar = heatmap.collections[0].colorbar
    cbar.set_label('Percentage (%)', labelpad=15)

    # --- Set labels, title, and ticks ---
    ax.set_xlabel('Order p', fontsize=12, labelpad=20)
    ax.set_ylabel('Order q', fontsize=12, labelpad=20)
    ax.set_title(title, fontsize=16, pad=60) # pad increased for title spacing

    ax.set_xticks(np.arange(data_matrix.shape[1]) + 0.5)
    ax.set_yticks(np.arange(data_matrix.shape[0]) + 0.5)
    ax.set_xticklabels(data_matrix.columns, rotation=0)
    ax.set_yticklabels(data_matrix.index, rotation=0)
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.tick_params(axis='x', pad=5)
    ax.tick_params(axis='y', pad=5)
    ax.tick_params(axis='both', which='minor', length=0) # Hide minor ticks

    # --- Invert Y-axis ---
    ax.invert_yaxis()

    # --- Add Marginal Totals (Row - Right Side) ---
    row_total_h_offset = 0.25
    for i, total in enumerate(row_sums):
        ax.text(data_matrix.shape[1] + row_total_h_offset, i + 0.5, f'{total}',
                va='center', ha='left', fontsize=10)

    # --- Add Marginal Totals (Column - Top Side - LOWERED POSITION) ---
    _, top_lim_axis = ax.get_ylim()
    col_total_v_offset = 0.35 
    col_total_y_pos = top_lim_axis + col_total_v_offset 

    for i, total in enumerate(col_sums):
        ax.text(i + 0.5, col_total_y_pos, f'{total}',
                ha='center', va='bottom', fontsize=10,
                clip_on=False) 

    # --- Add Overall Total (Top Right - Aligned Vertically with Column Totals) ---
    ax.text(data_matrix.shape[1] + row_total_h_offset, col_total_y_pos, f'Total(n): {total_for_percentage}',
            ha='left', va='bottom', fontsize=10, fontweight='bold',
            clip_on=False)

    # --- Adjust Layout for Spacing ---
    plt.subplots_adjust(left=0.15, right=0.85, top=0.80, bottom=0.18)

    # --- *** SAVE PLOT *** ---
    try:
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
        print(f"Plot saved successfully to: {save_path}")
    except Exception as e:
        print(f"Error saving plot to {save_path}: {e}")

    # --- Show Plot (Optional, can be commented out for batch processing) ---
    # plt.show() # Comment out if generating many plots to avoid display spam

    # --- *** CLOSE PLOT FIGURE *** ---
    plt.close(fig)


# %% Processing Function (MODIFIED TO GENERATE PLOT VARIANTS AND HANDLE SAVE PATH)
def process_and_plot(filepath, case_label, output_dir="ingarch_heatmaps"):
    """
    Loads data, processes each method, generates heatmaps with variations for zero-count cells,
    and saves them to a specified directory with descriptive names.

    Args:
        filepath (str): The path to the CSV file.
        case_label (str): Label for the dataset (e.g., "With Covariates"). Used in filename.
        output_dir (str): The directory where plots will be saved.
    """
    # --- Create Output Directory ---
    try:
        os.makedirs(output_dir, exist_ok=True) 
        print(f"Output directory set to: {output_dir}")
    except OSError as e:
        print(f"Error creating directory {output_dir}: {e}")
        return 

    # --- File Reading and Validation ---
    if not os.path.exists(filepath):
        print(f"Error: File not found at '{filepath}'. Please check the path.")
        return

    try:
        df = pd.read_csv(filepath)
        print(f"--- Processing file: {os.path.basename(filepath)} ({case_label}) ---")

        # --- Column Name Handling ---
        expected_cols = {'method': ['method', 'Method'],
                         'p': ['p', 'P', 'p_order'],
                         'q': ['q', 'Q', 'q_order'],
                         'count': ['count', 'Count', 'Freq', 'Frequency', 'freq']}
        rename_map = {}
        found_cols = {}

        for standard_name, possible_names in expected_cols.items():
            found = False
            for possible_name in possible_names:
                matching_cols = [col for col in df.columns if col.lower() == possible_name.lower()]
                if matching_cols:
                    actual_name = matching_cols[0]
                    if standard_name != actual_name and actual_name not in rename_map:
                         if actual_name not in rename_map.values(): 
                            rename_map[actual_name] = standard_name
                    found_cols[standard_name] = True
                    found = True
                    break
            if not found:
                print(f"Error: Required column '{standard_name}' (or variations like {possible_names}) not found in {filepath}.")
                print(f"   Available columns: {df.columns.tolist()}")
                return

        df.rename(columns=rename_map, inplace=True)
        if 'method' in df.columns:
            df['method'] = df['method'].astype(str).str.lower().str.strip()

        print(f"Using standardized columns: 'method', 'p', 'q', 'count'")

        total_models_expected = 1000 
        methods_in_file = df['method'].unique()
        print(f"Methods found (standardized): {methods_in_file}")

        # Define variations for zero-count cell coloring
        zero_marker_options = {
            "white_zeros": "white", # Filename suffix: actual color
            "black_zeros": "black"
        }

        for method in methods_in_file:
            print(f"\nProcessing Method: {method}")
            method_df = df[df['method'] == method].copy()

            # --- Data Cleaning ---
            for col in ['p', 'q', 'count']:
                method_df[col] = pd.to_numeric(method_df[col], errors='coerce')
            initial_rows = len(method_df)
            method_df.dropna(subset=['p', 'q', 'count'], inplace=True)
            if len(method_df) < initial_rows:
                print(f"Warning: Dropped {initial_rows - len(method_df)} rows with non-numeric/missing values for method '{method}'.")
            
            if method_df.empty: # Check if dataframe is empty after dropna
                print(f"Warning: No valid data remaining for method '{method}' after cleaning. Skipping plot.")
                continue
                
            method_df['p'] = method_df['p'].astype(int)
            method_df['q'] = method_df['q'].astype(int)
            method_df['count'] = method_df['count'].astype(int)

            # --- Pivot Data ---
            try:
                freq_matrix = method_df.pivot_table(index='q', columns='p', values='count',
                                                    fill_value=0, aggfunc='sum')
            except Exception as e:
                print(f"Error pivoting data for method '{method}': {e}")
                continue

            # --- Ensure 8x8 Matrix ---
            max_order = 7
            full_index = pd.Index(range(max_order + 1), name='q')
            full_columns = pd.Index(range(max_order + 1), name='p')
            freq_matrix = freq_matrix.reindex(index=full_index, columns=full_columns, fill_value=0).astype(int)

            # --- Calculate Totals ---
            row_totals = freq_matrix.sum(axis=1)
            col_totals = freq_matrix.sum(axis=0)
            calculated_total_from_data = freq_matrix.values.sum() # This is the sum of counts from the current method's data
            
            print(f"Data summary for method '{method}':")
            print(f"  (File: {os.path.basename(filepath)} - Case: {case_label} - Method: {method})")
            print(f"  - Row totals (q counts): {row_totals.tolist()}")
            print(f"  - Column totals (p counts): {col_totals.tolist()}")
            print(f"  - Calculated total count from data: {calculated_total_from_data}")

            # Use total_models_expected for percentage calculation, as in original script.
            # The warning for mismatch is good to keep.
            display_total_for_percentage = total_models_expected 
            if calculated_total_from_data != total_models_expected:
                print(f"  - Warning: Calculated total count from data ({calculated_total_from_data}) differs from the expected total for percentage calculation ({total_models_expected}). Percentages will be based on {total_models_expected}.")

            # --- Set Base Plot Title and Filename Method Part ---
            base_plot_title = ""
            filename_method_part = ""
            if method == 'stepwise':
                base_plot_title = 'Stepwise'
                filename_method_part = 'stepwise'
            elif method == 'grid_search' or method == 'gridsearch': # Accommodate 'gridsearch'
                base_plot_title = 'Grid Search'
                filename_method_part = 'grid_search'
            else:
                base_plot_title = method.capitalize()
                filename_method_part = method.replace(' ', '_').lower()
            
            case_suffix_for_filename = case_label.replace(' ', '_').lower()

            # Loop through zero marker options to create plot variants
            for zero_variant_name, actual_zero_color in zero_marker_options.items():
                # --- Construct Plot Title and Save Path for the Variant ---
                plot_title_variant = f"{base_plot_title}" # e.g., "Stepwise (white zeros)"
                filename = f"{filename_method_part}_{case_suffix_for_filename}_{zero_variant_name}.png"
                save_path = os.path.join(output_dir, filename)

                print(f"  Generating plot: {plot_title_variant} -> {filename}")
                # --- Generate and Save Heatmap for the Variant ---
                create_heatmap(freq_matrix, row_totals, col_totals,
                               display_total_for_percentage, # Total for percentage calculation
                               calculated_total_from_data,   # Actual sum from matrix (for reference, not directly used in % calc by default)
                               plot_title_variant,
                               save_path,
                               actual_zero_color) # Pass the specific color for zero cells

    # --- Error Handling ---
    except FileNotFoundError:
        print(f"Error: File not found at {filepath}")
    except pd.errors.EmptyDataError:
        print(f"Error: File {filepath} is empty or not valid CSV.")
    except KeyError as e:
        print(f"Error: Critical column {e} missing after attempting to standardize names in {filepath}. Cannot proceed.")
    except Exception as e:
        print(f"An unexpected error occurred while processing {filepath}: {e}")
        import traceback
        print(traceback.format_exc())


# %% --- Example Usage ---
# IMPORTANT: Replace these placeholder paths with the actual paths to your files!

# --- Define your file paths here ---
path_to_file_no_covariates = r'C:\Users\Rafael\Desktop\msc-thesis\Results\Simulation\no_covariates\ingarch_no_covariates_order_frequencies.csv'
path_to_file_with_covariates = r'C:\Users\Rafael\Desktop\msc-thesis\Results\Simulation\with_covariates\ingarch_with_covariates_order_frequencies.csv'

# --- Define Output Directory ---
plot_output_directory = "ingarch_order_heatmaps" # Define the folder name here (changed to distinguish)

# --- Run the processing and plotting ---
# Pass the output directory name to the function
if os.path.exists(path_to_file_with_covariates):
    process_and_plot(path_to_file_with_covariates, "With Covariates", output_dir=plot_output_directory)
else:
    print(f"File not found: {path_to_file_with_covariates}")

if os.path.exists(path_to_file_no_covariates):
    process_and_plot(path_to_file_no_covariates, "No Covariates", output_dir=plot_output_directory)
else:
    print(f"File not found: {path_to_file_no_covariates}")


print("\nProcessing complete.")