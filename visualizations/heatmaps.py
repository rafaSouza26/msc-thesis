# %% Import Libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os # Make sure os is imported

# %% Plotting Function (ADDED SAVE FUNCTIONALITY)
def create_heatmap(data_matrix, row_sums, col_sums, total_for_percentage, calculated_total, title, save_path): # Added save_path argument
    """
    Creates and saves a heatmap visualization with specified parameters.
    (Counts just above grid, large gap below title, saves plot to file)

    Args:
        data_matrix (pd.DataFrame): The 8x8 matrix of counts (p x q).
        row_sums (pd.Series): Marginal totals for rows (q orders).
        col_sums (pd.Series): Marginal totals for columns (p orders).
        total_for_percentage (int): The total number of models (e.g., 1000) to use for percentage calculation.
        calculated_total (int): The actual sum of counts from the data matrix.
        title (str): The title for the heatmap plot (e.g., "Stepwise", "Grid Search").
        save_path (str): The full path (directory + filename) where the plot should be saved.
    """
    # --- Figure Setup ---
    fig, ax = plt.subplots(figsize=(10, 10.5))

    # --- Percentage Calculation ---
    # (No changes here)
    if total_for_percentage == 0:
        percentage_matrix = data_matrix * 0
        print(f"Warning: Total sum for percentage calculation is zero for '{title}'. Percentages cannot be calculated.")
    else:
        percentage_matrix = (data_matrix.astype(float) / total_for_percentage) * 100

    # --- Create Heatmap ---
    # (No changes here)
    cbar_kws_params = {
        'orientation': 'horizontal',
        'location': 'bottom',
        'shrink': 0.7,
        'aspect': 30,
        'pad': 0.15
    }
    heatmap = sns.heatmap(percentage_matrix,
                          ax=ax,
                          annot=True,
                          fmt=".1f",
                          cmap='viridis_r',
                          linewidths=.5,
                          linecolor='lightgray',
                          cbar=True,
                          cbar_kws=cbar_kws_params,
                          vmin=0,
                          vmax=100)

    # --- Adjust Colorbar Label Spacing ---
    # (No changes here)
    cbar = heatmap.collections[0].colorbar
    cbar.set_label('Percentage (%)', labelpad=15)

    # --- Set labels, title, and ticks ---
    # (No changes here)
    ax.set_xlabel('Order p', fontsize=12, labelpad=20)
    ax.set_ylabel('Order q', fontsize=12, labelpad=20)
    ax.set_title(title, fontsize=16, pad=60)

    ax.set_xticks(np.arange(data_matrix.shape[1]) + 0.5)
    ax.set_yticks(np.arange(data_matrix.shape[0]) + 0.5)
    ax.set_xticklabels(data_matrix.columns, rotation=0)
    ax.set_yticklabels(data_matrix.index, rotation=0)
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.tick_params(axis='x', pad=5)
    ax.tick_params(axis='y', pad=5)
    ax.tick_params(axis='both', which='minor', length=0)

    # --- Invert Y-axis ---
    # (No changes here)
    ax.invert_yaxis()

    # --- Add Marginal Totals (Row - Right Side) ---
    # (No changes here)
    row_total_h_offset = 0.25
    for i, total in enumerate(row_sums):
        ax.text(data_matrix.shape[1] + row_total_h_offset, i + 0.5, f'{total}',
                va='center', ha='left', fontsize=10)

    # --- Add Marginal Totals (Column - Top Side - LOWERED POSITION) ---
    # (No changes here)
    _, top_lim_axis = ax.get_ylim()
    col_total_v_offset = 0.35
    col_total_y_pos = top_lim_axis + col_total_v_offset

    for i, total in enumerate(col_sums):
         ax.text(i + 0.5, col_total_y_pos, f'{total}',
                 ha='center', va='bottom', fontsize=10,
                 clip_on=False)

    # --- Add Overall Total (Top Right - Aligned Vertically with Column Totals) ---
    # (No changes here)
    ax.text(data_matrix.shape[1] + row_total_h_offset, col_total_y_pos, f'Total(n): {total_for_percentage}',
            ha='left', va='bottom', fontsize=10, fontweight='bold',
            clip_on=False)

    # --- Adjust Layout for Spacing ---
    # (No changes here)
    plt.subplots_adjust(left=0.15, right=0.85, top=0.80, bottom=0.18)

    # --- *** SAVE PLOT *** ---
    try:
        # Use bbox_inches='tight' to minimize whitespace and prevent cutoff
        # Use dpi=300 for higher resolution suitable for documents
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
        print(f"Plot saved successfully to: {save_path}")
    except Exception as e:
        print(f"Error saving plot to {save_path}: {e}")

    # --- Show Plot (Optional) ---
    plt.show()

    # --- *** CLOSE PLOT FIGURE *** ---
    # Close the figure to free up memory, important in loops
    plt.close(fig)


# %% Processing Function (MODIFIED TO HANDLE SAVE PATH AND DIRECTORY)
def process_and_plot(filepath, case_label, output_dir="ingarch_plots"): # Added output_dir argument
    """
    Loads data, processes each method, generates heatmaps, and saves them
    to a specified directory with descriptive names.

    Args:
        filepath (str): The path to the CSV file.
        case_label (str): Label for the dataset (e.g., "With Covariates"). Used in filename.
        output_dir (str): The directory where plots will be saved.
    """
    # --- Create Output Directory ---
    try:
        os.makedirs(output_dir, exist_ok=True) # exist_ok=True prevents error if dir exists
        print(f"Output directory set to: {output_dir}")
    except OSError as e:
        print(f"Error creating directory {output_dir}: {e}")
        return # Stop if we can't create the directory

    # --- File Reading and Validation ---
    # (No changes here)
    if not os.path.exists(filepath):
        print(f"Error: File not found at '{filepath}'. Please check the path.")
        return

    try:
        df = pd.read_csv(filepath)
        print(f"--- Processing file: {os.path.basename(filepath)} ({case_label}) ---")

        # --- Column Name Handling ---
        # (No changes here)
        expected_cols = {'method': ['method', 'Method'],
                         'p': ['p', 'P'],
                         'q': ['q', 'Q'],
                         'count': ['count', 'Count', 'Freq', 'Frequency']}
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


        # --- Process Each Method ---
        for method in methods_in_file:
            print(f"\nProcessing Method: {method}")
            method_df = df[df['method'] == method].copy()

            # --- Data Cleaning ---
            # (No changes here)
            for col in ['p', 'q', 'count']:
                method_df[col] = pd.to_numeric(method_df[col], errors='coerce')
            initial_rows = len(method_df)
            method_df.dropna(subset=['p', 'q', 'count'], inplace=True)
            if len(method_df) < initial_rows:
                print(f"Warning: Dropped {initial_rows - len(method_df)} rows with non-numeric/missing values for method '{method}'.")
            if not method_df.empty:
                 method_df['p'] = method_df['p'].astype(int)
                 method_df['q'] = method_df['q'].astype(int)
                 method_df['count'] = method_df['count'].astype(int)
            else:
                 print(f"Warning: No valid data remaining for method '{method}'. Skipping plot.")
                 continue

            # --- Pivot Data ---
            # (No changes here)
            try:
                freq_matrix = method_df.pivot_table(index='q', columns='p', values='count',
                                                    fill_value=0, aggfunc='sum')
            except Exception as e:
                print(f"Error pivoting data for method '{method}': {e}")
                continue

            # --- Ensure 8x8 Matrix ---
            # (No changes here)
            max_order = 7
            full_index = pd.Index(range(max_order + 1), name='q')
            full_columns = pd.Index(range(max_order + 1), name='p')
            freq_matrix = freq_matrix.reindex(index=full_index, columns=full_columns, fill_value=0).astype(int)

            # --- Calculate Totals ---
            # (No changes here)
            row_totals = freq_matrix.sum(axis=1)
            col_totals = freq_matrix.sum(axis=0)
            calculated_total = freq_matrix.values.sum()
            print(f"Data summary for method '{method}':")
            print(f"  (File: {os.path.basename(filepath)} - Case: {case_label} - Method: {method})")
            print(f"  - Row totals (q counts): {row_totals.tolist()}")
            print(f"  - Column totals (p counts): {col_totals.tolist()}")
            print(f"  - Calculated total count from data: {calculated_total}")
            if calculated_total != total_models_expected:
                print(f"  - Warning: Calculated total count ({calculated_total}) differs from the expected total ({total_models_expected}).")
            display_total = total_models_expected

            # --- Set Plot Title ---
            # (No changes here)
            if method == 'stepwise':
                plot_title = 'Stepwise'
                filename_method = 'stepwise' # Use consistent lowercase for filename part
            elif method == 'grid_search' or method == 'gridsearch':
                plot_title = 'Grid Search'
                filename_method = 'grid_search' # Use consistent lowercase for filename part
            else:
                plot_title = method.capitalize()
                filename_method = method.replace(' ', '_').lower() # Basic fallback for other methods

            # --- *** CONSTRUCT FILENAME AND SAVEPATH *** ---
            # Create a clean version of case_label for the filename
            case_suffix = case_label.replace(' ', '_').lower() # e.g., "with_covariates"
            filename = f"{filename_method}_{case_suffix}.png"
            save_path = os.path.join(output_dir, filename)

            # --- Generate and Save Heatmap ---
            # Pass the save_path to the plotting function
            create_heatmap(freq_matrix, row_totals, col_totals,
                           display_total,
                           calculated_total,
                           plot_title,
                           save_path) # Pass the full save path

    # --- Error Handling ---
    # (No changes here)
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
path_to_file_with_covariates = r'C:\Users\Rafael\Desktop\msc-thesis\results\ingarch_with_covariates_order_freq.csv'
path_to_file_no_covariates   = r'C:\Users\Rafael\Desktop\msc-thesis\results\ingarch_no_covariates_order_freq.csv'

# --- Define Output Directory ---
plot_output_directory = "ingarch_plots" # Define the folder name here

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