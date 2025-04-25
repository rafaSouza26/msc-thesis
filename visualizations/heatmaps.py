# %% Import Libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# %% Plotting Function (No Changes Here from last version)
def create_heatmap(data_matrix, row_sums, col_sums, total_for_percentage, calculated_total, title):
    """
    Creates a heatmap visualization with specified parameters.
    (Includes: y-axis inversion, vmax=100, top column totals, specific spacing, no marginal labels)

    Args:
        data_matrix (pd.DataFrame): The 8x8 matrix of counts (p x q).
        row_sums (pd.Series): Marginal totals for rows (q orders).
        col_sums (pd.Series): Marginal totals for columns (p orders).
        total_for_percentage (int): The total number of models (e.g., 1000) to use for percentage calculation.
        calculated_total (int): The actual sum of counts from the data matrix.
        title (str): The title for the heatmap plot.
    """
    # --- Figure Setup ---
    fig, ax = plt.subplots(figsize=(10.5, 9))

    # --- Percentage Calculation ---
    if total_for_percentage == 0:
        percentage_matrix = data_matrix * 0
        print(f"Warning: Total sum for percentage calculation is zero for '{title}'. Percentages cannot be calculated.")
    else:
        percentage_matrix = (data_matrix.astype(float) / total_for_percentage) * 100

    # --- Create Heatmap ---
    sns.heatmap(percentage_matrix,
                ax=ax,
                annot=True,
                fmt=".1f",
                cmap='viridis_r',
                linewidths=.5,
                linecolor='lightgray',
                cbar=True,
                cbar_kws={'label': 'Percentage (%)', 'shrink': 0.8},
                vmin=0,
                vmax=100)

    # --- Set labels, title, and ticks ---
    ax.set_xlabel('Order p', fontsize=12)
    ax.set_ylabel('Order q', fontsize=12)
    # Set the exact title passed to the function
    ax.set_title(title, fontsize=14, pad=50)

    ax.set_xticks(np.arange(data_matrix.shape[1]) + 0.5)
    ax.set_yticks(np.arange(data_matrix.shape[0]) + 0.5)
    ax.set_xticklabels(data_matrix.columns, rotation=0)
    ax.set_yticklabels(data_matrix.index, rotation=0)
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.tick_params(axis='both', which='minor', length=0)

    # --- Invert Y-axis ---
    ax.invert_yaxis()

    # --- Add Marginal Totals (Row - No Label) ---
    for i, total in enumerate(row_sums):
        ax.text(data_matrix.shape[1] + 0.1, i + 0.5, f'{total}', va='center', ha='left', fontsize=10)

    # --- Add Marginal Totals (Column - ON TOP - No Label) ---
    _, ymax_data = ax.get_ylim()
    offset_nums = 0.3

    for i, total in enumerate(col_sums):
        ax.text(i + 0.5, ymax_data + offset_nums, f'{total}',
                ha='center', va='bottom', fontsize=10,
                clip_on=False)

    # --- Add Overall Total (Top Right - Lowered Slightly) ---
    total_y_pos = ymax_data + offset_nums - 0.1
    ax.text(data_matrix.shape[1] + 0.1, total_y_pos, f'Total: {total_for_percentage}',
            ha='left', va='bottom', fontsize=10, fontweight='bold',
            clip_on=False)

    # --- Adjust Layout for Spacing ---
    plt.subplots_adjust(left=0.1, right=0.78, top=0.88, bottom=0.1)

    plt.show()


# %% Processing Function (TITLE UPDATED HERE)
def process_and_plot(filepath, case_label):
    """
    Loads data from a CSV, processes it for each method ('stepwise', 'grid_search'),
    and generates heatmaps according to user specifications.
    (Updated to use a static title)

    Args:
        filepath (str): The path to the CSV file.
        case_label (str): A label describing the dataset (e.g., "With Covariates", "No Covariates").
                         (This is no longer used in the title but kept for console messages).
    """
    # --- File Reading and Validation ---
    if not os.path.exists(filepath):
        print(f"Error: File not found at '{filepath}'. Please check the path.")
        return

    try:
        df = pd.read_csv(filepath)
        print(f"--- Processing file: {os.path.basename(filepath)} ({case_label}) ---")

        expected_cols = {'method': ['method', 'Method'],
                         'p': ['p', 'P'],
                         'q': ['q', 'Q'],
                         'count': ['count', 'Count', 'Freq', 'Frequency']}
        rename_map = {}
        found_cols = {}

        for standard_name, possible_names in expected_cols.items():
            found = False
            for possible_name in possible_names:
                if possible_name in df.columns:
                    if standard_name != possible_name:
                         rename_map[possible_name] = standard_name
                    found_cols[standard_name] = True
                    found = True
                    break
            if not found:
                 print(f"Error: Required column '{standard_name}' (or variations like {possible_names}) not found in {filepath}.")
                 return

        df.rename(columns=rename_map, inplace=True)
        print(f"Using columns: {df.columns.tolist()}")

        total_models_expected = 1000
        methods_in_file = df['method'].unique()
        print(f"Methods found: {methods_in_file}")

        # --- Process Each Method ---
        for method in methods_in_file:
            print(f"\nProcessing Method: {method}")
            method_df = df[df['method'] == method].copy()

            for col in ['p', 'q', 'count']:
                 method_df[col] = pd.to_numeric(method_df[col], errors='coerce')

            initial_rows = len(method_df)
            method_df.dropna(subset=['p', 'q', 'count'], inplace=True)
            if len(method_df) < initial_rows:
                print(f"Warning: Dropped {initial_rows - len(method_df)} rows with non-numeric p, q, or count values for method '{method}'.")

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

            # --- Ensure 8x8 Matrix (Orders 0-7) ---
            max_order = 7
            full_index = pd.Index(range(max_order + 1), name='q')
            full_columns = pd.Index(range(max_order + 1), name='p')
            freq_matrix = freq_matrix.reindex(index=full_index, columns=full_columns, fill_value=0).astype(int)

            # --- Calculate Totals ---
            row_totals = freq_matrix.sum(axis=1)
            col_totals = freq_matrix.sum(axis=0)
            calculated_total = freq_matrix.values.sum()

            print(f"Data summary for method '{method}':")
            # Keep console output specific to method and case for clarity
            print(f"  (File: {os.path.basename(filepath)} - Case: {case_label} - Method: {method})")
            print(f"  - Row totals (q counts): {row_totals.tolist()}")
            print(f"  - Column totals (p counts): {col_totals.tolist()}")
            print(f"  - Calculated total count: {calculated_total}")

            if calculated_total != total_models_expected:
                 print(f"  - Warning: Calculated total count ({calculated_total}) differs from the expected total ({total_models_expected}).")
                 print(f"    Percentages will be calculated based on {total_models_expected}, but marginal totals show actual counts.")

            # --- *** Set Static Plot Title *** ---
            plot_title = 'Order Selection Percentage' # Set the exact title here

            # --- Generate Heatmap ---
            create_heatmap(freq_matrix, row_totals, col_totals,
                           total_models_expected, calculated_total,
                           plot_title) # Pass the static title

    except FileNotFoundError:
        print(f"Error: File not found at {filepath}")
    except pd.errors.EmptyDataError:
        print(f"Error: File {filepath} is empty or not valid CSV.")
    except KeyError as e:
        print(f"Error: Column {e} not found in {filepath}. Please check CSV header.")
    except Exception as e:
        print(f"An unexpected error occurred while processing {filepath}: {e}")

# %% --- Example Usage ---
# IMPORTANT: Replace these placeholder paths with the actual paths to your files!

# --- Define your file paths here ---
path_to_file_with_covariates = r'C:\Users\Rafael\Desktop\msc-thesis\results\ingarch_with_covariates_order_freq_parallel.csv'
path_to_file_no_covariates   = r'C:\Users\Rafael\Desktop\msc-thesis\results\ingarch_no_covariates_order_freq.csv'

# --- Run the processing and plotting ---
process_and_plot(path_to_file_with_covariates, "With Covariates")
process_and_plot(path_to_file_no_covariates, "No Covariates")

print("\nProcessing complete.")