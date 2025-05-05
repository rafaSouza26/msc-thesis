import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

# --- Configuration ---
# Adjust these file paths if your CSV files are located elsewhere
GRID_SEARCH_FILE = 'ingarch_no_covariates_results_gs.csv'
STEPWISE_SEARCH_FILE = 'ingarch_no_covariates_results_sw.csv'
OUTPUT_DIR = 'ingarch_plots' # Directory to save plots

# Known true orders (since they are not in the files)
KNOWN_P_TRUE = 2
KNOWN_Q_TRUE = 6

# --- Column Names ---
# !!! IMPORTANT: Verify these column names match your CSV files !!!
# If your column names are different, update them here.
# Updated based on user feedback: estimated p/q columns are named 'p' and 'q'
P_ESTIMATED_COL = 'p'                 # Column for the estimated 'p' order
Q_ESTIMATED_COL = 'q'                 # Column for the estimated 'q' order
AIC_COL = 'aic'                 # Column for the AIC value
BIC_COL = 'bic'                 # Column for the BIC value
TIME_COL = 'time'               # Column for the execution time

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- Load Data ---
try:
    df_gs = pd.read_csv(GRID_SEARCH_FILE)
    df_sw = pd.read_csv(STEPWISE_SEARCH_FILE)
    print("--- Grid Search Data Info ---")
    df_gs.info()
    print("\n--- Stepwise Search Data Info ---")
    df_sw.info()
    print(f"\nSuccessfully loaded data.")
    print(f"Grid Search columns: {df_gs.columns.tolist()}")
    print(f"Stepwise Search columns: {df_sw.columns.tolist()}")
    print(f"\nUsing Known True Orders: p={KNOWN_P_TRUE}, q={KNOWN_Q_TRUE}")
    print("\n!!! Please verify the column names above match the configuration !!!\n")

except FileNotFoundError as e:
    print(f"Error loading file: {e}")
    print("Please ensure the CSV files are in the same directory as the script or provide the correct path.")
    exit()
except Exception as e:
    print(f"An error occurred loading data: {e}")
    exit()

# --- Data Validation ---
# Include BIC_COL in required columns
required_cols = [P_ESTIMATED_COL, Q_ESTIMATED_COL, AIC_COL, BIC_COL, TIME_COL]
missing_gs_cols = [col for col in required_cols if col not in df_gs.columns]
missing_sw_cols = [col for col in required_cols if col not in df_sw.columns]

# Check for missing columns in Grid Search file
if missing_gs_cols:
    print(f"Error: Missing required columns in Grid Search file ({GRID_SEARCH_FILE}): {missing_gs_cols}")
    print(f"Script expected columns based on configuration: {required_cols}")
    print(f"Columns found in file: {df_gs.columns.tolist()}")
    print("Please check the '--- Column Names ---' section in the script and your CSV file.")
    exit()
# Check for missing columns in Stepwise Search file
if missing_sw_cols:
    print(f"Error: Missing required columns in Stepwise Search file ({STEPWISE_SEARCH_FILE}): {missing_sw_cols}")
    print(f"Script expected columns based on configuration: {required_cols}")
    print(f"Columns found in file: {df_sw.columns.tolist()}")
    print("Please check the '--- Column Names ---' section in the script and your CSV file.")
    exit()

# --- Visualization 1: Percentage Hitting True Orders ---
def plot_order_hit_percentage(df_gs, df_sw, known_p_true, known_q_true, p_est, q_est, output_dir):
    """Calculates and plots the percentage of times p, q, and both are correctly estimated using known true orders."""
    print("Generating Plot 1: Order Hit Percentage...")
    results = {}
    for name, df in [('Grid Search', df_gs), ('Stepwise Search', df_sw)]:
        # Ensure p/q columns are numeric before comparison
        df[p_est] = pd.to_numeric(df[p_est], errors='coerce')
        df[q_est] = pd.to_numeric(df[q_est], errors='coerce')
        df_valid = df.dropna(subset=[p_est, q_est]) # Use only rows with valid p, q
        n_total = len(df_valid)

        if n_total == 0:
            print(f"Warning: No valid (p,q) pairs found for {name} in plot_order_hit_percentage.")
            results[name] = {'p_correct': 0, 'q_correct': 0, 'both_correct': 0}
            continue

        # Compare estimated p/q (now columns 'p', 'q') with known true p/q
        p_correct_count = (df_valid[p_est] == known_p_true).sum()
        q_correct_count = (df_valid[q_est] == known_q_true).sum()
        both_correct_count = ((df_valid[p_est] == known_p_true) & (df_valid[q_est] == known_q_true)).sum()

        results[name] = {
            'p_correct': (p_correct_count / n_total) * 100,
            'q_correct': (q_correct_count / n_total) * 100,
            'both_correct': (both_correct_count / n_total) * 100
        }

    df_plot = pd.DataFrame(results).T.reset_index()
    df_plot = df_plot.rename(columns={'index': 'Method',
                                      'p_correct': f'% Correct p (True={known_p_true})',
                                      'q_correct': f'% Correct q (True={known_q_true})',
                                      'both_correct': f'% Correct (p,q) (True={known_p_true},{known_q_true})'})

    df_melt = df_plot.melt(id_vars='Method', var_name='Metric', value_name='Percentage')

    plt.figure(figsize=(12, 6)) # Increased width for longer labels
    sns.barplot(data=df_melt, x='Metric', y='Percentage', hue='Method', palette='viridis')
    plt.title('Percentage of Correctly Estimated Orders (p, q)')
    plt.ylabel('Percentage (%)')
    plt.xlabel('Correct Order Metric')
    plt.ylim(0, 100)
    plt.xticks(rotation=10, ha='right') # Rotate labels slightly if needed
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '1_order_hit_percentage.png'))
    plt.close()
    print("Plot 1 saved.")

# --- Visualization 2: Heatmap of Selected Orders --- ## MODIFIED TO SEPARATE PLOTS ##
def plot_order_selection_heatmap(df_gs, df_sw, p_est_col, q_est_col, output_dir):
    """Generates separate heatmaps showing the frequency of selected (p, q) orders for each method."""
    print("Generating Plot 2: Order Selection Heatmaps...")

    plot_suffixes = {'Grid Search': 'gs', 'Stepwise Search': 'sw'}
    plot_letters = {'Grid Search': 'a', 'Stepwise Search': 'b'}

    for name, df in [('Grid Search', df_gs), ('Stepwise Search', df_sw)]:
        plot_suffix = plot_suffixes[name]
        plot_letter = plot_letters[name]
        print(f"  Generating heatmap for {name}...")

        # Ensure p/q columns are numeric and drop NAs for grouping
        df[p_est_col] = pd.to_numeric(df[p_est_col], errors='coerce')
        df[q_est_col] = pd.to_numeric(df[q_est_col], errors='coerce')
        df_valid = df.dropna(subset=[p_est_col, q_est_col])

        if df_valid.empty:
            print(f"  Warning: No valid (p,q) pairs found for {name}. Skipping heatmap.")
            continue

        # Convert p/q to integer for grouping if they are not already
        df_valid = df_valid.astype({p_est_col: int, q_est_col: int})

        # Calculate frequency of each (p, q) pair
        order_counts = df_valid.groupby([p_est_col, q_est_col]).size().unstack(fill_value=0)

        # Determine max p and q for consistent axes, ensure they are at least 0
        max_p = max(order_counts.columns.max(), 0)
        max_q = max(order_counts.index.max(), 0)

        # Reindex to ensure all p/q values up to the max are present, fill missing with 0
        p_range = range(int(order_counts.columns.min()), int(max_p) + 1)
        q_range = range(int(order_counts.index.min()), int(max_q) + 1)
        order_counts = order_counts.reindex(index=q_range, columns=p_range, fill_value=0)

        # Sort index (q) descending so higher q is at the top
        order_counts = order_counts.sort_index(ascending=False)

        # Create a new figure for this heatmap
        plt.figure(figsize=(7, 6))
        ax = plt.gca() # Get current axes

        # Generate heatmap
        sns.heatmap(order_counts, annot=True, fmt="d", cmap="viridis", ax=ax, cbar=True) # Show color bar
        ax.set_title(f'{name}: (p, q) Order Selection Frequency')
        ax.set_xlabel(f'Order p ({p_est_col})')
        ax.set_ylabel(f'Order q ({q_est_col})')

        plt.tight_layout()
        filename = f'2{plot_letter}_heatmap_{plot_suffix}.png'
        plt.savefig(os.path.join(output_dir, filename))
        plt.close()
        print(f"  Plot 2{plot_letter} ({filename}) saved.")


# --- Visualization 3: Difference between Original and Estimated Orders --- ## MODIFIED TO SEPARATE PLOTS ##
def plot_order_differences(df_gs, df_sw, known_p_true, known_q_true, p_est, q_est, output_dir):
    """Generates separate scatter plots of the differences (p_est - known_p_true, q_est - known_q_true) for each method."""
    print("Generating Plot 3: Order Differences Scatter Plots...")

    plot_suffixes = {'Grid Search': 'gs', 'Stepwise Search': 'sw'}
    plot_letters = {'Grid Search': 'a', 'Stepwise Search': 'b'}
    plot_colors = {'Grid Search': 'blue', 'Stepwise Search': 'orange'}

    for name, df in [('Grid Search', df_gs), ('Stepwise Search', df_sw)]:
        plot_suffix = plot_suffixes[name]
        plot_letter = plot_letters[name]
        plot_color = plot_colors[name]
        print(f"  Generating scatter plot for {name}...")

        # Calculate difference using known true values and estimated columns ('p', 'q')
        df['p_diff'] = pd.to_numeric(df[p_est], errors='coerce') - known_p_true
        df['q_diff'] = pd.to_numeric(df[q_est], errors='coerce') - known_q_true

        # Drop rows where differences couldn't be calculated
        df_diff_valid = df.dropna(subset=['p_diff', 'q_diff'])

        if df_diff_valid.empty:
            print(f"  Warning: No valid difference data found for {name}. Skipping scatter plot.")
            continue

        # Create a new figure for this scatter plot
        plt.figure(figsize=(7, 6))
        ax = plt.gca() # Get current axes

        # Generate scatter plot
        ax.scatter(df_diff_valid['p_diff'], df_diff_valid['q_diff'], alpha=0.6, s=50, c=plot_color)
        ax.set_title(f'{name}: Order Difference (Est - True[{known_p_true},{known_q_true}])')
        ax.set_xlabel(f'p_estimated - {known_p_true}')
        ax.set_ylabel(f'q_estimated - {known_q_true}')
        ax.axhline(0, color='grey', linestyle='--', lw=1)
        ax.axvline(0, color='grey', linestyle='--', lw=1)
        ax.grid(True, linestyle='--', alpha=0.5)

        # Ensure integer ticks if differences are integers
        max_abs_diff_p = df_diff_valid['p_diff'].abs().max()
        max_abs_diff_q = df_diff_valid['q_diff'].abs().max()
        min_diff_p = df_diff_valid['p_diff'].min()
        min_diff_q = df_diff_valid['q_diff'].min()

        if not pd.isna(max_abs_diff_p) and not pd.isna(min_diff_p):
            ax.set_xticks(np.arange(int(min_diff_p), int(max_abs_diff_p) + 2)) # Extend range slightly
        if not pd.isna(max_abs_diff_q) and not pd.isna(min_diff_q):
            ax.set_yticks(np.arange(int(min_diff_q), int(max_abs_diff_q) + 2)) # Extend range slightly

        plt.tight_layout()
        filename = f'3{plot_letter}_scatter_diff_{plot_suffix}.png'
        plt.savefig(os.path.join(output_dir, filename))
        plt.close()
        print(f"  Plot 3{plot_letter} ({filename}) saved.")


# --- Visualization 4: Compare AIC (Box plot of differences) --- ## RENUMBERED from 3 to 4 ##
def plot_aic_difference(df_gs, df_sw, aic_col, output_dir):
    """Merges data assuming rows correspond and plots the difference in AIC."""
    print("Generating Plot 4: AIC Difference...") # Renumbered
    # Assuming rows correspond to the same simulation runs.
    # If there's an ID column, merge on that instead.
    if len(df_gs) != len(df_sw):
        print("Warning: Grid Search and Stepwise data have different numbers of rows.")
        print("Cannot directly compare AIC row-by-row. Skipping AIC difference plot.")
        return

    # Convert AIC columns to numeric, coercing errors
    aic_gs_numeric = pd.to_numeric(df_gs[aic_col], errors='coerce')
    aic_sw_numeric = pd.to_numeric(df_sw[aic_col], errors='coerce')

    df_merged = pd.DataFrame({
        'AIC_gs': aic_gs_numeric,
        'AIC_sw': aic_sw_numeric
    })
    df_merged['AIC_diff (GS - SW)'] = df_merged['AIC_gs'] - df_merged['AIC_sw']

    # Drop rows where difference could not be calculated (due to NaNs in original AICs)
    df_merged.dropna(subset=['AIC_diff (GS - SW)'], inplace=True)

    if df_merged.empty:
        print(f"Warning: No valid comparable AIC data found in column '{aic_col}'. Skipping AIC difference plot.")
        return

    plt.figure(figsize=(8, 6))
    sns.boxplot(y=df_merged['AIC_diff (GS - SW)'], palette='viridis')
    plt.title('Difference in AIC (Grid Search - Stepwise Search)')
    plt.ylabel('AIC_gs - AIC_sw')
    plt.axhline(0, color='red', linestyle='--', lw=1, label='AIC_gs = AIC_sw')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '4_aic_difference_boxplot.png')) # Renamed file
    plt.close()
    print("Plot 4 saved.")

# --- Visualization 5: Compare BIC (Box plot of differences) --- ## RENUMBERED from 4 to 5 ##
def plot_bic_difference(df_gs, df_sw, bic_col, output_dir):
    """Merges data assuming rows correspond and plots the difference in BIC."""
    print("Generating Plot 5: BIC Difference...") # Renumbered
    # Assuming rows correspond to the same simulation runs.
    if len(df_gs) != len(df_sw):
        print("Warning: Grid Search and Stepwise data have different numbers of rows.")
        print("Cannot directly compare BIC row-by-row. Skipping BIC difference plot.")
        return

    # Convert BIC columns to numeric, coercing errors
    bic_gs_numeric = pd.to_numeric(df_gs[bic_col], errors='coerce')
    bic_sw_numeric = pd.to_numeric(df_sw[bic_col], errors='coerce')

    df_merged = pd.DataFrame({
        'BIC_gs': bic_gs_numeric,
        'BIC_sw': bic_sw_numeric
    })
    df_merged['BIC_diff (GS - SW)'] = df_merged['BIC_gs'] - df_merged['BIC_sw']

    # Drop rows where difference could not be calculated (due to NaNs in original BICs)
    df_merged.dropna(subset=['BIC_diff (GS - SW)'], inplace=True)

    if df_merged.empty:
        print(f"Warning: No valid comparable BIC data found in column '{bic_col}'. Skipping BIC difference plot.")
        return

    plt.figure(figsize=(8, 6))
    sns.boxplot(y=df_merged['BIC_diff (GS - SW)'], palette='plasma') # Different palette
    plt.title('Difference in BIC (Grid Search - Stepwise Search)')
    plt.ylabel('BIC_gs - BIC_sw')
    plt.axhline(0, color='red', linestyle='--', lw=1, label='BIC_gs = BIC_sw')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '5_bic_difference_boxplot.png')) # Renamed file
    plt.close()
    print("Plot 5 saved.")


# --- Visualization 6: Box plot with Times --- ## RENUMBERED from 5 to 6, UPDATED LABEL ##
def plot_times(df_gs, df_sw, time_col, output_dir):
    """Creates a box plot comparing execution times (in milliseconds)."""
    print("Generating Plot 6: Execution Times...") # Renumbered
    df_times = pd.DataFrame({
        'Time': pd.concat([df_gs[time_col], df_sw[time_col]], ignore_index=True),
        'Method': ['Grid Search'] * len(df_gs) + ['Stepwise Search'] * len(df_sw)
    })

    # Check for non-numeric time data
    if not pd.api.types.is_numeric_dtype(df_times['Time']):
         print(f"Warning: Time column ('{time_col}') is not numeric. Attempting conversion.")
         df_times['Time'] = pd.to_numeric(df_times['Time'], errors='coerce')
         # Drop rows where time could not be converted
         original_len = len(df_times)
         df_times = df_times.dropna(subset=['Time'])
         if len(df_times) < original_len:
             print(f"Warning: {original_len - len(df_times)} rows with non-numeric times were excluded from the plot.")


    if df_times.empty:
         print("Warning: No valid time data found. Skipping time plot.")
         return

    plt.figure(figsize=(8, 6))
    sns.boxplot(data=df_times, x='Method', y='Time', palette='viridis')
    plt.title('Comparison of Execution Times')
    plt.ylabel('Time (milliseconds)') # Updated label
    plt.xlabel('Search Method')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '6_time_comparison_boxplot.png')) # Renamed file
    plt.close()
    print("Plot 6 saved.")

# --- Generate Plots ---
# Updated calls to use known p/q and added heatmap plot call
# plot_order_hit_percentage(df_gs, df_sw, KNOWN_P_TRUE, KNOWN_Q_TRUE, P_ESTIMATED_COL, Q_ESTIMATED_COL, OUTPUT_DIR)
# plot_order_selection_heatmap(df_gs, df_sw, P_ESTIMATED_COL, Q_ESTIMATED_COL, OUTPUT_DIR) # Now generates 2a, 2b
# plot_order_differences(df_gs, df_sw, KNOWN_P_TRUE, KNOWN_Q_TRUE, P_ESTIMATED_COL, Q_ESTIMATED_COL, OUTPUT_DIR) # Now generates 3a, 3b
# plot_aic_difference(df_gs, df_sw, AIC_COL, OUTPUT_DIR)
# plot_bic_difference(df_gs, df_sw, BIC_COL, OUTPUT_DIR)
plot_times(df_gs, df_sw, TIME_COL, OUTPUT_DIR)

print(f"\nAll plots generated and saved in the '{OUTPUT_DIR}' directory.")
