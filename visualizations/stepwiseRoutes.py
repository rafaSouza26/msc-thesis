import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import json
import os

def load_json_summaries(json_filepath):
    """Loads the JSON summary file which is an array of run objects."""
    try:
        with open(json_filepath, 'r') as f:
            summaries_list = json.load(f)
        # Convert list of dicts to a dict keyed by run_id for easier lookup
        summaries_dict = {str(item['run_id']): item for item in summaries_list if 'run_id' in item}
        return summaries_dict
    except FileNotFoundError:
        print(f"Error: JSON summary file '{json_filepath}' not found.")
        return None
    except Exception as e:
        print(f"Error reading or parsing JSON summary file: {e}")
        return None

def get_final_orders_for_title(summary_data, run_id):
    """Extracts final p and q from the summary data for a given run_id."""
    if summary_data and run_id in summary_data:
        run_info = summary_data[run_id]
        if run_info and 'final_model_orders' in run_info and \
           isinstance(run_info['final_model_orders'], dict):
            p_final = run_info['final_model_orders'].get('p', 'N/A')
            q_final = run_info['final_model_orders'].get('q', 'N/A')
            if p_final is None: p_final = 'N/A' # Handle nulls from JSON
            if q_final is None: q_final = 'N/A' # Handle nulls from JSON
            return p_final, q_final
    return 'N/A', 'N/A'


def plot_single_stepwise_path(path_df_full, summary_data, target_run_id, 
                              p_max=7, q_max=7, output_filename="stepwise_path.png"):
    """
    Plots a single stepwise path for a specific run_id and saves it to a file.

    Args:
        path_df_full (pd.DataFrame): DataFrame containing all path data.
        summary_data (dict): Dictionary of run summaries keyed by run_id.
        target_run_id (str): The 'run_id' string to plot.
        p_max (int): Maximum order of p for the plot axes.
        q_max (int): Maximum order of q for the plot axes.
        output_filename (str): Filename to save the plot.
    """
    run_data = path_df_full[path_df_full['run_id'] == target_run_id].sort_values(by='step')

    if run_data.empty:
        print(f"No path data found in CSV for run_id: '{target_run_id}'. Skipping plot.")
        return

    final_p, final_q = get_final_orders_for_title(summary_data, target_run_id)
    plot_dynamic_title = f"Stepwise Search Path ({final_p},{final_q})"
    
    # Starting p,q for this run_id (from the CSV, added by R script)
    start_p_for_run = run_data['start_p_for_run'].iloc[0]
    start_q_for_run = run_data['start_q_for_run'].iloc[0]
    ic_used_for_run = run_data['ic_used_for_run'].iloc[0]

    fig, ax = plt.subplots(figsize=(7, 6)) # Create a new figure and axes for each plot

    p_values = run_data['p']
    q_values = run_data['q']
    steps = run_data['step']

    ax.plot(p_values, q_values, marker='o', linestyle='-', markersize=5, label='Path Taken', zorder=2)

    if not run_data.empty:
        # Start point of the path (first point in the CSV for this run_id)
        ax.plot(p_values.iloc[0], q_values.iloc[0], marker='^', color='green', markersize=10, 
                label=f"Path Start ({p_values.iloc[0]},{q_values.iloc[0]})", zorder=3)
        # End point of the path recorded in CSV
        ax.plot(p_values.iloc[-1], q_values.iloc[-1], marker='s', color='red', markersize=8, 
                label=f"Path End ({p_values.iloc[-1]},{q_values.iloc[-1]})", zorder=3)
        # Final selected model point (might be different from path end)
        if final_p != 'N/A' and final_q != 'N/A':
             ax.plot(final_p, final_q, marker='*', color='gold', markersize=12, 
                     label=f"Selected Model ({final_p},{final_q})", zorder=4, mec='black')


        for j, step_num in enumerate(steps):
            ax.text(p_values.iloc[j] + 0.15, q_values.iloc[j] + 0.15, str(int(step_num)), 
                    fontsize=7, ha='left', va='bottom', zorder=4)
    
    ax.set_title(plot_dynamic_title, fontsize=14)
    
    # Subtitle with more run info
    run_details_subtitle = f"Run ID: {target_run_id} (IC: {ic_used_for_run}, Start: p={start_p_for_run}, q={start_q_for_run})"
    fig.text(0.5, 0.90, run_details_subtitle, ha='center', fontsize=9, style='italic') # Add subtitle below main title

    ax.set_xlabel("Order p")
    ax.set_ylabel("Order q")
    
    ax.set_xlim(-0.5, p_max + 0.5)
    ax.set_ylim(-0.5, q_max + 0.5)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.grid(True, which='both', linestyle=':', linewidth=0.7, zorder=1)
    ax.set_aspect('equal', adjustable='box')
    ax.legend(fontsize=8, loc='best')

    plt.tight_layout(rect=[0, 0, 1, 0.88]) # Adjust layout to make space for title and subtitle
    
    try:
        plt.savefig(output_filename)
        print(f"Plot saved to '{output_filename}'")
    except Exception as e:
        print(f"Error saving plot '{output_filename}': {e}")
    plt.close(fig) # Close the figure to free memory, important for multiple plots

# --- Main part of the script ---
if __name__ == "__main__":
    # --- Configuration ---
    # Ensure these file paths are correct
    csv_paths_file = r"C:\Users\Rafael\Desktop\msc-thesis\all_stepwise_paths_simplified_ids.csv"
    json_summaries_file = r"C:\Users\Rafael\Desktop\msc-thesis\all_model_summaries_simplified_ids.json"
    output_plot_dir = "stepwise_plots" # Directory to save plots

    # Create output directory if it doesn't exist
    if not os.path.exists(output_plot_dir):
        os.makedirs(output_plot_dir)
        print(f"Created directory: {output_plot_dir}")

    # --- User Action Required ---
    # 1. Run your R script to generate/update the CSV and JSON files.
    # 2. Open 'all_model_summaries_simplified_ids.json'.
    # 3. For each run_id, find the 'final_model_orders' (p,q).
    # 4. Identify three 'run_id's (simple numbers as strings, e.g., "1", "5", "23") that correspond to:
    #    a) Your "true model" (e.g., final p=2, q=6)
    #    b) A model "close" to the true model
    #    c) A model "completely off" from the true model
    # 5. Update the 'run_ids_for_cases' dictionary below.
    #    The keys are descriptive names for your output files.
    #    The values are the 'run_id' strings from your JSON/CSV.
    # =========================================================================
    run_ids_for_cases = {
        "true_model_case": "PLEASE_REPLACE_WITH_RUN_ID_FOR_TRUE_MODEL",
        "close_model_case": "1",
        "off_model_case": "11"
    }
    # =========================================================================

    # Load all path data and summaries once
    all_paths_df = pd.read_csv(csv_paths_file)
    all_summaries = load_json_summaries(json_summaries_file)

    if all_paths_df is None or all_summaries is None:
        print("Could not load necessary data. Exiting.")
    else:
        plotted_any = False
        for case_name, run_id_to_plot in run_ids_for_cases.items():
            if "PLEASE_REPLACE" in run_id_to_plot:
                print(f"WARNING: Placeholder run_id found for '{case_name}'. Please update it with an actual run_id.")
                print(f"  To help, here are the first few available run_ids from your JSON summary: {list(all_summaries.keys())[:10]}")
                continue # Skip this case if placeholder is still there

            if run_id_to_plot not in all_summaries:
                print(f"WARNING: run_id '{run_id_to_plot}' for case '{case_name}' not found in JSON summaries. Skipping plot.")
                print(f"  Available run_ids: {list(all_summaries.keys())}")
                continue
            
            if run_id_to_plot not in all_paths_df['run_id'].unique():
                 print(f"WARNING: run_id '{run_id_to_plot}' for case '{case_name}' not found in CSV path data. Skipping plot.")
                 print(f"  Available run_ids in CSV: {all_paths_df['run_id'].unique().tolist()}")
                 continue


            output_file = os.path.join(output_plot_dir, f"path_{case_name}_run_{run_id_to_plot}.png")
            print(f"\nGenerating plot for: {case_name} (Run ID: {run_id_to_plot}) -> {output_file}")
            
            plot_single_stepwise_path(
                path_df_full=all_paths_df,
                summary_data=all_summaries,
                target_run_id=run_id_to_plot,
                p_max=7, 
                q_max=7,
                output_filename=output_file
            )
            plotted_any = True
        
        if not plotted_any and any("PLEASE_REPLACE" in rid for rid in run_ids_for_cases.values()):
             print("\nNo plots were generated because 'run_ids_for_cases' still contains placeholders.")
             print("Please edit the script to specify the correct run_ids based on your R script's JSON output.")