import csv
import re
import matplotlib.pyplot as plt
# import matplotlib.patches as patches # Rectangles no longer needed for path cells
import matplotlib.lines as mlines # For custom legend handles
import numpy as np
import os

# Helper function to parse path strings (remains the same)
def parse_path_string(path_str):
    coordinates = []
    steps = path_str.split(' -> ')
    for step in steps:
        match = re.search(r'\((\d+),\s*(\d+)\)', step)
        if match:
            p_val = int(match.group(1))
            q_val = int(match.group(2))
            coordinates.append((p_val, q_val))
    return coordinates

# MODIFIED Function to plot a single path on a grid
def plot_path_on_grid(path_coords, 
                        plot_title_custom, 
                        true_model_coords,   
                        max_order, 
                        output_filename):
    if not path_coords:
        print(f"Path data is empty. Skipping plot: {output_filename}")
        return

    fig, ax = plt.subplots(figsize=(8, 8.5)) 

    # Grid and axis setup
    ax.set_xlim(-0.5, max_order + 0.5)
    ax.set_ylim(-0.5, max_order + 0.5)
    ax.set_xticks(np.arange(max_order + 1))
    ax.set_yticks(np.arange(max_order + 1))
    ax.set_xlabel("p order", fontsize=12)
    ax.set_ylabel("q order", fontsize=12)
    ax.set_title(plot_title_custom, pad=20, fontsize=14) 
    ax.grid(True, linestyle='-', color='lightgrey', zorder=0)
    ax.set_aspect('equal', adjustable='box')

    # Plot dots for the path states
    for i, (p, q) in enumerate(path_coords):
        if 0 <= p <= max_order and 0 <= q <= max_order:
            dot_color = 'k' # Black for intermediate path points
            dot_size = 6
            dot_zorder = 2 # Path dots

            if i == 0: # Start point
                dot_color = 'g' # Green
                dot_size = 7
            elif i == len(path_coords) - 1: # End point
                dot_color = 'b' # Blue
                dot_size = 7
            
            ax.plot(p, q, marker='o', color=dot_color, markersize=dot_size, zorder=dot_zorder, linestyle='None')

    # Draw arrows between consecutive path states
    for i in range(len(path_coords) - 1):
        p_curr, q_curr = path_coords[i]
        p_next, q_next = path_coords[i+1]
        if (0 <= p_curr <= max_order and 0 <= q_curr <= max_order and
            0 <= p_next <= max_order and 0 <= q_next <= max_order):
            ax.annotate("",
                        xy=(p_next, q_next),
                        xytext=(p_curr, q_curr),
                        arrowprops=dict(arrowstyle="->", color="black", lw=1.5, mutation_scale=15), # Reduced mutation_scale for smaller arrowheads
                        zorder=1) # Arrows below dots

    # Plot "True Model" as a red dot if provided
    if true_model_coords:
        p_true, q_true = true_model_coords
        if 0 <= p_true <= max_order and 0 <= q_true <= max_order:
            ax.plot(p_true, q_true, 'ro', markersize=8, zorder=3, label='_nolegend_True Model Dot') # Red dot, higher zorder

    # Create legend handles
    legend_handles = [
        mlines.Line2D([], [], color='g', marker='o', linestyle='None', markersize=7, label='Start Point'),
        mlines.Line2D([], [], color='b', marker='o', linestyle='None', markersize=7, label='End Point'),
        mlines.Line2D([], [], color='k', marker='o', linestyle='None', markersize=6, label='Path Point')
    ]

    if true_model_coords: 
        legend_handles.append(
            mlines.Line2D([], [], color='r', marker='o', linestyle='None', markersize=8, label='True Model')
        )
    
    fig.legend(handles=legend_handles, loc='lower center', bbox_to_anchor=(0.5, 0.005), ncol=len(legend_handles), fancybox=True, fontsize=9) # Dynamic ncol
    plt.tight_layout(rect=[0, 0.08, 1, 0.95])

    try:
        plt.savefig(output_filename, dpi=300)
        print(f"Plot saved successfully to: {output_filename}")
    except Exception as e:
        print(f"Error saving plot to {output_filename}: {e}")
    plt.close(fig)

# MODIFIED Main function to generate grid path plots (structure remains same)
def generate_grid_path_plots_for_final_models(
    csv_filepath, 
    target_final_models_str_list, 
    plot_titles_list,                  
    true_model_coords_str,             
    output_dir="grid_path_plots" # Changed v3 to v4 for new version
):
    if not os.path.exists(csv_filepath):
        print(f"Error: CSV file not found at {csv_filepath}")
        return

    try:
        os.makedirs(output_dir, exist_ok=True)
        print(f"Output directory set to: {output_dir}")
    except OSError as e:
        print(f"Error creating directory {output_dir}: {e}")
        return

    parsed_true_model_coords = None
    if true_model_coords_str:
        try:
            p_true, q_true = map(int, true_model_coords_str.strip('()').split(','))
            parsed_true_model_coords = (p_true, q_true)
            print(f"True model coordinates set to: {parsed_true_model_coords}")
        except ValueError:
            print(f"Warning: Could not parse true_model_coords_str: '{true_model_coords_str}'. True model will not be plotted.")

    all_parsed_paths_data = []
    try:
        with open(csv_filepath, mode='r', encoding='utf-8-sig') as infile:
            reader = csv.DictReader(infile)
            if "Simulation_ID" not in reader.fieldnames or "Path_Taken" not in reader.fieldnames:
                print("Error: CSV must contain 'Simulation_ID' and 'Path_Taken' columns.")
                return
            for row in reader:
                sim_id = row.get("Simulation_ID")
                path_str = row.get("Path_Taken")
                if sim_id and path_str:
                    coordinates = parse_path_string(path_str)
                    if coordinates:
                        all_parsed_paths_data.append({'id': sim_id, 'path': coordinates})
    except Exception as e:
        print(f"Error reading or parsing CSV {csv_filepath}: {e}")
        return

    if not all_parsed_paths_data:
        print("No valid paths found in the CSV file.")
        return

    max_order_grid = 7
    plots_generated_for_model_type = {} 

    target_models_tuples = []
    for model_str in target_final_models_str_list:
        try:
            p, q = map(int, model_str.strip('()').split(','))
            target_models_tuples.append((p, q))
        except ValueError:
            print(f"Warning: Could not parse target model string: {model_str}.")
    
    if len(plot_titles_list) != len(target_models_tuples):
        print("Warning: Number of plot titles does not match number of target models. Default titles will be used where needed.")


    for i, target_model_tuple in enumerate(target_models_tuples):
        if target_model_tuple in plots_generated_for_model_type:
            continue

        current_plot_title = f"Path to Target Model {target_model_tuple}" 
        if i < len(plot_titles_list):
            current_plot_title = plot_titles_list[i]
        else:
            print(f"Using default title for target {target_model_tuple} as custom title was not provided for this index.")


        found_path_for_this_target = False
        for item in all_parsed_paths_data:
            path_coords = item['path']
            if not path_coords:
                continue
            
            actual_final_model_in_path = path_coords[-1]

            if actual_final_model_in_path == target_model_tuple:
                sim_id = item['id']
                print(f"\nGenerating plot for Simulation '{sim_id}' ending in target {target_model_tuple}.")
                print(f"  Using Title: {current_plot_title}")
                print(f"  Path: {path_coords}")
                
                safe_target_model_str = f"{target_model_tuple[0]}_{target_model_tuple[1]}"
                filename_sim_id_part = sim_id.replace("run_", "") 
                output_filename = os.path.join(output_dir, f"path_plot_target_{safe_target_model_str}_sim_{filename_sim_id_part}.png")
                
                plot_path_on_grid(
                    path_coords=path_coords, 
                    plot_title_custom=current_plot_title,
                    true_model_coords=parsed_true_model_coords,
                    max_order=max_order_grid, 
                    output_filename=output_filename
                )
                plots_generated_for_model_type[target_model_tuple] = True
                found_path_for_this_target = True
                break 

        if not found_path_for_this_target:
            print(f"\nNo paths found in the CSV ending specifically at {target_model_tuple} to generate a plot.")

# --- Main execution ---
if __name__ == "__main__":
    csv_file_path = r"C:\Users\Rafael\Desktop\msc-thesis\old\results\simulatedDataResults\auto_ingarch_stepwise_paths.csv" 
    
    target_models_list_str = ["(2,6)", "(2,0)", "(2,5)", "(4,7)"] 
    
    custom_titles = [
        "Stepwise Search Path to Model (2,6)",
        "Stepwise Search Path to Model (2,0)",
        "Stepwise Search Path to Model (2,5)",
        "Stepwise Search Path to Model (4,7)"
    ]

    true_model_specification_str = "(2,6)"
    
    output_plots_directory = "ingarch_sw_grid_path_plots"

    print(f"Attempting to read CSV file from: {csv_file_path}")
    print(f"Target final models for plotting: {target_models_list_str}")
    print(f"Custom titles for plots: {custom_titles}")
    print(f"True model specification: {true_model_specification_str}")
    print(f"Output directory for plots: {output_plots_directory}")

    generate_grid_path_plots_for_final_models(
        csv_filepath=csv_file_path, 
        target_final_models_str_list=target_models_list_str, 
        plot_titles_list=custom_titles,
        true_model_coords_str=true_model_specification_str,
        output_dir=output_plots_directory
    )
    
    print(f"\nProcessing complete. Grid path plots should be in the '{output_plots_directory}' directory.")