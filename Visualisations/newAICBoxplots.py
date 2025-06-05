# %% Import Libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import traceback # Kept for unexpected errors
from matplotlib.lines import Line2D # For custom legend

# %% Core Plotting Function (Simplified, No Mean Value Text)
def create_aic_difference_plot_for_scenario_simplified(
    best_models_filepath,
    scenario_name,
    df_ref_aic_full,
    output_dir="ingarch_aic_difference_boxplots",
    method_rename_map=None,
    custom_palette=None): # Removed show_mean_value_on_plot
    """
    Creates a box plot of AIC differences for a single scenario (simplified version).
    Only the mean dot is shown, not the numerical value.
    """
    REFERENCE_P_ORDER = 2
    REFERENCE_Q_ORDER = 6

    os.makedirs(output_dir, exist_ok=True)

    try:
        df_best_models = pd.read_csv(best_models_filepath)
    except Exception as e:
        print(f"Error reading best models file '{best_models_filepath}': {e}. Skipping scenario '{scenario_name}'.")
        return

    if df_best_models.empty:
        print(f"No data in '{best_models_filepath}'. Skipping scenario '{scenario_name}'.")
        return

    df_best_models.rename(columns={'aic': 'aic_best_method', 'AIC': 'aic_best_method'}, inplace=True)
    if not all(col in df_best_models.columns for col in ['method', 'sim_id', 'aic_best_method']):
        print(f"Missing essential columns in '{best_models_filepath}'. Need 'method', 'sim_id', 'aic_best_method'. Skipping.")
        return

    df_best_models['method'] = df_best_models['method'].astype(str).str.lower().str.strip()
    df_best_models['sim_id'] = pd.to_numeric(df_best_models['sim_id'], errors='coerce')
    df_best_models['aic_best_method'] = pd.to_numeric(df_best_models['aic_best_method'], errors='coerce')
    df_best_models.dropna(subset=['method', 'sim_id', 'aic_best_method'], inplace=True)
    if df_best_models.empty:
        print(f"No valid data after cleaning in '{best_models_filepath}'. Skipping.")
        return
    df_best_models['sim_id'] = df_best_models['sim_id'].astype(int)

    if method_rename_map:
        df_best_models['display_method'] = df_best_models['method'].map(method_rename_map).fillna(df_best_models['method'])
    else:
        df_best_models['display_method'] = df_best_models['method']
    
    df_ref_scenario = df_ref_aic_full[df_ref_aic_full['scenario'] == scenario_name].copy()
    if df_ref_scenario.empty:
        print(f"No reference AIC data for scenario '{scenario_name}'. Skipping.")
        return

    merged_df = pd.merge(
        df_best_models[['sim_id', 'display_method', 'aic_best_method']],
        df_ref_scenario[['sim_id', 'aic_ref_2_6']],
        on='sim_id',
        how='inner'
    )

    if merged_df.empty or not all(col in merged_df.columns for col in ['aic_ref_2_6', 'aic_best_method']):
        print(f"Merge resulted in empty data or missing AIC columns for '{scenario_name}'. Skipping.")
        return
    merged_df.dropna(subset=['aic_ref_2_6', 'aic_best_method'], inplace=True)
    if merged_df.empty :
        print(f"No valid AIC pairs after NaNs drop for '{scenario_name}'. Skipping.")
        return
        
    merged_df['aic_difference'] = merged_df['aic_ref_2_6'] - merged_df['aic_best_method']

    if merged_df['display_method'].nunique() == 0:
        print(f"No methods to plot for '{scenario_name}'. Skipping.")
        return

    # --- Plotting Aesthetics ---
    all_methods_in_plot_data = sorted(merged_df['display_method'].unique())
    preferred_order_list = ["Stepwise", "Grid"] 
    plot_order = [m for m in preferred_order_list if m in all_methods_in_plot_data]
    plot_order.extend([m for m in all_methods_in_plot_data if m not in plot_order])

    if not plot_order:
        print(f"No plot order for '{scenario_name}'. Skipping.")
        return

    active_palette = None
    if custom_palette:
        active_palette = {k: v for k, v in custom_palette.items() if k in plot_order}
    
    fig_width = max(6, len(plot_order) * 1.5 if len(plot_order) > 1 else 5)
    plt.figure(figsize=(fig_width, 6))

    ax = sns.boxplot(x='display_method', y='aic_difference', hue='display_method',
                     data=merged_df, palette=active_palette, order=plot_order,
                     width=0.6, dodge=False, legend=False)

    # Plot mean dot (without numerical annotation)
    for i, method_name_in_plot in enumerate(plot_order):
        method_data = merged_df[merged_df['display_method'] == method_name_in_plot]['aic_difference']
        if not method_data.empty:
            mean_val = method_data.mean()
            if pd.notnull(mean_val):
                ax.plot(i, mean_val, marker='o', color='red', markersize=7, zorder=5, linestyle='None')
                # Numerical annotation (ax.annotate call) has been removed

    plt.title(f'AIC Difference (Ref P={REFERENCE_P_ORDER}, Q={REFERENCE_Q_ORDER})\nScenario: {scenario_name}', fontsize=12)
    plt.xlabel('Method', fontsize=10)
    plt.ylabel(f'AIC(Ref) - AIC(Other)', fontsize=10)
    plt.axhline(0, color='grey', linestyle='--', linewidth=0.8)

    mean_dot_legend_entry = Line2D([0], [0], marker='o', color='w', 
                                   markerfacecolor='red', markersize=7, label='Mean')
    ax.legend(handles=[mean_dot_legend_entry], loc='best', fontsize=9)

    plt.grid(axis='y', linestyle=':', alpha=0.6)
    if len(plot_order) > 2:
        plt.xticks(rotation=20, ha='right', fontsize=9)
    else:
        plt.xticks(rotation=0, fontsize=9)
    plt.yticks(fontsize=9)
    plt.tight_layout()

    safe_scenario_name = "".join(c if c.isalnum() or c in ('_', '-') else '' for c in scenario_name).replace(' ', '_')
    filename = f"aic_diff_{safe_scenario_name.lower()}_ref_p{REFERENCE_P_ORDER}q{REFERENCE_Q_ORDER}.png"
    save_path = os.path.join(output_dir, filename)
    
    try:
        plt.savefig(save_path, dpi=250)
        print(f"Plot for '{scenario_name}' saved to: {save_path}")
    except Exception as e_save:
        print(f"Error saving plot {save_path}: {e_save}")
    plt.close()

# %% --- Main Execution ---
if __name__ == '__main__':
    print("--- AIC Difference Boxplot Generation (Simplified - Mean Dot Only) ---")

    # --- Configuration ---
    reference_aic_master_file = r"C:\Users\Rafael\Desktop\msc-thesis\Visualisations\Models_2_6_fitted_foreachscenario_by_simulation.csv"
    plot_output_main_directory = "ingarch_aic_difference_boxplots"

    user_method_display_names = {
        "stepwise": "Stepwise",
        "grid_search": "Grid",
        "grid": "Grid" 
    }
    user_plot_palette = {
        "Stepwise": "mediumseagreen", 
        "Grid": "#377E7F" 
    }
    # show_numerical_mean_on_plot variable removed as it's no longer used by the function

    # --- Load Master Reference AIC Data ---
    df_master_ref_aics = None
    if not os.path.exists(reference_aic_master_file):
        print(f"CRITICAL ERROR: Master reference AIC file not found: '{reference_aic_master_file}'.")
    else:
        try:
            df_master_ref_aics = pd.read_csv(reference_aic_master_file)
            df_master_ref_aics.rename(columns={'aic': 'aic_ref_2_6', 'AIC': 'aic_ref_2_6'}, inplace=True)
            if not all(col in df_master_ref_aics.columns for col in ['scenario', 'sim_id', 'aic_ref_2_6']):
                raise ValueError("Reference AIC file missing 'scenario', 'sim_id', or 'aic_ref_2_6' column.")
            
            df_master_ref_aics['sim_id'] = pd.to_numeric(df_master_ref_aics['sim_id'], errors='coerce')
            df_master_ref_aics.dropna(subset=['scenario', 'sim_id', 'aic_ref_2_6'], inplace=True)
            df_master_ref_aics['sim_id'] = df_master_ref_aics['sim_id'].astype(int)
            print(f"Loaded master reference AICs from: {reference_aic_master_file}")
        except Exception as e_load_ref:
            print(f"CRITICAL ERROR loading reference AIC file: {e_load_ref}")
            df_master_ref_aics = None

    # --- Define All Possible Scenarios and Their Data Files ---
    all_scenario_data_map = {
        "M1": r"C:\Users\Rafael\Desktop\msc-thesis\Results\Simulation\no_covariates\Model1\ingarch_no_covariates_results.csv",
        "M2": r"C:\Users\Rafael\Desktop\msc-thesis\Results\Simulation\no_covariates\Model2\ingarch_no_covariates_results.csv",
        "with_covariates": r"C:\Users\Rafael\Desktop\msc-thesis\Results\Simulation\with_covariates\Combined_cpu3\ingarch_with_covariates_results.csv"
    }

    # --- CHOOSE WHICH SCENARIO(S) TO RUN ---
    # scenarios_to_process = ["M1"] 
    # Or uncomment and set one of these:
    #scenarios_to_process = ["M2"]
    scenarios_to_process = ["with_covariates"]
    # scenarios_to_process = ["M1", "M2", "with_covariates"] # To run all

    if df_master_ref_aics is not None and not df_master_ref_aics.empty:
        for scenario_key in scenarios_to_process:
            if scenario_key in all_scenario_data_map:
                best_models_file = all_scenario_data_map[scenario_key]
                print(f"\nProcessing: Scenario '{scenario_key}', File '{best_models_file}'")
                try:
                    create_aic_difference_plot_for_scenario_simplified(
                        best_models_filepath=best_models_file,
                        scenario_name=scenario_key,
                        df_ref_aic_full=df_master_ref_aics,
                        output_dir=plot_output_main_directory,
                        method_rename_map=user_method_display_names,
                        custom_palette=user_plot_palette
                        # show_mean_value_on_plot argument removed from call
                    )
                except Exception as e_main_loop:
                    print(f"!! Error during processing of scenario '{scenario_key}' !!")
                    print(f"Details: {e_main_loop}")
            else:
                print(f"Warning: Scenario key '{scenario_key}' not found in 'all_scenario_data_map'. Skipping.")
    else:
        print("Master reference AIC data not loaded or empty. Cannot generate plots.")

    print("\n--- Plot generation attempt complete. ---")