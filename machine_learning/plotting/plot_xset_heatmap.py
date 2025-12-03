"""
Automatically plot the cross-dataset validation results heatmaps for all datasets.

- reads the cross-dataset validation results from a TSV file
- plots a heatmap of AUC scores where rows are training sets and columns are test sets
- writes a png and an svg file

"""

# Imports ------------------------------------------
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# Parameters ------------------------------------------
# If True, MODEL_NAME is ignored and best model is used
SHOW_BEST_MODEL = True
# Only used if SHOW_BEST_MODEL is False
MODEL_NAME = "xgboost"
# Output directory for heatmaps
OUTPUT_DIR_PATH = "/home/haig/Repos/meta-analysis/machine_learning/results/xset_heatmaps/"
# Dataset directories to plot heatmaps for
ALL_DATASET_DIR_PATHS = [
    "/home/haig/Repos/meta-analysis/machine_learning/datasets_main/duodenum_active_tss_after",
    "/home/haig/Repos/meta-analysis/machine_learning/datasets_main/stool_active_tss_after",
    "/home/haig/Repos/meta-analysis/machine_learning/datasets_main/stool_prospective_tss_after",
    "/home/haig/Repos/meta-analysis/machine_learning/datasets_main/stool_treated_tss_after",
]


# Main ------------------------------------------
for DATASET_DIR_PATH in ALL_DATASET_DIR_PATHS:

    DATASET_NAME = os.path.basename(DATASET_DIR_PATH)

    # If showing best model, need to identify it from best_models_results.tsv
    if SHOW_BEST_MODEL:
        PERFORMANCE_RESULTS_SUBDIR_FILE_NAME = "lodo_results/best_models_results.tsv"
        PERFORMANCE_RESULTS_PATH = os.path.join(DATASET_DIR_PATH, PERFORMANCE_RESULTS_SUBDIR_FILE_NAME)
        # Read performance results and get best model
        df = pd.read_csv(PERFORMANCE_RESULTS_PATH, sep='\t')
        MODEL_NAME = df[df['Replicate'] == 'Average'].loc[df[df['Replicate'] == 'Average']['Mean'].idxmax(), 'Model']

    # Path to cross-dataset validation results
    RESULTS_SUBDIR_FILE_NAME = f"lodo_results/xset_results_{MODEL_NAME}.tsv"
    RESULTS_PATH = os.path.join(DATASET_DIR_PATH, RESULTS_SUBDIR_FILE_NAME)


    # Outputs ------------------------------------------
    FILE_NAME = f"xset_heatmap_{DATASET_NAME}_{MODEL_NAME}"
    TITLE = f"Cross-Dataset Validation Results for {MODEL_NAME.upper()} Model\nin {DATASET_NAME} Dataset"


    # Plot ------------------------------------------
    # Create output directory if it doesn't exist
    if not os.path.exists(OUTPUT_DIR_PATH):
        os.makedirs(OUTPUT_DIR_PATH)

    # Read the cross-dataset validation results
    xset_df = pd.read_csv(RESULTS_PATH, sep='\t', index_col=0)

    # Create the heatmap
    plt.figure(figsize=(10, 8))
    sns.set_style("whitegrid")

    # Create heatmap with values shown in cells
    sns.heatmap(
        xset_df,
        annot=True,  # Show values in cells
        fmt='.3f',   # Format to 3 decimal places
        cmap='YlOrRd',  # Yellow-Orange-Red colormap
        vmin=0.5,    # Minimum value for colormap
        vmax=1.0,    # Maximum value for colormap
        square=True, # Make cells square
        cbar_kws={'label': 'AUC'}  # Add label to colorbar
    )

    # Customize the plot
    plt.title(TITLE)
    plt.xlabel('Test Set')
    plt.ylabel('Training Set')

    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)

    # Adjust layout to prevent label cutoff
    plt.tight_layout()


    # Export ------------------------------------------
    # Save as PNG
    plt.savefig(os.path.join(OUTPUT_DIR_PATH, FILE_NAME + '.png'), dpi=300, bbox_inches='tight')
    # Save as SVG
    # plt.savefig(os.path.join(OUTPUT_DIR_PATH, FILE_NAME + '.svg'), format='svg', bbox_inches='tight')

    # Show plot
    plt.show()

    plt.close()
