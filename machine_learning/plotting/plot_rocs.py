"""
Plot the ROC curve for a given model and dataset.

- reads the ROC curve data from a TSV file
- plots a ROC curve
- writes a png 
"""

# Imports ------------------------------------------
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import auc
import os

# Options ------------------------------------------
ANALYSIS_GROUPS = ["stool_prospective", "stool_active", "stool_treated", "duodenum_active"]
EVALUATION_METHODS = ["kfold", "lodo"]
BEST_MODEL = True
MODEL_NAME = "rf" # <- used only if BEST_MODEL is False

# Read summary AUCs --------------------------------
SUMMARY_PATH = '/home/haig/Repos/meta-analysis/machine_learning/results/main_results_summary.tsv'
summary_df = pd.read_csv(SUMMARY_PATH, sep='\t')
KFOLD_AUC_MAP = dict(zip(summary_df['Dataset'], summary_df['KFOLD AUC']))
LODO_AUC_MAP = dict(zip(summary_df['Dataset'], summary_df['LODO AUC']))

# Loop through all combinations
for EVALUATION_METHOD in EVALUATION_METHODS:
    for ANALYSIS_GROUP in ANALYSIS_GROUPS:
        print()

        if BEST_MODEL:
            # Read summary_best_model.tsv to get the model name
            summary_best_model_path = '/home/haig/Repos/meta-analysis/machine_learning/datasets_main/' + ANALYSIS_GROUP + '_tss_after/' + EVALUATION_METHOD + '_results/summary_best_model.tsv'
            summary_best_model_df = pd.read_csv(summary_best_model_path, sep='\t')
            # Get 2nd column name
            MODEL_NAME = summary_best_model_df.columns[1].lower()

        TSV_FILE = '/home/haig/Repos/meta-analysis/machine_learning/datasets_main/' + ANALYSIS_GROUP + '_tss_after/' + EVALUATION_METHOD + '_results/best_models/roc_' + MODEL_NAME + '_data.tsv'
        TITLE = ANALYSIS_GROUP + ' ' + MODEL_NAME.upper() + ' ' + EVALUATION_METHOD.upper() + ' ROC Curve'
        OUT_FILE = TITLE.lower().replace(' ', '_').replace('-', '_') + '.png'

        OUTPUT_DIR = '/home/haig/Repos/meta-analysis/machine_learning/results/roc_plots/'

        # Main ------------------------------------------
        # Read the data
        data = pd.read_csv(TSV_FILE, sep='	')
        fpr = data['fpr']
        tpr = data['tpr']

        # Get AUC from summary
        dataset_name = ANALYSIS_GROUP + '_tss_after'
        roc_auc = KFOLD_AUC_MAP[dataset_name] if EVALUATION_METHOD == 'kfold' else LODO_AUC_MAP[dataset_name]

        # Create plot
        plt.figure(figsize=(6, 6))

        # Plot the ROC curve
        plt.plot(fpr, tpr, color='red', lw=2)

        # Plot the diagonal line for reference
        plt.plot([0, 1], [0, 1], color='grey', lw=1, linestyle='--')

        # Set labels and title
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(TITLE)

        # Set axis limits
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.0])

        # Make it square
        plt.gca().set_aspect('equal', adjustable='box')

        # Display AUC value bottom right
        plt.text(0.95, 0.05, f'AUC = {roc_auc:.3f}',
                verticalalignment='bottom', horizontalalignment='right',
                transform=plt.gca().transAxes,
                fontsize=12)

        # Save the figure
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        output_path = os.path.join(OUTPUT_DIR, OUT_FILE)
        plt.savefig(output_path, bbox_inches='tight')

        print(f"Plot saved to {output_path}")

    # First set for this evaluation method processed
    # Now plot the ROC curves for all analysis groups of this evaluation method on one nice plot
    # Also write a TSV file with the data for the plot

    # Combined ------------------------------------------
    # Build a single ROC plot across all analysis groups for this evaluation method
    evaluation_method_name = "K-fold" if EVALUATION_METHOD == "kfold" else "LODO"
    COMBINED_TITLE = evaluation_method_name + ' ROC Curves for all analysis groups'
    COMBINED_OUT_FILE = COMBINED_TITLE.lower().replace(' ', '_').replace('-', '_') + '.png'
    COMBINED_TSV_FILE = COMBINED_TITLE.lower().replace(' ', '_').replace('-', '_') + '_data.tsv'

    plt.figure(figsize=(6, 6))

    combined_rows = []

    for ANALYSIS_GROUP in ANALYSIS_GROUPS:
        if BEST_MODEL:
            # Read best model per analysis group for this evaluation method
            summary_best_model_path = '/home/haig/Repos/meta-analysis/machine_learning/datasets_main/' + ANALYSIS_GROUP + '_tss_after/' + EVALUATION_METHOD + '_results/summary_best_model.tsv'
            summary_best_model_df = pd.read_csv(summary_best_model_path, sep='\t')
            model_name_for_group = summary_best_model_df.columns[1].lower()
        else:
            model_name_for_group = MODEL_NAME

        # Read ROC data for this analysis group and model
        tsv_file = '/home/haig/Repos/meta-analysis/machine_learning/datasets_main/' + ANALYSIS_GROUP + '_tss_after/' + EVALUATION_METHOD + '_results/best_models/roc_' + model_name_for_group + '_data.tsv'
        data = pd.read_csv(tsv_file, sep='\t')
        fpr = data['fpr']
        tpr = data['tpr']

        pretty_analysis_group_name = ANALYSIS_GROUP.replace('_', ' ').title()

        # Get AUC from summary and plot
        dataset_name = ANALYSIS_GROUP + '_tss_after'
        roc_auc = KFOLD_AUC_MAP[dataset_name] if EVALUATION_METHOD == 'kfold' else LODO_AUC_MAP[dataset_name]
        plt.plot(fpr, tpr, lw=2, label=pretty_analysis_group_name + ' (AUC=' + f'{roc_auc:.3f}' + ')')

        # Collect rows for combined TSV
        group_df = pd.DataFrame({
            'analysis_group': pretty_analysis_group_name,
            'model': model_name_for_group,
            'fpr': fpr,
            'tpr': tpr,
            'auc': round(roc_auc, 3)
        })
        combined_rows.append(group_df)

    # Plot the diagonal line for reference
    plt.plot([0, 1], [0, 1], color='grey', lw=1, linestyle='--')

    # Set labels and title
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(COMBINED_TITLE)

    # Set axis limits
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])

    # Make it square
    plt.gca().set_aspect('equal', adjustable='box')

    # Legend
    plt.legend(loc='lower right', fontsize=9)

    # Save the figure
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    combined_output_path = os.path.join(OUTPUT_DIR, COMBINED_OUT_FILE)
    plt.savefig(combined_output_path, bbox_inches='tight')
    print(f"Combined plot saved to {combined_output_path}")

    # Save the combined TSV data
    combined_df = pd.concat(combined_rows, ignore_index=True)
    combined_tsv_path = os.path.join(OUTPUT_DIR, COMBINED_TSV_FILE)
    combined_df.to_csv(combined_tsv_path, sep='\t', index=False)
    print(f"Combined ROC data written to {combined_tsv_path}")
