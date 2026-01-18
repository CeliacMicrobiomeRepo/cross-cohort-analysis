"""
Automatically plot the results of the best model between some selected datasets.

- plots a box and whisker plot of mean AUCs for each replicate of the best model for each dataset
- writes a png file
"""

# Imports ------------------------------------------
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns



# Inputs ------------------------------------------
# Analysis group name
ANALYSIS_GROUP_NAMES = ["stool_treated", "duodenum_active"]
# Eval method
EVAL_METHODS = ["lodo", "kfold"]

# Draws a red line at the highest median out of all datasets
TOP_MEDIAN_LINE_COLOR = True
# Maximum value for the y-axis
Y_AXIS_MAX = 1.0   # <- Make divisible by 0.05
# Minimum value for the y-axis
Y_AXIS_MIN = 0.5   # <- Make divisible by 0.05


# Outputs ------------------------------------------
OUTPUT_DIR_PATH = "/home/haig/Repos/cross-cohort-analysis/machine_learning/results/transformation_comparison/"

# Create a color mapping for transformations
TRANSFORMATION_COLORS = {
    'TSS': '#1f77b4',  # blue
    'Log10': '#2ca02c',  # green
    'Log10-sum': '#d62728',  # red
    'CLR': '#ff7f0e'  # orange
}


# Helpers ------------------------------------------
def ordinal(n: int) -> str:
    """Return ordinal string for an integer (1 -> 1st, 2 -> 2nd, ...)."""
    n_int = int(n)
    if 10 <= n_int % 100 <= 20:
        suffix = 'th'
    else:
        suffix = {1: 'st', 2: 'nd', 3: 'rd'}.get(n_int % 10, 'th')
    return f"{n_int}{suffix}"


def print_after_tss_ranks(df: pd.DataFrame, context_desc: str = "") -> None:
    """Print ranking stats for the 'After,TSS' transformation within the provided dataframe.

    Expects df with columns: 'Dataset' and 'AUC'.
    """
    target = 'After,TSS'
    if 'Dataset' not in df.columns or 'AUC' not in df.columns:
        return
    grouped = df.groupby('Dataset')['AUC']
    if target not in grouped.groups:
        print(f"{context_desc} 'After,TSS' not found among datasets; skipping rank summary.")
        return

    mins = grouped.min()
    q1s = grouped.quantile(0.25)
    medians = grouped.median()
    q3s = grouped.quantile(0.75)
    maxs = grouped.max()
    means = grouped.mean()

    rank_min = int(mins.rank(ascending=False, method='min').loc[target])
    rank_q1 = int(q1s.rank(ascending=False, method='min').loc[target])
    rank_med = int(medians.rank(ascending=False, method='min').loc[target])
    rank_q3 = int(q3s.rank(ascending=False, method='min').loc[target])
    rank_max = int(maxs.rank(ascending=False, method='min').loc[target])
    rank_mean = int(means.rank(ascending=False, method='min').loc[target])
    
    print(f"{context_desc}The After,TSS transformation had:")
    print(f" - The {ordinal(rank_min)} highest minimum: {mins.loc[target]:.3f}")
    print(f" - The {ordinal(rank_q1)} highest lower quartile: {q1s.loc[target]:.3f}")
    print(f" - The {ordinal(rank_med)} highest median: {medians.loc[target]:.3f}")
    print(f" - The {ordinal(rank_q3)} highest upper quartile: {q3s.loc[target]:.3f}")
    print(f" - The {ordinal(rank_max)} highest maximum: {maxs.loc[target]:.3f}")
    print(f" - The {ordinal(rank_mean)} highest mean: {means.loc[target]:.3f}")


# Accumulator to combine AUCs across all loops (all analysis groups and eval methods)
all_tests_results_list = []

# Loop through all tests
for ANALYSIS_GROUP_NAME in ANALYSIS_GROUP_NAMES:
    for EVAL_METHOD in EVAL_METHODS:
        print()

        # Dataset directory paths and names
        DATASET_DIR_PATHS_AND_NAMES = [
            ("/home/haig/Repos/cross-cohort-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_clr_after/", "After,CLR"),
            ("/home/haig/Repos/cross-cohort-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_clr_before/", "Before,CLR"),
            ("/home/haig/Repos/cross-cohort-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_clr_zscore_after/", "After,CLR,Z"),
            ("/home/haig/Repos/cross-cohort-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_clr_zscore_before/", "Before,CLR,Z"),

            ("/home/haig/Repos/cross-cohort-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_tss_after/", "After,TSS"),
            ("/home/haig/Repos/cross-cohort-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_tss_before/", "Before,TSS"),
            ("/home/haig/Repos/cross-cohort-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_tss_zscore_after/", "After,TSS,Z"),
            ("/home/haig/Repos/cross-cohort-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_tss_zscore_before/", "Before,TSS,Z"),

            ("/home/haig/Repos/cross-cohort-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_log10-sum_after/", "After,Log10-sum"),
            ("/home/haig/Repos/cross-cohort-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_log10-sum_before/", "Before,Log10-sum"),
            ("/home/haig/Repos/cross-cohort-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_log10-sum_zscore_after/", "After,Log10-sum,Z"),
            ("/home/haig/Repos/cross-cohort-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_log10-sum_zscore_before/", "Before,Log10-sum,Z"),

            ("/home/haig/Repos/cross-cohort-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_log10_after/", "After,Log10"),
            ("/home/haig/Repos/cross-cohort-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_log10_before/", "Before,Log10"),
            ("/home/haig/Repos/cross-cohort-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_log10_zscore_after/", "After,Log10,Z"),
            ("/home/haig/Repos/cross-cohort-analysis/machine_learning/datasets_transformation_comparison/" + ANALYSIS_GROUP_NAME + "_log10_zscore_before/", "Before,Log10,Z"),
        ]
        RESULTS_SUBDIR_FILE_NAME = EVAL_METHOD + "_results/best_models_results.tsv"


        # Plot ------------------------------------------
        # Create output directory if it doesn't exist
        if not os.path.exists(OUTPUT_DIR_PATH):
            os.makedirs(OUTPUT_DIR_PATH)

        # Read and process results from each dataset
        dataset_results = []
        for path, name in DATASET_DIR_PATHS_AND_NAMES:
            results_path = os.path.join(path, RESULTS_SUBDIR_FILE_NAME)
            if os.path.exists(results_path):
                df = pd.read_csv(results_path, sep='\t')
                # Identify the best model
                # Get only the average rows and find the best model based on mean AUC
                best_model = df[df['Replicate'] == 'Average'].loc[df[df['Replicate'] == 'Average']['Mean'].idxmax(), 'Model']
                df = df[df['Model'] == best_model]  # Keep only rows for best model
                # Remove the 'Average' rows
                df = df[df['Replicate'] != 'Average']
                # Create a DataFrame with the AUC values
                df_values = pd.DataFrame({
                    'AUC': df['Mean'].values,  # Use the Mean column for AUC values
                    'Dataset': [name] * len(df)  # Repeat dataset name for each row
                })
                dataset_results.append(df_values)
                # Accumulate for the final combined plot across all tests
                all_tests_results_list.append(df_values)
                # print(name, "Best Model:", best_model)

        # Combine all results
        all_results = pd.concat(dataset_results, ignore_index=True)

        # Create the plot
        plt.figure(figsize=(12, 9))
        sns.set_style("whitegrid")

        # Extract transformation type from dataset name (second part after comma)
        all_results['Transformation'] = all_results['Dataset'].apply(lambda x: x.split(',')[1])

        # Create boxplot with colors based on transformation
        ax = sns.boxplot(
            data=all_results,
            x='Dataset',
            y='AUC',
            hue='Transformation',
            width=0.7,
            showfliers=False,
            palette=TRANSFORMATION_COLORS  # Simply pass the color dictionary
        )
        # Remove legend to match previous behavior
        legend = ax.get_legend()
        if legend is not None:
            legend.remove()

        if TOP_MEDIAN_LINE_COLOR:
            medians = all_results.groupby('Dataset')['AUC'].median()
            top_median = medians.max()
            plt.axhline(y=top_median, color='red', linestyle='--', linewidth=1)

        # Customize the plot
        plt.xticks(rotation=45, ha='right')
        plt.xlabel('')
        plt.ylabel('AUC')
        plt.title('Best Model Performance Across Transformations')

        # Set y-axis limits and ticks
        plt.ylim(Y_AXIS_MIN, Y_AXIS_MAX)  # Set y-axis limits
        plt.yticks(np.arange(Y_AXIS_MIN, Y_AXIS_MAX + 0.001, 0.05))

        # Adjust layout to prevent label cutoff
        plt.tight_layout()


        # Export ------------------------------------------
        # Save as PNG
        plt.savefig(os.path.join(OUTPUT_DIR_PATH, EVAL_METHOD + '_' + ANALYSIS_GROUP_NAME + '.png'), dpi=300, bbox_inches='tight')

        # Show plot 
        # plt.show()

        plt.close()

        # Print ranking summary for 'After,TSS' within this specific test
        print_after_tss_ranks(all_results, context_desc=f"[{ANALYSIS_GROUP_NAME} | {EVAL_METHOD}] ")


# Create one more plot which compares all transformations across *all* tests
if len(all_tests_results_list) > 0:
    combined_results = pd.concat(all_tests_results_list, ignore_index=True)
    # Re-declare output and styling to ensure availability outside loops
    OUTPUT_DIR_PATH = "/home/haig/Repos/cross-cohort-analysis/machine_learning/results/transformation_comparison/"
    if not os.path.exists(OUTPUT_DIR_PATH):
        os.makedirs(OUTPUT_DIR_PATH)
    TRANSFORMATION_COLORS = {
        'TSS': '#1f77b4',
        'Log10': '#2ca02c',
        'Log10-sum': '#d62728',
        'CLR': '#ff7f0e'
    }
    Y_AXIS_MAX = 1.0
    Y_AXIS_MIN = 0.5

    plt.figure(figsize=(12, 9))
    sns.set_style("whitegrid")

    combined_results['Transformation'] = combined_results['Dataset'].apply(lambda x: x.split(',')[1])
    ax = sns.boxplot(
        data=combined_results,
        x='Dataset',
        y='AUC',
        hue='Transformation',
        width=0.7,
        showfliers=False,
        palette=TRANSFORMATION_COLORS
    )
    legend = ax.get_legend()
    if legend is not None:
        legend.remove()

    medians = combined_results.groupby('Dataset')['AUC'].median()
    top_median = medians.max()
    plt.axhline(y=top_median, color='red', linestyle='--', linewidth=1)

    plt.xticks(rotation=45, ha='right')
    plt.xlabel('')
    plt.ylabel('AUC')
    plt.title('Best Model Performance Across Transformations (All Tests Combined)')

    plt.ylim(Y_AXIS_MIN, Y_AXIS_MAX)
    plt.yticks(np.arange(Y_AXIS_MIN, Y_AXIS_MAX + 0.001, 0.05))

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR_PATH, 'all_tests_combined.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # Print ranking summary for 'After,TSS' across all tests combined
    print()
    print_after_tss_ranks(combined_results, context_desc='[ALL TESTS COMBINED] ')
