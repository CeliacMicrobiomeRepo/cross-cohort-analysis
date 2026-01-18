"""
Automatically plot the feature importance results of the best models in multiple datasets.

- the columns of the input file are: Feature (e.g. ASV sequence), followed by each model name (e.g. rf, lr, svm, mlp)
- the values under each model name are the feature importances for the ASV of that row
- plots a box and whisker plot of feature importances for the top N features across all best models in one dataset
- y-axis runs from 0 to 1
- writes a 3:4 png and an svg file

"""

# Imports ------------------------------------------
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Parameters ------------------------------------------
DOTS = True
ALL_FEATURES = False
TOP_N_FEATURES = 30  # Only if ALL_FEATURES is False
BEST_MODEL_ONLY = True

# Using taxonomic labels for ASVs
USE_TAX_LABELS = True
LABEL_RANKS = [4, 5, 6]              # 0=Domain, 1=Phylum, 2=Class, 3=Order, 4=Family, 5=Genus, 6=Species
RANK_NAMES = ["d", "p", "c", "o", "f", "g", "sp"]

# Inputs ------------------------------------------
ALL_DATASET_DIR_PATHS = [
    "/home/haig/Repos/cross-cohort-analysis/machine_learning/datasets_main/stool_prospective_tss_after",
    "/home/haig/Repos/cross-cohort-analysis/machine_learning/datasets_main/stool_active_tss_after",
    "/home/haig/Repos/cross-cohort-analysis/machine_learning/datasets_main/stool_treated_tss_after",
    "/home/haig/Repos/cross-cohort-analysis/machine_learning/datasets_main/duodenum_active_tss_after"
]
# Subdirectory file names
RESULTS_SUBDIR_FILE_NAME = "kfold_results/feature_importance_best_models.tsv"
PERFORMANCE_RESULTS_SUBDIR_FILE_NAME = "kfold_results/best_models_results.tsv"


# Outputs ------------------------------------------
OUTPUT_DIR_PATH = "/home/haig/Repos/cross-cohort-analysis/machine_learning/results/feature_importance/"



def get_taxonomy_label(taxonomy_string):
    """Given a full taxonomy string (from the "Taxonomy" column of the feature_importance_df), 
    return the label at the specified rank"""
    # Split the taxonomy string into a list of taxonomic levels
    taxonomy_levels = taxonomy_string.split(";")

    # Get the label at each rank
    label = ''
    for rank in LABEL_RANKS:
        # Check if we have the taxonomic resolution at the specified rank
        if len(taxonomy_levels) > rank:
            # Return the taxonomic classification at the specified rank
            label += taxonomy_levels[rank].strip() + " | "
        else:
            # Insufficient taxonomic resolution, "No ID at {RANK}"
            label += "no_ID_rank_" + RANK_NAMES[rank] + " | "
            break
    # Remove the trailing " | "
    label = label.rstrip(" | ")
    return label




# Main ------------------------------------------
for DATASET_DIR_PATH in ALL_DATASET_DIR_PATHS:

    # Inputs ------------------------------------------
    RESULTS_PATH = os.path.join(DATASET_DIR_PATH, RESULTS_SUBDIR_FILE_NAME)
    PERFORMANCE_RESULTS_PATH = os.path.join(DATASET_DIR_PATH, PERFORMANCE_RESULTS_SUBDIR_FILE_NAME)
    DATASET_NAME = os.path.basename(DATASET_DIR_PATH)

    # Outputs ------------------------------------------
    FILE_NAME = ("all" if ALL_FEATURES else "top" + str(TOP_N_FEATURES)) + "_features_" + DATASET_NAME

    PART_A = 'All' if ALL_FEATURES else ('Top ' + str(TOP_N_FEATURES))
    PART_B = "for The Best Model" if BEST_MODEL_ONLY else "Across Best Models"
    TITLE = DATASET_NAME + ": " + PART_A + ' Feature Importances ' + PART_B



    # Plot ------------------------------------------
    # Create output directory if it doesn't exist
    if not os.path.exists(OUTPUT_DIR_PATH):
        os.makedirs(OUTPUT_DIR_PATH)

    # Read the feature importance data
    feature_importance_df = pd.read_csv(RESULTS_PATH, sep='\t')

    # Fix labels of ASVs ---
    # If using taxonomic labels for ASVs
    if USE_TAX_LABELS:
        # Replace ASV column with taxonomic labels at the specified rank
        feature_importance_df['ASV'] = feature_importance_df['Taxonomy'].apply(get_taxonomy_label)
        # If there are duplicates, add a suffix to their ASV labels (e.g. Turicibacter_1, Turicibacter_2, etc.)
        duplicate_mask = feature_importance_df['ASV'].duplicated(keep=False)
        if duplicate_mask.any():
            for name in feature_importance_df.loc[duplicate_mask, 'ASV'].unique():
                # Get indices of duplicates for this name
                dup_indices = feature_importance_df[feature_importance_df['ASV'] == name].index
                # Add numbered suffix to duplicates
                feature_importance_df.loc[dup_indices, 'ASV'] = [
                    f'{name}_({i+1})' for i in range(len(dup_indices))
                ]
    # Else using made up ASV IDs
    else:
        # Replace ASV sequences with IDs
        feature_importance_df['ASV'] = [f'ASV_{i+1}' for i in range(len(feature_importance_df))]

    if BEST_MODEL_ONLY:
        # Identify the best model
        df = pd.read_csv(PERFORMANCE_RESULTS_PATH, sep='\t')
        best_model = df[df['Replicate'] == 'Average'].loc[df[df['Replicate'] == 'Average']['Mean'].idxmax(), 'Model']
        # Filter feature importance data to only include the best model
        feature_importance_df = feature_importance_df[['ASV', best_model]]

    # Normalize importances within each model (if not already normalized)
    feature_importance_df.iloc[:, 1:] = feature_importance_df.iloc[:, 1:].div(feature_importance_df.iloc[:, 1:].sum(axis=0), axis=1)

    # Get the top N features across all models (by getting mean importance across all models)
    top_n_features = int(1e6) if ALL_FEATURES else TOP_N_FEATURES
    top_features = feature_importance_df.iloc[:, 1:].mean(axis=1).nlargest(top_n_features).index

    # Filter the dataframe to include only top features
    feature_importance_df = feature_importance_df.loc[top_features]

    # Melt the dataframe to create a long format suitable for plotting
    feature_importance_df = pd.melt(
        feature_importance_df, 
        id_vars=['ASV'],  # The Feature column stays as is
        var_name='Model',     # The model names become a new column
        value_name='Importance'  # The importance values go here
    )

    # Set the plotting style
    sns.set(style="whitegrid")

    # Initialize the matplotlib figure with 3:4 aspect ratio
    plt.figure(figsize=(12, 12))

    if DOTS:
        # Create a stripplot (dots) instead of boxplot
        sns.stripplot(
            x='Importance',
            y='ASV',
            data=feature_importance_df,
            orient='h',
            size=10,  # Size of the dots
            jitter=False  # No jittering since we want exact positions
        )
    else:
        # Create a boxplot
        sns.boxplot(
            x='Importance',
            y='ASV',
            data=feature_importance_df,
            orient='h'
        )

    # Set the title and labels
    plt.title(TITLE)
    plt.xlabel('Normalized Importance')
    plt.ylabel('Features')

    # Tight layout for better spacing
    plt.tight_layout()

    # Export ------------------------------------------
    # Save the plot as PNG with 3:4 aspect ratio
    png_output_path = os.path.join(OUTPUT_DIR_PATH, FILE_NAME + '.png')
    plt.savefig(png_output_path, dpi=300, bbox_inches='tight', format='png')

    # Save the plot as SVG
    # svg_output_path = os.path.join(OUTPUT_DIR_PATH, FILE_NAME + '.svg')
    # plt.savefig(svg_output_path, dpi=300, bbox_inches='tight', format='svg')

    # Show the plot
    plt.show()

    plt.close()
