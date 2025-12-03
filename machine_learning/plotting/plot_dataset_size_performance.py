"""
Plot machine learning performance against dataset size metrics.

Reads the summary of model performance from main_results_summary.tsv
and generates scatter plots to visualize the relationship between
performance (AUC) and various dataset size metrics (# samples, 
# positive samples, # datasets).

Points in the plots are colored by the cross-validation method used.
"""

# Imports ------------------------------------------
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# File Paths ------------------------------------------
BASE_DIR = "/home/haig/Repos/meta-analysis/machine_learning"
INPUT_FILE = os.path.join(BASE_DIR, "results", "main_results_summary.tsv")
OUTPUT_DIR = os.path.join(BASE_DIR, "results", "dataset_sizes")

# Option to not plot the XSET AUC
PLOT_XSET_AUC = False

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Main ------------------------------------------
# Read the summary data
df = pd.read_csv(INPUT_FILE, sep='\t')

# Reshape the DataFrame from wide to long format for easier plotting
id_vars = ['Dataset', 'Body Site', 'Stage', '# Datasets', '# Samples', '# Positive Samples', '# Negative Samples']
value_vars = ['KFOLD AUC', 'LODO AUC']
if PLOT_XSET_AUC:
    value_vars.append('Mean XSET AUC')
df_long = pd.melt(df, id_vars=id_vars, value_vars=value_vars, var_name='Cross-Validation Type', value_name='AUC')

# Clean up the 'Cross-Validation Type' names for the legend
df_long['Cross-Validation Type'] = df_long['Cross-Validation Type'].str.replace(' AUC', '').replace({'Mean XSET': 'XSET'})

# Define the x-axes for the plots
plot_definitions = {
    '# Samples': 'Number of Samples',
    '# Positive Samples': 'Number of Positive Samples',
    '# Datasets': 'Number of Datasets'
}

# Generate and save each plot
for x_col, x_label in plot_definitions.items():
    plt.figure(figsize=(10, 8))
    
    sns.scatterplot(
        data=df_long,
        x=x_col,
        y='AUC',
        hue='Cross-Validation Type',
        palette='viridis',
        s=150,  # size of points
        alpha=0.8,
        edgecolor='k'
    )
    
    plt.title(f'Model Performance (AUC) vs. {x_label}', fontsize=16)
    plt.xlabel(x_label, fontsize=12)
    plt.ylabel('AUC Score', fontsize=12)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.legend(title='CV Type', fontsize=10)
    plt.ylim(0.5, 1.0) # Set y-axis limit for better comparison
    
    plt.tight_layout()
    
    # Construct output file path
    sanitized_xlabel = x_label.lower().replace(' ', '_').replace('#', 'n')
    output_filename = f'auc_vs_{sanitized_xlabel}.png'
    output_path = os.path.join(OUTPUT_DIR, output_filename)
    
    # Save the plot
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    
    print(f"Plot saved to {output_path}")
