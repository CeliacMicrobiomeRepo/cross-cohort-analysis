"""
This script visualises the ASVs for each analysis group after running the `combine_trunc_asvs_to_phyloseq.R` script.

It plots:
 - The distribution of ASV abundances across datasets
 - The number of ASVs before and after filtering
 - The number of shared ASVs between datasets
 - The proportion of shared ASVs between datasets
 - The effect of a global prevalence filter on the mean proportion of total abundance removed
"""


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import re


ROOT_PHYLOSEQ_OBJECTS_DIR = "/home/haig/Repos/meta-analysis/preprocessing/phyloseq_objects"

PHYLOSEQ_OBJECT_DIRS = [
    ROOT_PHYLOSEQ_OBJECTS_DIR + "/stool_active_phyloseq_objects",
    ROOT_PHYLOSEQ_OBJECTS_DIR + "/stool_treated_phyloseq_objects",
    ROOT_PHYLOSEQ_OBJECTS_DIR + "/prospective_phyloseq_objects",
    ROOT_PHYLOSEQ_OBJECTS_DIR + "/duodenum_phyloseq_objects",
]

def plot_abundance_distribution(abundance_df, output_dir):
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=abundance_df, x="Dataset", y="Core_ASV_Abundance")
    plt.xticks(rotation=45, ha='right')
    plt.title("Distribution of Core ASV Abundance Across Datasets")
    plt.ylabel("Percentage Abundance of Core ASVs")
    plt.xlabel("Dataset")
    plt.tight_layout()
    output_path = os.path.join(output_dir, "core_asv_abundance_distribution.png")
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    print(f"Saved abundance distribution plot to {output_path}")

def plot_asv_counts(counts_df, output_dir):
    individual_counts = counts_df[~counts_df['set'].str.contains('&&')].copy()
    individual_counts.set_index('set', inplace=True)
    
    individual_counts.plot(kind='bar', figsize=(12, 7))
    plt.title('ASV Counts Before and After Filtering per Dataset')
    plt.ylabel('Number of ASVs')
    plt.xlabel('Dataset')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    output_path = os.path.join(output_dir, "asv_counts_before_after_filtering.png")
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    print(f"Saved ASV counts plot to {output_path}")

def plot_shared_asv_heatmap(counts_df, output_dir, count_column, title_suffix):
    shared_counts = counts_df[counts_df['set'].str.contains('&&')].copy()
    if shared_counts.empty:
        print("No shared dataset counts to plot.")
        return

    datasets = sorted(list(set(sum([s.split('&&') for s in shared_counts['set']], []))))
    heatmap_data = pd.DataFrame(index=datasets, columns=datasets, dtype=float)

    for _, row in shared_counts.iterrows():
        ds1, ds2 = row['set'].split('&&')
        heatmap_data.loc[ds1, ds2] = row[count_column]
        heatmap_data.loc[ds2, ds1] = row[count_column]
    
    individual_counts = counts_df[~counts_df['set'].str.contains('&&')].copy()
    for _, row in individual_counts.iterrows():
        ds = row['set']
        if ds in heatmap_data.index:
            heatmap_data.loc[ds, ds] = row[count_column]

    plt.figure(figsize=(10, 8))
    sns.heatmap(heatmap_data.astype(float), annot=True, fmt=".0f", cmap="viridis")
    plt.title(f"Shared ASVs Between Datasets ({title_suffix})")
    plt.tight_layout()
    output_path = os.path.join(output_dir, f"shared_asvs_heatmap_{count_column}.png")
    plt.savefig(output_path)
    plt.close()
    print(f"Saved shared ASV heatmap to {output_path}")

def plot_proportional_heatmap(counts_df, output_dir, count_column, title_suffix):
    individual_counts = counts_df[~counts_df['set'].str.contains('&&')].set_index('set')
    shared_counts = counts_df[counts_df['set'].str.contains('&&')].copy()
    
    datasets = sorted(list(individual_counts.index))
    proportions = pd.DataFrame(index=datasets, columns=datasets, dtype=float)

    for row_ds in datasets:
        for col_ds in datasets:
            if row_ds == col_ds:
                proportions.loc[row_ds, col_ds] = 1.0
                continue
            
            total_asvs_in_row = individual_counts.loc[row_ds, count_column]
            if total_asvs_in_row == 0:
                proportions.loc[row_ds, col_ds] = 0.0
                continue

            pair_set1 = f"{row_ds}&&{col_ds}"
            pair_set2 = f"{col_ds}&&{row_ds}"
            
            shared_row = shared_counts[(shared_counts['set'] == pair_set1) | (shared_counts['set'] == pair_set2)]
            
            if not shared_row.empty:
                shared_asvs = shared_row.iloc[0][count_column]
                proportions.loc[row_ds, col_ds] = shared_asvs / total_asvs_in_row
            else:
                proportions.loc[row_ds, col_ds] = 0.0

    plt.figure(figsize=(11, 9))
    sns.heatmap(proportions, annot=True, fmt=".2f", cmap="viridis")
    plt.title(f"Proportion of Shared ASVs ({title_suffix})")
    plt.xlabel("in Dataset")
    plt.ylabel("ASVs from Dataset")
    plt.tight_layout()
    
    filename = f"proportional_shared_asvs_{count_column}.png"
    output_path = os.path.join(output_dir, filename)
    plt.savefig(output_path)
    plt.close()
    print(f"Saved proportional heatmap to {output_path}")


def plot_prevalence_filtering_effect(sample_df, abundance_df, output_dir):
    """
    For each dataset, plots the effect of prevalence filtering on the mean proportion of total abundance removed.
    """
    dataset_ids = sample_df['Dataset_ID'].unique()
    
    abundance_df = abundance_df.set_index('Sample.ID')

    for dataset_id in dataset_ids:
        dataset_sample_ids = sample_df[sample_df['Dataset_ID'] == dataset_id]['Sample.ID']
        
        dataset_abundance_df = abundance_df.loc[dataset_sample_ids]
        
        n_samples = len(dataset_sample_ids)
        if n_samples == 0:
            print(f"Skipping {dataset_id} as it has no samples.")
            continue
            
        print(f"Processing {dataset_id} with {n_samples} samples...")

        prevalence_thresholds_percent = []
        mean_filtered_proportions = []

        presence_df = (dataset_abundance_df > 0).astype(int)
        asv_prevalence = presence_df.sum(axis=0)

        for x in range(n_samples + 1):
            asvs_to_filter = asv_prevalence[asv_prevalence < x].index
            
            if len(asvs_to_filter) > 0:
                filtered_abundance_per_sample = dataset_abundance_df[asvs_to_filter].sum(axis=1)
                mean_proportion_filtered = filtered_abundance_per_sample.mean()
            else:
                mean_proportion_filtered = 0.0

            prevalence_thresholds_percent.append((x / n_samples) * 100)
            mean_filtered_proportions.append(mean_proportion_filtered)

        plt.figure(figsize=(10, 6))
        plt.plot(prevalence_thresholds_percent, mean_filtered_proportions, marker='.')
        plt.title(f'Proportion of Abundance Filtered by Prevalence for {dataset_id}')
        plt.xlabel('Prevalence Filter (% of samples ASV must be present in)')
        plt.ylabel('Mean Proportion of Abundance Filtered Out')
        plt.grid(True)
        plt.tight_layout()
        
        output_filename = f"proportion_filtered_{dataset_id}.png"
        output_path = os.path.join(output_dir, output_filename)
        plt.savefig(output_path)
        plt.close()
        print(f"Saved prevalence filtering plot for {dataset_id} to {output_path}")


def plot_prevalence_filtering_effect_all_datasets(sample_df, abundance_df, output_dir):
    """
    Plots the effect of a global prevalence filter on the mean proportion of total abundance removed for each dataset.
    """
    dataset_ids = sample_df['Dataset_ID'].unique()
    abundance_df = abundance_df.set_index('Sample.ID')
    
    n_total_samples = len(sample_df)
    if n_total_samples == 0:
        print("No samples found. Skipping combined plot.")
        return
        
    print(f"Processing combined plot with {n_total_samples} total samples across {len(dataset_ids)} datasets...")

    presence_df = (abundance_df > 0).astype(int)
    asv_prevalence_all = presence_df.sum(axis=0)

    prevalence_thresholds_percent = []
    results_per_dataset = {ds_id: [] for ds_id in dataset_ids}

    for x in range(n_total_samples + 1):
        prevalence_thresholds_percent.append((x / n_total_samples) * 100)
        asvs_to_filter = asv_prevalence_all[asv_prevalence_all < x].index

        if len(asvs_to_filter) == 0:
            for ds_id in dataset_ids:
                results_per_dataset[ds_id].append(0.0)
            continue

        for ds_id in dataset_ids:
            dataset_sample_ids = sample_df[sample_df['Dataset_ID'] == ds_id]['Sample.ID']
            if len(dataset_sample_ids) == 0:
                results_per_dataset[ds_id].append(0.0)
                continue
            
            dataset_abundance_df = abundance_df.loc[dataset_sample_ids]
            
            filtered_abundance_per_sample = dataset_abundance_df[asvs_to_filter].sum(axis=1)
            mean_proportion_filtered = filtered_abundance_per_sample.mean()
            results_per_dataset[ds_id].append(mean_proportion_filtered)

    plt.figure(figsize=(12, 8))
    for ds_id, mean_proportions in results_per_dataset.items():
        plt.plot(prevalence_thresholds_percent, mean_proportions, marker='.', linestyle='-', label=ds_id)
        
    plt.title('Proportion of Abundance Filtered by Global Prevalence')
    plt.xlabel('Global Prevalence Filter (% of all samples ASV must be present in)')
    plt.ylabel('Mean Proportion of Abundance Filtered Out')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    output_filename = "proportion_filtered_all_datasets.png"
    output_path = os.path.join(output_dir, output_filename)
    plt.savefig(output_path)
    plt.close()
    print(f"Saved combined prevalence filtering plot to {output_path}")


def generate_asv_counts_summary(phyloseq_dirs, output_dir):
    asv_counts = {}
    
    dir_to_group_name = {
        "stool_active_phyloseq_objects": "Stool Active",
        "stool_treated_phyloseq_objects": "Stool Treated",
        "prospective_phyloseq_objects": "Stool Prospective",
        "duodenum_phyloseq_objects": "Duodenal Active"
    }
    
    for dir_path in phyloseq_dirs:
        log_file = os.path.join(dir_path, "combination_console_output.log")
        if not os.path.exists(log_file):
            print(f"Log file not found: {log_file}")
            continue

        dir_name = os.path.basename(dir_path)
        group_name = dir_to_group_name.get(dir_name)
        if not group_name:
            print(f"Skipping directory with no group name mapping: {dir_name}")
            continue

        before_count = None
        after_count = None

        with open(log_file, 'r') as f:
            for line in f:
                if "Total ASVs before filtering:" in line:
                    match = re.search(r'(\d+)', line)
                    if match:
                        before_count = int(match.group(1))
                elif "ASVs to keep based on prevalence:" in line:
                    match = re.search(r'(\d+)', line)
                    if match:
                        after_count = int(match.group(1))
        
        if before_count is not None and after_count is not None:
            asv_counts[group_name] = {'before': before_count, 'after': after_count}
        else:
            print(f"Could not find ASV counts in {log_file}")

    output_file_path = os.path.join(output_dir, "asv_counts.txt")
    with open(output_file_path, 'w') as f:
        f.write("Number of ASVs Per Group:\n\n")
        
        group_order = ["Stool Prospective", "Stool Active", "Stool Treated", "Duodenal Active"]
        for group_name in group_order:
            if group_name in asv_counts:
                counts = asv_counts[group_name]
                f.write(f"{group_name}:\n")
                f.write(f" - # ASVs before filtering = {counts['before']}\n")
                f.write(f" - # ASVs after filtering = {counts['after']}\n\n")

    print(f"ASV counts summary saved to {output_file_path}")


if __name__ == "__main__":
    generate_asv_counts_summary(PHYLOSEQ_OBJECT_DIRS, ROOT_PHYLOSEQ_OBJECTS_DIR)

    for PHYLOSEQ_OBJECT_DIR in PHYLOSEQ_OBJECT_DIRS:

        abundance_file = os.path.join(PHYLOSEQ_OBJECT_DIR, "core_asv_sample_abundance.csv")
        counts_file = os.path.join(PHYLOSEQ_OBJECT_DIR, "core_asv_dataset_counts.csv")
        sample_data_file = os.path.join(PHYLOSEQ_OBJECT_DIR, "ps2_sample_data.tsv")
        asv_table_file = os.path.join(PHYLOSEQ_OBJECT_DIR, "ps2_asv_table.tsv")

        plots_dir = os.path.join(PHYLOSEQ_OBJECT_DIR, "plots")
        os.makedirs(plots_dir, exist_ok=True)

        if os.path.exists(abundance_file) and os.path.exists(counts_file):
            abundance_df = pd.read_csv(abundance_file)
            counts_df = pd.read_csv(counts_file)

            plot_abundance_distribution(abundance_df, plots_dir)
            plot_asv_counts(counts_df, plots_dir)
            plot_shared_asv_heatmap(counts_df, plots_dir, 'num_asvs_before_filter', 'Before Filtering')
            plot_shared_asv_heatmap(counts_df, plots_dir, 'num_asvs_after_filter', 'After Filtering')
            plot_proportional_heatmap(counts_df, plots_dir, 'num_asvs_before_filter', 'Before Filtering')
            plot_proportional_heatmap(counts_df, plots_dir, 'num_asvs_after_filter', 'After Filtering')
            print("\nFinished generating original plots.")
        else:
            print("Skipping original plots as one or more input files were not found.")

        if os.path.exists(sample_data_file) and os.path.exists(asv_table_file):
            sample_df = pd.read_csv(sample_data_file, sep='\t')
            abundance_df = pd.read_csv(asv_table_file, sep='\t')
            
            plot_prevalence_filtering_effect(sample_df, abundance_df, plots_dir)
            plot_prevalence_filtering_effect_all_datasets(sample_df, abundance_df, plots_dir)
            print("\nPrevalence filtering plots generated successfully.")
        else:
            print(f"Skipping prevalence filtering plots due to missing files:")
            if not os.path.exists(sample_data_file):
                print(f"  - Not found: {sample_data_file}")
            if not os.path.exists(asv_table_file):
                print(f"  - Not found: {asv_table_file}")
