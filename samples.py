"""
samples.py

Breaks down the sample/dataset filtering processes seen in the paper, 
displaying details of samples and datasets at each step.

This was not the code used for the cross-cohort analysis, but simulates the process accurately.
"""


# SET UP --------------------------------------------------------------

# Panada !
import pandas as pd


# Print per-dataset breakdowns?
PER_DATASET_BREAKDOWNS = True


# These files were obtained from the Celiac Microbiome Repository version 1.0
# https://github.com/CeliacMicrobiomeRepo/celiac-repository
# They may be updated in the future. Version 1.0 is used here.
DATASET_FILE = 'all_samples.tsv'
LOW_READ_FILE = 'low_read_samples.tsv'


# Functions
def print_sample_summary(filtered_df):
    """
    Prints summary statistics and per-dataset breakdown for the given filtered DataFrame.
    """
    unique_datasets = filtered_df['Original_Dataset_ID'].nunique()
    matching_count = len(filtered_df)
    print(f"- Number of datasets: {unique_datasets}")
    print(f"- Number of samples: {matching_count}")

    prospective_celiac = (filtered_df['Will_Develop_Celiac'] == True).sum()
    diagnosed_celiac = (filtered_df['Diagnosed_Celiac'] == True).sum()
    print(f"- Diagnosed celiac samples: {diagnosed_celiac}")
    print(f"- Prospective CeD samples: {prospective_celiac}")
    print(f"- Celiac samples (either): {diagnosed_celiac + prospective_celiac}")

    if PER_DATASET_BREAKDOWNS:
        dataset_counts = filtered_df['Dataset_ID'].value_counts()
        print(f"\nPer dataset breakdown:")
        for dataset, count in sorted(dataset_counts.items()):
            dataset_samples = filtered_df[filtered_df['Dataset_ID'] == dataset]
            celiac_mask = (dataset_samples['Will_Develop_Celiac'] == True) | \
                        (dataset_samples['Diagnosed_Celiac'] == True)
            celiac_samples = celiac_mask.sum()
            n_treated_celiac = (dataset_samples[celiac_mask]['Gluten_Free_Diet'] == True).sum()
            n_active_celiac = celiac_samples - n_treated_celiac
            healthy_controls = count - celiac_samples
            gap = ' ' * (25 - len(dataset))
            print(f"  {dataset}: {gap} {count} samples [HC={healthy_controls}, ACeD={n_active_celiac}, TCeD={n_treated_celiac}]")




# ------------------------------------------------------------------------------------
# MAIN 
# ------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------
# 0. We start with all_samples.tsv, which contains all data in the Celiac Microbiome Repository.
# Read the TSV file
df = pd.read_csv(DATASET_FILE, sep='\t')

# Add a column for the "Original_Dataset_ID" (altered later on)
df['Original_Dataset_ID'] = df['Dataset_ID']
print("-" * 50)
print(f"\n---\n0. All samples in all_samples.tsv:")

# Print summary
print_sample_summary(df)



# ------------------------------------------------------------------------------------------------
# 1. The meta-analysis is focused on stool/duodenal samples that are 16S that captured V4 or shotgun sequencing.
# Apply the 1st filtering criteria
site_filter = df['Sample_Site'].isin(['stool', 'duodenal'])

# Create filters for 16S V4 and Shotgun
s16_v4_filter = (df['Sequencing_Type'] == '16S') & df['Amplicon_Region'].str.contains('V4', na=False)
sg_filter = df['Sequencing_Type'] == 'SG'

# For section 1.5, find 16S samples excluded for not being V4
excluded_by_amplicon_df = df[site_filter & (df['Sequencing_Type'] == '16S') & ~df['Amplicon_Region'].str.contains('V4', na=False)]

# Combine filters and apply to the DataFrame
combined_filter = site_filter & (s16_v4_filter | sg_filter)
matching_count = combined_filter.sum()
df = df[combined_filter]
print(f"\n---\n1. After filtering for stool/duodenal samples that are 16S that captured V4 or Shotgun:")

# Print summary
print_sample_summary(df)



# ------------------------------------------------------------------------------------------------
# 1.5. Print the samples and datasets that were excluded specifically from the amplicon region filter
dataset_counts = excluded_by_amplicon_df['Dataset_ID'].value_counts()
print(f"\n---\n1.5. Samples and datasets excluded specifically due to the V4 region filter:")
print(f"\n{len(dataset_counts)} Datasets were excluded:")
for dataset, count in sorted(dataset_counts.items()):
    dataset_samples = excluded_by_amplicon_df[excluded_by_amplicon_df['Dataset_ID'] == dataset]
    celiac_mask = (dataset_samples['Will_Develop_Celiac'] == True) | \
                (dataset_samples['Diagnosed_Celiac'] == True)
    celiac_samples = celiac_mask.sum()
    n_treated_celiac = (dataset_samples[celiac_mask]['Gluten_Free_Diet'] == True).sum()
    n_active_celiac = celiac_samples - n_treated_celiac
    healthy_controls = count - celiac_samples
    gap = ' ' * (25 - len(dataset))
    print(f"  {dataset}: {gap} {count} samples [HC={healthy_controls}, ACeD={n_active_celiac}, TCeD={n_treated_celiac}]")



# ------------------------------------------------------------------------------------------------
# 2. Remove samples that have significant factors, have <1000 reads, 
# Apply the 2nd filtering criteria

# Load list of samples that did not meet the read-depth threshold (<1000 non-chimera reads)
low_read_samples_df = pd.read_csv(LOW_READ_FILE, sep='\t')
low_read_ids = set(low_read_samples_df['Sample_ID'])

# Build boolean filters
significant_factor_filter = ~(df['Any_Significant_Factor'].astype(str).str.upper() == 'TRUE')
adequate_reads_filter = ~df['Sample_ID'].isin(low_read_ids)
second_filter = significant_factor_filter & adequate_reads_filter

df = df[second_filter]
print(f"\n---\n2. After removing samples with significant factors or fewer than 1000 reads:")
print_sample_summary(df)



# ------------------------------------------------------------------------------------------------
# 3. Remove samples with more targeted filters: 
#  - Down-sampling of Milletich's HCs to 26
#  - "Illumina" not in Seq_Tech (this excludes 16S_27_Fornasaro)
#  - Duodenal samples from SG_80_Mouzan (low quality after host read removal)
# Apply the 3rd filtering criteria

# Remove duodenal samples from SG_80_Mouzan
mouzan_duodenal_mask = (df['Dataset_ID'] == 'SG_80_Mouzan') & (df['Sample_Site'] == 'duodenal')
df = df[~mouzan_duodenal_mask]

# Filter Milletich's HCs
# (Discard all but a random 26 of the healthy controls in the 16S_1211_Milletich dataset)
# Identify healthy-control samples within the Milletich dataset
milletich_mask = df['Dataset_ID'] == '16S_1211_Milletich'

# Select healthy controls
is_celiac = (df['Will_Develop_Celiac'].astype(str).str.upper() == 'TRUE')
milletich_hc_indices = df[milletich_mask & ~is_celiac].index

# If there are more than the desired number of HCs, randomly select those to keep
HC_KEEP_COUNT = 26
if len(milletich_hc_indices) > HC_KEEP_COUNT:
    keep_indices = df.loc[milletich_hc_indices].sample(
        n=HC_KEEP_COUNT).index
    drop_indices = milletich_hc_indices.difference(keep_indices)
    df = df.drop(index=drop_indices)

# Filter out samples from non-Illumina sequencers
illumina_mask = df['Seq_Tech'].str.contains('Illumina', na=False)
df = df[illumina_mask]

print(f"\n---\n3. After applying targeted filters (Milletich HC down-sampling, Non-Illumina exclusion, Mouzan duodenal filter):")
print_sample_summary(df)



# ------------------------------------------------------------------------------------------------
# 3.5. Split datasets by sample site
# Datasets with samples from multiple sites are split into new, site-specific datasets.
# e.g. "16S_102_Bodkhe" -> "16S_102_Bodkhe_(Stool)" + "16S_102_Bodkhe_(Duodenal)"

# Identify datasets with more than one unique sample site
multi_site_mask = df.groupby('Dataset_ID')['Sample_Site'].transform('nunique') > 1

# For rows in these datasets, append the sample site to the Dataset_ID
df.loc[multi_site_mask, 'Dataset_ID'] = df.loc[multi_site_mask, 'Dataset_ID'] + '_(' + df.loc[multi_site_mask, 'Sample_Site'].str.capitalize() + ')'

print(f"\n---\n3.5. After splitting multi-site datasets:")
print_sample_summary(df)



# ------------------------------------------------------------------------------------------------
# 4. Remove datasets with an imbalance of case-control samples
# Apply the 4th filtering criteria

def dataset_is_balanced(n_treated_celiac, n_active_celiac, n_healthy):
    """Checks if a dataset has an acceptable case-control balance."""

    is_balanced_for_treated = n_healthy >= 7 and n_treated_celiac >= 7
    is_balanced_for_active = n_healthy >= 7 and n_active_celiac >= 7
    return is_balanced_for_treated, is_balanced_for_active

# Determine which datasets have a balanced number of celiac and healthy control samples
balanced_dataset_ids = []

# Group the DataFrame by Dataset_ID and evaluate balance for each dataset
for dataset_id, dataset_df in df.groupby('Dataset_ID'):
    # Identify celiac samples (diagnosed or prospective)
    celiac_mask = (dataset_df['Will_Develop_Celiac'] == True) | (dataset_df['Diagnosed_Celiac'] == True)
    n_celiac = celiac_mask.sum()

    # Split celiac samples into treated (on a gluten-free diet) and active
    n_treated_celiac = (dataset_df[celiac_mask]['Gluten_Free_Diet'] == True).sum()
    n_active_celiac = n_celiac - n_treated_celiac

    # Healthy controls are the remaining samples in the dataset
    n_healthy = len(dataset_df) - n_celiac
    # Assess balance
    is_balanced_for_treated, is_balanced_for_active = dataset_is_balanced(n_treated_celiac, n_active_celiac, n_healthy)
    # Keep dataset if it satisfies the balance criterion in either group
    if is_balanced_for_treated or is_balanced_for_active:
        balanced_dataset_ids.append(dataset_id)
        # If it is not balanced for treated, remove treated samples (if any)
        if not is_balanced_for_treated:
            df = df[~((df['Dataset_ID'] == dataset_id) & (df['Diagnosed_Celiac'] == True) & (df['Gluten_Free_Diet'] == True))]
        # If it is not balanced for active, remove active samples (if any)
        if not is_balanced_for_active:
            df = df[~((df['Dataset_ID'] == dataset_id) & (df['Diagnosed_Celiac'] == True) & (df['Gluten_Free_Diet'] == False))]

# Filter the main DataFrame to include only balanced datasets
df = df[df['Dataset_ID'].isin(balanced_dataset_ids)]

print(f"\n---\n4. After removing datasets with imbalanced celiac vs healthy controls:")
print_sample_summary(df)




# ------------------------------------------------------------------------------------------------
# 5. Split into groups for analysis:
#  - Stool Prospective CeD
#  - Stool Active CeD
#  - Stool Treated CeD
#  - Duodenal Active CeD

def run_and_print_analysis_groups(dataframe, title):
    """
    Takes a dataframe and a title, and prints the analysis group breakdowns.
    """
    print(f"\n---\n{title}")
    print(f"Groups for Analysis:")

    def get_analysis_group_stats(df, case_mask, control_mask):
        """
        Filters a DataFrame to find datasets containing both cases and controls,
        and returns summary statistics and the filtered DataFrame for that group.
        """
        # Find datasets that have both cases and controls
        datasets_with_cases = set(df.loc[case_mask, 'Dataset_ID'].unique())
        datasets_with_controls = set(df.loc[control_mask, 'Dataset_ID'].unique())
        datasets_with_both = list(datasets_with_cases.intersection(datasets_with_controls))
        # Filter for samples belonging to these datasets AND being either a case or a control
        group_df = df[
            df['Dataset_ID'].isin(datasets_with_both) &
            (case_mask | control_mask)
        ]
        # Calculate stats from the final filtered group
        n_cases = case_mask.loc[group_df.index].sum()
        n_controls = control_mask.loc[group_df.index].sum()
        n_total = n_cases + n_controls
        n_datasets = len(datasets_with_both)
        return n_total, n_controls, n_cases, n_datasets, group_df

    def print_analysis_group_breakdown(analysis_df, case_mask, control_mask):
        """Prints a per-dataset breakdown for an analysis group."""
        # Group by dataset and iterate
        for dataset_id, dataset_df in sorted(analysis_df.groupby('Dataset_ID')):
            # The masks are for the larger dataframe, so we use the index
            num_cases = case_mask.loc[dataset_df.index].sum()
            num_controls = control_mask.loc[dataset_df.index].sum()
            total_samples = num_cases + num_controls
            dataset_id_print = dataset_id.replace('_(Stool)', '').replace('_(Duodenal)', '')
            gap = ' ' * (30 - len(dataset_id_print))
            print(f"     - {dataset_id_print}:{gap} {total_samples} samples [HC={num_controls}, CeD={num_cases}]")

    # Identify each body site and control group
    stool_df = dataframe[dataframe['Sample_Site'] == 'stool']
    control_mask_stool_non_prospective = (stool_df['Diagnosed_Celiac'] == False) & (stool_df['Gluten_Free_Diet'] == False)
    control_mask_stool_prospective = stool_df['Will_Develop_Celiac'] == False
    duodenal_df = dataframe[dataframe['Sample_Site'] == 'duodenal']
    control_mask_duodenal = (duodenal_df['Diagnosed_Celiac'] == False) & (duodenal_df['Gluten_Free_Diet'] == False)

    # 1. Stool Prospective CeD
    case_mask_prosp = stool_df['Will_Develop_Celiac'] == True
    n_total, n_hc, n_ced, n_datasets, group_df = get_analysis_group_stats(stool_df, case_mask_prosp, control_mask_stool_prospective)
    if n_total > 0:
        print(f"   Stool Prospective CeD   (N={n_total}, HC={n_hc}, CeD={n_ced})   [Across {n_datasets} Datasets]")
        print_analysis_group_breakdown(group_df, case_mask_prosp, control_mask_stool_prospective)

    # 2. Stool Active CeD
    case_mask_active = (stool_df['Diagnosed_Celiac'] == True) & (stool_df['Gluten_Free_Diet'] == False)
    n_total, n_hc, n_ced, n_datasets, group_df = get_analysis_group_stats(stool_df, case_mask_active, control_mask_stool_non_prospective)
    if n_total > 0:
        print(f"   Stool Active CeD        (N={n_total}, HC={n_hc}, CeD={n_ced})   [Across {n_datasets} Datasets]")
        print_analysis_group_breakdown(group_df, case_mask_active, control_mask_stool_non_prospective)

    # 3. Stool Treated CeD
    case_mask_treated = (stool_df['Diagnosed_Celiac'] == True) & (stool_df['Gluten_Free_Diet'] == True)
    n_total, n_hc, n_ced, n_datasets, group_df = get_analysis_group_stats(stool_df, case_mask_treated, control_mask_stool_non_prospective)
    if n_total > 0:
        print(f"   Stool Treated CeD       (N={n_total}, HC={n_hc}, CeD={n_ced})   [Across {n_datasets} Datasets]")
        print_analysis_group_breakdown(group_df, case_mask_treated, control_mask_stool_non_prospective)

    # 4. Duodenal Active CeD
    case_mask_duodenal_active = (duodenal_df['Diagnosed_Celiac'] == True) & (duodenal_df['Gluten_Free_Diet'] == False)
    n_total, n_hc, n_ced, n_datasets, group_df = get_analysis_group_stats(duodenal_df, case_mask_duodenal_active, control_mask_duodenal)
    if n_total > 0:
        print(f"   Duodenal Active CeD     (N={n_total}, HC={n_hc}, CeD={n_ced})   [Across {n_datasets} Datasets]")
        print_analysis_group_breakdown(group_df, case_mask_duodenal_active, control_mask_duodenal)


# Run the analysis for SG, 16S, and combined datasets
sg_df = df[df['Sequencing_Type'] == 'SG']
s16_df = df[df['Sequencing_Type'] == '16S']
run_and_print_analysis_groups(sg_df, "5a. Analysis Groups (Shotgun Datasets)")
run_and_print_analysis_groups(s16_df, "5b. Analysis Groups (16S Datasets)")
run_and_print_analysis_groups(df, "5c. Analysis Groups (Shotgun + 16S Datasets)")


# ------------------------------------------------------------------------------------------------
print("\n" + "-" * 50 + "\n")

