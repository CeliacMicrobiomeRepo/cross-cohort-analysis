# This script is for checking the read counts of samples in all analysis groups
# after combination of datasets and filtering of ASVs with combine_trunc_asvs_to_phyloseq.R
#
# It reads the ps0 and ps1 phyloseq objects from each subdirectory
# of PHYLOSEQ_DIR_PATH and reports the min, median, and max sample read counts
# for each dataset within each analysis group and phyloseq object stage.
#
# Expected directory structure:
#
# PHYLOSEQ_DIR_PATH/
# │
# ├── analysis_group_1/
# │   ├── ps0.rds
# │   └── ps1.rds
# │
# ├── analysis_group_2/
# │   ├── ps0.rds
# │   └── ps1.rds
# │
# ...
#

# SET UP ====================================
library(phyloseq)

# Path to the directory containing the phyloseq objects
PHYLOSEQ_DIR_PATH <- "/home/haig/Repos/cross-cohort-analysis/preprocessing/phyloseq_objects"

# Get a list of all analysis group directories
analysis_dirs <- list.dirs(PHYLOSEQ_DIR_PATH, full.names = TRUE, recursive = FALSE)




# FUNCTION ====================================

# Function to get and print read count summaries
get_read_count_summary <- function(phyloseq_obj, stage_name, analysis_group_name) {
    if (is.null(phyloseq_obj)) {
        cat(paste("    -", stage_name, "is NULL, skipping.\n"))
        return()
    }

    cat(paste("  - Stage:", stage_name, "\n"))

    # Get sample data and read counts
    sample_data_df <- as(sample_data(phyloseq_obj), "data.frame")
    read_counts <- sample_sums(phyloseq_obj)

    # Combine into a data frame
    read_counts_df <- data.frame(
        Dataset = sample_data_df$Dataset_ID,
        Read_Count = read_counts
    )

    # Summarize by dataset
    summary_df <- aggregate(Read_Count ~ Dataset, data = read_counts_df, FUN = function(x) {
        c(
            min = min(x),
            median = median(x),
            max = max(x)
        )
    })

    # Print summary
    for (i in 1:nrow(summary_df)) {
        dataset_name <- summary_df$Dataset[i]
        counts <- summary_df$Read_Count[i, ]
        cat(sprintf("    - Dataset: %-25s | Min: %-7.0f | Median: %-7.0f | Max: %-7.0f\n",
                    dataset_name, counts['min'], counts['median'], counts['max']))
    }
}




# MAIN ====================================

# Record results in a text file
sink(file.path(PHYLOSEQ_DIR_PATH, "read_counts.txt"))

# Iterate over each analysis group directory
for (dir_path in analysis_dirs) {
    analysis_group_name <- basename(dir_path)
    cat(paste0("\n--- Analysis Group: ", analysis_group_name, " ---\n"))

    # Load phyloseq objects
    ps0_path <- file.path(dir_path, "ps0.rds")
    ps1_path <- file.path(dir_path, "ps1.rds")

    ps0 <- if (file.exists(ps0_path)) readRDS(ps0_path) else NULL
    ps1 <- if (file.exists(ps1_path)) readRDS(ps1_path) else NULL

    # Get and print summaries
    get_read_count_summary(ps0, "ps0 (Unfiltered)", analysis_group_name)
    get_read_count_summary(ps1, "ps1 (ASV Filtered)", analysis_group_name)
}

sink()


