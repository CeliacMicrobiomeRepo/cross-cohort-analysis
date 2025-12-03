
# This script takes phyloseq object (.rds file) for a single analysis group originating from the combine_trunc_asvs_phyloseq.R script
# It extracts the ASV abundance table and sample metadata, and exports them as separate TSV files ready for machine learning


# SET UP ====================================
# Load required library
library(phyloseq)

# Set the working directory
setwd("/home/haig/Repos/meta-analysis/machine_learning/")

# Input file
# [!!!] Change groups each time
# e.g.   duodenum_phyloseq_objects   prospective_phyloseq_objects   stool_active_phyloseq_objects   stool_treated_phyloseq_objects
PHYLOSEQ_OBJECT_PATH <- "/home/haig/Repos/meta-analysis/preprocessing/phyloseq_objects/duodenum_phyloseq_objects/ps0.rds"

# Output files
# [!!!] Change groups each time
# e.g.   duodenum_active   stool_prospective   stool_active   stool_treated
OUTPUT_ABUNDANCES_PATH <- "/home/haig/Repos/meta-analysis/machine_learning/all_data/duodenum_active/unfiltered_asv_table.tsv"
OUTPUT_METADATA_PATH <- "/home/haig/Repos/meta-analysis/machine_learning/all_data/duodenum_active/unfiltered_sample_data.tsv"
OUTPUT_TAXONOMIES_PATH <- "/home/haig/Repos/meta-analysis/machine_learning/all_data/duodenum_active/unfiltered_taxonomies.tsv"

# Create output directories if they don't exist
dir.create(dirname(OUTPUT_ABUNDANCES_PATH), recursive = TRUE)
dir.create(dirname(OUTPUT_METADATA_PATH), recursive = TRUE)
dir.create(dirname(OUTPUT_TAXONOMIES_PATH), recursive = TRUE)


# PROCESSING ====================================

# Read the phyloseq object
ps <- readRDS(PHYLOSEQ_OBJECT_PATH)

# Extract ASV table and convert to matrix
asv_table <- data.frame(Sample_ID = rownames(otu_table(ps)), otu_table(ps))

# Print overview of ASV table (rows=samples, cols=ASVs)
cat("ASV table dimensions:", dim(asv_table)[1], "samples x", dim(asv_table)[2], "ASVs\n")
cat("First 3 samples:\n")
print(asv_table$Sample_ID[1:3])

# Extract sample metadata
sample_data <- data.frame(Sample_ID = rownames(sample_data(ps)), sample_data(ps))

# Print overview of sample metadata
cat("Sample metadata dimensions:", dim(sample_data)[1], "samples x", dim(sample_data)[2], "variables\n")
cat("Variable names in sample metadata:\n")
print(names(sample_data))

# Extract taxonomy table
tax_table <- data.frame(
    ASV = rownames(tax_table(ps)),
    Taxonomy = apply(tax_table(ps), 1, function(x) paste(na.omit(x), collapse="; "))
)

# Print overview of taxonomy table
cat("\nTaxonomy table dimensions:", dim(tax_table)[1], "ASVs x", dim(tax_table)[2], "columns\n")
cat("First 3 taxonomic assignments:\n")
print(head(tax_table, 3))

# Write transposed ASV table to file
write.table(data.frame(asv_table), 
            file = OUTPUT_ABUNDANCES_PATH, 
            sep = "\t", 
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

# Write metadata file
write.table(sample_data, 
            file = OUTPUT_METADATA_PATH, 
            sep = "\t", 
            row.names = FALSE,
            quote = FALSE)

# Write taxonomy file
write.table(tax_table,
            file = OUTPUT_TAXONOMIES_PATH,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)
