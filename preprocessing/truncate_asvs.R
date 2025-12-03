
# Script for combination of 16S datasets according to their ASV alignments

# Inputs:
#  - This script runs on the contents of the Celiac Microbiome Repository (https://github.com/CeliacMicrobiomeRepo/celiac-repository/)
#  - Specifically, the contents of the 16S_datasets directory
#  - Each dataset is a subdirectory with the following files required for this script to run:
#    - seqs.fna: FASTA file with all ASV sequences for the dataset
#    - asv_abundances.tsv: tsv file with ASV abundances, with samples as rows and sequences as columns
#     (both of these are output from DADA2 or a similar process)

# Steps:
# 1. Writes a single FASTA file with all ASV sequences combined from all datasets
#   - ASV IDs are updated by adding their dataset ID (e.g. "ASV_0028" -> "16S_24_Giacomin;ASV_0028")
# 2. Run mothur on the FASTA file to align to the reference alignment
#   - An existing alignment of full length 16S genes is used...
#   - Like silva.nr_v138_1.align (from: https://mothur.org/wiki/silva_reference_files/)
# 3. Truncates the alignment so that the overlap is the best
#   - Ensuring that all ASVs overlap with no overhangs the best you can
# 4. For each dataset, replaces the ASVs with truncated versions and export to a subdirectory
#   - ASVs which failed to align well are discarded
#   - Truncated ASVs which are identical are aggregated, summing abundances
#   - ASVs IDs are changed and are still specific to each dataset


# Outputs:
#  A) For each dataset directory processed, a subdirectory is created named OUTPUT_SUBDIR_NAME
#    - This contains the truncated ASVs specific to the dataset
#    - Inside this subdirectory, the following files are created:
#      - asv_abundances.tsv             -  Abundance table of samples (rows) by truncated ASVs (columns)
#      - seqs.fna                       -  FASTA file with truncated ASV sequences, with generic ASV_# IDs
#      - seqs_table.rds                 -  DADA2-style sequence table with truncated ASV sequences as columns
#      - asv_abundances_transposed.tsv  -  Transposed abundance table of ASVs (rows) by samples (columns)
#      - about.txt                      -  Information about the truncation parameters and dataset
#      - taxonomy.tsv                   -  Table of taxonomy for each ASV (if SAVE_PHYLOSEQ_AND_TAXONOMY is TRUE)
#      - phyloseq_obj.rds               -  Phyloseq object with all of the above (if SAVE_PHYLOSEQ_AND_TAXONOMY is TRUE)
#  B) In the current working directory (WORKING_DIR), a subdirectory is created also named OUTPUT_SUBDIR_NAME
#    - This contains summary information about the alignment and truncation across all datasets:
#      - start_freq.csv                         -  Table of start position frequencies from the alignment
#      - end_freq.csv                           -  Table of end position frequencies from the alignment
#      - truncated_asv_lengths_boxplot.png      -  Boxplot of truncated ASV lengths per dataset
#      - truncated_asv_lengths_density_plot.png -  Density plot of truncated ASV lengths per dataset
#      - discarded_asvs.csv                     -  Table of discarded ASVs per dataset
#    - It also contains intermediate and final sequence files:
#      - all_v4_asvs.fasta                     -  FASTA file with all ASV sequences from all datasets combined
#      - all_v4_asvs.align                     -  Alignment of all ASV sequences from all datasets combined
#      - all_v4_asvs.filter.fasta              -  Alignment of all ASV sequences from all datasets combined, filtered for empty columns
#      - all_v4_asvs_trunc.fasta               -  FASTA file with all truncated ASV sequences from all datasets combined


# Requirements:
#  - R 4.5.1
#  - R packages: dada2, phyloseq, DECIPHER, Biostrings, ggplot2
#  - mothur 1.49.0
#  - SILVA database (mothur): silva.nr_v138_1
#  - GTDB databases (DADA2): sbdi-gtdb-sativa.r10rs226.1genome.assignTaxonomy_formatted.fna.gz, sbdi-gtdb-sativa.r10rs226.20genomes.addSpecies_formatted.fna.gz


# SET UP ====================================
# Load necessary packages
library(dada2)
library(phyloseq)
library(DECIPHER)
library(Biostrings)
library(ggplot2)


# OPTIONS ===================================
# Directories of datasets to include
# Assumed to each include the following files:
#   - seqs.fna: FASTA file with all ASV sequences (output from DADA2)
#   - asv_abundances.tsv: TSV file with ASV abundances (samples as rows, sequences as columns)
DATASET_DIRS <- c(
    # V4 datasets ------------
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_27_Fornasaro/",
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_49_Turjeman/",         # Stool Active & Stool Treated
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_102_Bodkhe/",          # Stool Active & Duodenal Active
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_50_Bibbo/",
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_80_Garcia/",           # Duodenal Active
    "/home/haig/Repos/celiac-repository/16S_datasets/16S_68_Girdhar/",          # Stool Prospective

    # V3-V4 datasets ------------
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_136_Nobel/",           # Stool Treated
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_60_Shi/",              # Stool Active
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_96_Quagliariello/",    # Stool Treated
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_227_Oscarsson/",
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_99_Lionetti/",
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_5_Senicar/",
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_27_Federica/",
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_20_Rawson/",           # Stool Treated
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_24_Giacomin/",
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_62_Tian/",
    "/home/haig/Repos/celiac-repository/16S_datasets/16S_1211_Milletich/"      # Stool Prospective
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_119_Salamon/"          # Duodenal Active

    # V4-V6 datasets ------------
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_56_Laffaldano/",

    # non-V4 datasets ------------
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_20_Caminero/",
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_179_Verdu/",
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_25_Francavilla/",
    # "/home/haig/Repos/celiac-repository/16S_datasets/16S_39_Olivares/",
    )

# Subdirectory to output combination into for each dataset
OUTPUT_SUBDIR_NAME <- "v4_truncation_stool_prospective"

# All original ASVs (for mothur alignment)
ALL_ASVS_FASTA_PATH <- paste0(OUTPUT_SUBDIR_NAME, "/all_v4_asvs.fasta")

# All ASVs immediately after truncated (not cleaned up yet)
ALL_ASVS_TRUNC_FASTA_PATH <- paste0(OUTPUT_SUBDIR_NAME, "/all_v4_asvs_trunc.fasta")

# All ASVs after truncation and length filtering, poor alignment filtering and canonicalisation
ALL_ASVS_TRUNC_FILTERED_FASTA_PATH <- paste0(OUTPUT_SUBDIR_NAME, "/all_v4_asvs_trunc_filtered.fasta")

# About the output data
ABOUT_TXT_STR <- "This truncation is of all V4 datasets included in the meta-analysis. It's truncation positions are conservative, keeping maximum number of ASVs, while sacrificing some sequence length."

# Option to skip taxonomic and phyloseq stuff (slow)
SAVE_PHYLOSEQ_AND_TAXONOMY <- TRUE

# Set the working directory (for mothur)
WORKING_DIR <- "/home/haig/Repos/meta-analysis/preprocessing"
setwd(WORKING_DIR)

# Create the subdirectory
dir.create(OUTPUT_SUBDIR_NAME, recursive = TRUE)

# Paths to databases (These taxonomic assignments are actually discarded and replaced in the combine_trunc_asvs_to_phyloseq.R script anyway)
# TAXONOMY_TRAIN_SET <- "/mnt/secondary/16S_databases/silva_nr99_v138.1_train_set.fa.gz"    #  <--- SILVA 138.1 (not used in meta-analysis)
# SPECIES_ASSIGNMENT_SET <- "/mnt/secondary/16S_databases/silva_species_assignment_v138.1.fa.gz"    #  <--- SILVA 138.1 (not used in meta-analysis)
TAXONOMY_TRAIN_SET <- "/mnt/secondary/16S_databases/sbdi-gtdb-sativa.r10rs226.1genome.assignTaxonomy_formatted.fna.gz"    #  <--- GTDB r226 (https://figshare.scilifelab.se/articles/dataset/SBDI_Sativa_curated_16S_GTDB_database/14869077/9)
SPECIES_ASSIGNMENT_SET <- "/mnt/secondary/16S_databases/sbdi-gtdb-sativa.r10rs226.20genomes.addSpecies_formatted.fna.gz"    #  <--- GTDB r226

# Record all console output to a log file
sink(file.path(OUTPUT_SUBDIR_NAME, "truncation_console_output.log"), split = TRUE)


# Combine all ASVs into one FASTA file ===========================

# Open connection to FASTA file in "append" mode
con <- file(ALL_ASVS_FASTA_PATH, open = "a")

# For all datasets
for (dataset_dir in DATASET_DIRS) {
    # Read all ASVs from the dataset
    seqs_path <- file.path(dataset_dir, "seqs.fna")
    seqs <- readLines(seqs_path)
    
    # Extract the directory name
    dataset_name <- basename(dataset_dir)

    # Filter out sequences with non-ACTG characters and log them
    headers <- seqs[c(TRUE, FALSE)]
    sequences <- seqs[c(FALSE, TRUE)]
    valid <- grepl("^[ACGTacgt]+$", sequences)
    if (any(!valid)) {
        cat(paste("INFO: Discarding", sum(!valid), "sequences from", dataset_name, "due to non-ACTG characters.\n"))
    }
    seqs <- as.vector(rbind(headers[valid], sequences[valid]))

    # Give ASVs new IDs (descriptions)
    # (e.g. "ASV_0028" -> "16S_24_Giacomin;ASV_0028")
    seqs <- gsub("^>(ASV_\\d+)", paste0(">", dataset_name, ";\\1"), seqs)
    
    # Write to the output FASTA file
    writeLines(seqs, con)
}

# Close the connection
close(con)






# Run mothur to align ASVs ==========================

# Run mothur (in terminal)
#   mothur

# Run alignment
#   e.g. 
#   align.seqs(fasta=v4_truncation_stool_prospective/all_v4_asvs.fasta, reference=silva.nr_v138_1/silva.nr_v138_1.align)
#   which outputs -> all_v4_asvs.align
alignment_path <- sub("\\.fasta$", ".align", ALL_ASVS_FASTA_PATH)

# Filter empty columns in alignment (optional - reduces RAM usage)
#   e.g. 
#   filter.seqs(fasta=v4_truncation_stool_prospective/all_v4_asvs.align)
#   which outputs -> all_v4_asvs.filter.fasta
alignment_path <- sub("\\.fasta$", ".filter.fasta", ALL_ASVS_FASTA_PATH)

# Read output alignment file
alignment <- readDNAStringSet(alignment_path, format="fasta")
alignment_names <- names(alignment)

# Convert to a matrix
seqs <- as.character(alignment)
seqs_split <- strsplit(seqs, '')
alignment_mat <- do.call(rbind, seqs_split)
rownames(alignment_mat) <- alignment_names

# Remove objects from memory
rm(alignment)
rm(seqs_split)
gc()

# Print before
cat("Dimensions of whole alignment matrix:\n")
print(dim(alignment_mat))





# Decide truncation positions ============================

# For every sequence, get starting and ending positions relative to the alignment
start_positions <- apply(alignment_mat, 1, function(seq) min(which(!seq %in% c(".", "-"))))
end_positions <- apply(alignment_mat, 1, function(seq) max(which(!seq %in% c(".", "-"))))

# Print min, max, and median for starting and ending positions
print("Start Positions:")
print(summary(start_positions))
print("End Positions:")
print(summary(end_positions))

# Print dataframe of frequencies of starting positions
start_freq <- as.data.frame(table(start_positions))
colnames(start_freq) <- c("Position", "Frequency")
print("Frequencies of starting positions:")
print(start_freq)
# 
# Print dataframe of frequencies of ending positions
end_freq <- as.data.frame(table(end_positions))
colnames(end_freq) <- c("Position", "Frequency")
print("Frequencies of ending positions:")
print(end_freq)

# Define truncation positions 
trunc_start <- 321
trunc_end <- 463






# Truncate ASVs ============================
cat("Truncating sequences from position", trunc_start, "to", trunc_end, "\n")

# Truncate the whole alignment from trunc_start to trunc_end
alignment_mat_truncated <- alignment_mat[, trunc_start:trunc_end]

# Reconstruct sequences
truncated_sequences <- apply(alignment_mat_truncated, 1, paste0, collapse='')
names(truncated_sequences) <- rownames(alignment_mat_truncated)

# Remove gaps from truncated sequences for mapping
truncated_sequences_no_gaps <- gsub("[.-]", "", truncated_sequences)

# Write to file (ALL_ASVS_TRUNC_FASTA_PATH)
trunc_seqs_no_gaps_dna <- DNAStringSet(truncated_sequences_no_gaps)
names(trunc_seqs_no_gaps_dna) <- alignment_names
writeXStringSet(trunc_seqs_no_gaps_dna, filepath=ALL_ASVS_TRUNC_FASTA_PATH, format="fasta")

# Create mapping of all original ASV sequences to truncated sequences
sequence_mapping <- truncated_sequences_no_gaps
names(sequence_mapping) <- names(truncated_sequences)  # names are 'dataset_name;ASV_ID'

# Set poorly aligned sequences in the mapping to NA
# (sequences which do not cover this truncation range - i.e. they start after trunc_start or end before trunc_end)
sequence_mapping <- as.list(sequence_mapping)
sequence_mapping[start_positions > trunc_start | end_positions < trunc_end] <- NA

# Clean up all subsequences ----
# (canonicalise sequences to remove residual overhangs)
valid_idx <- which(!sapply(sequence_mapping, is.null) & !is.na(sequence_mapping))
seq_vec   <- unlist(sequence_mapping[valid_idx], use.names = FALSE)
unique_seqs <- unique(seq_vec)
unique_seqs <- unique_seqs[order(nchar(unique_seqs))]
canon_map <- setNames(unique_seqs, unique_seqs)
for (i in seq_along(unique_seqs)) {
  pattern <- unique_seqs[i]
  if (canon_map[[pattern]] != pattern) next
  if (i == length(unique_seqs)) break
  longer  <- unique_seqs[(i + 1):length(unique_seqs)]
  hits    <- longer[vcountPattern(pattern,
                                  DNAStringSet(longer),
                                  fixed = TRUE) > 0]
  if (length(hits))
    canon_map[hits] <- pattern
}
# Flatten any accidental chains
get_root <- function(s, map) { while (map[[s]] != s) s <- map[[s]]; s }
canon_map <- setNames(vapply(names(canon_map),
                             get_root,
                             FUN.VALUE = character(1),
                             map = canon_map),
                      names(canon_map))
# Replace sequences with their canonical forms
sequence_mapping[valid_idx] <- lapply(seq_vec, function(s) canon_map[[s]])
cat("Subsequence clean‑up complete –", length(unique(unlist(sequence_mapping[valid_idx]))),
    "unique canonical sequences remain.\n")
# ----

# Print min, max and median lengths of successfully truncated sequences
print("Truncated sequences:")
sequence_lengths_with_na <- sapply(sequence_mapping, function(x) if(!is.null(x)) nchar(x) else NA)
sequence_lengths <- sequence_lengths_with_na[!is.na(sequence_lengths_with_na)]
cat("Min length:", min(sequence_lengths), "\n")
cat("Median length:", median(sequence_lengths), "\n")
cat("Max length:", max(sequence_lengths), "\n")
cat("Std length:", sd(sequence_lengths), "\n")

# Apply minimum and maximum length filters
min_length <- 90
max_length <- 100
cat("Applying minimum and maximum length filters: ", min_length, " to ", max_length, "\n")
sequences_before_len_filter <- sum(!sapply(sequence_mapping, is.null) & !is.na(sequence_mapping))
cat("Number of sequences before length filtering:", sequences_before_len_filter, "\n")
filter_mask <- sequence_lengths_with_na >= min_length & sequence_lengths_with_na <= max_length
filter_mask[is.na(filter_mask)] <- FALSE
sequence_mapping <- sequence_mapping[filter_mask]
cat("Number of sequences after length filtering:", length(sequence_mapping), "\n")

# Write to file
filtered_seqs_dna <- DNAStringSet(unlist(sequence_mapping))
names(filtered_seqs_dna) <- names(sequence_mapping)
writeXStringSet(filtered_seqs_dna, filepath=ALL_ASVS_TRUNC_FILTERED_FASTA_PATH, format="fasta")


# Export truncated ASVs ========================

# Initialize a list to store truncated lengths for each dataset
truncated_lengths_list <- list()

# Initialize a list to store discard statistics
discard_stats_list <- list()

# For all datasets
for (dataset_dir in DATASET_DIRS) {
    
    # Extract the directory name (last part of the path)
    dataset_name <- basename(dataset_dir)
    cat("\n\n---\nStarting dataset:", dataset_name, "\n")
    
    # Load the ASV abundance table (matrix where column headers are original sequences)
    abundance_path <- file.path(dataset_dir, "asv_abundances.tsv")
    if (!file.exists(abundance_path)) {
        cat("(!!!) Skipping dataset", dataset_name, "because asv_abundances.tsv not found.\n")
        next
    }
    dataset_sequence_table <- read.table(abundance_path, sep="\t", header=TRUE, row.names=1, check.names=FALSE, stringsAsFactors=FALSE)
    dataset_sequence_table <- as.matrix(dataset_sequence_table)
    
    # Load the sequences and get mapping between ASV IDs and sequences
    seqs_fna_path <- file.path(dataset_dir, "seqs.fna")
    asv_seqs <- readDNAStringSet(seqs_fna_path, format="fasta")
    
    # Create mapping from ASV IDs to sequences
    asv_seq_strings <- as.character(asv_seqs)
    names(asv_seq_strings) <- names(asv_seqs)  # ASV IDs
    
    # Create full IDs used in the alignment
    full_ids <- paste0(dataset_name, ";", names(asv_seqs))
    
    # Get truncated sequences corresponding to each ASV
    truncated_seqs_for_dataset <- sequence_mapping[full_ids]
    
    # Store discard statistics and calculate percentage of poorly aligned ASVs
    num_discarded <- sum(sapply(truncated_seqs_for_dataset, function(x) is.null(x) || is.na(x)))
    num_total <- length(truncated_seqs_for_dataset)
    num_kept <- num_total - num_discarded
    null_percentage <- (num_discarded / num_total) * 100
    
    discard_stats_list[[dataset_name]] <- data.frame(
        Num_ASVs_Discarded = num_discarded,
        Num_ASVs_Kept = num_kept,
        Percent_ASVs_Discarded = null_percentage,
        Percent_ASVs_Kept = (num_kept / num_total) * 100
    )
    
    # Calculate and print the percentage of ASVs with NULL headers (poorly aligned)
    cat("Percentage of ASVs with NULL truncated sequences due to poor alignment:", round(null_percentage, 4), "%  (", num_discarded, "/", num_total, ")\n")
    
    # If no sequences worked
    if (sum(is.na(truncated_seqs_for_dataset)) == length(truncated_seqs_for_dataset)) {
        cat("WARNING: 0 sequences were able to be truncated for dataset ", dataset_name, "\n")
        cat("Skipping processing ", dataset_name, "\n")
        next
    }
    
    # Create mapping from original sequences to truncated sequences
    original_to_truncated <- truncated_seqs_for_dataset
    names(original_to_truncated) <- asv_seq_strings  # original sequences
    
    # Update column headers with truncated sequences using the sequence mapping
    # Replace original column names with truncated sequences if they exist in the mapping
    colnames(dataset_sequence_table) <- sapply(colnames(dataset_sequence_table), function(seq) {
        if (!is.null(original_to_truncated[[seq]])) {
            # cat("Original:", seq, "-> Truncated:", original_to_truncated[[seq]], "\n")
            original_to_truncated[[seq]]
        } else {
            #cat("Original:", seq, "-> Truncated: NA (no mapping available)\n")
            NA  # Mark as NA if there's no truncated sequence available
        }
    })
    # Drop columns with NA values (where truncated sequences were not available)
    dataset_sequence_table <- dataset_sequence_table[, !is.na(colnames(dataset_sequence_table))]
    
    # If no valid sequences remain after filtering, skip to the next dataset
    if (ncol(dataset_sequence_table) == 0) {
      cat("WARNING: All sequences for dataset", dataset_name, "were discarded after filtering. Skipping.\n")
      next
    }
    
    # Aggregate by duplicate column names, summing abundances of identical truncated sequences
    truncated_aggregated_sequence_table <- t(rowsum(t(dataset_sequence_table), group=colnames(dataset_sequence_table)))
    
    # Print before duplicate aggregation
    cat("Number of ASVs before aggregation of duplicates after truncation:", ncol(dataset_sequence_table), "\n")
    cat("Number of ASVs after aggregation of duplicates after truncation:", ncol(truncated_aggregated_sequence_table), "\n")
    
    # Get truncated sequences for length calculation
    truncated_sequences_dataset <- colnames(truncated_aggregated_sequence_table)
    truncated_lengths <- nchar(truncated_sequences_dataset)
    
    # Store the truncated lengths in the list
    truncated_lengths_list[[dataset_name]] <- truncated_lengths
    
    # Print to user
    min_length <- min(truncated_lengths)
    median_length <- median(truncated_lengths)
    max_length <- max(truncated_lengths)
    cat("\nMinimum length of the truncated ASVs:", min_length, "\n")
    cat("Median length of the truncated ASVs:", median_length, "\n")
    cat("Maximum length of the truncated ASVs:", max_length, "\n\n")
    
    if (SAVE_PHYLOSEQ_AND_TAXONOMY) {
        # Assign taxonomy and species
        cat("Assigning taxonomy to truncated sequences...\n")
        truncated_sequence_strings <- colnames(truncated_aggregated_sequence_table)
        truncated_sequence_set <- DNAStringSet(truncated_sequence_strings)
        taxonomy <- assignTaxonomy(truncated_sequence_strings, TAXONOMY_TRAIN_SET, multithread=TRUE)
        taxonomy <- addSpecies(taxonomy, SPECIES_ASSIGNMENT_SET)
        
        # Sample names
        sample_names <- rownames(truncated_aggregated_sequence_table)
        
        # Create a metadata data frame (with no metadata)
        sample_metadata <- data.frame(sample_names = sample_names)
        rownames(sample_metadata) <- sample_metadata$sample_names
        sample_metadata <- sample_data(sample_metadata)
        
        # Create phyloseq object
        phyloseq_obj <- phyloseq(
            otu_table(truncated_aggregated_sequence_table, taxa_are_rows=FALSE),
            sample_metadata,
            tax_table(taxonomy)
        )
        
        # Add DNA sequences to phyloseq object
        dna_sequences <- Biostrings::DNAStringSet(taxa_names(phyloseq_obj))
        names(dna_sequences) <- taxa_names(phyloseq_obj)
        phyloseq_obj <- merge_phyloseq(phyloseq_obj, dna_sequences)
    }
    
    
    # Write outputs ---------------------------
    
    cat("Writing outputs for this dataset\n")
    output_dir <- file.path(dataset_dir, OUTPUT_SUBDIR_NAME)
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    # Export ASV table as a tsv file (asv_abundances.tsv into dataset_dir)
    asv_abundances_path <- file.path(output_dir, "asv_abundances.tsv")
    asv_abundances <- as.data.frame(otu_table(phyloseq_obj))
    write.table(asv_abundances, file=asv_abundances_path, sep="\t", quote=FALSE, col.names=NA)
    
    # Export ASV sequences as a fasta file (seqs.fna into dataset_dir)
    asv_seqs_path <- file.path(output_dir, "seqs.fna")
    asv_seqs <- colnames(truncated_aggregated_sequence_table)
    names(asv_seqs) <- paste0(">ASV_", seq_along(asv_seqs))
    asv_seqs_fasta <- paste(names(asv_seqs), asv_seqs, sep="\n")
    write(asv_seqs_fasta, file=asv_seqs_path)
    
    if (SAVE_PHYLOSEQ_AND_TAXONOMY) {
        # Export taxonomy as a tsv file (taxonomy.tsv into dataset_dir)
        taxonomy_path <- file.path(output_dir, "taxonomy.tsv")
        taxonomy_df <- as.data.frame(tax_table(phyloseq_obj))
        write.table(taxonomy_df, file=taxonomy_path, sep="\t", quote=FALSE, col.names=NA)
        
        # Export the phyoseq object as a RDS file (phyloseq_obj.rds into dataset_dir)
        phyloseq_obj_path <- file.path(output_dir, "phyloseq_obj.rds")
        saveRDS(phyloseq_obj, file = phyloseq_obj_path)
    }
    
    # Export truncated_aggregated_sequence_table as a RDS file
    truncated_aggregated_sequence_table_path <- file.path(output_dir, "seqs_table.rds")
    saveRDS(truncated_aggregated_sequence_table, file = truncated_aggregated_sequence_table_path)
    
    # Transpose ASV abundances table and rename sequences with ASV IDs and "ASV_ID" in the first corner cell
    asv_abundances_transposed <- t(asv_abundances)
    asv_ids <- paste0("ASV_", seq_len(nrow(asv_abundances_transposed)))
    rownames(asv_abundances_transposed) <- asv_ids
    asv_abundances_transposed <- as.data.frame(asv_abundances_transposed)
    asv_abundances_transposed <- cbind(ASV_ID = rownames(asv_abundances_transposed), asv_abundances_transposed)
    rownames(asv_abundances_transposed) <- NULL
    asv_abundances_transposed_path <- file.path(output_dir, "asv_abundances_transposed.tsv")
    write.table(asv_abundances_transposed, file=asv_abundances_transposed_path, sep="\t", quote=FALSE, row.names=FALSE)
    
    # Write about text
    about_txt_path <- file.path(output_dir, "about.txt")
    about_txt_content <- paste0(
        ABOUT_TXT_STR, "\n\n",
        "Date of Combination: ", Sys.Date(), "\n",
        "Dataset: ", dataset_name, "\n",
        "Taxonomy Training Set Path: ", TAXONOMY_TRAIN_SET, "\n",
        "Species Assignment Set Path: ", SPECIES_ASSIGNMENT_SET, "\n",
        "Sequences were truncated based on positions from ", trunc_start, " to ", trunc_end, " in the alignment."
    )
    write(about_txt_content, file = about_txt_path)
    
    
    cat(paste("Dataset done:", dataset_dir, "\n"))
}
cat("Done!\n")




# Convert the list to a data frame for ggplot
truncated_lengths_df <- do.call(rbind, lapply(names(truncated_lengths_list), function(name) {
    data.frame(Dataset = name, Truncated_Length = truncated_lengths_list[[name]])
}))


# Box-and-whisker plot for all datasets' truncated lengths
boxplot <- ggplot(truncated_lengths_df, aes(x = Dataset, y = Truncated_Length)) +
    geom_boxplot() +
    labs(title = "Truncated ASV Lengths per Dataset", x = "Dataset", y = "Truncated Length") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(boxplot)


# Prepare a data frame with normalized densities for each dataset
density_data <- do.call(rbind, lapply(names(truncated_lengths_list), function(dataset_name) {
    lengths <- truncated_lengths_list[[dataset_name]]
    density_obj <- density(lengths)  # Calculate density
    
    # Normalize the density by dividing by the maximum value
    data.frame(
        Truncated_Length = density_obj$x,
        Density = density_obj$y / max(density_obj$y),  # Normalize to max 1
        Dataset = dataset_name
    )
}))

# Frequencies of lengths between datasets (density plot)
density_plot <- ggplot(density_data, aes(x = Truncated_Length, y = Density, color = Dataset)) +
    geom_line() +
    labs(title = "Normalized Density of Truncated ASV Lengths Across Datasets",
         x = "Truncated Length",
         y = "Normalized Density") +
    theme_minimal() +
    theme(legend.position = "right")
print(density_plot)


# Make a directory in the working directory called OUTPUT_SUBDIR_NAME and put in the plots and also write the dataframes end_freq and start_freq
# Define the output directory name
output_main_dir <- file.path(getwd(), OUTPUT_SUBDIR_NAME)
if (!dir.exists(output_main_dir)) dir.create(output_main_dir, recursive = TRUE)

# Save start_freq and end_freq data frames as CSV files
write.csv(start_freq, file = file.path(output_main_dir, "start_freq.csv"), row.names = FALSE)
write.csv(end_freq, file = file.path(output_main_dir, "end_freq.csv"), row.names = FALSE)

# Convert discard_stats_list to a data frame and save it as a CSV file
discard_stats_df <- do.call(rbind, discard_stats_list)
rownames(discard_stats_df) <- names(discard_stats_list)
write.csv(discard_stats_df, file = file.path(output_main_dir, "discarded_asvs.csv"), row.names = TRUE)

# Save plots to the output directory
# Save the boxplot as a PNG
boxplot_path <- file.path(output_main_dir, "truncated_asv_lengths_boxplot.png")
ggsave(boxplot_path, plot = boxplot, width = 8, height = 6, dpi = 300, bg = "white")

# Save the density plot as a PNG
density_plot_path <- file.path(output_main_dir, "truncated_asv_lengths_density_plot.png")
ggsave(density_plot_path, plot = density_plot, width = 8, height = 6, dpi = 300, bg = "white")

cat("Data frames and plots have been saved to the output directory:", output_main_dir, "\n")


# Stop logging
sink()



# Optional extra section: Verify truncation ============================
# This code:
# 1. Visualises all of the unique truncated sequences in a web browser
# 2. Verifies that there are no subsequences in the truncated sequences


# 0) Read the sequences ----------
# Read the sequences from the FASTA file
sequences <- readDNAStringSet(ALL_ASVS_TRUNC_FILTERED_FASTA_PATH)

# Remove exact duplicate sequences to get a set of unique sequences
unique_sequences <- unique(sequences)

# Print the number of sequences and empty sequences
cat("Number of sequences:", length(sequences), "\n")
cat("Number of empty sequences:", sum(width(sequences) == 0), "\n")

# Filter out empty sequences that may have been in the input file
unique_sequences <- unique_sequences[width(unique_sequences) > 0]
cat("Number of unique sequences:", length(unique_sequences), "\n")

# 1) Visualise ----------

# Perform multiple sequence alignment on the unique sequences
alignment <- AlignSeqs(unique_sequences)
# Open the alignment in a web browser for visual inspection
BrowseSeqs(alignment)

# 2) Check for subsequences ----------
cat("\nChecking for subsequences...\n")
# Sort unique sequences by length in ascending order to optimize the search
sorted_sequences <- unique_sequences[order(width(unique_sequences))]
n_seqs <- length(sorted_sequences)
found_subsequence <- FALSE
# Loop through each sequence and check if it's a subsequence of any *longer* sequence
for (i in 1:(n_seqs - 1)) {
  # Progress statement
  cat(sprintf("Processing sequence %d of %d...\r", i, n_seqs - 1))
  pattern_seq <- sorted_sequences[[i]]
  # Only compare against sequences that are longer
  # This is the main optimization
  subject_seqs <- sorted_sequences[(i + 1):n_seqs]
  # Use vcountPattern to find occurrences. It's vectorized and efficient.
  match_counts <- vcountPattern(pattern_seq, subject_seqs, fixed = TRUE)
  if (any(match_counts > 0)) {
    match_indices <- which(match_counts > 0)
    for (match_idx in match_indices) {
      subject_seq <- subject_seqs[[match_idx]]
      # We only care about proper subsequences, so lengths must be different
      if (length(pattern_seq) < length(subject_seq)) {
        if (!found_subsequence) {
          cat("\n") # Print a newline before the first finding
        }
        found_subsequence <- TRUE
        cat("--- Found a subsequence ---\n")
        cat("Subsequence (", length(pattern_seq), "bp):\n", sep="")
        print(pattern_seq)
        cat("Is a subsequence of (", length(subject_seq), "bp):\n", sep="")
        print(subject_seq)
        cat("---------------------------\n")
      }
    }
  }
}
cat(rep(" ", 40), "\r")
if (!found_subsequence) {
  cat("No subsequences found.\n")
} else {
  cat("Subsequence check complete.\n")
}

