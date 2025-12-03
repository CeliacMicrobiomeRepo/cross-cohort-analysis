# Shotgun group processing  (script 2)

# Reads in shotgun dataset phyloseq objects produced using "shotgun_processing.R".
# Both MetaPhlan4 relative abundance objects and deterministic absolute counts.
# M = MetaPhlAn (10% prevalence filter)
# D = Deterministic (15% prevalence filter)

# Requirements:
#  - R 4.5.1
#  - R packages: phyloseq


# ── 0) Set Up ────────────────────────────

# Set working directory
setwd ("...shotgun_objects")

# Imports
library(phyloseq)



# ── 1) Group Prep ────────────────────────────────────────────

# Read in objects
SG_132_Francavilla_D <- readRDS("132_Francavilla_deterministic.rds")
SG_132_Francavilla_M <- readRDS("132_Francavilla_relative_abundance.rds")
SG_95_Costigan_D <- readRDS("95_Costigan_deterministic.rds")
SG_95_Costigan_M <- readRDS("95_Costigan_relative_abundance.rds")
SG_80_Mouzan_D <- readRDS("80_Mouzan_deterministic.rds")
SG_80_Mouzan_M <- readRDS("80_Mouzan_relative_abundance.rds")
SG_118_Leonard_D <- readRDS("118_Leonard_deterministic.rds")
SG_118_Leonard_M <- readRDS("118_Leonard_relative_abundance.rds")

# Combine and remove samples < 1000 reads
SG_All_D <- merge_phyloseq(SG_132_Francavilla_D,SG_80_Mouzan_D, SG_95_Costigan_D)
SG_All_D <- subset_samples(SG_All_D, Sample_Site != "duodenal")
SG_All_D <- prune_samples(sample_sums(SG_All_D) >= 1000, SG_All_D)

# Seperate by Group

SG_active_CeD_D <- subset_samples(SG_All_D, Group %in% c("HC", "ACD"))
SG_active_CeD_D <- subset_samples(SG_active_CeD_D, Dataset_ID != "SG_132_Francavilla")
SG_active_CeD_D <- prune_taxa(taxa_sums(SG_active_CeD_D) > 0, SG_active_CeD_D)

SG_treated_CeD_D <- subset_samples(SG_All_D, Group %in% c("HC", "TCD"))
SG_treated_CeD_D <- subset_samples(SG_treated_CeD_D, Dataset_ID != "SG_80_Mouzan")
SG_treated_CeD_D <- prune_taxa(taxa_sums(SG_treated_CeD_D) > 0, SG_treated_CeD_D)


# ── 2) Prevalence filter by dataset ────────────────────────────────────────────
apply_prevalence_by_dataset <- function(ps,
                                        dataset_col = "Dataset_ID",
                                        min_prevalence = 0.15,
                                        min_abund = 0) {   
    stopifnot(phyloseq::nsamples(ps) > 0)
    sdat <- as(phyloseq::sample_data(ps), "data.frame")
    if (!dataset_col %in% colnames(sdat)) {
        stop(sprintf("Column '%s' not found in sample_data()", dataset_col))
    }
    
    # OTU matrix with samples as rows
    M <- as(phyloseq::otu_table(ps), "matrix")
    if (phyloseq::taxa_are_rows(ps)) M <- t(M)
    
    # Align to samples present in both places
    keep_samps <- intersect(rownames(sdat), rownames(M))
    sdat <- sdat[keep_samps, , drop = FALSE]
    M    <- M[keep_samps, , drop = FALSE]
    
    # Group by dataset (drop NA)
    datasets <- unique(sdat[[dataset_col]])
    datasets <- datasets[!is.na(datasets)]
    
    asv_keep_union <- character(0)
    for (ds in datasets) {
        ds_ids <- rownames(sdat)[sdat[[dataset_col]] == ds]
        if (length(ds_ids) == 0) next
        
        # ceil(10% of samples in this dataset)
        min_n <- ceiling(min_prevalence * length(ds_ids))
        
        # ASVs present in >= min_n samples of this dataset at > min_abund
        present_counts <- colSums(M[ds_ids, , drop = FALSE] > min_abund, na.rm = TRUE)
        keep_asvs_ds   <- names(present_counts)[present_counts >= min_n]
        asv_keep_union <- union(asv_keep_union, keep_asvs_ds)
    }
    
    # Prune taxa to union across datasets
    if (length(asv_keep_union) == 0) {
        message("No ASVs passed the prevalence threshold; returning original object.")
        return(ps)
    }
    phyloseq::prune_taxa(phyloseq::taxa_names(ps) %in% asv_keep_union, ps)
}

# ---- Apply to your two objects ----------------------------------------------
MIN_PREVALENCE_PER_DATASET <- 0.15

SG_active_CeD_ps1_D  <- apply_prevalence_by_dataset(
    SG_active_CeD_D,
    dataset_col   = "Dataset_ID",
    min_prevalence = MIN_PREVALENCE_PER_DATASET,
    min_abund      = 0    # strictly >0 counts as present
)

SG_treated_CeD_ps1_D <- apply_prevalence_by_dataset(
    SG_treated_CeD_D,
    dataset_col   = "Dataset_ID",
    min_prevalence = MIN_PREVALENCE_PER_DATASET,
    min_abund      = 0
)
SG_pro_ps1_D <- apply_prevalence_by_dataset(
    SG_118_Leonard_D,
    dataset_col    = "Dataset_ID",
    min_prevalence = MIN_PREVALENCE_PER_DATASET,
    min_abund      = 0
)

# (Optional) drop zero-sum taxa after filtering
SG_active_CeD_ps1_D  <- prune_taxa(taxa_sums(SG_active_CeD_ps1_D)  > 0, SG_active_CeD_ps1_D)
SG_treated_CeD_ps1_D <- prune_taxa(taxa_sums(SG_treated_CeD_ps1_D) > 0, SG_treated_CeD_ps1_D)
SG_pro_ps1_D <- prune_taxa(taxa_sums(SG_pro_ps1_D) > 0, SG_pro_ps1_D)

saveRDS(SG_active_CeD_ps1_D, file = "SG_active_CeD_ps1_D")
saveRDS(SG_treated_CeD_ps1_D, file = "SG_treated_CeD_ps1_D")
saveRDS(SG_pro_ps1_D, file = "SG_pro_ps1_D")



# ── 2) Group Prep (MetaPhlAn) ────────────────────────────────────────────

# Combine and remove samples < 1000 reads

SG_All_M <- merge_phyloseq(SG_132_Francavilla_M,SG_80_Mouzan_M, SG_95_Costigan_M)

SG_All_M <- subset_samples(SG_All_M, Sample_Site != "duodenal")


# Seperate by Group

SG_active_CeD_M <- subset_samples(SG_All_M, Group %in% c("HC", "ACD"))
SG_active_CeD_M <- subset_samples(SG_active_CeD_M, Dataset_ID != "SG_132_Francavilla")
SG_active_CeD_M <- prune_taxa(taxa_sums(SG_active_CeD_M) > 0, SG_active_CeD_M)

SG_treated_CeD_M <- subset_samples(SG_All_M, Group %in% c("HC", "TCD"))
SG_treated_CeD_M <- subset_samples(SG_treated_CeD_M, Dataset_ID != "SG_80_Mouzan")
SG_treated_CeD_M <- prune_taxa(taxa_sums(SG_treated_CeD_M) > 0, SG_treated_CeD_M)

# Save

saveRDS(SG_active_CeD_M, file = "SG_active_CeD_ps0_M")
saveRDS(SG_treated_CeD_M, file = "SG_treated_CeD_ps0_M")


# ── 3) Prevalence filter (keep taxa passing in ANY dataset) ────────────────────
prevalence_filter_any_dataset <- function(ps,
                                          dataset_col   = "Dataset_ID",
                                          min_prevalence = 0.10, 
                                          min_abund      = 0) {  
    stopifnot(phyloseq::nsamples(ps) > 0)
    
    # Sample data
    sdat <- as(phyloseq::sample_data(ps), "data.frame")
    if (!dataset_col %in% colnames(sdat)) {
        stop(sprintf("Column '%s' not found in sample_data()", dataset_col))
    }
    
    # OTU/ASV matrix with samples as rows
    M <- as(phyloseq::otu_table(ps), "matrix")
    if (phyloseq::taxa_are_rows(ps)) M <- t(M)
    
    # Align samples present in both
    keep_samps <- intersect(rownames(sdat), rownames(M))
    sdat <- sdat[keep_samps, , drop = FALSE]
    M    <- M[keep_samps, , drop = FALSE]
    
    # Datasets (drop NAs)
    ds_vals <- sdat[[dataset_col]]
    ds_levels <- unique(ds_vals[!is.na(ds_vals)])
    
    # Union of features that pass >= min_prevalence in >=1 dataset
    keep_union <- character(0)
    for (ds in ds_levels) {
        idx <- which(ds_vals == ds)
        if (!length(idx)) next
        min_n <- ceiling(min_prevalence * length(idx))
        present_counts <- colSums(M[idx, , drop = FALSE] > min_abund, na.rm = TRUE)
        keep_ds <- names(present_counts)[present_counts >= min_n]
        keep_union <- union(keep_union, keep_ds)
    }
    
    if (!length(keep_union)) {
        message("No taxa passed the prevalence threshold in any dataset; returning original object.")
        return(ps)
    }
    
    phyloseq::prune_taxa(phyloseq::taxa_names(ps) %in% keep_union, ps)
}

# ── 4) Apply to your two relative-abundance objects ────────────────────────────
MIN_PREV  <- 0.10

SG_active_CeD_M_ps1  <- prevalence_filter_any_dataset(
    SG_active_CeD_M, dataset_col = "Dataset_ID",
    min_prevalence = MIN_PREV
)

SG_treated_CeD_M_ps1 <- prevalence_filter_any_dataset(
    SG_treated_CeD_M, dataset_col = "Dataset_ID",
    min_prevalence = MIN_PREV
)


SG_pro_CeD_M_ps1 <- prevalence_filter_any_dataset(
    SG_118_Leonard_M, dataset_col = "Dataset_ID",
    min_prevalence = MIN_PREV
)

saveRDS(SG_active_CeD_M_ps1, file = "SG_active_CeD_ps1_M")
saveRDS(SG_treated_CeD_M_ps1, file = "SG_treated_CeD_ps1_M")
saveRDS(SG_pro_CeD_M_ps1, file = "SG_pro_CeD_ps1_M")
