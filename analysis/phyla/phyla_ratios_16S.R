# Abundance ratios (16S)

# Computes per‑sample log₂(Bacillota/Bacteroidota) ratios.

# Requirements:
#  - R 4.5.1
#  - R packages: phyloseq, microbiome, dplyr, tibble, metafor



# ── 0) Set Up ────────────────────────────

# Imports
library(phyloseq)
library(microbiome)
library(dplyr)
library(tibble)
library(metafor)

# Seed
set.seed(7)

# Primary biological factor
#   ("Diagnosed_Celiac" or "Will_Develop_Celiac")
GROUPING_VAR    <- "Diagnosed_Celiac"

# Create output folder
# e.g:
#   stool_prospective
#   stool_active
#   stool_treated
#   duodenum_active
out_dir <- "16S/abundance_ratios_results/stool_active"
dir.create(out_dir, recursive = TRUE)

# Avoids log2(0)
PSEUDOCOUNT <- 1e-6  



# ── 1) Load Dataset for Analysis ────────────────────────────

# Load the phyloseq object for the selected analysis
# The ps1 object should be used (filtered ASVs, but non-transformed counts)
# e.g:
#   prospective_phyloseq_objects
#   stool_active_phyloseq_objects
#   stool_treated_phyloseq_objects
#   duodenum_phyloseq_objects
ps <- readRDS("/home/haig/Repos/cross-cohort-analysis/preprocessing/phyloseq_objects/stool_active_phyloseq_objects/ps1.rds")

# Set factor level
sd <- sample_data(ps)
sd$Diagnosed_Celiac <- factor(sd$Diagnosed_Celiac, levels = c(FALSE,TRUE))
sample_data(ps) <- sd

# 1. Coerce your sample_data to a data.frame
sd_df <- as(sample_data(ps), "data.frame")

# 2. Find all Dataset_IDs that have at least one sample
dataset_ids <- names(which(table(sd_df$Dataset_ID) > 0))

# 3. Split into a named list, pruning both samples and zero‐count taxa
ps_list <- lapply(dataset_ids, function(d) {
    # get the sample names for this dataset
    samps <- rownames(sd_df)[sd_df$Dataset_ID == d]
    # prune to only those samples
    ps_sub <- prune_samples(samps,ps)                        
    # drop any taxa that now have zero counts
    prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
})
names(ps_list) <- dataset_ids

# 4. Quick check
str(ps_list)

# Note some ASVs will have zero variance and will need to be removed before ANCOM
bad_asvs <- c("")

# Iterate over each element of ps_list, removing the unwanted ASVs
ps_list_no_bad <- lapply(
    ps_list,
    function(ps) {
        keep <- !taxa_names(ps) %in% bad_asvs   # logical vector: TRUE = keep this ASV
        prune_taxa(keep, ps)
    })
ps_list <- ps_list_no_bad



# ── 2) Functions ────────────────────────────

# ───── compute log₂(numerator / denominator) at a given rank ──────
compute_log2_ratio <- function(ps, numerator, denominator, rank, pseudocount = 1e-6) {
    # 3a) relative abundance & collapse
    ps_rel  <- microbiome::transform(ps, "compositional")
    ps_rank <- microbiome::aggregate_taxa(ps_rel, rank)
    otu_mat <- as.matrix(otu_table(ps_rank))
    
    # 3b) taxonomy vector & genus harmonisation
    tax_vec <- as.character(tax_table(ps_rank)[, rank])
    if (rank == "Genus") tax_vec <- sub("_.*", "", tax_vec)
    
    # 3c) locate taxa (may span multiple rows)
    num_idx   <- which(tax_vec == numerator)
    denom_idx <- which(tax_vec == denominator)
    if (length(num_idx) == 0 || length(denom_idx) == 0)
        return(rep(NA_real_, nsamples(ps_rank)))  # taxon missing ⇒ all NA
    
    # 3d) aggregate & compute ratio
    num_abund   <- colSums(otu_mat[num_idx, , drop = FALSE])
    denom_abund <- colSums(otu_mat[denom_idx, , drop = FALSE])
    ratio <- (num_abund + pseudocount) / (denom_abund + pseudocount)
    log2(ratio)
}

# ───── Per‑study analysis ─────────────────────────────────────────────────
analyze_study <- function(ps, study_id, outdir = ".") {
    # 4a) sample metadata
    meta_df <- data.frame(sample_data(ps)) %>%
        mutate(Sample.ID = rownames(.),
               Read_depth = sample_sums(ps))
    
    # 4b) per‑sample ratios
    ratio_df <- tibble(
        Sample.ID = sample_names(ps),
        log2_Bacillota_Bacteroidota = compute_log2_ratio(ps, "Bacillota", "Bacteroidota", rank = "Phylum", pseudocount = PSEUDOCOUNT)
        
    )
    
    # 4c) merge
    merged_df <- left_join(meta_df, ratio_df, by = "Sample.ID")
    if (!GROUPING_VAR %in% colnames(merged_df))
        stop("Column '", GROUPING_VAR, "' missing for study ", study_id)
    merged_df[[GROUPING_VAR]] <- as.factor(merged_df[[GROUPING_VAR]])
    
    # 4d) summary stats
    summ <- ratio_df %>%
        summarise(across(starts_with("log2_"),
                         list(mean = ~mean(., na.rm = TRUE),
                              sd   = ~sd(.,   na.rm = TRUE),
                              n    = ~sum(!is.na(.))),
                         .names = "{.col}_{.fn}")) %>%
        mutate(Study = study_id)
    
    # 4e) linear models
    lm_df <- lapply(c("log2_Bacillota_Bacteroidota"), function(metric) {
        if (all(is.na(merged_df[[metric]]))) return(NULL)
        mod <- lm(as.formula(paste0(metric, " ~ ", GROUPING_VAR, " + Read_depth")), data = merged_df)
        coef_tab <- summary(mod)$coefficients
        diag_row <- grep(GROUPING_VAR, rownames(coef_tab))
        if (length(diag_row) == 0) return(NULL)
        tibble(metric    = metric,
               estimate  = coef_tab[diag_row, 1],
               se        = coef_tab[diag_row, 2],
               t_value   = coef_tab[diag_row, 3],
               p_value   = coef_tab[diag_row, 4],
               n         = sum(!is.na(merged_df[[metric]])),
               Study     = study_id)
    }) %>% bind_rows()
    
    # 4f) write outputs
    write.csv(summ, file = file.path(outdir, paste0(study_id, "_core_ratio_summary.csv")), row.names = FALSE)
    write.csv(lm_df, file = file.path(outdir, paste0(study_id, "_core_ratio_lm.csv")),      row.names = FALSE)
    message("Wrote ", study_id, "_core_ratio_summary.csv & _lm.csv")
    
    invisible(list(summary = summ, lm = lm_df, merged = merged_df))
}






# ── 3) Main ────────────────────────────


# Per‑study loop 
results    <- lapply(names(ps_list), function(id) analyze_study(ps_list[[id]], id, outdir = out_dir))
study_stats <- bind_rows(lapply(results, `[[`, "summary"))
lm_stats    <- bind_rows(lapply(results, `[[`, "lm"))



# Random‑effects meta‑analysis on Diagnosed_Celiac β 
meta_out <- lapply(c("log2_Bacillota_Bacteroidota"), function(metric) {
    df <- filter(lm_stats, metric == !!metric)
    if (nrow(df) < 2) {
        warning(metric, " skipped: <2 studies with data.")
        return(NULL)
    }
    res <- rma(yi = df$estimate, sei = df$se, method = "REML")
    tibble(metric = metric, k = res$k, tau2 = res$tau2,
           est = as.numeric(res$beta), ci.lb = res$ci.lb, ci.ub = res$ci.ub, pval = res$pval)
}) %>% bind_rows()
write.csv(meta_out, file.path(out_dir, "core_ratio_meta_analysis.csv"), row.names = FALSE)
message("Wrote core_ratio_meta_analysis.csv")



# Combined (pooled) per‑sample table & dataset‑adjusted LM 
ps_comb <- do.call(merge_phyloseq, ps_list)
combined_ratio_df <- tibble(
    Sample.ID = sample_names(ps_comb),
    log2_Bacillota_Bacteroidota = compute_log2_ratio(ps_comb, "Bacillota", "Bacteroidota", rank = "Phylum", pseudocount = PSEUDOCOUNT)
  
)
combined_meta_df <- bind_rows(lapply(names(ps_list), function(id) {
    ps <- ps_list[[id]]
    data.frame(sample_data(ps)) %>%
        mutate(Sample.ID = rownames(.), Read_depth = sample_sums(ps), Dataset_ID = id)
}))
combined_df <- left_join(combined_meta_df, combined_ratio_df, by = "Sample.ID")
combined_df[[GROUPING_VAR]] <- as.factor(combined_df[[GROUPING_VAR]])
combined_df$Dataset_ID <- as.factor(combined_df$Dataset_ID)
write.csv(combined_df, file.path(out_dir, "combined_core_ratio_per_sample.csv"), row.names = FALSE)
message("Wrote combined_core_ratio_per_sample.csv")


combined_lm <- lapply(c("log2_Bacillota_Bacteroidota"), function(metric) {
    if (all(is.na(combined_df[[metric]]))) return(NULL)
    mod <- lm(as.formula(paste0(metric, " ~ ", GROUPING_VAR, " + Read_depth + Dataset_ID")), data = combined_df)
    coef_tab <- summary(mod)$coefficients
    diag_row <- grep(GROUPING_VAR, rownames(coef_tab))
    tibble(metric = metric,
           estimate = coef_tab[diag_row,1],
           se       = coef_tab[diag_row,2],
           t_value  = coef_tab[diag_row,3],
           p_value  = coef_tab[diag_row,4])
}) %>% bind_rows()
write.csv(combined_lm, file.path(out_dir, "combined_core_ratio_lm.csv"), row.names = FALSE)
message("Wrote combined_core_ratio_lm.csv")