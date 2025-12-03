# Phylum Composition Testing (16S)

# Requirements:
#  - R 4.5.1
#  - R packages: phyloseq, microbiome, dplyr, purrr, tibble, broom, metafor



# ── 0) Set Up ────────────────────────────

# Imports
library(phyloseq)
library(microbiome)
library(dplyr)
library(purrr)
library(broom)
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
out_dir <- "16S/phylum_compositions_results/stool_active"
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
ps <- readRDS("/home/haig/Repos/meta-analysis/preprocessing/phyloseq_objects/stool_active_phyloseq_objects/ps1.rds")

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

# ---- Helper: collapse suffixes-----------------
collapse_suffix <- function(x) sub("_.*", "", x)

# ---- Helper: log2 abundance of a phylum ------------------------------------
compute_phylum_log2 <- function(ps, phylum_target, synonyms = NULL, pseudocount = 1e-6) {
    # Relative abundance & collapse to phylum level
    ps_rel   <- microbiome::transform(ps, "compositional")
    ps_phyla <- phyloseq::tax_glom(ps_rel, taxrank = "Phylum", NArm = FALSE)
    otu_mat  <- as.matrix(phyloseq::otu_table(ps_phyla))
    tax_raw  <- as.character(phyloseq::tax_table(ps_phyla)[, "Phylum"])
    tax_norm <- collapse_suffix(tax_raw)
    
    # Indices for target phylum
    idx <- which(tax_norm %in% c(phylum_target, synonyms))
    
    # Orientation‑safe abundance per sample
    if (length(idx) == 0) {
        abund <- rep(0, phyloseq::nsamples(ps_phyla))
    } else if (phyloseq::taxa_are_rows(ps_phyla)) {
        abund <- colSums(otu_mat[idx, , drop = FALSE])
    } else {
        abund <- rowSums(otu_mat[, idx, drop = FALSE])
    }
    
    log2(abund + pseudocount)
}

# ---- Per‑study analysis -----------------------------------------------------
analyze_study <- function(ps, study_id) {
    # Sample metadata
    meta_df <- data.frame(phyloseq::sample_data(ps), check.names = FALSE) %>%
        mutate(Sample.ID  = rownames(.),
               Read_depth = phyloseq::sample_sums(ps),
               Dataset_ID = study_id)
    
    # Add phylum log2 abundances
    # Targets with GTDB synonyms
  meta_df <- meta_df %>%
    mutate(
      # Firmicutes target should include Bacillota + split clades
      log2_Bacillota     = compute_phylum_log2(ps, "Firmicutes",
                                 synonyms = c("Bacillota", "Bacillota_A", "Bacillota_B", "Bacillota_C",
                                              "Firmicutes_A", "Firmicutes_B"),
                                 pseudocount = PSEUDOCOUNT),
      log2_Pseudomonadota = compute_phylum_log2(ps, "Proteobacteria",
                                 synonyms = c("Pseudomonadota"),  # GTDB synonym
                                 pseudocount = PSEUDOCOUNT),
      log2_Bacteroidota   = compute_phylum_log2(ps, "Bacteroidota",
                                 synonyms = c("Bacteroidetes"),   # legacy synonym
                                 pseudocount = PSEUDOCOUNT)
    )

  phylum_cols <- c("log2_Bacillota", "log2_Pseudomonadota", "log2_Bacteroidota")
    
    # Fit LMs per phylum
    lm_tbl <- map_dfr(phylum_cols, function(col) {
        fit <- lm(as.formula(paste(col, "~", GROUPING_VAR, "+ Read_depth")), data = meta_df)
        coef_row <- broom::tidy(fit) %>% filter(grepl(paste0("^", GROUPING_VAR), term))
        if (nrow(coef_row) == 0) {
            return(tibble(metric = col, beta = NA_real_, se = NA_real_, p = NA_real_, study = study_id, n = nrow(meta_df)))
        }
        tibble(metric = col,
               beta   = coef_row$estimate,
               se     = coef_row$std.error,
               p      = coef_row$p.value,
               study  = study_id,
               n      = nrow(meta_df))
    })
    
    # Write per‑study CSV
    write.csv(lm_tbl, paste0(out_dir, "/", study_id, "_phylum_lm.csv"), row.names = FALSE)
    message("✓ Wrote ", study_id, "_phylum_lm.csv")
    
    list(lm_tbl = lm_tbl, sample_df = meta_df)
}





# ── 3) Main ────────────────────────────

# ---- Run per‑study loop -----------------------------------------------------
study_res <- map2(ps_list, names(ps_list), analyze_study)

# Gather per‑study LM results -------------------------------------------------
per_study_lm <- bind_rows(map(study_res, "lm_tbl"))
write.csv(per_study_lm, paste0(out_dir, "/phylum_per_study_summary.csv"), row.names = FALSE)

# ---- Random‑effects meta‑analysis ------------------------------------------
meta_res <- per_study_lm %>%
    group_by(metric) %>%
    group_modify(~ {
        if (nrow(.x) < 2 || anyNA(.x$se))
            return(tibble(beta = NA, se = NA, z = NA, p = NA, ci_lb = NA, ci_ub = NA, k = nrow(.x)))
        m <- metafor::rma.uni(yi = .x$beta, sei = .x$se, method = "REML")
        tibble(beta  = as.numeric(m$b),
               se    = m$se,
               z     = m$zval,
               p     = m$pval,
               ci_lb = m$ci.lb,
               ci_ub = m$ci.ub,
               k     = m$k)
    }) %>%
    ungroup()

write.csv(meta_res, paste0(out_dir, "/phylum_meta_analysis.csv"), row.names = FALSE)

# ---- Combined dataset analysis ---------------------------------------------
combined_df <- bind_rows(map(study_res, "sample_df"))

combined_lm_tbl <- map_dfr(c("log2_Bacillota", "log2_Pseudomonadota", "log2_Bacteroidota"),
                           function(col) {
                               fit <- lm(as.formula(paste(col, "~", GROUPING_VAR, "+ Read_depth + Dataset_ID")),
                                         data = combined_df)
                               coef_row <- broom::tidy(fit) %>% filter(grepl(paste0("^", GROUPING_VAR), term))
                               tibble(metric = col,
                                      beta   = coef_row$estimate,
                                      se     = coef_row$std.error,
                                      p      = coef_row$p.value,
                                      n      = nrow(combined_df))
                           })

write.csv(combined_lm_tbl, paste0(out_dir, "/phylum_combined_lm.csv"), row.names = FALSE)

message("Pipeline complete: per‑study, meta‑analysis, and combined LM CSVs written.")

