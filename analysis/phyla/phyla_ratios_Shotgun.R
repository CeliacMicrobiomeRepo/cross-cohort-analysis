# Abundance ratios (Shotgun)

# Computes per‑sample log₂(Bacillota/Bacteroidota) ratios.

# Requirements:
#  - R 4.5.1
#  - R packages: phyloseq, microbiome, dplyr, tibble, metafor, purrr, broom



# ── 0) Set Up ────────────────────────────────────────────────────────────────
library(phyloseq)
library(microbiome)
library(dplyr)
library(tibble)
library(metafor)
library(purrr)
library(broom)

# Seed
set.seed(7)

# Primary biological factor
#   ("Diagnosed_Celiac" or "Will_Develop_Celiac")
GROUPING_VAR <- "Diagnosed_Celiac"  # or "Will_Develop_Celiac"


# Create output folder
# e.g:
#   stool_prospective
#   stool_active
#   stool_treated
#   duodenum_active
out_dir <- "Shotgun/abundance_ratios_results/stool_active"
dir.create(out_dir, recursive = TRUE)

PSEUDOCOUNT <- 1e-6  # avoid log2(0)


# ── 1) Load Dataset for Analysis ────────────────────────────

# Load the phyloseq object for the selected analysis
# The ps1_M object should be used (filtered MetaPhlAn features)

# ── 1) Load dataset (RA) & standardize to proportions ────────────────────────
ps <- readRDS("...active_shotgun_ps1_M")

# Set factor level
sd <- sample_data(ps)
sd$Diagnosed_Celiac <- factor(sd$Diagnosed_Celiac, levels = c(FALSE,TRUE))
sample_data(ps) <- sd

# Normalize OTU table to proportions (if 0–100, convert; replace NA with 0)
otu <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) otu <- t(otu)
otu[is.na(otu)] <- 0
medval <- suppressWarnings(median(otu, na.rm = TRUE))
if (is.finite(medval) && medval > 1) otu <- otu / 100
rs <- rowSums(otu, na.rm = TRUE)
nz <- rs > 0
if (any(nz)) otu[nz, ] <- otu[nz, , drop = FALSE] / rs[nz]
otu_table(ps) <- otu_table(otu, taxa_are_rows = FALSE)

# Split by dataset; drop taxa with zero total abundance in each dataset
sd_df <- as(sample_data(ps), "data.frame")
dataset_ids <- names(which(table(sd_df$Dataset_ID) > 0))
ps_list <- lapply(dataset_ids, function(d) {
  samps  <- rownames(sd_df)[sd_df$Dataset_ID == d]
  ps_sub <- prune_samples(samps, ps)
  prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
})
names(ps_list) <- dataset_ids

# Optional manual removals
bad_asvs <- c("")  # fill if needed
if (length(bad_asvs)) {
  ps_list <- lapply(ps_list, function(x) prune_taxa(!(taxa_names(x) %in% bad_asvs), x))
}


# ── 2) Helpers ───────────────────────────────────────────────────────────────
collapse_suffix <- function(x) sub("_.*", "", x)

# phylum matcher to cover GTDB naming (Firmicutes ≈ Bacillota + split clades)
is_target_phylum <- function(phylum_vec, target, synonyms = NULL) {
  lab <- collapse_suffix(as.character(phylum_vec))
  lab %in% c(target, synonyms)
}

# Compute log2(numerator/denominator) at a given rank on RA
compute_log2_ratio <- function(ps, numerator, denominator, rank,
                               pseudocount = 1e-6,
                               num_synonyms = NULL, denom_synonyms = NULL) {
  # Ensure compositional
  ps_rel  <- microbiome::transform(ps, "compositional")

  # Collapse to rank
  ps_rank <- suppressWarnings(phyloseq::tax_glom(ps_rel, taxrank = rank, NArm = FALSE))
  otu_mat <- as.matrix(phyloseq::otu_table(ps_rank))

  # Rank labels
  tax_vec <- as.character(phyloseq::tax_table(ps_rank)[, rank])
  # Harmonize genus names like Prevotella_9 -> Prevotella
  if (rank == "Genus") tax_vec <- collapse_suffix(tax_vec)

  # Locate taxa (may span multiple rows)
  if (rank == "Phylum") {
    num_idx   <- which(is_target_phylum(tax_vec, numerator,  num_synonyms))
    denom_idx <- which(is_target_phylum(tax_vec, denominator, denom_synonyms))
  } else {
    num_idx   <- which(tax_vec %in% c(numerator,  num_synonyms))
    denom_idx <- which(tax_vec %in% c(denominator, denom_synonyms))
  }

  if (length(num_idx) == 0 || length(denom_idx) == 0) {
    return(rep(NA_real_, phyloseq::nsamples(ps_rank)))
  }

  # Orientation-safe sums
  if (phyloseq::taxa_are_rows(ps_rank)) {
    num_abund   <- colSums(otu_mat[num_idx,   , drop = FALSE])
    denom_abund <- colSums(otu_mat[denom_idx, , drop = FALSE])
  } else {
    num_abund   <- rowSums(otu_mat[, num_idx,   drop = FALSE])
    denom_abund <- rowSums(otu_mat[, denom_idx, drop = FALSE])
  }

  log2((num_abund + pseudocount) / (denom_abund + pseudocount))
}

# Per-study analysis on RA (no Read_depth term)
analyze_study <- function(ps, study_id, outdir = ".") {
  meta_df <- data.frame(sample_data(ps), check.names = FALSE) %>%
    mutate(Sample.ID = rownames(.), Dataset_ID = study_id)

  ratio_df <- tibble(
    Sample.ID = sample_names(ps),
    # Firmicutes/Bacteroidota @ Phylum
    log2_Bacillota_Bacteroidota =
      compute_log2_ratio(
        ps, numerator = "Firmicutes", denominator = "Bacteroidota", rank = "Phylum",
        pseudocount = PSEUDOCOUNT,
        num_synonyms   = c("Bacillota", "Bacillota_A", "Bacillota_B", "Bacillota_C", "Firmicutes_A", "Firmicutes_B"),
        denom_synonyms = c("Bacteroidetes")  # legacy synonym
      ))

  merged_df <- left_join(meta_df, ratio_df, by = "Sample.ID")
  if (!GROUPING_VAR %in% colnames(merged_df)) {
    stop("Column '", GROUPING_VAR, "' missing for study ", study_id)
  }
  merged_df[[GROUPING_VAR]] <- as.factor(merged_df[[GROUPING_VAR]])

  # Per-study summaries
  summ <- ratio_df %>%
    summarise(across(starts_with("log2_"),
                     list(mean = ~mean(., na.rm = TRUE),
                          sd   = ~sd(.,   na.rm = TRUE),
                          n    = ~sum(!is.na(.))),
                     .names = "{.col}_{.fn}")) %>%
    mutate(Study = study_id)

  # Per-study LMs (RA: drop Read_depth)
  lm_df <- lapply(c("log2_Bacillota_Bacteroidota"), function(metric) {
    if (all(is.na(merged_df[[metric]]))) return(NULL)
    mod <- lm(as.formula(paste0(metric, " ~ ", GROUPING_VAR)), data = merged_df)
    coef_tab <- summary(mod)$coefficients
    diag_row <- grep(paste0("^", GROUPING_VAR), rownames(coef_tab))
    if (!length(diag_row)) return(NULL)
    tibble(metric    = metric,
           estimate  = coef_tab[diag_row, 1],
           se        = coef_tab[diag_row, 2],
           t_value   = coef_tab[diag_row, 3],
           p_value   = coef_tab[diag_row, 4],
           n         = sum(!is.na(merged_df[[metric]])),
           Study     = study_id)
  }) %>% bind_rows()

  write.csv(summ, file = file.path(outdir, paste0(study_id, "_core_ratio_summary.csv")), row.names = FALSE)
  write.csv(lm_df, file = file.path(outdir, paste0(study_id, "_core_ratio_lm.csv")),      row.names = FALSE)
  message("✓ Wrote ", study_id, "_core_ratio_summary.csv & _lm.csv")

  invisible(list(summary = summ, lm = lm_df, merged = merged_df))
}


# ── 3) Main ──────────────────────────────────────────────────────────────────
# Per-study loop
results     <- lapply(names(ps_list), function(id) analyze_study(ps_list[[id]], id, outdir = out_dir))
study_stats <- bind_rows(lapply(results, `[[`, "summary"))
lm_stats    <- bind_rows(lapply(results, `[[`, "lm"))

# Random-effects meta-analysis of GROUPING_VAR beta
meta_out <- lapply(c("log2_Bacillota_Bacteroidota"), function(metric) {
  df <- dplyr::filter(lm_stats, metric == !!metric)
  if (nrow(df) < 2 || anyNA(df$se)) {
    warning(metric, " skipped: insufficient studies or missing SE.")
    return(NULL)
  }
  m <- metafor::rma.uni(yi = df$estimate, sei = df$se, method = "REML")
  tibble(metric = metric, k = m$k, tau2 = m$tau2,
         est = as.numeric(m$b), se = m$se, z = m$zval, pval = m$pval,
         ci.lb = m$ci.lb, ci.ub = m$ci.ub)
}) %>% bind_rows()
write.csv(meta_out, file.path(out_dir, "core_ratio_meta_analysis.csv"), row.names = FALSE)
message("Wrote core_ratio_meta_analysis.csv")

# Combined dataset (pooled) with dataset adjustment (no Read_depth)
ps_comb <- do.call(merge_phyloseq, ps_list)
combined_ratio_df <- tibble(
  Sample.ID = sample_names(ps_comb),
  log2_Bacillota_Bacteroidota =
    compute_log2_ratio(
      ps_comb, "Firmicutes", "Bacteroidota", "Phylum", PSEUDOCOUNT,
      num_synonyms   = c("Bacillota", "Bacillota_A", "Bacillota_B", "Bacillota_C", "Firmicutes_A", "Firmicutes_B"),
      denom_synonyms = c("Bacteroidetes")
    ))

combined_meta_df <- bind_rows(lapply(names(ps_list), function(id) {
  x <- ps_list[[id]]
  data.frame(sample_data(x)) %>%
    mutate(Sample.ID = rownames(.), Dataset_ID = id)
}))
combined_df <- left_join(combined_meta_df, combined_ratio_df, by = "Sample.ID")
combined_df[[GROUPING_VAR]] <- as.factor(combined_df[[GROUPING_VAR]])
combined_df$Dataset_ID      <- as.factor(combined_df$Dataset_ID)
write.csv(combined_df, file.path(out_dir, "combined_core_ratio_per_sample.csv"), row.names = FALSE)
message("Wrote combined_core_ratio_per_sample.csv")

combined_lm <- lapply(c("log2_Bacillota_Bacteroidota"), function(metric) {
  if (all(is.na(combined_df[[metric]]))) return(NULL)
  mod <- lm(as.formula(paste0(metric, " ~ ", GROUPING_VAR, " + Dataset_ID")), data = combined_df)
  ct  <- summary(mod)$coefficients
  row <- grep(paste0("^", GROUPING_VAR), rownames(ct))
  tibble(metric = metric,
         estimate = ct[row, 1],
         se       = ct[row, 2],
         t_value  = ct[row, 3],
         p_value  = ct[row, 4])
}) %>% bind_rows()
write.csv(combined_lm, file.path(out_dir, "combined_core_ratio_lm.csv"), row.names = FALSE)
message("Wrote combined_core_ratio_lm.csv")
