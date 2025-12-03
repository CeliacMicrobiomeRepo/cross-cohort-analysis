# Phylum Composition Testing (Shotgun)

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

# Primary biological factor ("Diagnosed_Celiac" or "Will_Develop_Celiac")
GROUPING_VAR <- "Diagnosed_Celiac"

# Create output folder
# e.g:
#   stool_prospective
#   stool_active
#   stool_treated
#   duodenum_active
out_dir <- "Shotgun/phylum_compositions_results/stool_active"
dir.create(out_dir, recursive = TRUE)

# Avoid log2(0)
PSEUDOCOUNT <- 1e-6


# ── 1) Load RA Dataset ───────────────────────────────────────────────────────
# ps1 should contain filtered features at relative abundance (MetaPhlAn)
ps <- readRDS("......active_shotgun_ps1_M")

# Set factor level
sd <- sample_data(ps)
sd$Diagnosed_Celiac <- factor(sd$Diagnosed_Celiac, levels = c(FALSE,TRUE))
sample_data(ps) <- sd

# Standardize to proportions (handles 0–100% too)
otu <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) otu <- t(otu)
otu[is.na(otu)] <- 0
medval <- median(otu, na.rm = TRUE)
if (is.finite(medval) && medval > 1) otu <- otu / 100
row_sums <- rowSums(otu, na.rm = TRUE)
nz <- row_sums > 0
if (any(nz)) otu[nz, ] <- otu[nz, , drop = FALSE] / row_sums[nz]
otu_table(ps) <- otu_table(otu, taxa_are_rows = FALSE)

# Sample data
sd_df <- as(sample_data(ps), "data.frame")
dataset_ids <- names(which(table(sd_df$Dataset_ID) > 0))

# Split by dataset; drop taxa with zero total RA within each dataset
ps_list <- lapply(dataset_ids, function(d) {
  samps  <- rownames(sd_df)[sd_df$Dataset_ID == d]
  ps_sub <- prune_samples(samps, ps)
  prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
})
names(ps_list) <- dataset_ids

# Optional hand curation
bad_asvs <- c("")  # fill if needed
if (length(bad_asvs)) {
  ps_list <- lapply(ps_list, function(x) prune_taxa(!(taxa_names(x) %in% bad_asvs), x))
}


# ── 2) Functions ─────────────────────────────────────────────────────────────
collapse_suffix <- function(x) sub("_.*", "", x)

# Map GTDB → target phylum label (so “Firmicutes” captures Bacillota variants)
is_target_phylum <- function(phylum_vec, target, synonyms = NULL) {
  lab <- collapse_suffix(as.character(phylum_vec))
  lab %in% c(target, synonyms)
}

compute_phylum_log2 <- function(ps, phylum_target, synonyms = NULL, pseudocount = 1e-6) {
  # Ensure compositional scale (sums to 1 per sample)
  ps_rel   <- microbiome::transform(ps, "compositional")
  ps_phyla <- phyloseq::tax_glom(ps_rel, taxrank = "Phylum", NArm = FALSE)

  otu_mat <- as.matrix(phyloseq::otu_table(ps_phyla))
  phy_col <- phyloseq::tax_table(ps_phyla)[, "Phylum"]
  phy_lab <- if (!is.null(phy_col)) as.vector(phy_col) else rep(NA_character_, ntaxa(ps_phyla))

  idx <- which(is_target_phylum(phy_lab, phylum_target, synonyms))

  if (length(idx) == 0) {
    abund <- rep(0, phyloseq::nsamples(ps_phyla))
  } else if (phyloseq::taxa_are_rows(ps_phyla)) {
    abund <- colSums(otu_mat[idx, , drop = FALSE])
  } else {
    abund <- rowSums(otu_mat[, idx, drop = FALSE])
  }
  log2(abund + pseudocount)
}

# Per-study analysis on RA (NO Read_depth term)
analyze_study <- function(ps, study_id) {
  meta_df <- data.frame(phyloseq::sample_data(ps), check.names = FALSE) %>%
    mutate(Sample.ID = rownames(.), Dataset_ID = study_id)

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

  # Fit LMs per phylum on RA (drop Read_depth)
  lm_tbl <- map_dfr(phylum_cols, function(col) {
    fit <- lm(as.formula(paste(col, "~", GROUPING_VAR)), data = meta_df)
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

  write.csv(lm_tbl, file.path(out_dir, paste0(study_id, "_phylum_lm.csv")), row.names = FALSE)
  message("✓ Wrote ", study_id, "_phylum_lm.csv")

  list(lm_tbl = lm_tbl, sample_df = meta_df)
}


# ── 3) Main ──────────────────────────────────────────────────────────────────
# Per-study loop
study_res <- map2(ps_list, names(ps_list), analyze_study)

# Collect per-study results
per_study_lm <- bind_rows(map(study_res, "lm_tbl"))
write.csv(per_study_lm, file.path(out_dir, "phylum_per_study_summary.csv"), row.names = FALSE)

# Random-effects meta-analysis
meta_res <- per_study_lm %>%
  group_by(metric) %>%
  group_modify(~ {
    if (nrow(.x) < 2 || anyNA(.x$se)) {
      return(tibble(beta = NA, se = NA, z = NA, p = NA, ci_lb = NA, ci_ub = NA, k = nrow(.x)))
    }
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

write.csv(meta_res, file.path(out_dir, "phylum_meta_analysis.csv"), row.names = FALSE)

# Combined dataset analysis (RA; control for Dataset_ID only)
combined_df <- bind_rows(map(study_res, "sample_df"))
combined_lm_tbl <- map_dfr(
  c("log2_Bacillota", "log2_Pseudomonadota", "log2_Bacteroidota"),
  function(col) {
    fit <- lm(as.formula(paste(col, "~", GROUPING_VAR, "+ Dataset_ID")), data = combined_df)
    coef_row <- broom::tidy(fit) %>% filter(grepl(paste0("^", GROUPING_VAR), term))
    tibble(metric = col,
           beta   = coef_row$estimate,
           se     = coef_row$std.error,
           p      = coef_row$p.value,
           n      = nrow(combined_df))
  }
)
write.csv(combined_lm_tbl, file.path(out_dir, "phylum_combined_lm.csv"), row.names = FALSE)

message("Pipeline complete: per-study, meta-analysis, and combined LM CSVs written (RA-based).")
