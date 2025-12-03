# Differential Abundance Code (Species)

# Runs differential abundance analysis using ANCOM-BC2 at the species level for 16S/shotgun phyloseq objects

# Requirements:
#  - R 4.5.1
#  - R packages: phyloseq, ANCOMBC, dplyr, tidyr, metafor, tibble




# ── 0) Set Up ────────────────────────────

# Imports
library(phyloseq)
library(ANCOMBC)
library(dplyr)
library(tidyr)
library(metafor)
library(tibble)


# Seed
set.seed(7)



# ── 1) Load Dataset for Analysis ────────────────────────────

# Load the phyloseq object for the selected analysis
# The ps1 object should be used (filtered ASVs/features, but non-transformed counts)
# e.g:
#   prospective_phyloseq_objects
#   stool_active_phyloseq_objects
#   stool_treated_phyloseq_objects
#   duodenum_phyloseq_objects
ps = readRDS("/home/haig/Repos/meta-analysis/preprocessing/phyloseq_objects/stool_active_phyloseq_objects/ps1.rds")

# Use 1a for 16S or 1b for shotgun due to different species column formatting

# ── 2)  16S species glom ────────────────────────────
# Species tax_glom
tax <- as.data.frame(tax_table(ps))

tax$Genus_Species <- ifelse(
    is.na(tax$Species) | tax$Species == "",
    NA_character_,                              # keep as NA if Species missing
    paste(tax$Genus, tax$Species, sep = "_")    # Genus_Species when Species exists
)

tax_table(ps) <- as.matrix(tax)
ps <- tax_glom(ps, taxrank = "Genus_Species")

taxa_names(ps) <- as.character(tax_table(ps)[, "Genus_Species"])


# ── 3)  Shotgun species glom ────────────────────────────
# Species tax_glom
# ps <- tax_glom(ps, taxrank = "Species")

#tax <- as.data.frame(tax_table(ps))
#taxa_names(ps) <- as.character(tax_table(ps)[, "Species"])


#Set factor level
sd <- sample_data(ps)
sd$Diagnosed_Celiac <- factor(sd$Diagnosed_Celiac, levels = c(FALSE, TRUE))
sample_data(ps) <- sd


# ── 4) Create output folder ─────────────────────────────────────────────────────
outdir <- "16S/differential_abundance/stool_active/Species"
dir.create(outdir, showWarnings = FALSE)

# ── 5)  Create ps_list ────────────────────────────
#ps_list should contain each groups datasets


# 1. Coerce your sample_data to a data.frame
sd_df <- as(sample_data(ps), "data.frame")

# 2. Find all Dataset_IDs that have at least one sample
dataset_ids <- names(which(table(sd_df$Dataset_ID) > 0))

# 3. Split into a named list, pruning both samples and zero‐count taxa
ps_list <- lapply(dataset_ids, function(d) {
  # get the sample names for this dataset
  samps <- rownames(sd_df)[sd_df$Dataset_ID == d]
  # prune to only those samples
  ps_sub <- prune_samples(samps, ps)                        
  # drop any taxa that now have zero counts
  prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
})
names(ps_list) <- dataset_ids


# ── 6)  Remove zero-variance Species ────────────────────────────
#Note some species will have zero variance and will need to be removed before ANCOM c("Species","Species")
bad_asvs <- c("")

# Iterate over each element of ps_list, removing the unwanted ASVs
ps_list_no_bad <- lapply(
    ps_list,
    function(ps) {
        keep <- !taxa_names(ps) %in% bad_asvs   # logical vector: TRUE = keep this ASV
        prune_taxa(keep, ps)
    })

ps_list <- ps_list_no_bad


# ── 7) Drop zero‑variance ASVs ────────────────────────────────────────────────
ps_list_clean <- lapply(ps_list, function(ps) filter_taxa(ps, function(x) var(x) > 0, prune = TRUE))

# ── 8) Helper functions ──────────────────────────────────────────────────────
get_slope_cols <- function(df, prefix) grep(paste0("^", prefix), names(df), value = TRUE)
get_pq_col     <- function(df, type = c("p", "q")) {
  type <- match.arg(type)
  pat  <- paste0("^", type, "(_val)?_.*Diagnosed_Celiac")
  grep(pat, names(df), value = TRUE)[1]
}

add_depth_offset <- function(ps, col = "log_libsize") {
  sample_data(ps)[[col]] <- log(sample_sums(ps))
  ps
}

fix_formula_core <- function(ps) {
  sd <- sample_data(ps)
  if ("Sex" %in% colnames(sd) && nlevels(factor(sd$Sex)) > 1)
    "Diagnosed_Celiac + Sex"           # return character string (no leading ~)
  else
    "Diagnosed_Celiac"
}


# ── 9) Per‑dataset ANCOM‑BC² runs ─────────────────
ancom_res_list  <- list()
per_dataset_sig <- list()

for (d in names(ps_list_clean)) {
  ps_d   <- add_depth_offset(ps_list_clean[[d]])
  fix_fm <- fix_formula_core(ps_d)        # formula object

  res <- ancombc2(
    data          = ps_d,
    fix_formula   = fix_fm,
    rand_formula  = NULL,                 # single study ⇒ no random terms
    group         = "Diagnosed_Celiac",
    p_adj_method  = "fdr",
    struc_zero    = TRUE,
    pseudo_sens   = TRUE,
    prv_cut       = 0,
    lib_cut       = 1000,
    alpha         = 0.05,
    n_cl          = 4,
    verbose       = TRUE
  )$res

  # Flag ASVs significant at p & q < .05
  p_col <- get_pq_col(res, "p");  q_col <- get_pq_col(res, "q")
  res$sig_pq <- if (!is.na(p_col) && !is.na(q_col))
                  !is.na(res[[p_col]]) & !is.na(res[[q_col]]) &
                  res[[p_col]] < 0.05 & res[[q_col]] < 0.05
                else FALSE

  write.csv(res, file = file.path(outdir, sprintf("ancombc2_%s.csv", d)), row.names = FALSE)

  ancom_res_list[[d]]  <- res
  per_dataset_sig[[d]] <- tibble(taxon = res$taxon, sig = res$sig_pq)
}

## ── 10) Helpers ──────────────────────────────────────────────────────────────
get_beta_col <- function(df) {
  # first column for the contrast that is *not* an se_
  setdiff(grep("Diagnosed_Celiac", names(df), value = TRUE),
          grep("^se_", names(df), value = TRUE))[1]
}
get_se_col <- function(df) {
  # first SE column for the contrast
  grep("^se_.*Diagnosed_Celiac", names(df), value = TRUE)[1]
}

## ── 11) Assemble long table of per-study effect sizes ─────────────────────
meta_tbl <- lapply(names(ancom_res_list), function(d) {
  res <- ancom_res_list[[d]]

  beta_col <- get_beta_col(res)
  se_col   <- get_se_col(res)

  if (is.na(beta_col) || is.na(se_col))
    stop("Could not locate slope/SE columns in dataset ", d,
         ".  β-col = ", beta_col, ", SE-col = ", se_col)

  tibble(
    Dataset_ID = d,
    taxon      = res$taxon,
    yi         = res[[beta_col]],    # log-fold-change
    vi         = res[[se_col]]^2     # variance = SE²
  )
}) %>% bind_rows()

## ── 12) Random-effects REML meta-analysis per taxon ───────────────────────
meta_results <- meta_tbl %>%
  group_by(taxon) %>%
  filter(!is.na(yi) & !is.na(vi)) %>%            # drop incomplete rows
  group_modify(~{
    if (nrow(.x) < 2) {
      tibble(estimate = .x$yi,
             se       = sqrt(.x$vi),
             zval     = NA_real_,
             pval     = NA_real_,
             ci.lb    = NA_real_,
             ci.ub    = NA_real_,
             I2       = NA_real_,
             tau2     = NA_real_)
    } else {
      mod <- rma.uni(yi, vi, data = .x, method = "REML")
      tibble(estimate = mod$b[1,1],
             se       = mod$se,
             zval     = mod$zval,
             pval     = mod$pval,
             ci.lb    = mod$ci.lb,
             ci.ub    = mod$ci.ub,
             I2       = mod$I2,
             tau2     = mod$tau2)
    }
  }) %>%
  ungroup() %>%
  mutate(qval = p.adjust(pval, method = "fdr"),
         sig  = pval < 0.05 & qval < 0.05)

## ── 13) Write full & significant-only results ─────────────────────────────
write.csv(meta_results,
          file = file.path(outdir, "meta_reml_results.csv"),
          row.names = FALSE)

meta_sig <- filter(meta_results, sig)
write.csv(meta_sig,
          file = file.path(outdir, "meta_reml_sig_only.csv"),
          row.names = FALSE)

message("✓ REML meta-analysis complete:",
        "\n  • Full table    → ", file.path(outdir, "meta_reml_results.csv"),
        "\n  • Sig-only taxa → ", file.path(outdir, "meta_reml_sig_only.csv"))


## ── 14) Pooled-sample ANCOM-BC² (combined analysis)  ───────────────────────────────


# 1) Make sure library-size offset is present
ps_comb <- add_depth_offset(ps)


# 2) Run ANCOM-BC² on the pooled data
comb_res <- ancombc2(
  data         = ps_comb,
  fix_formula  = "Diagnosed_Celiac + Dataset_ID",      # study indicator as fixed effect
  rand_formula = NULL,
  group        = "Diagnosed_Celiac",
  p_adj_method = "fdr",
  struc_zero   = TRUE,
  pseudo_sens  = TRUE,
  prv_cut      = 0,
  lib_cut      = 1000,
  alpha        = 0.05,
  n_cl         = 4,
  verbose      = TRUE
)$res


## ── 15) Tidy & flag significance ────────────────────────────────────────────
lfc_col <- grep("^lfc_.*Diagnosed_Celiac", names(comb_res), value = TRUE)[1]
se_col  <- grep("^se_.*Diagnosed_Celiac",  names(comb_res), value = TRUE)[1]
p_col   <- sub("^lfc_", "p_",  lfc_col)
q_col   <- sub("^lfc_", "q_",  lfc_col)

if (any(is.na(c(lfc_col, se_col, p_col, q_col))))
  stop("Could not locate the expected lfc / se / p / q columns in pooled ANCOM-BC² output.")

lfc_thresh <- 0.50   # ≈ 1.65-fold change (natural-log scale)

comb_tidy <- comb_res %>%
  transmute(
    taxon        = taxon,
    beta_ancom   = .data[[lfc_col]],
    se_ancom     = .data[[se_col]],
    p_ancom      = .data[[p_col]],
    q_ancom      = .data[[q_col]],
    sig_ancom    = p_ancom < 0.05 & q_ancom < 0.05 & abs(beta_ancom) >= lfc_thresh,
    sign_ancom   = sign(beta_ancom)
  )

comb_sig <- filter(comb_tidy, sig_ancom)


## ── 16) Write output ───────────────────────────────────────────────────────
write.csv(comb_tidy,
          file = file.path(outdir, "combined_ancombc2_results.csv"),
          row.names = FALSE)

write.csv(comb_sig,
          file = file.path(outdir, "combined_ancombc2_sig_only.csv"),
          row.names = FALSE)

message("✓ Pooled ANCOM-BC² complete:",
        "\n  • Full table    → ", file.path(outdir, "combined_ancombc2_results.csv"),
        "\n  • Sig-only taxa → ", file.path(outdir, "combined_ancombc2_sig_only.csv"))


##  ── CSV with *all* significant taxa (meta + pooled) + taxonomy ──────────────────────────


# Pick out the columns we want to keep from each table
meta_keep  <- meta_sig  %>%                              # from REML meta-analysis
  transmute(
    taxon,
    estimate_meta = estimate,
    se_meta       = se,
    p_meta        = pval,
    q_meta        = qval
  )

pooled_keep <- comb_sig %>%                              # from pooled ANCOM-BC²
  transmute(
    taxon,
    lfc_pooled = beta_ancom,
    se_pooled  = se_ancom,
    p_pooled   = p_ancom,
    q_pooled   = q_ancom
  )

# Union of significant taxa, keeping all columns from both methods
sig_tbl <- full_join(meta_keep, pooled_keep, by = "taxon")

# Add taxonomy from the original phyloseq object
tax_tbl <- as.data.frame(tax_table(ps))          # matrix → data.frame
tax_tbl$taxon <- rownames(tax_tbl)

sig_tbl <- left_join(sig_tbl, tax_tbl, by = "taxon")

##  ── 17) Add mean relative-abundance columns (overall & by Celiac status) ─────────────────────


## -- Relative-abundance matrix ------------------------------------------
rel_abund <- transform_sample_counts(ps, function(x) x / sum(x))

abund_mat <- as(otu_table(rel_abund), "matrix")
if (!taxa_are_rows(rel_abund))
  abund_mat <- t(abund_mat)

## Grouping
diag_vec     <- sample_data(rel_abund)$Diagnosed_Celiac
diag_logical <- as.logical(diag_vec)

is_true  <- diag_logical == TRUE
is_false <- diag_logical == FALSE

## Means (with na.rm = TRUE just to be safe)
overall_mean <- rowMeans(abund_mat, na.rm = TRUE)

mean_true <- if (any(is_true)) {
  rowMeans(abund_mat[, is_true, drop = FALSE], na.rm = TRUE)
} else {
  rep(NA_real_, nrow(abund_mat))
}

mean_false <- if (any(is_false)) {
  rowMeans(abund_mat[, is_false, drop = FALSE], na.rm = TRUE)
} else {
  rep(NA_real_, nrow(abund_mat))
}

## Prevalence
overall_prev <- rowMeans(abund_mat > 0, na.rm = TRUE)

prev_true <- if (any(is_true)) {
  rowMeans(abund_mat[, is_true, drop = FALSE] > 0, na.rm = TRUE)
} else {
  rep(NA_real_, nrow(abund_mat))
}

prev_false <- if (any(is_false)) {
  rowMeans(abund_mat[, is_false, drop = FALSE] > 0, na.rm = TRUE)
} else {
  rep(NA_real_, nrow(abund_mat))
}

abund_tbl <- tibble(
  taxon                           = rownames(abund_mat),
  mean_rel_abund_overall          = overall_mean,
  mean_rel_abund_celiac_TRUE      = mean_true,
  mean_rel_abund_celiac_FALSE     = mean_false,
  prev_overall                    = overall_prev,
  prev_celiac_TRUE                = prev_true,
  prev_celiac_FALSE               = prev_false
)

##  -- Merge abundances into sig_tbl, safely ----------------------------
sig_tbl <- left_join(
  sig_tbl,
  abund_tbl,
  by     = "taxon",
  suffix = c("", ".abund")   # ← everything from abund_tbl gets ".abund"
)

##  -- (e) Coalesce & remove duplicates ------------------------------------
for (base in c("mean_rel_abund_overall",
               "mean_rel_abund_celiac_TRUE",
               "mean_rel_abund_celiac_FALSE")) {

  dup <- paste0(base, ".abund")

  if (dup %in% names(sig_tbl)) {              # only if the duplicate exists
    sig_tbl[[base]] <- dplyr::coalesce(sig_tbl[[base]], sig_tbl[[dup]])
    sig_tbl[[dup]]  <- NULL                   # drop the suffixed column
  }
}


## -- 18) Write the final CSV -------------------------------------------------
out_csv <- file.path(outdir, "sig_taxa_combined.csv")
write.csv(sig_tbl, file = out_csv, row.names = FALSE)

message("✓ Combined significant-taxa table (with relative abundances) written to ",
        out_csv)

