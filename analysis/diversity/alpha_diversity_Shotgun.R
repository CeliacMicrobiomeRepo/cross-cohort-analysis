
# Alpha Diversity Code (shotgun)

# Takes a phyloseq object with transformed counts and computes four alpha diversity metrics 

# Requirements:
#  - R 4.5.1
#  - R packages: microbiome, meta


# ── 0) Set Up ────────────────────────────
library(microbiome)
library(meta)

# Seed
set.seed(7)

# Primary biological factor
#   ("Diagnosed_Celiac" or "Will_Develop_Celiac")
GROUPING_VAR <- "Diagnosed_Celiac"

# Create output folder
# e.g:
#   stool_prospective
#   stool_active
#   stool_treated
#   duodenum_active
out_dir <- "Shotgun/alpha_diversity_results/stool_active"
dir.create(out_dir, recursive = TRUE)


# Output file name
out_file_name <- "meta_summary_metrics.csv"


# ── 1) Load Dataset for Analysis ────────────────────────────

# Load the phyloseq object for the selected analysis
# The ps0 object should be used (unfiltered MetaPhlAn features)

# ── Load Dataset (MetaPhlAn Relative Abundance) ───────────────────────────
ps <- readRDS("...active_shotgun_ps0_M")

# Set factor level
sd <- sample_data(ps)
sd$Diagnosed_Celiac <- factor(sd$Diagnosed_Celiac, levels = c(FALSE,TRUE))
sample_data(ps) <- sd

# Ensure orientation and convert to proportions (MetaPhlAn can be 0–1 or 0–100)
otu <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) otu <- t(otu)

# Convert percent → proportion if median > 1
medval <- median(otu, na.rm = TRUE)
if (is.finite(medval) && medval > 1) {
  otu <- otu / 100
}
# Replace NA with 0
otu[is.na(otu)] <- 0

# TSS-normalize each sample to sum = 1 (safe if already normalized)
cs <- rowSums(otu, na.rm = TRUE)
nz <- cs > 0
if (any(nz)) {
  otu[nz, ] <- otu[nz, , drop = FALSE] / cs[nz]
}

# Put back into phyloseq (taxa_are_rows = FALSE here, we’ll keep it that way)
otu_table(ps) <- otu_table(otu, taxa_are_rows = FALSE)

# 1. Coerce your sample_data to a data.frame
sd_df <- as(sample_data(ps), "data.frame")

# 2. Datasets that actually have samples
dataset_ids <- names(which(table(sd_df$Dataset_ID) > 0))

# 3. Split list, prune zero-sum taxa (across the subset) so metrics behave
ps_list <- lapply(dataset_ids, function(d) {
  samps  <- rownames(sd_df)[sd_df$Dataset_ID == d]
  ps_sub <- prune_samples(samps, ps)
  # drop taxa with zero total abundance across this dataset
  prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
})
names(ps_list) <- dataset_ids

# Optionally drop hand-curated “bad ASVs”
bad_asvs <- c("")   # fill if needed
if (length(bad_asvs)) {
  ps_list <- lapply(ps_list, function(x) prune_taxa(!(taxa_names(x) %in% bad_asvs), x))
}

# ── 2) Metrics & helpers (RA-friendly) ───────────────────────────────────────
metric_mapping <- list(
  observed       = "observed",
  shannon        = "diversity_shannon",
  dominance      = "dominance_simpson",
  rare_abundance = "rarity_low_abundance"
)

calculate_stats_for_metric <- function(metric, merged_df, mapping, pseudocount = 1e-6) {
  actual_metric <- mapping[[metric]]
  if (is.null(actual_metric) || !(actual_metric %in% colnames(merged_df))) {
    stop("Metric ", metric, " not found in merged data.")
  }
  df <- merged_df
  df$log2_metric <- log2(df[[actual_metric]] + pseudocount)

  # NOTE: For RA we DROP Read_depth from the model
  lm_formula <- as.formula(paste("log2_metric ~", GROUPING_VAR))
  lmcoefs <- summary(lm(lm_formula, data = df))$coefficients

  mean_all <- mean(df$log2_metric, na.rm = TRUE)
  se_all   <- sd(df$log2_metric,   na.rm = TRUE)/sqrt(nrow(df))

  get_stats <- function(flag) {
    v <- df$log2_metric[df[[GROUPING_VAR]] == flag]
    c(N  = length(v),
      M  = if (length(v)) mean(v, na.rm = TRUE) else NA,
      SD = if (length(v)) sd(v,   na.rm = TRUE) else NA)
  }
  s_false <- get_stats("FALSE")
  s_true  <- get_stats("TRUE")

  coef_name <- paste0(GROUPING_VAR, "TRUE")
  if (coef_name %in% rownames(lmcoefs)) {
    est  <- lmcoefs[coef_name, "Estimate"]
    se   <- lmcoefs[coef_name, "Std. Error"]
    pval <- lmcoefs[coef_name, "Pr(>|t|)"]
  } else {
    est <- se <- pval <- NA
  }
  sig <- ifelse(!is.na(pval) && pval < .05, "Significant", "ns")

  res <- data.frame(
    Metric           = metric,
    Mean_Log2_All    = mean_all,
    SE_Log2_All      = se_all,
    N_FALSE          = s_false["N"],
    Mean_Log2_FALSE  = s_false["M"],
    SD_Log2_FALSE    = s_false["SD"],
    N_TRUE           = s_true["N"],
    Mean_Log2_TRUE   = s_true["M"],
    SD_Log2_TRUE     = s_true["SD"],
    Coef_Group       = est,
    SE_Group         = se,
    P_Group          = pval,
    Significance     = sig,
    stringsAsFactors = FALSE
  )
  # keep original naming pattern but w/o “Read_depth”
  names(res)[names(res) == 'Coef_Group'] <- paste0('Coef_', GROUPING_VAR)
  names(res)[names(res) == 'SE_Group']   <- paste0('SE_', GROUPING_VAR)
  names(res)[names(res) == 'P_Group']    <- paste0('P_', GROUPING_VAR)
  res
}

run_alpha_analysis <- function(physeq_obj, pseudocount = 1e-6) {
  # microbiome::alpha works on proportions for these indices
  alpha_div <- microbiome::alpha(
    physeq_obj,
    index = c("observed","shannon","fisher","bulla",
              "simpson","dominance","low_abundance","rare_abundance")
  )
  alpha_div$Sample.ID <- row.names(alpha_div)

  meta_df <- data.frame(sample_data(physeq_obj))
  meta_df$Sample.ID <- row.names(meta_df)

  merged_df <- merge(alpha_div, meta_df, by = "Sample.ID")

  valid_map <- metric_mapping[
    vapply(metric_mapping, `%in%`, logical(1), colnames(merged_df))
  ]
  valid_names <- names(valid_map)
  if (length(valid_names) == 0) stop("No metrics found in merged data!")

  do.call(rbind, lapply(valid_names, function(m) {
    calculate_stats_for_metric(m, merged_df, valid_map, pseudocount)
  }))
}

# ── 3) Run per-dataset & write CSVs ──────────────────────────────────────────
for (id in names(ps_list)) {
  tab <- run_alpha_analysis(ps_list[[id]])
  tab$Study <- id
  write.csv(tab, file.path(out_dir, paste0(id, "_alpha_analysis.csv")), row.names = FALSE)
  message("Wrote ", id, "_alpha_analysis.csv")
}

# ── 4) Meta-analysis across datasets ─────────────────────────────────────────
files  <- list.files(path = out_dir, pattern = "_alpha_analysis\\.csv$", full.names = TRUE)
allres <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))

meta_out <- lapply(unique(allres$Metric), function(metric_name) {
  dat <- subset(allres, Metric == metric_name)
  m <- metacont(
    n.e       = dat$N_TRUE,
    mean.e    = dat$Mean_Log2_TRUE,
    sd.e      = dat$SD_Log2_TRUE,
    n.c       = dat$N_FALSE,
    mean.c    = dat$Mean_Log2_FALSE,
    sd.c      = dat$SD_Log2_FALSE,
    studlab   = dat$Study,
    data      = dat,
    sm        = "SMD",
    method.tau = "REML",
    common    = FALSE,
    random    = TRUE
  )
  list(metric = metric_name, meta = m)
})

# ── 5) Print summaries ───────────────────────────────────────────────────────
for (res in meta_out) {
  cat("\n===== Metric:", res$metric, "=====\n")
  print(summary(res$meta))
}

# ── 6) Export meta-analysis summary ──────────────────────────────────────────
meta_summary_list <- lapply(meta_out, function(x) {
  m <- x$meta
  get1 <- function(obj, name) {
    v <- obj[[name]]
    if (!is.null(v) && length(v) == 1) v else NA_real_
  }
  data.frame(
    Metric      = x$metric,
    TE_random   = get1(m, "TE.random"),
    seTE_random = get1(m, "seTE.random"),
    z_random    = get1(m, "zval.random"),
    p_random    = get1(m, "pval.random"),
    ci_lower    = get1(m, "lower.random"),
    ci_upper    = get1(m, "upper.random"),
    tau2        = if (!is.null(m[["tau2"]]) && length(m[["tau2"]]) == 1) m[["tau2"]] else
                  if (!is.null(m[["tau^2"]]) && length(m[["tau^2"]]) == 1) m[["tau^2"]] else NA_real_,
    I2          = get1(m, "I2"),
    stringsAsFactors = FALSE
  )
})
meta_summary_df <- do.call(rbind, meta_summary_list)
write.csv(meta_summary_df, file.path(out_dir, out_file_name), row.names = FALSE)
