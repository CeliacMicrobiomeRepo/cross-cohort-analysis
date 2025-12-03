# Constrained Beta-Diversity Analysis (CAP) (shotgun)

# Takes a phyloseq object, performs a TSS transformation, and then runs a 
# constrained analysis of principal coordinates (CAP) to identify how much
# of the variation in community composition can be explained by a given variable,
# while conditioning on a confounding variable.

# Requirements:
#  - R 4.5.1
#  - R packages: phyloseq, vegan, ggplot2, dplyr


# ── 0) Set Up ────────────────────────────────────────────────────────────────
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)

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
out_dir <- "Shotgun/constrained_beta_diversity_results/stool_active"
dir.create(out_dir, recursive = TRUE)

out_sp_scores_file_name     <- "cap_loading_data.csv"
out_sample_scores_file_name <- "cap_sample_data.csv"
out_anova_df_file_name      <- "cap_anova_results.csv"
out_var_explained_file_name <- "var_explained.csv"

# Small epsilon for numerical stability in distances
JSD_EPS <- 0

# ── 1) Load Dataset for Analysis ────────────────────────────

# Load the phyloseq object for the selected analysis
# The ps1_M object should be used (filtered MetaPhlAn features)

# ── Load Dataset (MetaPhlAn Relative Abundance) ───────────────────────────
ps <- readRDS("...active_shotgun_ps1_M")

# Normalize to per-sample proportions; convert %→proportion; replace NA with 0
otu <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) otu <- t(otu)
otu[is.na(otu)] <- 0
medval <- suppressWarnings(median(otu, na.rm = TRUE))
if (is.finite(medval) && medval > 1) otu <- otu / 100
rs <- rowSums(otu, na.rm = TRUE)
nz <- rs > 0
if (any(nz)) otu[nz, ] <- otu[nz, , drop = FALSE] / rs[nz]
otu_table(ps) <- otu_table(otu, taxa_are_rows = FALSE)

# (Optional) drop samples with zero total abundance after filtering
keep_samples <- rowSums(otu_table(ps)) > 0
ps <- prune_samples(keep_samples, ps)

# ── 2) Prepare distance & metadata ───────────────────────────────────────────
# JSD on compositional data
# phyloseq::distance("jsd") expects sample rows in otu matrix internally (handled)
dist_mat <- phyloseq::distance(ps, method = "jsd", type = "samples", pseudocount = JSD_EPS)

# Community matrix for species scores (sites = samples, species = taxa)
comm_mat <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) comm_mat <- t(comm_mat)

# Metadata
metadata <- as(sample_data(ps), "data.frame")
metadata[[GROUPING_VAR]] <- as.factor(metadata[[GROUPING_VAR]])
# Condition covariate
if (!"Dataset_ID" %in% names(metadata)) stop("Missing 'Dataset_ID' in sample_data.")
metadata$Dataset <- as.factor(metadata$Dataset_ID)

# Ensure metadata rows align to dist labels
samp <- intersect(labels(dist_mat), rownames(metadata))
dist_mat <- as.dist(as.matrix(dist_mat)[samp, samp, drop = FALSE])
metadata <- metadata[samp, , drop = FALSE]
comm_mat <- comm_mat[samp, , drop = FALSE]

# ── 3) Run CAP Analysis (condition on Dataset) ───────────────────────────────
cap_model <- capscale(
  as.formula(paste0("dist_mat ~ ", GROUPING_VAR, " + Condition(Dataset)")),
  data = metadata,
  add = TRUE
)

# Sample scores (sites)
sample_scores_df <- as.data.frame(scores(cap_model, display = "sites"))
sample_scores_df <- cbind(sample_scores_df, metadata)
# Optional renaming CAP2→MDS1 for consistency with your previous outputs
colnames(sample_scores_df) <- sub("CAP2", "MDS1", colnames(sample_scores_df))

# Species (taxa) scores as weighted averages
sp_scores_df <- as.data.frame(wascores(scores(cap_model, display = "sites"), comm_mat))
sp_scores_df$ASV <- rownames(sp_scores_df)
colnames(sp_scores_df) <- sub("CAP2", "MDS1", colnames(sp_scores_df))

# Add taxonomy & a readable label
tax_df <- as.data.frame(tax_table(ps))
tax_df$ASV <- rownames(tax_df)

sp_scores_df <- left_join(sp_scores_df, tax_df, by = "ASV") %>%
  mutate(
    Label = dplyr::coalesce(Genus, Family, Order, Class, Phylum, Domain, ASV)
  )

# (Optional) Top 20 taxa by loading magnitude
if (all(c("CAP1","MDS1") %in% names(sp_scores_df))) {
  sp_scores_df <- sp_scores_df %>%
    mutate(loading = sqrt(CAP1^2 + MDS1^2)) %>%
    arrange(desc(loading)) %>%
    slice_head(n = 20)
}

# % variance explained for first two constrained axes
eig <- eigenvals(cap_model)
var_explained <- (eig / sum(eig)) * 100
var_explained <- var_explained[1:2]
var_explained <- data.frame(Axis = c("CAP1","CAP2"), Percent = as.numeric(var_explained))

# ANOVA (terms) with BH-adjusted p-values
anova_results <- anova.cca(cap_model, by = "terms", permutations = 999)
anova_df <- as.data.frame(anova_results)
if ("Pr(>F)" %in% names(anova_df)) {
  anova_df$P_adj_BH <- p.adjust(anova_df$`Pr(>F)`, method = "BH")
}

# ── 4) Save Results ──────────────────────────────────────────────────────────
write.csv(sp_scores_df,    file.path(out_dir, out_sp_scores_file_name),     row.names = FALSE)
write.csv(sample_scores_df,file.path(out_dir, out_sample_scores_file_name), row.names = FALSE)
write.csv(anova_df,        file.path(out_dir, out_anova_df_file_name))
write.csv(var_explained,   file.path(out_dir, out_var_explained_file_name), row.names = FALSE)
