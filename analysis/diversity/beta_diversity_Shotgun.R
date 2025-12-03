# Beta‑diversity PERMANOVA (Shotgun)

# Takes a phyloseq object with transformed counts and computes unconstrained beta-diversity
# with the following parameters.

# For each phyloseq object in `ps_list`:
#   1. Apply Total‑Sum Scaling (TSS) → Jensen–Shannon (JSD), Bray–Curtis, Euclidean.
#   2. Apply centred‑log‑ratio (CLR) → Euclidean.
#   3. Run PERMANOVA (adonis2) on a primary variable (`group_var`).
#      When analysing the pooled object, also test `Dataset_ID`.
# Results: tidy tibble — Dataset | Transformation_Distance | Variable | R2 | p_value

# Requirements:
#  - R 4.5.1
#  - R packages: phyloseq, vegan, philentropy, compositions, tibble, dplyr, purrr


# ── 0) Set Up ────────────────────────────

library(phyloseq)
library(vegan)
library(philentropy)
library(compositions)
library(tibble)
library(dplyr)
library(purrr)

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
out_dir <- "Shotgun/beta_diversity_results/stool_active"
dir.create(out_dir, recursive = TRUE)


# Output file name
out_file_name <- "beta_diversity_results_stool_active.csv"


# ── 1) Load Dataset for Analysis ────────────────────────────

# Load the phyloseq object for the selected analysis
# The ps1_M object should be used (filtered MetaPhlAn features)

# ── Load Dataset (MetaPhlAn Relative Abundance) ───────────────────────────
ps <- readRDS("...active_shotgun_ps1_M")

# Small epsilon for CLR zero-replacement
CLR_EPS <- 1e-6

# Coerce to proportions (accepts 0–1 or 0–100); replace NA with 0
otu <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) otu <- t(otu)
otu[is.na(otu)] <- 0
medval <- suppressWarnings(median(otu, na.rm = TRUE))
if (is.finite(medval) && medval > 1) otu <- otu / 100
rs <- rowSums(otu, na.rm = TRUE)
nz <- rs > 0
if (any(nz)) otu[nz, ] <- otu[nz, , drop = FALSE] / rs[nz]
otu_table(ps) <- otu_table(otu, taxa_are_rows = FALSE)

# Split by dataset; drop taxa with zero total RA within each dataset
sd_df <- as(sample_data(ps), "data.frame")
dataset_ids <- names(which(table(sd_df$Dataset_ID) > 0))
ps_list <- lapply(dataset_ids, function(d) {
    samps  <- rownames(sd_df)[sd_df$Dataset_ID == d]
    ps_sub <- prune_samples(samps, ps)
    prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
})
names(ps_list) <- dataset_ids


# ── 2) Functions ─────────────────────────────────────────────────────────────
as_matrix <- function(x) {
    m <- as.matrix(unclass(x))
    dimnames(m) <- dimnames(x)
    m
}

otu_matrix <- function(ps) {
    mat <- as(otu_table(ps), "matrix")
    if (taxa_are_rows(ps)) mat <- t(mat) # rows = samples
    mat
}

# TSS is identity on RA (kept for clarity/safety)
tss_transform <- function(ps) {
    transform_sample_counts(ps, function(x) {
        s <- sum(x, na.rm = TRUE); if (s == 0) return(x * 0)
        x / s
    })
}

# CLR with small multiplicative replacement (works for RA)
clr_transform <- function(ps, eps = CLR_EPS) {
    X <- otu_matrix(ps)
    # replace zeros, renormalize
    X <- sweep(X + eps, 1, rowSums(X + eps), "/")
    clrX <- compositions::clr(X)               # rows = samples
    otu_table(ps) <- otu_table(t(clrX), taxa_are_rows = TRUE)
    ps
}

distance_matrix <- function(mat, method) {
    mat <- as_matrix(mat)
    # philentropy expects rows = samples, non-negative, rows sum to 1 for JSD
    safe_jsd <- function(m) {
        tryCatch({
            d <- philentropy::distance(m, method = "jensen-shannon")
            dimnames(d) <- list(rownames(m), rownames(m))
            as.dist(d)
        }, error = function(e) {
            warning("JSD failed: ", conditionMessage(e)); NULL
        })
    }
    switch(method,
           JSD       = safe_jsd(mat),
           bray      = vegan::vegdist(mat, method = "bray"),
           euclidean = vegan::vegdist(mat, method = "euclidean"),
           stop("Unsupported distance: ", method))
}

run_permanova <- function(dist_mat, meta, var, permutations = 999) {
    meta_df <- as.data.frame(meta)
    
    # 0) column exists?
    if (!is.character(var) || length(var) != 1 || !var %in% colnames(meta_df)) {
        return(tibble(Variable = var, R2 = NA_real_, p_value = NA_real_))
    }
    
    # 1) shared sample IDs (keep order of dist)
    samp <- intersect(labels(dist_mat), rownames(meta_df))
    if (length(samp) < 2) {
        return(tibble(Variable = var, R2 = NA_real_, p_value = NA_real_))
    }
    
    # 2) robust column extraction → atomic vector
    col_raw <- meta_df[samp, var, drop = FALSE][[1]]  # [[1]] ensures a vector even if df/matrix/list
    if (is.list(col_raw)) col_raw <- unlist(col_raw, use.names = FALSE)
    col_raw <- as.vector(col_raw)
    
    # 3) drop rows with NA in the grouping var (and mirror in the distance)
    keep <- !is.na(col_raw)
    if (!any(keep)) {
        return(tibble(Variable = var, R2 = NA_real_, p_value = NA_real_))
    }
    samp2   <- samp[keep]
    grp_vec <- col_raw[keep]
    
    # 4) need >= 2 levels
    grp <- factor(grp_vec)
    if (nlevels(grp) < 2) {
        return(tibble(Variable = var, R2 = NA_real_, p_value = NA_real_))
    }
    
    # 5) subset the distance matrix to matching samples
    D <- as.matrix(dist_mat)
    D <- as.dist(D[samp2, samp2, drop = FALSE])
    
    perm <- vegan::adonis2(D ~ grp, permutations = permutations)
    tibble(Variable = var, R2 = perm$R2[1], p_value = perm$`Pr(>F)`[1])
}


analyse_dataset <- function(ps, dataset_name, group_var,
                            dataset_id_var = NULL, permutations = 999) {
    run_all <- function(d, m, vars) map_dfr(vars, run_permanova, dist_mat = d, meta = m, permutations = permutations)
    
    res <- list()
    
    # TSS (on RA it's a no-op normalization, safe for JSD/Bray/Euclid)
    ps_tss   <- tss_transform(ps)
    mat_tss  <- otu_matrix(ps_tss)
    meta_tss <- as.data.frame(sample_data(ps_tss))
    vars_tss <- c(group_var, if (!is.null(dataset_id_var) && dataset_name == "Pooled") dataset_id_var)
    
    dist_jsd <- distance_matrix(mat_tss, "JSD")
    if (!is.null(dist_jsd)) res[["TSS_JSD"]] <- run_all(dist_jsd, meta_tss, vars_tss)
    res[["TSS_Bray"]]      <- run_all(distance_matrix(mat_tss, "bray"),      meta_tss, vars_tss)
    res[["TSS_Euclidean"]] <- run_all(distance_matrix(mat_tss, "euclidean"), meta_tss, vars_tss)
    
    # CLR (epsilon replacement) → Euclidean
    ps_clr   <- clr_transform(ps, eps = CLR_EPS)
    mat_clr  <- otu_matrix(ps_clr)
    meta_clr <- as.data.frame(sample_data(ps_clr))
    vars_clr <- c(group_var, if (!is.null(dataset_id_var) && dataset_name == "Pooled") dataset_id_var)
    res[["CLR_Euclidean"]] <- run_all(distance_matrix(mat_clr, "euclidean"), meta_clr, vars_clr)
    
    bind_rows(res, .id = "Transformation_Distance") %>%
        mutate(Dataset = dataset_name, .before = 1)
}

gather_results <- function(ps_list, pooled_ps, group_var,
                           dataset_id_var = "Dataset_ID", permutations = 999) {
    per_ds <- imap(ps_list, analyse_dataset, group_var = group_var,
                   dataset_id_var = NULL, permutations = permutations)
    pooled <- analyse_dataset(pooled_ps, "Pooled", group_var,
                              dataset_id_var = dataset_id_var,
                              permutations = permutations)
    bind_rows(per_ds, pooled)
}


# ── 3) Main ──────────────────────────────────────────────────────────────────
results_tbl <- gather_results(
    ps_list,
    pooled_ps      = ps,          # already standardized to RA proportions
    group_var      = GROUPING_VAR,
    dataset_id_var = "Dataset_ID"
)

print(results_tbl)
write.csv(results_tbl, file.path(out_dir, out_file_name), row.names = FALSE)