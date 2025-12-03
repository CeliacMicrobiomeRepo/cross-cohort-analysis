# Beta‑diversity PERMANOVA (16S)

# Takes a phyloseq object with non-transformed counts and filters ASVs to remove zero variance ASVs
# This phyloseq object likely originated from the combine_trunc_asvs_to_phyloseq.R script

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

# Imports
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
GROUPING_VAR    <- "Diagnosed_Celiac"

# Create output folder
# e.g:
#   stool_prospective
#   stool_active
#   stool_treated
#   duodenum_active
out_dir <- "beta_diversity_results/stool_active"
dir.create(out_dir, recursive = TRUE)

# Output file name
out_file_name <- "beta_diversity_results_stool_active.csv"



# ── 1) Load Dataset for Analysis ────────────────────────────

# Load the phyloseq object for the selected analysis
# The ps1 object should be used (filtered ASVs, but non-transformed counts)
# e.g:
#   prospective_phyloseq_objects
#   stool_active_phyloseq_objects
#   stool_treated_phyloseq_objects
#   duodenum_phyloseq_objects
ps <- readRDS("/home/haig/Repos/meta-analysis/preprocessing/phyloseq_objects/stool_active_phyloseq_objects/ps1.rds")

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

as_matrix <- function(x) {
  m <- as.matrix(unclass(x))
  dimnames(m) <- dimnames(x)
  m
}

otu_matrix <- function(ps) {
  mat <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) mat <- t(mat)
  mat
}

tss_transform <- function(ps, pseudo = 0) {
  transform_sample_counts(ps, function(x) {
    tot <- sum(x)
    if (tot == 0) return(rep(0, length(x)))
    (x + pseudo) / (tot + pseudo * length(x))
  })
}

clr_transform <- function(ps) {
  otu <- otu_matrix(ps) + 1
  clr_otu <- compositions::clr(otu)
  otu_table(ps) <- otu_table(t(clr_otu), taxa_are_rows = TRUE)
  ps
}

distance_matrix <- function(mat, method) {
  mat <- as_matrix(mat)
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
         stop("Unsupported distance method: ", method))
}

run_permanova <- function(dist_mat, meta, var, permutations = 999) {
  meta_df <- as.data.frame(meta)
  if (!var %in% colnames(meta_df))
    return(tibble(Variable = var, R2 = NA_real_, p_value = NA_real_))
  samp <- intersect(labels(dist_mat), rownames(meta_df))
  if (length(samp) < 2)
    return(tibble(Variable = var, R2 = NA_real_, p_value = NA_real_))
  grp <- factor(meta_df[samp, var, drop = FALSE][[1]])
  if (nlevels(grp) < 2)
    return(tibble(Variable = var, R2 = NA_real_, p_value = NA_real_))
  perm <- vegan::adonis2(as.dist(as.matrix(dist_mat)[samp, samp]) ~ grp,
                         permutations = permutations)
  tibble(Variable = var, R2 = perm$R2[1], p_value = perm$`Pr(>F)`[1])
}

analyse_dataset <- function(ps, dataset_name, group_var,
                            dataset_id_var = NULL, permutations = 999) {
  run_all <- function(d, m, vars) map_dfr(vars, run_permanova, dist_mat = d, meta = m, permutations = permutations)

  res <- list()

  # TSS branch
  ps_tss  <- tss_transform(ps)
  mat_tss <- otu_matrix(ps_tss)
  meta_tss <- sample_data(ps_tss) |> as.data.frame()
  vars_tss <- c(group_var, if (!is.null(dataset_id_var) && dataset_name == "Pooled") dataset_id_var)
  dist_jsd <- distance_matrix(mat_tss, "JSD"); if (!is.null(dist_jsd)) res[["TSS_JSD"]] <- run_all(dist_jsd, meta_tss, vars_tss)
  res[["TSS_Bray"]]      <- run_all(distance_matrix(mat_tss, "bray"),      meta_tss, vars_tss)
  res[["TSS_Euclidean"]] <- run_all(distance_matrix(mat_tss, "euclidean"), meta_tss, vars_tss)

  # CLR Euclidean
  ps_clr <- clr_transform(ps)
  mat_clr <- otu_matrix(ps_clr)
  meta_clr <- sample_data(ps_clr) |> as.data.frame()
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



# ── 3) Main ────────────────────────────

results_tbl <- gather_results(
    ps_list,
    pooled_ps      = ps,           
    group_var      = GROUPING_VAR,
    dataset_id_var = "Dataset_ID"
)
print(results_tbl)

write.csv(results_tbl, file.path(out_dir, out_file_name), row.names = FALSE)


