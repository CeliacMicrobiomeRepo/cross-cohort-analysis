# Alpha Diversity Code (16S)

# Takes a phyloseq object with non-transformed counts and filters ASVs to remove zero variance ASVs
# This phyloseq object likely originated from the combine_trunc_asvs_to_phyloseq.R script

# Requirements:
#  - R 4.5.1
#  - R packages: microbiome, meta



# ── 0) Set Up ────────────────────────────

# Imports
library(microbiome)
library(meta)

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
out_dir <- "16S/alpha_diversity_results/stool_active_yuck"
dir.create(out_dir, recursive = TRUE)

# Output file name
out_file_name <- "meta_summary_metrics.csv"



# ── 1) Load Dataset for Analysis ────────────────────────────

# Load the phyloseq object for the selected analysis
# The ps0 object should be used (unfiltered ASVs)
# e.g:
#   prospective_phyloseq_objects
#   stool_active_phyloseq_objects
#   stool_treated_phyloseq_objects
#   duodenum_phyloseq_objects
ps <- readRDS("/home/haig/Repos/meta-analysis/preprocessing/phyloseq_objects/stool_active_phyloseq_objects/ps0.rds")

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



# ── 2) Define Analysis Functions and Metrics ──────────────────────────

# Define your metric mapping
metric_mapping <- list(
    observed       = "observed",
    shannon        = "diversity_shannon",
    dominance      = "dominance_simpson",
    rare_abundance = "rarity_low_abundance"
)

# function to get per‐metric stats + counts
calculate_stats_for_metric <- function(metric, merged_df, mapping, pseudocount = 1e-6) {
    actual_metric <- mapping[[metric]]
    if (is.null(actual_metric) || !(actual_metric %in% colnames(merged_df))) {
        stop("Metric ", metric, " not found in merged data.")
    }
    
    df <- merged_df
    df$log2_metric <- log2(df[[actual_metric]] + pseudocount)
    
    # linear model for GROUPING_VAR effect
    lm_formula <- as.formula(paste("log2_metric ~", GROUPING_VAR, "+ Read_depth"))
    lmcoefs <- summary(lm(lm_formula, data = df))$coefficients
    
    # overall
    mean_all <- mean(df$log2_metric, na.rm=TRUE)
    se_all   <- sd(df$log2_metric,   na.rm=TRUE)/sqrt(nrow(df))
    
    # group‐wise mean, sd & N
    get_stats <- function(flag) {
        v <- df$log2_metric[df[[GROUPING_VAR]] == flag]
        c(N   = length(v),
          M   = if(length(v)) mean(v, na.rm=TRUE) else NA,
          SD  = if(length(v)) sd(v,   na.rm=TRUE) else NA)
    }
    s_false <- get_stats("FALSE")
    s_true  <- get_stats("TRUE")
    
    # extract effect
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
        Metric             = metric,
        Mean_Log2_All      = mean_all,
        SE_Log2_All        = se_all,
        N_FALSE            = s_false["N"],
        Mean_Log2_FALSE    = s_false["M"],
        SD_Log2_FALSE      = s_false["SD"],
        N_TRUE             = s_true["N"],
        Mean_Log2_TRUE     = s_true["M"],
        SD_Log2_TRUE       = s_true["SD"],
        Coef_Celiac        = est,
        SE_Celiac          = se,
        P_Celiac           = pval,
        Significance       = sig,
        stringsAsFactors   = FALSE
    )
    names(res)[names(res) == 'Coef_Celiac'] <- paste0('Coef_', GROUPING_VAR)
    names(res)[names(res) == 'SE_Celiac']   <- paste0('SE_', GROUPING_VAR)
    names(res)[names(res) == 'P_Celiac']    <- paste0('P_', GROUPING_VAR)
    res
}

# wrapper to run on one phyloseq object
run_alpha_analysis <- function(physeq_obj, pseudocount=1e-6) {
    alpha_div <- microbiome::alpha(
        physeq_obj,
        index = c("observed","shannon","fisher","bulla",
                  "simpson","dominance","low_abundance","rare_abundance")
    )
    alpha_div$Sample.ID <- row.names(alpha_div)
    
    meta_df <- data.frame(sample_data(physeq_obj))
    meta_df$Sample.ID   <- row.names(meta_df)
    meta_df$Read_depth  <- sample_sums(physeq_obj)
    
    merged_df <- merge(alpha_div, meta_df, by="Sample.ID")
    
    # only keep metrics we actually have
    valid_map   <- metric_mapping[
        vapply(metric_mapping, `%in%`, logical(1), colnames(merged_df))
    ]
    valid_names <- names(valid_map)
    if(length(valid_names)==0) stop("No metrics found in merged data!")
    
    do.call(rbind,
            lapply(valid_names, function(m) {
                calculate_stats_for_metric(m, merged_df, valid_map, pseudocount)
            })
    )
}



# ── 3) Run Per-Dataset Analysis ─────────────────────────────────────

# run across your list and write CSVs
#    assume ps_list is named, e.g. list(Study1=phy1, Study2=phy2, ...)
for(id in names(ps_list)) {
    tab <- run_alpha_analysis(ps_list[[id]])
    tab$Study <- id
    write.csv(tab, file.path(out_dir, paste0(id, "_alpha_analysis.csv")), row.names=FALSE)
    message("Wrote ", id, "_alpha_analysis.csv")
}



# ── 4) Perform Meta-Analysis Across Datasets ───────────────────────────

# load them back & do random‐effects meta‐analysis per metric
files <- list.files(path = out_dir, pattern = "_alpha_analysis\\.csv$", full.names = TRUE)
allres <- do.call(rbind, lapply(files, read.csv, stringsAsFactors=FALSE))

meta_out <- lapply(unique(allres$Metric), function(metric_name) {
    dat <- subset(allres, Metric==metric_name)
    m <- metacont(
        n.e      = dat$N_TRUE,
        mean.e   = dat$Mean_Log2_TRUE,
        sd.e     = dat$SD_Log2_TRUE,
        n.c      = dat$N_FALSE,
        mean.c   = dat$Mean_Log2_FALSE,
        sd.c     = dat$SD_Log2_FALSE,
        studlab  = dat$Study,
        data     = dat,
        sm       = "SMD",       # standardized mean difference
        method.tau = "REML",    # random‐effects τ² by REML
        common = FALSE,
        random= TRUE
    )
    list(metric=metric_name, meta=m)
})



# ── 5) Display Meta-Analysis Results ────────────────────────────────

# print summaries
for(res in meta_out) {
    cat("\n===== Metric:", res$metric, "=====\n")
    print(summary(res$meta))
}



# ── 6) Summarize and Export Meta-Analysis Results ───────────────────

# Tidy up the main output stats and write them to CSV, now grabbing zval.random & pval.random 
meta_summary_list <- lapply(meta_out, function(x) {
    m <- x$meta
    # helper to extract a scalar or return NA
    get1 <- function(obj, name) {
        v <- obj[[name]]
        if (!is.null(v) && length(v) == 1) v else NA_real_
    }
    
    data.frame(
        Metric       = x$metric,
        TE_random    = get1(m, "TE.random"),
        seTE_random  = get1(m, "seTE.random"),
        z_random     = get1(m, "zval.random"),    # <— changed here
        p_random     = get1(m, "pval.random"),    # <— and here
        ci_lower     = get1(m, "lower.random"),
        ci_upper     = get1(m, "upper.random"),
        tau2         = if (!is.null(m[["tau2"]]) && length(m[["tau2"]])==1) {
            m[["tau2"]]
        } else if (!is.null(m[["tau^2"]]) && length(m[["tau^2"]])==1) {
            m[["tau^2"]]
        } else NA_real_,
        I2           = get1(m, "I2"),
        stringsAsFactors = FALSE
    )
})

meta_summary_df <- do.call(rbind, meta_summary_list)
write.csv(
    meta_summary_df,
    file      = file.path(out_dir, out_file_name),
    row.names = FALSE
)

