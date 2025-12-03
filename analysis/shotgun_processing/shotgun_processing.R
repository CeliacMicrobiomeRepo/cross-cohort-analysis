# Shotgun Processing (script 1)

# Uses the celiac repository to produce both MetaPhlAn and deterministic phyloseq objects
# for each shotgun dataset. Outputs are used in the shotgun_group_seperation.R script.

# MetaPhlAn (relative abundances) are used for all analysis except differential abundance.
# Determinstic abslute counts are used for differential abundance only.

# Requirements:
#  - R 4.5.1
#  - R packages: readr, dplyr, stringr, tibble, purrr, phyloseq, gert



# ── 0) Set Up ────────────────────────────

# Imports
library(readr)
library(dplyr)
library(stringr)
library(tibble)
library(purrr)
library(phyloseq)
library(gert)



# ── 1) Clone Celiac Microbiome Repository ────────────────────────────────────
if (!dir.exists("celiac-repository")) {
    git_clone(
        url  = "https://github.com/CeliacMicrobiomeRepo/celiac-repository.git",
        path = "celiac-repository"
    )
}

# ── 2) CONFIG (edit here) ────────────────────────────────────────────────────
repo_dir        <- "celiac-repository"
meta_dir        <- file.path(repo_dir, "SG_datasets", "SG_80_Mouzan")
profiles_dir    <- file.path(meta_dir, "metaphlan_outputs", "gtdb")
sra_table_file  <- file.path(meta_dir, "SraRunTable.csv")
sra_result_file <- file.path(meta_dir, "sra_result.csv")
metadata_file   <- file.path(repo_dir, "all_samples.tsv")
out_prefix      <- "80_Mouzan"

paired_end <- TRUE # <--- set FALSE for single-end; TRUE for paired-end
# =============================================================================



# ── 3) Read & merge MetaPhlAn profiles (SRR) ─────────────────────────────────
profile_files <- list.files(
    profiles_dir,
    pattern = "_profile\\.(txt|tsv)$",
    full.names = TRUE,
    recursive = FALSE
)
stopifnot(length(profile_files) > 0)

read_profile <- function(f) {
    sample_id <- basename(f) |> str_remove("_profile\\.(txt|tsv)$")
    df <- readr::read_tsv(f, comment = "#", col_names = FALSE, show_col_types = FALSE)
    if (ncol(df) < 2) stop("Unexpected MetaPhlAn format: ", f)
    clade <- df[[1]]
    ra    <- if (ncol(df) >= 3 && is.numeric(df[[3]])) df[[3]] else df[[2]]
    tibble(clade_name = as.character(clade), !!sample_id := as.numeric(ra))
}

merged <- Reduce(function(a, b) dplyr::full_join(a, b, by = "clade_name"),
                 lapply(profile_files, read_profile)) %>%
    distinct(clade_name, .keep_all = TRUE) %>%
    filter(clade_name != "UNCLASSIFIED")

# keep species/strain only
subset_ab <- merged %>% filter(str_detect(clade_name, "(^|[|;])(s__|t__)"))

# ── 4) Parse taxonomy ────────────────────────────────────────────────────────
get_rank <- function(x, tag) {
    hit <- str_extract(x, paste0("(^|[|;])", tag, "[^|;]*"))
    ifelse(is.na(hit), NA_character_, str_replace(hit, paste0("(^|[|;])", tag), ""))
}

taxonomy_df <- subset_ab %>%
    transmute(
        taxon   = clade_name,
        Domain  = get_rank(clade_name, "d__"),
        Phylum  = get_rank(clade_name, "p__"),
        Class   = get_rank(clade_name, "c__"),
        Order   = get_rank(clade_name, "o__"),
        Family  = get_rank(clade_name, "f__"),
        Genus   = get_rank(clade_name, "g__"),
        Species = get_rank(clade_name, "s__") |> str_replace_all("_", " "),
        Strain  = get_rank(clade_name, "t__")
    ) %>%
    mutate(across(-taxon, ~ na_if(.x, ""))) %>%
    as.data.frame()
rownames(taxonomy_df) <- taxonomy_df$taxon
taxonomy_df$taxon <- NULL
taxonomy <- as.matrix(taxonomy_df)

# ── 5) Build OTU (RA). Convert % -> proportions; TSS per sample ───────────────
sample_cols <- setdiff(names(subset_ab), "clade_name")
otu_prop_df <- subset_ab %>%
    select(clade_name, all_of(sample_cols)) %>%
    as.data.frame()
rownames(otu_prop_df) <- otu_prop_df$clade_name
otu_prop_df$clade_name <- NULL

median_val <- median(unlist(otu_prop_df), na.rm = TRUE)
if (is.finite(median_val) && median_val > 1) {
    otu_prop_df[] <- lapply(otu_prop_df, function(x) as.numeric(x) / 100)
} else {
    otu_prop_df[] <- lapply(otu_prop_df, as.numeric)
}
otu_prop_df[is.na(otu_prop_df)] <- 0

otu_prop <- as.matrix(otu_prop_df)
# TSS per sample to guard against drift in GTDB outputs
col_sums <- colSums(otu_prop, na.rm = TRUE)
nz <- which(col_sums > 0)
if (length(nz) > 0) {
    otu_prop[, nz] <- sweep(otu_prop[, nz, drop = FALSE], 2, col_sums[nz], "/")
}

ps_ra <- phyloseq(otu_table(otu_prop, taxa_are_rows = TRUE), tax_table(taxonomy))


# ── 6) Remove samples with no species/strain mass ────────────────────────────
# mat_prop is the species/strain RA matrix (after %→prop, NA→0, TSS on nonzero cols)
mat_prop <- as.matrix(otu_prop)
if (!is.null(dimnames(otu_prop))) dimnames(mat_prop) <- dimnames(otu_prop)

zero_cols <- which(colSums(mat_prop) == 0 | is.na(colSums(mat_prop)))
drop_ids  <- colnames(mat_prop)[zero_cols]

if (length(drop_ids) > 0) {
  message("Dropping ", length(drop_ids), " zero-signal runs: ",
          paste(drop_ids, collapse = ", "))

  # Remove from RA matrix
  mat_prop <- mat_prop[, setdiff(colnames(mat_prop), drop_ids), drop = FALSE]

  # Keep SRA mapping and table aligned to remaining runs
  sra_tbl   <- sra_tbl   %>% filter(Run %in% colnames(mat_prop))
  run_to_exp <- setNames(str_trim(sra_tbl$Experiment), str_trim(sra_tbl$Run))

  # Rebuild ps_ra on the filtered set
  ps_ra <- phyloseq(
    otu_table(mat_prop, taxa_are_rows = TRUE),
    tax_table(taxonomy)
  )
} else {
  message("No zero-signal runs to drop.")
}

# ── 7) SRA tables & depths ───────────────────────────────────────────────────
sra_tbl <- read_csv(sra_table_file,  show_col_types = FALSE)
sra_res <- read_csv(sra_result_file, show_col_types = FALSE)

multiplier <- if (paired_end) 2 else 1
reads_srx <- sra_res %>%
    transmute(
        Experiment      = str_trim(`Experiment Accession`),
        total_spots_num = as.numeric(gsub("[^0-9]", "", `Total Spots`)),
        reads           = round(total_spots_num * multiplier)
    ) %>%
    transmute(Experiment, reads) %>%
    deframe()

run_to_exp <- setNames(str_trim(sra_tbl$Experiment), str_trim(sra_tbl$Run))

# ── 8) Deterministic counts via Largest-Remainder ────────────────────────────
largest_remainder_counts <- function(p, depth) {
    if (depth <= 0 || all(p == 0)) return(integer(length(p)))
    raw <- p * depth
    flo <- floor(raw)
    need <- depth - sum(flo)
    if (need > 0) {
        rem <- raw - flo
        idx <- order(rem, decreasing = TRUE)[seq_len(need)]
        flo[idx] <- flo[idx] + 1L
    }
    as.integer(flo)
}

mat_prop <- as(otu_table(ps_ra), "matrix")
if (!otu_table(ps_ra)@taxa_are_rows) mat_prop <- t(mat_prop)

runs <- colnames(mat_prop)
reads_run <- reads_srx[ run_to_exp[runs] ]
names(reads_run) <- runs
stopifnot(!any(is.na(reads_run)))  # fail fast if any run is missing a depth

count_mat <- vapply(
    runs,
    FUN.VALUE = integer(nrow(mat_prop)),
    FUN = function(run_id) largest_remainder_counts(mat_prop[, run_id], as.integer(reads_run[[run_id]]))
)
dimnames(count_mat) <- dimnames(mat_prop)

# ── 9) Deduplicate taxa & sync taxonomy ──────────────────────────────────────
if (anyDuplicated(rownames(count_mat))) {
    count_mat <- rowsum(count_mat, group = rownames(count_mat)) |> as.matrix()
}
uniq_tax <- taxonomy[!duplicated(rownames(taxonomy)), , drop = FALSE]
common   <- intersect(rownames(count_mat), rownames(uniq_tax))
count_mat <- count_mat[common, , drop = FALSE]
uniq_tax  <- uniq_tax[common,  , drop = FALSE]

# ── 10) Sample metadata ─────────────────────────────────────────────────────
sample_md <- sra_tbl %>% select(Run, Experiment) %>% column_to_rownames("Run")

ps_counts <- phyloseq(
    otu_table(count_mat, taxa_are_rows = TRUE),
    tax_table(uniq_tax),
    sample_data(sample_md[colnames(count_mat), , drop = FALSE])
)

# ── 11) Attach metadata to BOTH ps_counts and ps_ra ──────────────────────────
if (file.exists(metadata_file)) {
    meta_df_raw <- read_tsv(metadata_file, show_col_types = FALSE)
    
    # We'll allow several ways to match:
    #  1) rownames(meta_df) == Run
    #  2) rownames(meta_df) == Experiment
    #  3) meta_df has explicit Run / Experiment columns to join on
    meta_df <- meta_df_raw
    
    # Standardize: keep both a rowname index (Sample_ID) and any explicit columns if present
    if ("Sample_ID" %in% names(meta_df)) {
        meta_df <- meta_df %>% mutate(Sample_ID = as.character(Sample_ID))
    } else {
        # create a synthetic Sample_ID if not present
        meta_df <- meta_df %>% mutate(Sample_ID = dplyr::coalesce(.data$Run, .data$Experiment))
    }
    
    # Helper: build a run-indexed sample_data frame for a given phyloseq object
    make_md_for_ps <- function(ps_obj) {
        runs_here <- sample_names(ps_obj)
        
        # Start from SRA table (Run + Experiment) so we always have those two
        md <- sra_tbl %>%
            transmute(Run = str_trim(Run), Experiment = str_trim(Experiment)) %>%
            filter(Run %in% runs_here)
        
        # Strategy A: Sample_ID == Run
        md <- md %>%
            left_join(meta_df %>% select(where(~!is.list(.x))) %>% # avoid list-cols
                          distinct(Sample_ID, .keep_all = TRUE),
                      by = c("Run" = "Sample_ID"))
        
        # Strategy B: if we still have many NAs, try Sample_ID == Experiment
        still_sparse <- function(df) {
            # consider “sparse” if >50% rows have all metadata cols NA (beyond Run/Experiment)
            meta_cols <- setdiff(names(df), c("Run","Experiment"))
            if (length(meta_cols) == 0) return(TRUE)
            mean(rowSums(!is.na(df[meta_cols])) == 0) > 0.5
        }
        if (still_sparse(md)) {
            md <- sra_tbl %>%
                transmute(Run = str_trim(Run), Experiment = str_trim(Experiment)) %>%
                filter(Run %in% runs_here) %>%
                left_join(meta_df %>% distinct(Sample_ID, .keep_all = TRUE),
                          by = c("Experiment" = "Sample_ID"))
        }
        
        # Strategy C: explicit columns in metadata (Run / Experiment) if available
        if (still_sparse(md)) {
            tmp <- sra_tbl %>%
                transmute(Run = str_trim(Run), Experiment = str_trim(Experiment)) %>%
                filter(Run %in% runs_here)
            if ("Run" %in% names(meta_df)) {
                tmp <- tmp %>%
                    left_join(meta_df %>% mutate(Run = str_trim(Run)), by = "Run")
            }
            if (still_sparse(tmp) && "Experiment" %in% names(meta_df)) {
                tmp <- sra_tbl %>%
                    transmute(Run = str_trim(Run), Experiment = str_trim(Experiment)) %>%
                    filter(Run %in% runs_here) %>%
                    left_join(meta_df %>% mutate(Experiment = str_trim(Experiment)), by = "Experiment")
            }
            # Use tmp if it improved coverage
            md <- tmp
        }
        
        # Finalize rownames = Run, keep order = sample_names(ps_obj)
        md <- md %>% tibble::column_to_rownames("Run")
        md <- md[runs_here, , drop = FALSE]
        
        # Drop duplicate columns that arose from joins
        md <- md[, !duplicated(names(md)), drop = FALSE]
        
        md
    }
    
    md_counts <- make_md_for_ps(ps_counts)
    md_ra     <- make_md_for_ps(ps_ra)
    
    sample_data(ps_counts) <- sample_data(md_counts)
    sample_data(ps_ra)     <- sample_data(md_ra)
}


# ── 12) Sanity checks ─────────────────────────────────────────────────────────
stopifnot(all(colSums(count_mat) == as.integer(reads_run[colnames(count_mat)])))
stopifnot(all(count_mat >= 0))

# ── 13) Save ──────────────────────────────────────────────────────────────────
saveRDS(ps_counts, paste0(out_prefix, "_deterministic.rds"))
saveRDS(ps_ra,     paste0(out_prefix, "_relative_abundance.rds"))
