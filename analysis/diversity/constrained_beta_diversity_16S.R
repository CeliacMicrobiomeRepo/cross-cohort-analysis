# Constrained Beta-Diversity Analysis (CAP)

# Takes a phyloseq object, performs a TSS transformation, and then runs a 
# constrained analysis of principal coordinates (CAP) to identify how much
# of the variation in community composition can be explained by a given variable,
# while conditioning on a confounding variable.

# Requirements:
#  - R 4.5.1
#  - R packages: phyloseq, vegan, ggplot2, dplyr



# ── 0) Set Up ────────────────────────────

# Imports
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)

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
out_dir <- "16S/constrained_beta_diversity_results/stool_active"
dir.create(out_dir, recursive = TRUE)

# Output file names
out_sp_scores_file_name <- "cap_loading_data.csv"
out_sample_scores_file_name <- "cap_sample_data.csv"
out_anova_df_file_name <- "cap_anova_results.csv"
out_var_explained_file_name <- "var_explained.csv"



# ── 1) Load Dataset for Analysis ────────────────────────────

# Load the phyloseq object for the selected analysis
# The ps2 object should be used (filtered ASVs with TSS-transformed counts)
# e.g:
#   prospective_phyloseq_objects
#   stool_active_phyloseq_objects
#   stool_treated_phyloseq_objects
#   duodenum_phyloseq_objects
ps <- readRDS("/home/haig/Repos/meta-analysis/preprocessing/phyloseq_objects/stool_active_phyloseq_objects/ps2.rds")



# ── 2) Prepare Data ────────────────────────────

# Produce distance object this example is JSD
dist_mat <- phyloseq::distance(ps, method = "jsd", type = "samples")

# Prepare the community matrix and metadata (as before)
comm_mat <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) {
  comm_mat <- t(comm_mat)
}

metadata <- as(sample_data(ps), "data.frame")
metadata[[GROUPING_VAR]] <- as.factor(metadata[[GROUPING_VAR]])
metadata$Dataset <- as.factor(metadata$Dataset_ID)



# ── 3) Run CAP Analysis ────────────────────────────

# Run CAP analysis (condition on Dataset.ID)
cap_model  <- capscale(
  as.formula(paste0("dist_mat ~ ", GROUPING_VAR, " + Condition(Dataset)")),
  data = metadata, add = TRUE
)


# Extract sample scores and rename CAP2 -> MDS1 if you prefer
sample_scores_df <- as.data.frame(scores(cap_model, display = "sites"))
sample_scores_df <- cbind(sample_scores_df, metadata)
colnames(sample_scores_df) <- sub("CAP2", "MDS1", colnames(sample_scores_df))

# Extract species (ASV) scores, rename CAP2 -> MDS1
sp_scores_df <- as.data.frame(wascores(scores(cap_model, display = "sites"), comm_mat))
sp_scores_df$ASV <- rownames(sp_scores_df)
colnames(sp_scores_df) <- sub("CAP2", "MDS1", colnames(sp_scores_df))

# Merge taxonomy and create a label for each ASV
tax_df <- as.data.frame(tax_table(ps))
tax_df$ASV <- rownames(tax_df)

sp_scores_df <- left_join(sp_scores_df, tax_df, by = "ASV") %>%
  mutate(
    Label = case_when(
      !is.na(Genus)   ~ Genus,
      !is.na(Family)  ~ Family,
      !is.na(Order)   ~ Order,
      !is.na(Class)   ~ Class,
      !is.na(Phylum)  ~ Phylum,
      !is.na(Kingdom) ~ Kingdom,
      TRUE            ~ ASV  # fallback if all else is NA
    )
  )

# (Optional) Filter top 20 taxa by loading
sp_scores_df <- sp_scores_df %>%
  mutate(loading = sqrt(CAP1^2 + MDS1^2)) %>%
  arrange(desc(loading)) %>%
  slice_head(n = 20)

# Calculate variation for CAP1 and MDS1
var_explained <- eigenvals(cap_model) / sum(eigenvals(cap_model)) * 100
var_explained <- var_explained[1:2]  # for CAP1 and MDS1

anova_results <- anova.cca(cap_model, by = "terms", permutations = 999)
anova_df <- as.data.frame(anova_results)
anova_df$P_adj_BH <- p.adjust(anova_df$`Pr(>F)`, method = "BH")



# ── 4) Save Results ────────────────────────────────────────────────────

write.csv(sp_scores_df,   file = file.path(out_dir, out_sp_scores_file_name))
write.csv(sample_scores_df,file = file.path(out_dir, out_sample_scores_file_name))
write.csv(anova_df,file = file.path(out_dir, out_anova_df_file_name)) 
write.csv(var_explained,  file = file.path(out_dir, out_var_explained_file_name))


