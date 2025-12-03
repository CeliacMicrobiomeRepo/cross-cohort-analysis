"""
pearson_agreement_within_groups.py

Computes Pearson correlation of effect sizes within groups (e.g. 16S vs shotgun)
Please change file names below as needed
Uses ANCOM-BC2 sig_taxa_combined files at both genus and species levels
"""

# Import libraries
import os
import pandas as pd
import numpy as np
from scipy.stats import pearsonr


# ── 1) Set working directory ────────────────────────────────────────────
os.chdir(r"C:/PhD/Meta_analysis/FINAl/Agreement")  # your folder


# ── 2) Helper: pick effect size (lfc_pooled, else meta) ────────────────────
def add_effect_column(df,
                      lfc_col="lfc_pooled",
                      meta_candidates=("meta_pooled", "estimate_meta"),
                      out_col="effect"):
    """
    Creates df[out_col] using:
      1) df[lfc_col]
      2) if NA, fallback to first existing column in meta_candidates
    """
    df = df.copy()

    if lfc_col in df.columns:
        effect = pd.to_numeric(df[lfc_col], errors="coerce")
    else:
        effect = pd.Series([np.nan] * len(df), index=df.index)

    meta_col = None
    for c in meta_candidates:
        if c in df.columns:
            meta_col = c
            break

    if meta_col is not None:
        meta_vals = pd.to_numeric(df[meta_col], errors="coerce")
        effect = effect.where(~effect.isna(), meta_vals)

    df[out_col] = effect
    return df



# ── 3) Helpers for taxon columns & Pearson agreement ───────────────────────
def choose_taxon_col(df: pd.DataFrame, candidates):
    """
    Pick the first column name from `candidates` that exists in df.
    e.g. candidates = ("genus", "Genus")
    """
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(f"None of {candidates} found in columns: {list(df.columns)}")

def pearson_agreement_16s_vs_shotgun(
    df_16s: pd.DataFrame,
    df_shotgun: pd.DataFrame,
    taxon_candidates,
    min_n: int = 2
):
    """
    Aligns 16S and shotgun on a taxon column (e.g. 'genus'/ 'Genus'),
    computes Pearson r.

    Returns:
        r (float or np.nan),
        p (float or np.nan),
        merged (DataFrame),
        n_overlap (int)
    """
    df_16s = add_effect_column(df_16s, out_col="effect_16s")
    df_shotgun = add_effect_column(df_shotgun, out_col="effect_shotgun")

    col_16s = choose_taxon_col(df_16s, taxon_candidates)
    col_shotgun = choose_taxon_col(df_shotgun, taxon_candidates)

    # Keep taxa + effect, rename taxon to a common name for the merge
    df_16s_small = (
        df_16s[[col_16s, "effect_16s"]]
        .dropna(subset=["effect_16s"])
        .rename(columns={col_16s: "taxon"})
    )
    df_shotgun_small = (
        df_shotgun[[col_shotgun, "effect_shotgun"]]
        .dropna(subset=["effect_shotgun"])
        .rename(columns={col_shotgun: "taxon"})
    )

    merged = pd.merge(df_16s_small, df_shotgun_small, on="taxon", how="inner")
    merged = merged.dropna(subset=["effect_16s", "effect_shotgun"])
    n_overlap = merged.shape[0]

    # Give back a nice taxon column name in the output
    out_taxon_name = col_16s if col_16s is not None else col_shotgun
    merged = merged.rename(columns={"taxon": out_taxon_name})

    if n_overlap < min_n:
        return np.nan, np.nan, merged, n_overlap

    r, p = pearsonr(merged["effect_16s"], merged["effect_shotgun"])
    return r, p, merged, n_overlap


# ── 4) Genus level ───────────────────────────────────────────────────────────
# Change file names as needed
genus_16s_path = "sig_taxa_combined_genus_active_16S.csv"
genus_shotgun_path = "sig_taxa_combined_genus_active_shotgun.csv"

genus_16s = pd.read_csv(genus_16s_path)
genus_shotgun = pd.read_csv(genus_shotgun_path)

r_genus, p_genus, genus_merged, n_genus = pearson_agreement_16s_vs_shotgun(
    genus_16s,
    genus_shotgun,
    taxon_candidates=("genus", "Genus")
)

if np.isnan(r_genus):
    print(f"Genus-level: not enough overlapping taxa for Pearson r (n_shared={n_genus})")
else:
    print("Genus-level Pearson r:", r_genus, "p-value:", p_genus, f"(n_shared={n_genus})")

genus_merged.to_csv("genus_16s_shotgun_aligned_for_pearson.csv", index=False)


# ── 5) Species level ──────────────────────────────────────────────────────────
# Change file names as needed
species_16s_path = "sig_taxa_combined_species_active_16S.csv"
species_shotgun_path = "sig_taxa_combined_species_active_shotgun.csv"

species_16s = pd.read_csv(species_16s_path)
species_shotgun = pd.read_csv(species_shotgun_path)

r_species, p_species, species_merged, n_species = pearson_agreement_16s_vs_shotgun(
    species_16s,
    species_shotgun,
    taxon_candidates=("species", "Species")
)

if np.isnan(r_species):
    print(f"Species-level: not enough overlapping taxa for Pearson r (n_shared={n_species})")
else:
    print("Species-level Pearson r:", r_species, "p-value:", p_species, f"(n_shared={n_species})")

species_merged.to_csv("species_16s_shotgun_aligned_for_pearson.csv", index=False)


# ── 6) Summary file (incl. shared counts) ────────────────────────────────────
summary = pd.DataFrame({
    "level": ["genus", "species"],
    "n_shared_taxa": [n_genus, n_species],
    "pearson_r": [r_genus, r_species],
    "p_value": [p_genus, p_species],
})

summary.to_csv("pearson_agreement_summary.csv", index=False)
print("Summary saved to pearson_agreement_summary.csv")
