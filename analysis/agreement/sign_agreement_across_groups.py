"""
sign_agreement_across_groups.py

Sign agreement across groups at genus or species level
Ran for all 16S or shotgun groups
Uses sig_taxa_combined_genus_*.csv files as input
Please ensure the sig_taxa_combined_genus_*.csv files are in the specified DATA_DIR
and are named correctly- e.g sig_taxa_combined_genus_active_16S.csv
"""

# Import libraries
import pandas as pd
import numpy as np
from pathlib import Path
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns


# ── 1) Parameters ────────────────────────────────────────────────────────────
DATA_DIR = Path("C:/PhD/Meta_analysis/FINAl/Agreement/across group")

# Choose taxonomic level: "species" or "genus"
LEVEL = "species"   # change to "species" if needed

if LEVEL == "genus":
    pattern = "sig_taxa_combined_genus_*.csv"
    id_col = "Genus"
else:
    pattern = "sig_taxa_combined_species_*.csv"
    id_col = "Species"


# ── 2) Load all files ─────────────────────────────────────────────────────────
files = sorted(DATA_DIR.glob(pattern))
if not files:
    raise FileNotFoundError(f"No files found for {pattern}")

def get_col(df, name):
    for c in df.columns:
        if c.lower() == name.lower():
            return c
    return None

def group_name(p: Path) -> str:
    return p.stem.replace(f"sig_taxa_combined_{LEVEL}_", "").replace(" - Copy", "").strip("_")

series_by_group = {}
for p in files:
    df = pd.read_csv(p)
    tax_col = get_col(df, id_col)
    lfc_col = get_col(df, "lfc_pooled")
    est_col = get_col(df, "estimate_meta")
    if tax_col is None:
        raise KeyError(f"{id_col} column not found in {p.name}")

    if lfc_col is not None and est_col is not None:
        effect = df[lfc_col].combine_first(df[est_col])
    elif lfc_col is not None:
        effect = df[lfc_col]
    else:
        effect = df[est_col]

    s = pd.Series(effect.values, index=df[tax_col].astype(str)).dropna()
    s = s.groupby(level=0).mean()
    series_by_group[group_name(p)] = s

if len(series_by_group) < 2:
    raise RuntimeError("Need at least two datasets for comparison")


# ── 3) Compute sign agreement ─────────────────────────────────────────────────
records = []
groups = list(series_by_group.keys())

for a, b in combinations(groups, 2):
    joined = pd.concat([series_by_group[a], series_by_group[b]], axis=1, join="inner").dropna()
    n = len(joined)
    if n == 0:
        continue
    sign_same = np.mean(np.sign(joined.iloc[:, 0]) == np.sign(joined.iloc[:, 1]))
    records.append(dict(pair=f"{a} vs {b}", a=a, b=b, n=n, sign_agreement=sign_same))


# ── 4) Save summary ───────────────────────────────────────────────────────────
df_sign = pd.DataFrame(records)
out_csv = DATA_DIR / f"sign_agreement_{LEVEL}.csv"
df_sign.to_csv(out_csv, index=False)
print(f"Saved: {out_csv}")


# ── 5) Make symmetric matrix for heatmap ───────────────────────────────────────
matrix = pd.DataFrame(index=groups, columns=groups, dtype=float)
for _, row in df_sign.iterrows():
    matrix.loc[row["a"], row["b"]] = row["sign_agreement"]
    matrix.loc[row["b"], row["a"]] = row["sign_agreement"]

# Fill diagonal with 1 (perfect self-agreement)
np.fill_diagonal(matrix.values, 1.0)


# ── 6) Plot heatmap ───────────────────────────────────────────────────────────
plt.figure(figsize=(6,5))
sns.heatmap(matrix.astype(float), annot=True, cmap="RdYlGn", vmin=0, vmax=1)
plt.title(f"Pairwise sign agreement across groups ({LEVEL.capitalize()} level)")
plt.tight_layout()

out_png = DATA_DIR / f"sign_agreement_{LEVEL}.png"
plt.savefig(out_png, dpi=300)
print(f"Saved: {out_png}")
