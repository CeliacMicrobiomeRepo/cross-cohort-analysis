"""
pearson_agreement_across_groups_species.py

Pearson agreement across groups at species level
Ran for all 16S or shotgun groups
Uses sig_taxa_combined_species_*.csv files as input
Please ensure the sig_taxa_combined_species_*.csv files are in the specified DATA_DIR
and are named correctly- e.g sig_taxa_combined_species_active_16S.csv
"""

# Import libraries
import pandas as pd
import numpy as np
from pathlib import Path
from itertools import combinations
from scipy.stats import pearsonr
import matplotlib.pyplot as plt


# ── 1) Set working directory ────────────────────────────────────────────
DATA_DIR = Path("C:/PhD/Meta_analysis/FINAl/Agreement/across group")


# ── 2) Load species-level CSVs ────────────────────────────────────────────
files = sorted(DATA_DIR.glob("sig_taxa_combined_species_*.csv"))
if not files:
    raise FileNotFoundError("No files matching 'sig_taxa_combined_species_*.csv' found.")

def group_name(p: Path) -> str:
    # Make short labels
    name = p.stem.replace("sig_taxa_combined_species_", "").replace(" - Copy", "").strip("_")
    return name

def get_col(df, name):
    for c in df.columns:
        if c.lower() == name.lower():
            return c
    return None


# ── 3) Load each file as a Series (species: Effect) ───────────────────────
series_by_group = {}
for p in files:
    df = pd.read_csv(p)

    species_col = get_col(df, "species")
    lfc_col = get_col(df, "lfc_pooled")
    est_col = get_col(df, "estimate_meta")

    if species_col is None:
        raise KeyError(f"'species' column not found in {p.name}")
    if lfc_col is None and est_col is None:
        raise KeyError(f"Neither 'lfc_pooled' nor 'estimate_meta' found in {p.name}")

    # Prefer lfc_pooled, fallback to estimate_meta
    if lfc_col is not None and est_col is not None:
        effect = df[lfc_col].combine_first(df[est_col])
    elif lfc_col is not None:
        effect = df[lfc_col]
    else:
        effect = df[est_col]

    s = pd.Series(effect.values, index=df[species_col].astype(str)).dropna()
    s = s.groupby(level=0).mean()  # average duplicates
    series_by_group[group_name(p)] = s

if len(series_by_group) < 2:
    raise RuntimeError("Need at least two species files to compute pairwise correlations.")


# ── 4) Pairwise Pearson correlations ────────────────────────────────────────
records = []
global_max_abs = 0.0

for a, b in combinations(series_by_group.keys(), 2):
    joined = pd.concat([series_by_group[a], series_by_group[b]], axis=1, join="inner")
    joined.columns = [a, b]
    joined = joined.dropna()
    n = len(joined)
    if n < 2:
        continue

    x = joined[a].to_numpy()
    y = joined[b].to_numpy()
    r, pval = pearsonr(x, y)
    m, c = np.polyfit(x, y, 1)

    records.append(dict(pair=f"{a} vs {b}", a=a, b=b, n=n, r=r, p=pval, slope=m, intercept=c))
    global_max_abs = max(global_max_abs, np.max(np.abs(np.concatenate([x, y]))))

if not records:
    raise RuntimeError("No overlapping species found between any pair of groups.")


# ── 5) Save correlation summary ─────────────────────────────────────────────
summary = pd.DataFrame.from_records(records).sort_values("r", ascending=False)
summary_path = DATA_DIR / "pairwise_agreement_species_summary.csv"
summary.to_csv(summary_path, index=False)
print(f"Saved: {summary_path}")


# ── 6) Plot agreement lines ─────────────────────────────────────────────────
xmax = max(global_max_abs, 1e-6)
xs = np.linspace(-xmax, xmax, 300)

plt.figure(figsize=(8, 8))
for rec in records:
    ys = rec["slope"] * xs + rec["intercept"]
    plt.plot(xs, ys, label=f"{rec['a']} vs {rec['b']} (r={rec['r']:.2f}, n={rec['n']})")

plt.plot(xs, xs, "--", linewidth=1, color="black")
plt.title("Pairwise Pearson agreement across groups (species level)")
plt.xlabel("Effect (group A scale)")
plt.ylabel("Effect (group B scale)")
plt.axhline(0, linewidth=0.5)
plt.axvline(0, linewidth=0.5)
plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
plt.tight_layout()

plot_path = DATA_DIR / "pairwise_agreement_species_lines.png"
plt.savefig(plot_path, dpi=300)
print(f"Saved: {plot_path}")
