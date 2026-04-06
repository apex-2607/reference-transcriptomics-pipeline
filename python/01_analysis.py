#!/usr/bin/env python3
"""
============================================================
REFERENCE-BASED TRANSCRIPTOMICS — Python Analysis
Tools: pydeseq2, pandas, seaborn, matplotlib, scipy
============================================================
PURPOSE:
    Alternative to R-based DESeq2 using the PyDESeq2 library,
    plus additional Python-specific analyses:
      1. DESeq2-equivalent analysis with PyDESeq2
      2. TPM/RPKM normalization and visualization
      3. Gene expression heatmaps
      4. Co-expression network analysis
      5. Export for downstream tools
============================================================
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.preprocessing import StandardScaler

# Try to import PyDESeq2; fall back to custom analysis if not installed
try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    PYDESEQ2_AVAILABLE = True
except ImportError:
    print("PyDESeq2 not installed. Run: pip install pydeseq2")
    print("Proceeding with normalization and visualization only.")
    PYDESEQ2_AVAILABLE = False

# ============================================================
# SECTION 1: Load count matrix and metadata
# ============================================================

COUNT_PATH  = "results/counts/gene_count_matrix.tsv"
OUTPUT_DIR  = "results/python_analysis"
os.makedirs(OUTPUT_DIR, exist_ok=True)

counts_df = pd.read_csv(COUNT_PATH, sep="\t", index_col="GeneID")
print(f"Count matrix: {counts_df.shape[0]} genes x {counts_df.shape[1]} samples")

# Sample metadata
metadata = pd.DataFrame({
    "condition": ["Control", "Control", "Control",
                  "Treatment", "Treatment", "Treatment"],
    "batch":     ["A", "A", "B", "A", "A", "B"]
}, index=counts_df.columns)

# Pre-filter: remove genes with fewer than 10 total counts
counts_filt = counts_df[counts_df.sum(axis=1) >= 10]
print(f"After filtering: {counts_filt.shape[0]} genes")

# ============================================================
# SECTION 2: TPM normalization (from gene length)
# ============================================================
# TPM = (count / gene_length_kb) / (sum of all RPK) * 1e6
# We need gene lengths from the GTF or annotation file.
# If not available, we use CPM as an approximation.

def compute_cpm(df: pd.DataFrame) -> pd.DataFrame:
    """Counts Per Million — simple library-size normalization."""
    return df.divide(df.sum(axis=0), axis=1) * 1e6

def compute_tpm(df: pd.DataFrame, gene_lengths: pd.Series) -> pd.DataFrame:
    """Transcripts Per Million — accounts for gene length bias."""
    rpk = df.divide(gene_lengths / 1000, axis=0)
    return rpk.divide(rpk.sum(axis=0), axis=1) * 1e6

cpm = compute_cpm(counts_filt)
log2_cpm = np.log2(cpm + 1)

print("CPM normalization complete.")

# ============================================================
# SECTION 3: PyDESeq2 Differential Expression (if available)
# ============================================================

if PYDESEQ2_AVAILABLE:
    print("\nRunning PyDESeq2...")

    dds = DeseqDataSet(
        counts   = counts_filt.T,          # PyDESeq2 expects samples x genes
        metadata = metadata,
        design_factors = "condition",
        refit_cooks = True,
    )
    dds.deseq2()

    stat_res = DeseqStats(dds, contrast=["condition", "Treatment", "Control"])
    stat_res.summary()
    stat_res.lfc_shrink(coeff="condition_Treatment_vs_Control")

    results_df = stat_res.results_df.reset_index()
    results_df.columns = ["gene_id"] + list(results_df.columns[1:])
    results_df = results_df.sort_values("padj")
    results_df.to_csv(f"{OUTPUT_DIR}/PyDESeq2_results.csv", index=False)

    # Significant DEGs
    sig_df = results_df[(results_df["padj"] < 0.05) &
                         (results_df["log2FoldChange"].abs() > 1)]
    sig_df.to_csv(f"{OUTPUT_DIR}/PyDESeq2_significant.csv", index=False)
    print(f"Significant DEGs: {len(sig_df)}")

    # Volcano
    results_df["group"] = "NS"
    results_df.loc[(results_df["padj"] < 0.05) & (results_df["log2FoldChange"] > 1),
                   "group"] = "Up"
    results_df.loc[(results_df["padj"] < 0.05) & (results_df["log2FoldChange"] < -1),
                   "group"] = "Down"

    palette = {"Up": "#E41A1C", "Down": "#377EB8", "NS": "#CCCCCC"}
    fig, ax = plt.subplots(figsize=(8, 6))
    for grp, color in palette.items():
        sub = results_df[results_df["group"] == grp]
        ax.scatter(sub["log2FoldChange"], -np.log10(sub["padj"] + 1e-300),
                   c=color, alpha=0.5, s=8, label=grp)
    ax.axvline(x=1,  linestyle="--", color="grey", linewidth=0.8)
    ax.axvline(x=-1, linestyle="--", color="grey", linewidth=0.8)
    ax.axhline(y=-np.log10(0.05), linestyle="--", color="grey", linewidth=0.8)
    ax.set_xlabel("log2 Fold Change", fontsize=12)
    ax.set_ylabel("-log10(Adj. p-value)", fontsize=12)
    ax.set_title("Volcano Plot: Treatment vs Control (PyDESeq2)", fontsize=13)
    ax.legend(title="Regulation")
    fig.tight_layout()
    fig.savefig(f"{OUTPUT_DIR}/volcano_pydeseq2.png", dpi=150)
    plt.close()

# ============================================================
# SECTION 4: Correlation Matrix and Hierarchical Clustering
# ============================================================

print("\nComputing sample correlation matrix...")

corr_matrix = log2_cpm.corr(method="pearson")

# Hierarchical clustering on samples
Z = linkage(corr_matrix.values, method="ward")

fig, axes = plt.subplots(1, 2, figsize=(14, 6),
                          gridspec_kw={"width_ratios": [1, 3]})

# Dendrogram
dendrogram(Z, labels=list(corr_matrix.columns),
           orientation="right", ax=axes[0], leaf_font_size=9)
axes[0].set_title("Sample Dendrogram")
axes[0].invert_yaxis()

# Heatmap
im = axes[1].imshow(corr_matrix.values, cmap="RdBu_r", vmin=0.8, vmax=1.0,
                    aspect="auto")
axes[1].set_xticks(range(len(corr_matrix.columns)))
axes[1].set_yticks(range(len(corr_matrix.columns)))
axes[1].set_xticklabels(corr_matrix.columns, rotation=45, ha="right", fontsize=9)
axes[1].set_yticklabels(corr_matrix.columns, fontsize=9)
plt.colorbar(im, ax=axes[1], label="Pearson r", fraction=0.046, pad=0.04)
axes[1].set_title("Sample Correlation Heatmap (log2 CPM)")

fig.tight_layout()
fig.savefig(f"{OUTPUT_DIR}/sample_correlation.png", dpi=150)
plt.close()

print("Correlation analysis done.")

# ============================================================
# SECTION 5: Gene expression heatmap — top variable genes
# ============================================================

# Identify top 100 most variable genes
gene_var   = log2_cpm.var(axis=1).sort_values(ascending=False)
top100     = gene_var.head(100).index
mat_heat   = log2_cpm.loc[top100]

# Z-score per gene
mat_z = mat_heat.apply(lambda row: (row - row.mean()) / row.std(), axis=1)

fig, ax = plt.subplots(figsize=(10, 14))
sns.heatmap(mat_z,
            cmap        = "RdBu_r",
            center      = 0,
            xticklabels = list(log2_cpm.columns),
            yticklabels = False,
            linewidths  = 0,
            vmin        = -2, vmax = 2,
            ax          = ax,
            cbar_kws    = {"label": "Z-score"})

# Annotate conditions
condition_colors = metadata["condition"].map(
    {"Control": "#1B7837", "Treatment": "#762A83"}).values
for i, color in enumerate(condition_colors):
    ax.add_patch(plt.Rectangle((i, -2), 1, 2, color=color,
                                transform=ax.get_xaxis_transform(), clip_on=False))

ax.set_title("Top 100 Variable Genes (Z-scored log2 CPM)", fontsize=13, pad=15)
ax.set_xlabel("Sample")
ax.set_ylabel("Gene")
fig.tight_layout()
fig.savefig(f"{OUTPUT_DIR}/heatmap_top100_variable.png", dpi=150)
plt.close()

# ============================================================
# SECTION 6: Export normalized expression tables
# ============================================================

log2_cpm.to_csv(f"{OUTPUT_DIR}/log2_CPM_normalized.csv")
cpm.to_csv(f"{OUTPUT_DIR}/CPM_normalized.csv")

print(f"\nAll outputs saved to: {OUTPUT_DIR}/")
print("Python reference-based transcriptomics analysis COMPLETE.")
