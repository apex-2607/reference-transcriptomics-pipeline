# 🧬 Reference-Based Transcriptomics Analysis Pipeline

> **Complete RNA-seq workflow using a reference genome** — from raw reads to differential expression, pathway analysis, and publication-ready figures.

---

## 📌 Why Reference-Based Transcriptomics?

When a **high-quality reference genome and annotation** are available for your organism, aligning reads to the reference is far more accurate and interpretable than de novo assembly. Reference-based approaches offer:

- **Higher sensitivity** — known transcripts are detected even at low expression
- **Splice-aware alignment** — reads crossing exon junctions are correctly mapped
- **Standardized gene IDs** — results are directly compatible with GO/KEGG databases
- **Isoform-level analysis** — can detect alternative splicing events
- **Faster analysis** — no computationally expensive assembly step

**When to choose reference-based over de novo:**
| Feature | Reference-Based | De Novo |
|---------|----------------|---------|
| Genome available? | ✅ Yes | ❌ No |
| Speed | Faster | Slow |
| Novel transcript discovery | Limited | ✅ Yes |
| Annotation dependency | High | None |
| Recommended for | Model organisms | Non-model, fungi, environmental |

---

## 🔬 How the Pipeline Works

```
Raw FASTQ reads
     │
     ▼
[1] QC + Trimming (FastQC → Trimmomatic → MultiQC)
     │  Remove adapters, low-quality bases; assess read quality
     ▼
[2] STAR Alignment (Splice-Aware)
     │  Map reads to reference genome using STAR 2-pass mode
     │  Output: sorted BAM files with splice junction information
     ▼
[3] Read Counting (featureCounts)
     │  Count fragments mapping to each annotated gene in GTF
     │  Output: raw count matrix (genes × samples)
     ▼
[4] Differential Expression (DESeq2 in R)
     │  Negative binomial GLM, LFC shrinkage, FDR correction
     │  Optional: batch correction via design formula
     ▼
[5] Pathway Analysis (clusterProfiler in R)
     │  GO ORA / KEGG enrichment for up- and down-regulated genes
     ▼
[6] Python Analysis (PyDESeq2 / visualization)
     └  Alternative DEA, correlation analysis, variable gene heatmaps
```

---

## 📁 Repository Structure

```
repo_reference_transcriptomics/
├── linux/
│   ├── 01_qc_and_alignment.sh    # FastQC, Trimmomatic, STAR 2-pass
│   └── 02_read_counting.sh       # featureCounts, count matrix cleanup
├── rstudio/
│   └── 01_deseq2_full_pipeline.R # DESeq2, LFC shrinkage, all plots, GO/KEGG
├── python/
│   └── 01_analysis.py            # PyDESeq2, CPM normalization, heatmaps
└── README.md
```

---

## 🛠️ Software Requirements

### Linux
| Tool | Purpose |
|------|---------|
| FastQC ≥0.12 | Read quality assessment |
| Trimmomatic ≥0.39 | Adapter trimming |
| MultiQC ≥1.14 | Aggregate QC reports |
| STAR ≥2.7.10 | Splice-aware alignment |
| Samtools ≥1.17 | BAM indexing |
| Subread/featureCounts ≥2.0 | Gene-level counting |

### R
```r
BiocManager::install(c("DESeq2", "apeglm", "clusterProfiler",
                       "enrichplot", "ComplexHeatmap",
                       "org.Hs.eg.db"))     # or appropriate OrgDb
install.packages(c("ggplot2", "dplyr", "readr", "pheatmap",
                   "RColorBrewer", "ggrepel", "circlize", "tibble"))
```

### Python
```bash
pip install pydeseq2 pandas numpy matplotlib seaborn scipy scikit-learn
```

---

## 📂 Data Directory Layout

```
data/raw/                ← Raw FASTQ (*_R1_001.fastq.gz)
genome/
├── genome.fa            ← Reference genome FASTA
├── annotation.gtf       ← Gene annotation (Ensembl/NCBI)
└── star_index/          ← Auto-built by 01_qc_and_alignment.sh
results/
├── qc/
├── aligned/             ← BAM files per sample
├── counts/              ← gene_count_matrix.tsv
├── DEG/                 ← All DEA outputs, plots
└── python_analysis/
```

---

## 🚀 Quick Start

```bash
# 1. Obtain reference genome and annotation
# Example: Fungal genome from NCBI/Ensembl Fungi
wget https://ftp.ncbi.nlm.nih.gov/genomes/.../genome.fa.gz -O genome/genome.fa.gz
gunzip genome/genome.fa.gz

# 2. Run Linux pipeline
bash linux/01_qc_and_alignment.sh
bash linux/02_read_counting.sh

# 3. Open RStudio → run rstudio/01_deseq2_full_pipeline.R

# 4. Run Python analysis
python python/01_analysis.py
```

---

## ⚠️ Key Parameters to Check

| Script | Parameter | Notes |
|--------|-----------|-------|
| `01_qc_and_alignment.sh` | `--sjdbOverhang` | Set to read_length - 1 |
| `01_qc_and_alignment.sh` | `--SS_lib_type` | Strand-specific: RF (dUTP), FR, or remove |
| `02_read_counting.sh` | `-s 2` | Strandedness: 0=unstranded, 1=forward, 2=reverse |
| `01_deseq2_full_pipeline.R` | `design` | Add `+ batch` if batch effects present |
| `01_deseq2_full_pipeline.R` | `organism =` | Change KEGG organism code: hsa, mmu, sce |

---

## 📊 Expected Key Outputs

| File | Description |
|------|-------------|
| `results/counts/gene_count_matrix.tsv` | Raw integer count matrix |
| `results/DEG/DESeq2_all_genes.csv` | Full results with shrunk LFC |
| `results/DEG/significant_DEGs.csv` | FDR<0.05, |LFC|>1 |
| `results/DEG/volcano.pdf` | Volcano plot |
| `results/DEG/DEG_heatmap.pdf` | ComplexHeatmap of top DEGs |
| `results/DEG/GO_enrichment_up.csv` | Enriched GO terms |
| `results/DEG/KEGG_dotplot.pdf` | KEGG pathway dotplot |

---

## 📚 References

1. Dobin A et al. (2013) *STAR: ultrafast universal RNA-seq aligner.* **Bioinformatics.**
2. Liao Y et al. (2014) *featureCounts: an efficient general purpose program for assigning sequence reads to genomic features.* **Bioinformatics.**
3. Love MI et al. (2014) *DESeq2.* **Genome Biology.**
4. Zhu A et al. (2019) *Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences.* **Bioinformatics.** (apeglm LFC shrinkage)

---

## 👤 Author

**Apex** | M.Tech Biotechnology | MS University of Baroda | DBT-funded
