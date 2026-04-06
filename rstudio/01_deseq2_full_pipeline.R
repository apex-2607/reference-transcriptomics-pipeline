# ============================================================
# REFERENCE-BASED TRANSCRIPTOMICS — DESEQ2 ANALYSIS (R)
# Tools: DESeq2, ggplot2, ComplexHeatmap, clusterProfiler
# ============================================================
# PURPOSE:
#   featureCounts produces a count matrix (genes x samples).
#   DESeq2 normalizes these counts, models biological variability,
#   and identifies statistically significant expression changes
#   between conditions. Critically, DESeq2 works with raw counts
#   (NOT RPKM/TPM) — it does its own internal normalization.
# ============================================================

# ---- Install / Load packages ----
pkgs_bioc <- c("DESeq2", "apeglm", "clusterProfiler",
               "org.Hs.eg.db", "enrichplot", "ComplexHeatmap")
pkgs_cran <- c("ggplot2", "dplyr", "readr", "pheatmap",
               "RColorBrewer", "ggrepel", "circlize", "tibble")

for (p in pkgs_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p)
}
for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(DESeq2); library(apeglm); library(clusterProfiler)
  library(ggplot2); library(dplyr); library(readr)
  library(pheatmap); library(RColorBrewer); library(ggrepel)
  library(tibble); library(ComplexHeatmap); library(circlize)
})

# ============================================================
# SECTION 1: Load count matrix
# ============================================================

count_mat <- read_tsv("results/counts/gene_count_matrix.tsv") %>%
  column_to_rownames("GeneID") %>%
  as.matrix()

# Round to integers (featureCounts with --fraction may produce decimals)
count_mat <- round(count_mat)
mode(count_mat) <- "integer"

message("Count matrix: ", nrow(count_mat), " genes x ", ncol(count_mat), " samples")

# Sample metadata — EDIT THIS to match your experimental design
coldata <- data.frame(
  sample    = colnames(count_mat),
  condition = c("Control", "Control", "Control",
                "Treatment", "Treatment", "Treatment"),
  batch     = c("A", "A", "B", "A", "A", "B"),   # Remove if no batch effect
  row.names = colnames(count_mat)
)

# Verify order matches
stopifnot(all(rownames(coldata) == colnames(count_mat)))

# ============================================================
# SECTION 2: DESeq2 — construct, filter, run
# ============================================================

# Include batch in design formula if batch effects are present.
# If no batch: design = ~ condition
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = coldata,
  design    = ~ batch + condition
)

# Pre-filter: retain genes with at least 10 counts in ≥ 3 samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
message("After pre-filtering: ", nrow(dds), " genes")

# Set reference condition
dds$condition <- relevel(dds$condition, ref = "Control")

# Run DESeq2 (estimates size factors, dispersions, GLM fit)
dds <- DESeq(dds, parallel = FALSE)

# Dispersion plot — check model fit quality
pdf("results/DEG/dispersion_plot.pdf", width = 6, height = 5)
plotDispEsts(dds, main = "DESeq2 Dispersion Estimates")
dev.off()

# ============================================================
# SECTION 3: Extract results with LFC shrinkage
# ============================================================

res <- results(dds,
               contrast  = c("condition", "Treatment", "Control"),
               alpha     = 0.05)

# LFC shrinkage using apeglm (for plotting and ranking)
res_shrunk <- lfcShrink(dds,
                        coef = "condition_Treatment_vs_Control",
                        type = "apeglm")

summary(res_shrunk)

# Save full results
res_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene_id") %>%
  arrange(padj)

dir.create("results/DEG", showWarnings = FALSE)
write_csv(res_df, "results/DEG/DESeq2_all_genes.csv")

# Significant DEGs
sig <- res_df %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)
up  <- sig %>% filter(log2FoldChange > 1)
dn  <- sig %>% filter(log2FoldChange < -1)

write_csv(sig, "results/DEG/significant_DEGs.csv")
write_csv(up,  "results/DEG/upregulated_DEGs.csv")
write_csv(dn,  "results/DEG/downregulated_DEGs.csv")

message("Total DEGs: ", nrow(sig), " (Up: ", nrow(up), ", Down: ", nrow(dn), ")")

# ============================================================
# SECTION 4: Visualization Suite
# ============================================================

# ---- 4A: PCA ----
vsd <- vst(dds, blind = FALSE)
pca_data    <- plotPCA(vsd, intgroup = c("condition", "batch"), returnData = TRUE)
percentVar  <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2,
                               color = condition, shape = batch, label = name)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  scale_color_manual(values = c("Control" = "#1B7837", "Treatment" = "#762A83")) +
  theme_classic(base_size = 13) +
  ggtitle("PCA — Sample Clustering")

ggsave("results/DEG/PCA.pdf", p_pca, width = 7, height = 5)

# ---- 4B: Volcano ----
volcano_df <- res_df %>%
  mutate(
    label = ifelse(rank(padj) <= 20 & abs(log2FoldChange) > 1, gene_id, ""),
    group = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

p_vol <- ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj + 1e-300),
                                 color = group, label = label)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_text_repel(size = 2.5, color = "black", max.overlaps = 30) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c("Up" = "#D7191C", "Down" = "#2C7BB6", "NS" = "#AAAAAA")) +
  labs(title = "Volcano: Treatment vs Control",
       x = "log2 Fold Change", y = "-log10(Adj. p-value)") +
  theme_bw(base_size = 13)

ggsave("results/DEG/volcano.pdf", p_vol, width = 8, height = 6)

# ---- 4C: Sample-Sample Correlation Heatmap ----
sampleDist <- dist(t(assay(vsd)))
sampleDistMat <- as.matrix(sampleDist)

color_pal <- colorRampPalette(rev(brewer.pal(9, "Blues")))(100)
annot_col <- data.frame(Condition = coldata$condition,
                        Batch = coldata$batch,
                        row.names = rownames(coldata))

pdf("results/DEG/sample_correlation_heatmap.pdf", width = 7, height = 6)
pheatmap(sampleDistMat,
         clustering_distance_rows = sampleDist,
         clustering_distance_cols = sampleDist,
         annotation_col = annot_col,
         col = color_pal,
         main = "Sample-to-Sample Distance")
dev.off()

# ---- 4D: DEG Heatmap ----
top_degs <- head(sig$gene_id, 50)
mat_deg  <- assay(vsd)[top_degs, ]
mat_deg  <- mat_deg - rowMeans(mat_deg)

col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#D6604D"))

ha <- HeatmapAnnotation(
  Condition = coldata$condition,
  Batch     = coldata$batch,
  col = list(
    Condition = c("Control" = "#1B7837", "Treatment" = "#762A83"),
    Batch     = c("A" = "#FDB863", "B" = "#B2ABD2")
  )
)

pdf("results/DEG/DEG_heatmap.pdf", width = 8, height = 10)
Heatmap(mat_deg,
        top_annotation    = ha,
        col               = col_fun,
        cluster_rows      = TRUE,
        cluster_columns   = TRUE,
        show_row_names    = TRUE,
        show_column_names = TRUE,
        row_names_gp      = gpar(fontsize = 6),
        column_names_gp   = gpar(fontsize = 9),
        heatmap_legend_param = list(title = "Z-score"),
        column_title      = "Top 50 DEGs")
dev.off()

# ============================================================
# SECTION 5: Pathway Enrichment (clusterProfiler)
# ============================================================
# Using org.Hs.eg.db for human; change to org.Mm.eg.db for mouse,
# or provide a custom TERM2GENE for non-model organisms.

if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  library(org.Hs.eg.db)

  # Convert gene symbols to Entrez IDs
  convert_ids <- function(gene_ids, from = "SYMBOL") {
    bitr(gene_ids, fromType = from, toType = "ENTREZID",
         OrgDb = org.Hs.eg.db, drop = TRUE)
  }

  up_entrez  <- convert_ids(up$gene_id)$ENTREZID
  dn_entrez  <- convert_ids(dn$gene_id)$ENTREZID
  all_entrez <- convert_ids(res_df$gene_id)$ENTREZID

  # GO ORA (upregulated)
  go_up <- enrichGO(gene         = up_entrez,
                    universe     = all_entrez,
                    OrgDb        = org.Hs.eg.db,
                    ont          = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    readable     = TRUE)

  if (!is.null(go_up) && nrow(go_up) > 0) {
    write_csv(as.data.frame(go_up), "results/DEG/GO_enrichment_up.csv")
    p_dot <- dotplot(go_up, split = "ONTOLOGY", showCategory = 10) +
      facet_grid(ONTOLOGY ~ ., scale = "free") +
      ggtitle("GO Enrichment: Upregulated Genes")
    ggsave("results/DEG/GO_dotplot_up.pdf", p_dot, width = 9, height = 10)
  }

  # KEGG pathway enrichment
  kegg_res <- enrichKEGG(gene         = up_entrez,
                          universe     = all_entrez,
                          organism     = "hsa",  # Change: mmu=mouse, eco=ecoli
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05)

  if (!is.null(kegg_res) && nrow(kegg_res) > 0) {
    write_csv(as.data.frame(kegg_res), "results/DEG/KEGG_enrichment.csv")
    p_kegg <- dotplot(kegg_res, showCategory = 20) +
      ggtitle("KEGG Pathway Enrichment")
    ggsave("results/DEG/KEGG_dotplot.pdf", p_kegg, width = 9, height = 8)
  }

} else {
  message("NOTE: org.Hs.eg.db not found. Skipping built-in enrichment.")
  message("For non-model organisms, use custom TERM2GENE from Trinotate output.")
}

message("\nAll results saved to: results/DEG/")
message("Reference-based transcriptomics R analysis COMPLETE.")
