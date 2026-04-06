#!/bin/bash
# ============================================================
# STEP 2: READ COUNTING WITH featureCounts
# Tool: Subread/featureCounts
# ============================================================
# PURPOSE:
#   After alignment, we count how many reads map to each gene
#   in the GTF annotation. featureCounts is faster and more
#   memory-efficient than HTSeq for large datasets.
#
# COUNTING MODES:
#   -p         : Paired-end reads (count fragments, not individual reads)
#   -s 2       : Reverse-stranded library (dUTP method); use -s 1 for
#                forward-stranded, -s 0 for unstranded
#   -t exon    : Count reads overlapping "exon" features in GTF
#   -g gene_id : Aggregate exon counts to gene level using "gene_id"
#   --fracOverlap 0.2 : Require 20% overlap to count a read to a feature
# ============================================================

set -euo pipefail

ALIGNED_DIR="results/aligned"
COUNTS_DIR="results/counts"
GENOME_GTF="genome/annotation.gtf"
THREADS=8

mkdir -p "$COUNTS_DIR" logs

echo "=============================================="
echo " Reference Transcriptomics: featureCounts"
echo "=============================================="

# Collect all BAM files
BAM_FILES=$(find "$ALIGNED_DIR" -name "Aligned.sortedByCoord.out.bam" | sort | tr '\n' ' ')
echo "BAM files found: $(echo $BAM_FILES | wc -w)"

# ---------- Run featureCounts ----------
featureCounts \
    -T "$THREADS" \
    -p \
    --countReadPairs \
    -s 2 \
    -t exon \
    -g gene_id \
    -a "$GENOME_GTF" \
    -o "$COUNTS_DIR/raw_counts.txt" \
    --fracOverlap 0.2 \
    --minOverlap 10 \
    -M \
    --fraction \
    $BAM_FILES \
    2> logs/featurecounts.log

echo "  -> featureCounts complete."

# ---------- Clean the output to a simple count matrix ----------
# featureCounts output has 6 metadata columns (Chr, Start, End, Strand, Length, ...)
# We extract only the gene ID column + count columns

echo "Cleaning count matrix..."

awk 'NR>1 {
    printf $1
    for (i=7; i<=NF; i++) printf "\t"$i
    printf "\n"
}' "$COUNTS_DIR/raw_counts.txt" > "$COUNTS_DIR/count_matrix.txt"

# Fix column headers (strip full BAM paths, keep sample names only)
head -2 "$COUNTS_DIR/raw_counts.txt" | tail -1 | \
    awk '{
        printf "GeneID"
        for (i=7; i<=NF; i++) {
            split($i, a, "/")
            printf "\t"a[length(a)-1]
        }
        printf "\n"
    }' > "$COUNTS_DIR/count_matrix_header.txt"

cat "$COUNTS_DIR/count_matrix_header.txt" \
    "$COUNTS_DIR/count_matrix.txt" > "$COUNTS_DIR/gene_count_matrix.tsv"

rm "$COUNTS_DIR/count_matrix.txt" "$COUNTS_DIR/count_matrix_header.txt"

# ---------- Count matrix summary ----------
NGENES=$(wc -l < "$COUNTS_DIR/gene_count_matrix.tsv")
echo "Total genes in count matrix: $((NGENES - 1))"

# Print assignment stats from log
echo ""
echo "--- Assignment Summary ---"
grep "Successfully assigned" logs/featurecounts.log || true
grep "Unassigned" logs/featurecounts.log || true

echo ""
echo "=============================================="
echo " COUNTING COMPLETE"
echo " Count matrix: $COUNTS_DIR/gene_count_matrix.tsv"
echo "=============================================="
