#!/bin/bash
# ============================================================
# REFERENCE-BASED TRANSCRIPTOMICS PIPELINE
# STEP 1+2: QC → ALIGNMENT TO REFERENCE GENOME
# Tools: FastQC, Trimmomatic, STAR (or HISAT2)
# ============================================================
# PURPOSE:
#   When a reference genome is available, reads are aligned
#   directly to the genome (STAR) rather than assembled de novo.
#   STAR is splice-aware, meaning it handles reads that span
#   exon-exon junctions — critical for eukaryotic RNA-seq.
#   HISAT2 is an alternative that uses less RAM.
# ============================================================

set -euo pipefail

# ---------- USER CONFIGURATION ----------
RAW_DIR="data/raw"
TRIMMED_DIR="data/trimmed"
ALIGNED_DIR="results/aligned"
QC_DIR="results/qc"
GENOME_DIR="genome"
GENOME_FASTA="$GENOME_DIR/genome.fa"
GENOME_GTF="$GENOME_DIR/annotation.gtf"
STAR_INDEX="$GENOME_DIR/star_index"
THREADS=16

mkdir -p "$TRIMMED_DIR" "$ALIGNED_DIR" "$QC_DIR" "$STAR_INDEX" logs

echo "=============================================="
echo " Reference-Based Transcriptomics: Alignment"
echo "=============================================="
echo "Start: $(date)"

# ---------- STEP 1: FastQC + Trimmomatic ----------
echo "[1/4] Quality control and trimming..."

fastqc --threads "$THREADS" --outdir "$QC_DIR/raw" "$RAW_DIR"/*.fastq.gz \
    2> logs/fastqc_raw.log

for R1 in "$RAW_DIR"/*_R1_001.fastq.gz; do
    SAMPLE=$(basename "$R1" _R1_001.fastq.gz)
    R2="$RAW_DIR/${SAMPLE}_R2_001.fastq.gz"

    trimmomatic PE -threads "$THREADS" -phred33 \
        "$R1" "$R2" \
        "$TRIMMED_DIR/${SAMPLE}_R1_paired.fastq.gz"   \
        "$TRIMMED_DIR/${SAMPLE}_R1_unpaired.fastq.gz" \
        "$TRIMMED_DIR/${SAMPLE}_R2_paired.fastq.gz"   \
        "$TRIMMED_DIR/${SAMPLE}_R2_unpaired.fastq.gz" \
        ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
        2> logs/trim_"$SAMPLE".log
    echo "  Trimmed: $SAMPLE"
done

multiqc "$QC_DIR/raw" --outdir "$QC_DIR/raw/multiqc" 2>/dev/null

# ---------- STEP 2: Build STAR Genome Index ----------
# STAR requires a genome index built from the reference FASTA and GTF.
# This step only needs to be run ONCE per genome/annotation combination.
# --sjdbOverhang: read length - 1 (for 150 bp reads, use 149)

echo "[2/4] Building STAR genome index (run once)..."

if [ ! -f "$STAR_INDEX/Genome" ]; then
    STAR \
        --runMode genomeGenerate \
        --genomeDir "$STAR_INDEX" \
        --genomeFastaFiles "$GENOME_FASTA" \
        --sjdbGTFfile "$GENOME_GTF" \
        --sjdbOverhang 149 \
        --runThreadN "$THREADS" \
        --genomeSAindexNbases 11 \
        2> logs/star_index.log
    echo "  -> STAR index built."
else
    echo "  -> STAR index already exists, skipping."
fi

# ---------- STEP 3: Align each sample with STAR ----------
# Key flags:
#   --outSAMtype BAM SortedByCoordinate : output sorted BAM directly
#   --quantMode GeneCounts              : also output raw gene counts
#   --outFilterMismatchNmax 2           : strict mismatch filter
#   --alignSJDBoverhangMin 1            : allow 1 nt overhang at splice sites
#   --twopassMode Basic                 : 2-pass mapping for novel splice sites

echo "[3/4] Aligning samples with STAR..."

for R1 in "$TRIMMED_DIR"/*_R1_paired.fastq.gz; do
    SAMPLE=$(basename "$R1" _R1_paired.fastq.gz)
    R2="$TRIMMED_DIR/${SAMPLE}_R2_paired.fastq.gz"

    mkdir -p "$ALIGNED_DIR/$SAMPLE"
    echo "  Aligning: $SAMPLE"

    STAR \
        --runMode alignReads \
        --genomeDir "$STAR_INDEX" \
        --readFilesIn "$R1" "$R2" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS NM MD \
        --outFileNamePrefix "$ALIGNED_DIR/$SAMPLE/" \
        --runThreadN "$THREADS" \
        --quantMode GeneCounts \
        --outFilterMismatchNmax 2 \
        --alignSJDBoverhangMin 1 \
        --outSJtype Standard \
        --twopassMode Basic \
        2> logs/star_"$SAMPLE".log

    # Index BAM for downstream tools
    samtools index "$ALIGNED_DIR/$SAMPLE/Aligned.sortedByCoord.out.bam"
    echo "    -> $SAMPLE aligned and indexed."
done

# ---------- STEP 4: Alignment QC with RSeQC / MultiQC ----------
echo "[4/4] Running alignment QC..."

multiqc \
    "$ALIGNED_DIR" \
    logs/ \
    --outdir "$QC_DIR/aligned/multiqc" \
    --filename "multiqc_alignment_report" \
    2> logs/multiqc_aligned.log

echo ""
echo "=============================================="
echo " ALIGNMENT COMPLETE"
echo " BAM files: $ALIGNED_DIR/<sample>/Aligned.sortedByCoord.out.bam"
echo " Gene counts: $ALIGNED_DIR/<sample>/ReadsPerGene.out.tab"
echo " QC report: $QC_DIR/aligned/multiqc/"
echo "=============================================="
echo "End: $(date)"
