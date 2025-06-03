#!/bin/bash

# Quality trimming and adapter removal with fastp
# Input: paired-end FASTQ files
# Output: trimmed FASTQ files and QC reports

# Usage: ./fastp_trim.sh SAMPLE_ID

SAMPLE="$1"

# Paths
IN_DIR="results/fastqscreen"
TRIM_DIR="results/trimmed_data"
FASTQC_DIR="results/fastqc_trimmed"
ADAPTERS="short_reads_adapters.fa"

# Create output directories if they don't exist
mkdir -p "$TRIM_DIR" "$FASTQC_DIR"

# Run fastp
fastp \
  --in1 "${IN_DIR}/${SAMPLE}_1.tagged_filter.fastq.gz" \
  --in2 "${IN_DIR}/${SAMPLE}_2.tagged_filter.fastq.gz" \
  --out1 "${TRIM_DIR}/${SAMPLE}_1_trimmed.fq.gz" \
  --out2 "${TRIM_DIR}/${SAMPLE}_2_trimmed.fq.gz" \
  --length_required 50 \
  --correction \
  --adapter_fasta "$ADAPTERS" \
  --cut_right \
  --cut_right_window_size 4 \
  --cut_right_mean_quality 15 \
  --json "${TRIM_DIR}/${SAMPLE}_trimmed.json" \
  --html "${TRIM_DIR}/${SAMPLE}_trimmed.html" \
  -w 8

# Run FastQC on trimmed reads
fastqc --outdir "$FASTQC_DIR" --threads 2 \
  "${TRIM_DIR}/${SAMPLE}_1_trimmed.fq.gz" \
  "${TRIM_DIR}/${SAMPLE}_2_trimmed.fq.gz"