#!/bin/bash

# ---------------------------
# FastQ Screen + FastQC Script
# ---------------------------

# INPUT:  Paired-end FASTQ (.fq.gz)
# TOOL:   fastq_screen with bowtie2 aligner
# CONFIG: User-defined fastqscreen configuration file
# OUTPUT: Tagged FASTQs + FastQC reports

# Usage: ./fastqscreen.sh SAMPLE_ID

# Set paths
IN_DIR="data/rawdata"
OUT_DIR="results/fastqscreen"
FASTQC_OUT="results/fastqc_screened"
CONFIG="fastqscreen.conf"

# Sample name passed as argument
SAMPLE="$1"

# Temporary and output FASTQ paths
TEMP_FASTQ_F="${OUT_DIR}/${SAMPLE}_1.fastq"
TEMP_FASTQ_R="${OUT_DIR}/${SAMPLE}_2.fastq"
FILTERED_FASTQ_F="${OUT_DIR}/${SAMPLE}_1.tagged_filter.fastq"
FILTERED_FASTQ_R="${OUT_DIR}/${SAMPLE}_2.tagged_filter.fastq"

# Create output directories
mkdir -p "$OUT_DIR" "$FASTQC_OUT"

# Check input files
if [[ ! -f "${IN_DIR}/${SAMPLE}_1.fq.gz" || ! -f "${IN_DIR}/${SAMPLE}_2.fq.gz" ]]; then
  echo "ERROR: Input files not found for $SAMPLE"
  exit 1
fi

# Unzip raw FASTQs temporarily
gunzip -c "${IN_DIR}/${SAMPLE}_1.fq.gz" > "$TEMP_FASTQ_F"
gunzip -c "${IN_DIR}/${SAMPLE}_2.fq.gz" > "$TEMP_FASTQ_R"

# Run FastQ Screen
fastq_screen --conf "$CONFIG" --aligner bowtie2 --nohits "$TEMP_FASTQ_F" --threads 8 --outdir "$OUT_DIR" 2> "$OUT_DIR/${SAMPLE}_1_fastqscreen.log"
fastq_screen --conf "$CONFIG" --aligner bowtie2 --nohits "$TEMP_FASTQ_R" --threads 8 --outdir "$OUT_DIR" 2> "$OUT_DIR/${SAMPLE}_2_fastqscreen.log"

# Gzip filtered output
[ -f "$FILTERED_FASTQ_F" ] && gzip "$FILTERED_FASTQ_F"
[ -f "$FILTERED_FASTQ_R" ] && gzip "$FILTERED_FASTQ_R"

# Run FastQC if gzipped filtered files exist
if [[ -f "${FILTERED_FASTQ_F}.gz" && -f "${FILTERED_FASTQ_R}.gz" ]]; then
  fastqc --outdir "$FASTQC_OUT" --threads 8 "${FILTERED_FASTQ_F}.gz" "${FILTERED_FASTQ_R}.gz"
else
  echo "WARNING: Filtered gzipped files not found for FastQC"
fi

# Clean up temporary uncompressed FASTQs
rm -f "$TEMP_FASTQ_F" "$TEMP_FASTQ_R"

echo " FastQ Screen completed for $SAMPLE"