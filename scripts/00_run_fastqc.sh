#!/bin/bash

# Run FastQC on trimmed paired-end FASTQ files
# Usage: ./fastqc_trimmed.sh SAMPLE_ID

SAMPLE="$1"

# Paths
OUT_DIR=""

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Run FastQC
fastqc --outdir "$OUT_DIR" --threads 2 ${SAMPLE}