#!/bin/bash

# Align trimmed paired-end reads to the N. furzeri genome using STAR

# Usage: ./star_align.sh SAMPLE_ID

# Input/output setup
SAMPLE="$1"
IN_DIR="results/trimmed_data"
OUT_DIR="results/STAR_alignment"
GENOME_DIR="data/genome/STAR_index"

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Run STAR
STAR --genomeDir "$GENOME_DIR" \
     --readFilesIn "${IN_DIR}/${SAMPLE}_1_trimmed.fq.gz" "${IN_DIR}/${SAMPLE}_2_trimmed.fq.gz" \
     --readFilesCommand zcat \
     --outFileNamePrefix "${OUT_DIR}/${SAMPLE}_" \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within \
     --outSAMattributes Standard \
     --runThreadN 4