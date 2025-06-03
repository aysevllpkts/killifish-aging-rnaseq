#!/bin/bash

# ----------------------------------------
# FeatureCounts: Gene-level read counting
# Genome: Nothobranchius furzeri (GCF_043380555.1_NfurGRZ-RIMD1)
# Data: Paired-end RNA-seq
# Output: Gene-level counts
# ----------------------------------------


# Define paths
BAM_DIR="results/STAR_alignment_clean/sortedByReadName"
OUT_DIR="results/featureCounts"
GTF_FILE="data/genome/GCF_043380555.1_NfurGRZ-RIMD1_genomic.gtf"

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Run featureCounts
echo "Running featureCounts..."
featureCounts -T 8 \
  -p --countReadPairs \
  -t exon -g gene_id \
  -a "$GTF_FILE" \
  -o "$OUT_FILE" \
  "$BAM_DIR"/*_Aligned.sortedByReadName.out.bam