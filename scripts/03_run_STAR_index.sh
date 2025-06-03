#!/bin/bash

# STAR genome index generation for Nothobranchius furzeri (GCF_043380555.1_NfurGRZ-RIMD1)

# Usage: ./star_index.sh

# Paths
GENOME_DIR="data/genome/STAR_index"
FASTA="data/genome/GCF_043380555.1_NfurGRZ-RIMD1_genomic.fasta"
GTF="data/genome/GCF_043380555.1_NfurGRZ-RIMD1_genomic.gtf"

# Create output directory
mkdir -p "$GENOME_DIR"

# Run STAR genomeGenerate
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir "$GENOME_DIR" \
     --genomeFastaFiles "$FASTA" \
     --sjdbGTFfile "$GTF" \
     --genomeSAindexNbases 14