# Killifish Aging RNA-seq

This repository contains code and resources for analyzing RNA-seq data from aging brain, heart, muscle, and spleen tissues of male and female African turquoise killifish (*Nothobranchius furzeri*) between two age groups (6 weeks vs. 16 weeks). The project includes quality control, trimming, alignment, and differential gene expression analysis using DESeq2. Results are visualized through an interactive Shiny app designed to explore gene expression changes across tissues, sexes, and age groups.

## Folder Structure

- `scripts/` – Shell and R scripts for trimming, alignment, quantification, and DE analysis
- `shiny_app/` – Shiny application files
- `data/` – Only metadata and count files

## Dataset

Raw RNA-seq data is from:  
**Xu et al. 2023, Sci Data**  
DOI: [10.1038/s41597-023-02609-x](https://doi.org/10.1038/s41597-023-02609-x)  
SRA accession: [SRP430823](https://www.ncbi.nlm.nih.gov/sra/SRP430823)

## Quick Start

Clone the repository and explore the Shiny app:

```bash
git clone https://github.com/aysevllpkts/killifish-aging-rnaseq.git
cd killifish-aging-rnaseq/shiny_app
R -e "shiny::runApp()"


