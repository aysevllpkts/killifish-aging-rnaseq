
# Killifish RNA-seq Explorer

This Shiny app lets you explore differential gene expression across tissues, age, and sex in *Nothobranchius furzeri* using interactive plots and tables.

## Folder Structure

```
killifish_rnaseq_app/
├── app.R
└── data/
    ├── normalized_counts.rds
    ├── colData.rds
    ├── res_age.rds
    ├── res_sex.rds
    ├── res_age_brain.rds
    ├── res_sex_brain.rds
    ├── res_age_heart.rds
    ├── res_sex_heart.rds
    └── ... (other tissues)
```

## Features

- Select multiple genes and visualize normalized expression by age and sex
- See differential expression summary across all tissues and full dataset
- Log2 fold changes colored by direction (red = up, blue = down)
- Significant genes (padj < 0.05) are marked with *
- Download plots and DE results as PDF or CSV

## How to Run

### 1. Install R and RStudio (if not already)

- R: https://cran.r-project.org
- RStudio: https://posit.co/download/rstudio-desktop

### 2. Install required R packages

Open R or RStudio and run:

```r
install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages(c("shiny", "ggplot2", "dplyr", "reshape2", "bslib", "DT", "tibble", "RColorBrewer"))
```

### 3. Run the app

Open RStudio, set your working directory to the folder containing `app.R`, and run:

```r
shiny::runApp()
```

Or:

```r
shiny::runApp("path/to/killifish_rnaseq_app")
```

---

## Notes

- If you're running on a new system, make sure your `data/` folder contains all `.rds` files for the app to load correctly.
- DESeq2 requires a working installation of Bioconductor.
- If you face issues with fonts or themes, try restarting RStudio or updating `bslib`.

---

Created by Aysevil Pektas
