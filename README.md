# Selection Index Analysis App

This repository contains an R Shiny application for an integrated plant-breeding analysis pipeline. The current pipeline incorporates the models and analysis updates maintained in the `Update` folder.

## Features

- Upload Excel files and preview the raw data.
- Shared input validation and pipeline-safety checks.
- **Module 1 - Breeding:** heritability, genetic gain, response per year, realized gain, and generation summaries.
- **Module 2 - Genetic Diversity Analysis:** Mahalanobis D2 distances, clustering, PCA, genotype superiority, and trait correlations.
- **Module 3 - Mating:**
  - Griffing Method I
  - Griffing Method II
  - Griffing Method III
  - Griffing Method IV
  - Partial diallel
  - Line x Tester
- **Module 4 - Selection Index:** single-trait selection and multi-trait LPSI analysis.
- **Module 5 - Multi-Environment Trial:** BLUP, Finlay-Wilkinson stability, AMMI, GGE, quality control, and integrated ranking outputs.
- View interactive result tables in the app.
- Download analysis results as Excel files.
- Download selected charts as image files.

## Project Files

- `app.R` - Main Shiny application, user interface, server orchestration, analysis pipelines, plotting, and export functions.
- `modules/module_1_breeding.R` - Standalone breeding and genetic-gain functions.
- `modules/module_3_mating.R` - Standalone mating-design functions.
- `modules/shared_data_validation.R` - Validation and cross-module pipeline-safety helpers.
- `modules/advanced_analysis_extensions.R` - Advanced extensions used by genetic diversity, selection-index, and MET workflows.
- `Update/` - Supplied update snapshot retained for comparison and future reference.
- `.gitignore` - Files and folders excluded from Git tracking.

## Required R Packages

The app uses these R packages:

```r
install.packages(c(
  "shiny",
  "bslib",
  "readxl",
  "tidyverse",
  "emmeans",
  "multcomp",
  "multcompView",
  "pheatmap",
  "ggplot2",
  "DT",
  "writexl",
  "lme4",
  "lmerTest",
  "patchwork"
))
```

## How to Run

Open the project folder in RStudio or VS Code, then run:

```r
shiny::runApp()
```

or open `app.R` and run the app from your R environment.

## Input Data Notes

The app expects an Excel file with appropriate columns for the selected analysis.

For general selection analysis, the default identifier and replication columns are:

- `Variety`
- `Rep`

For MET analysis, the app expects:

- `Genotype` or `Variety`
- `Environment` or `Location`
- One or more numeric trait columns

For mating-design analysis, select the correct parent, line, tester, replication, type, and trait columns in the app after uploading the file.

## Output

Depending on the selected analysis, the app can produce:

- ANOVA tables
- GCA and SCA tables
- Variance component tables
- Trait summaries
- Mean comparison tables
- Superiority index tables
- Selection ranking tables
- BLUP tables
- Stability and biplot outputs
- Excel workbooks and chart downloads

## Git Backup Workflow

After editing the project, save a new backup to GitHub with:

```powershell
git status
git add .
git commit -m "Describe the update"
git push
```
