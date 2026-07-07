# Selection Index Analysis App

This repository contains an R Shiny application for plant breeding selection analysis. The app supports mating-design analysis, linear phenotypic selection index analysis, and multi-environment trial analysis from uploaded Excel data.

## Features

- Upload Excel files and preview the raw data.
- Run mating-design analysis for:
  - Griffing Method I
  - Griffing Method II
  - Griffing Method III
  - Griffing Method IV
  - Partial diallel
  - Line x Tester
- Run LPSI analysis to combine multiple traits into one selection score.
- Run MET analysis across environments, including BLUP, Finlay-Wilkinson stability, AMMI, GGE, and integrated ranking outputs.
- View interactive result tables in the app.
- Download analysis results as Excel files.
- Download selected charts as image files.

## Project Files

- `app.R` - Main Shiny application, user interface, server logic, LPSI analysis, MET analysis, plotting, and export functions.
- `mating_design_module.R` - Standalone statistical functions for mating-design analysis.
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

