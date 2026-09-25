# PhenoSelect Analysis (Shiny)

A Shiny app for plant-breeding trial analysis. Breeders upload one Excel
workbook per trial and get statistically-sound trait analysis, a weighted
selection index, and advance/retest/discard recommendations — for both a
**single environment (Single-Location Trial / LPSI)** and **multiple
environments (Multi-Environment Trial / MET)** — without touching R code.

## What it does

| Analysis | Trigger | Core question it answers |
|---|---|---|
| Pre-Breeding | `analysis_method = BREEDING` | How has the breeding population changed across generations/years? |
| Genetic Diversity | `analysis_method = DIVERSITY` | How distinct/related are the genotypes? |
| Mating Design | `analysis_method = MATING` | GCA/SCA from Griffing I–IV, partial diallel, or Line × Tester |
| **Multi-factor** | `analysis_method = MULTIFACTOR` | Compare two or three factors and their interactions in factorial or split-plot trials, with mean separation and check-relative superiority. |
| **Single-Location Trial (LPSI)** | `analysis_method = LPSI` | Within one environment: which genotypes rank best on a weighted, multi-trait index, and should each be **advanced / retested / discarded**? |
| **Multi-Environment Trial (MET)** | `analysis_method = MET` | Across environments: which genotypes are high-performing **and** stable (AMMI, GGE, Finlay–Wickham, ASV), and what's the integrated multi-trait ranking (MGIDI-style)? |

This repo's docs focus on the **LPSI** and **MET** pipelines, since that's the
selection-index core of the app; the other modules share the same app shell,
upload flow, and export system.

## Project layout

```
project/
├── app.R                              # UI + server, wires every module together
├── modules/
│   ├── shared_data_validation.R       # upload parsing and shared helpers
│   ├── advanced_analysis_extensions.R # selection-index extensions
│   ├── module_1_breeding.R
│   ├── module_2_genetic_diversity.R
│   ├── module_3_mating.R
│   ├── module_multi_factor.R          # factorial and split-plot analysis
│   ├── module_4_selection_index.R     # thin wrapper -> advanced_analysis_extensions.R
│   └── module_5_met.R                 # MET engine: models, AMMI/GGE, FW, decision board
├── AGENTS.md
├── README.md
└── docs/
    ├── PRD.md
    └── ARCHITECTURE.md
```

## Input data format

One Excel file, one trial (sheet 1 is read). Expected columns:

- `Variety` — genotype/entry ID (configurable via `id_col`)
- `Rep` — replication (configurable via `rep_col`)
- An optional entry-type column (`Type`, `Entry_Type`, …) flagging `CHECK` /
  `CONTROL` / `COMMERCIAL` entries used as benchmarks
- One column per trait
- For MET: an `Environment` column (location/season/site)
- Two special rows layered into the trait columns to drive the index:
  - a **weight** row (`WEIGHT`/`IMPORTANCE`/…) — relative importance per trait
  - a **direction** row (`DIRECTION`/…) — whether higher or lower is better
    per trait

## Running locally

```r
# from the project root
install.packages(c(
  "shiny", "bslib", "readxl", "tidyverse", "emmeans", "multcomp",
  "multcompView", "pheatmap", "ggplot2", "DT", "writexl", "lme4",
  "lmerTest", "sommer", "patchwork", "zip"
))
shiny::runApp(".")
```

Then in the app: **Data** tab → upload workbook → **Analyze** tab → pick a
trait, pick "Single-Location Trial" or "Multi-Environment Trial" → **Run
analysis** → **Results** / **Charts** tabs → **Export** to download.

## Output

Every analysis produces a set of result tables (visible in **Results**,
charted in **Charts**) that export to a single `.xlsx` workbook per analysis,
or all of them zipped together from **Export → Download everything**.

## Docs

- [`docs/PRD.md`](docs/PRD.md) — who this is for, what it must do, what's out of scope
- [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) — data flow, statistical
  methods, module responsibilities, state management
- [`AGENTS.md`](AGENTS.md) — conventions for anyone (human or AI) editing this code
