# AGENTS.md

Instructions for any coding agent (or human) working in this repo. Read this
before editing `app.R` or anything in `modules/`.

## What this project is

A single-file-ish Shiny app (`app.R`) for plant-breeding trial analysis, plus
numbered `modules/*.R` files it sources. Two of those modules are the
selection-index core:

- `module_4_selection_index.R` — thin wrapper/name registry; the LPSI
  (Single-Location Trial) engine and its helpers live in `app.R`,
  `shared_data_validation.R`, and `advanced_analysis_extensions.R`.
- `module_multi_factor.R` — factorial and split-plot analysis for two or
  three treatment factors.
- `module_5_met.R` — the full Multi-Environment Trial engine: outlier
  cleaning, LMM/ANOVA modeling, AMMI, GGE, Finlay–Wickham, ASV, integrated
  multi-trait ranking, and the ADVANCE/RETEST/DISCARD decision board.

See `docs/ARCHITECTURE.md` for the full data flow and `docs/PRD.md` for scope.

## Ground rules

1. **Don't break the module contract.** `app.R` calls specific function
   names and expects specific list fields back (e.g. `result$anova_full`,
   `result$met_by_trait[[trait]]$genotype_summary`). If you rename or
   restructure a function's return value, grep `app.R` for every call site
   and every `results$<field>` / `result$<field>` access before finishing.
2. **Keep `id_col`, `rep_col`, and the weight/direction row labels in one
   place.** They're defined once near the top of `app.R`
   (`id_col`, `rep_col`, `weight_row_labels`, `direction_row_labels`,
   `type_col_candidates`, `check_type_labels`). Don't hardcode `"Variety"`
   or `"Rep"` elsewhere — reference these constants.
3. **New analysis types follow the existing five-part pattern:**
   - a `run_*` function that takes the uploaded data.frame + settings and
     returns a named list of result tables/plots
   - an entry in `saved_results` (reactiveValues in `server()`)
   - a branch in `build_export_tables()` for the xlsx export
   - `nav_panel` / `conditionalPanel` wiring in the UI for "Analyze",
     "Results", and "Charts"
   - an entry in `analysis_method` choices and `download_all`'s
     `export_names`
   Match this pattern rather than inventing a parallel mechanism.
4. **MET-specific:** trait columns can carry per-trait `weight`, `direction`,
   and `target_value` metadata (read via `make_met_data()` /
   `prepare_met_trait_settings()`). Respect `trait_direction` ("Higher
   better" vs "Lower better"/target-based) in any new scoring or ranking
   logic — don't assume higher-is-always-better.
5. **MET model type is a user choice**, not a hardcoded model. Any new MET
   output should work under `"LMM"`, `"SOMMER"`, and `"ANOVA_RCBD"`, and
   label estimates correctly (`BLUP` vs adjusted mean/BLUE) — see
   `estimate_label` in `run_met_pipeline()`.
6. **Never silently drop genotypes/environments.** Outlier filtering,
   missing-check-genotype handling, and single-environment genotypes all
   push human-readable strings into a `notes` vector that surfaces in the
   pipeline-review output. If you add filtering, add a matching note.
7. **Plots use `theme_bw()` and a small fixed palette** (e.g. `#3498DB` kept,
   `#E74C3C` outlier/flag, `#2ECC71` clean/good). Match this rather than
   introducing new ad hoc colors, so Charts stays visually consistent.
8. **Downloads are the contract with the breeder.** Every result field that
   should be exportable must be added to `build_export_tables()` — a table
   only visible in the Results tab but missing from exports is a bug.

## Module placement

Shared upload and selection helpers live in `modules/shared_data_validation.R`.
Selection-index extensions live in `modules/advanced_analysis_extensions.R`.
Keep the naming convention (`si_*` for shared/selection-index helpers,
`mf_*` for Multi-factor, and `met_*`/`si_met_*` for MET) and update
`docs/ARCHITECTURE.md` when adding modules.

## Before you open a PR / hand back code

- [ ] Grep `app.R` for every reference to any function/field you changed
- [ ] If you touched `run_met_pipeline()` or `build_export_tables()`, confirm
      both the `LMM` and `ANOVA_RCBD` code paths still run
- [ ] If you added a result table, confirm it's wired into: Results tab
      (`conditionalPanel` + `DTOutput`), Export (`build_export_tables`), and
      (if relevant) Charts
- [ ] Don't commit real breeding-trial data — use synthetic Variety/Rep/trait
      data for testing
- [ ] Match existing style: base R + tidyverse pipes (`%>%`), snake_case
      function names, `si_`/`met_` prefixes for module-scoped helpers
