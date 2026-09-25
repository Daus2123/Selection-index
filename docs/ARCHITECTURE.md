# Architecture — Selection Index Pipeline

## 1. Shape of the app

A single Shiny app (`bslib` page with `nav_panel`s: **Data → Analyze →
Results → Charts → Export**), backed by numbered `modules/*.R` files sourced
once at startup:

```
app.R
 ├─ modules/shared_data_validation.R        upload parsing + shared helpers
 ├─ modules/advanced_analysis_extensions.R  LPSI / Smith-Hazel / MTSI engine
 ├─ modules/module_1_breeding.R             Pre-Breeding analysis
 ├─ modules/module_2_genetic_diversity.R    Genetic Diversity analysis
 ├─ modules/module_3_mating.R               Griffing / diallel / Line×Tester
 ├─ modules/module_multi_factor.R           factorial and split-plot analysis
 ├─ modules/module_4_selection_index.R      wrapper -> extensions engine
 └─ modules/module_5_met.R                  MET engine (this is the big one)
```

`module_4_selection_index.R` is intentionally thin — it just names the
module (`MODULE_4_NAME`) and documents the primary entry point
(`si_ext_run_multitrait_selection()`), which actually lives in
`advanced_analysis_extensions.R`. `module_5_met.R` (~4,100 lines) is
self-contained: data prep, modeling, stability analysis, ranking, and
plotting all live in that one file.

The Multi-factor module (`module_multi_factor.R`) accepts two or three
treatment factors, factorial or split-plot design, CRD or RCBD base
randomization, a response trait, and check levels. It returns settings, raw
cell summaries, Type III ANOVA, estimated marginal means with comparison
letters, pairwise contrasts, and superiority against checks within levels of
the other factors. Results compare every main effect and interaction, showing
raw mean ± SD and model-based letters when the matching ANOVA term is
significant. CV and per-pair LSD/Tukey thresholds are reported. Heritability
is left unestimated until a genotype variance model is specified. Charts use
bars for mean comparisons and interactions and a heatmap for superiority.
Tables and charts follow the existing saved-result and workbook-export flow.

## 2. Request/response flow (one analysis run)

```
Excel upload
   │  si_read_excel_upload() + si_validate_uploaded_file()
   ▼
uploaded_data()            reactive, raw data.frame (+ validation report attr)
   │
   │  user picks trait(s), analysis_method, and per-analysis settings
   │  clicks "Run analysis"
   ▼
run_*_pipeline(df, ...)    e.g. run_selection_pipeline() [LPSI] or
                                 run_met_pipeline()       [MET]
   │  returns a big named list of tables / plots / settings
   ▼
analysis_results()         reactiveVal holding the latest run's output
analysis_used()            reactiveVal: "LPSI" | "MET" | "BREEDING" | ...
saved_results$<TYPE>       reactiveValues cache, one slot per analysis type,
                            so switching tabs doesn't lose earlier runs
   │
   ├──▶ Results tab   — conditionalPanel + DTOutput per result_view
   ├──▶ Charts tab    — chart_sidebar_menu + plotOutput / downloadHandler
   └──▶ Export tab    — build_export_tables(analysis_type, results)
                          → writexl::write_xlsx() → .xlsx / zipped .zip
```

Key point: **`saved_results` is a cache, `analysis_results()` is "what's
active right now."** Downstream reactives (e.g. `lpsi_result()`,
`met_result_for_table()`) check `analysis_used()` and fall back to
`saved_results$LPSI` / `saved_results$MET` so a breeder can, say, run MET,
switch to look at an old LPSI export, and still get the right data — this is
why almost every result-reading reactive has an `if (identical(analysis_used(),
"X")) analysis_results() else saved_results$X` branch.

## 3. Data model

Input workbook (first sheet), long format, one row per plot/observation:

| Column | Role | Configured by |
|---|---|---|
| `Variety` | genotype/entry ID | `id_col` |
| `Rep` | replication | `rep_col` |
| `Type` / `Entry_Type` / … | optional; flags CHECK/CONTROL entries | `type_col_candidates`, `check_type_labels` |
| `Environment` | required for MET only | fixed name, read in `make_met_data()` |
| one column per trait | trait values | — |
| a `WEIGHT`/`IMPORTANCE` row | per-trait importance, layered into the sheet | `weight_row_labels` |
| a `DIRECTION` row | `Higher better` / `Lower better` (+ optional target) per trait | `direction_row_labels` |

These weight/direction/target values are extracted once (LPSI: inside the
`advanced_analysis_extensions.R` engine; MET: `make_met_data()` →
`prepare_met_trait_settings()`) and threaded through every downstream
scoring function as `trait_weight`, `trait_direction`, `target_value` — no
analysis function should assume "higher is better."

## 4. Single-Location Trial (LPSI) — statistical pipeline

Engine: `run_selection_pipeline()` in `app.R`, with helpers in
`shared_data_validation.R` and `advanced_analysis_extensions.R`. The
pipeline is:

1. **Per-trait ANOVA** (RCBD; `run_simple_anova`) and, optionally, **LSD or
   Tukey mean comparison** (`run_lsd_test`, `lsd_significance_alpha`) with
   compact-letter-display grouping (`multcomp`/`multcompView`).
2. **Standardization** of trait means (so traits on different scales combine
   fairly) → `standardized_scores`.
3. **Weighted combination** into a single linear phenotypic selection index
   (LPSI / Smith-Hazel style) → `weighted_contributions`, `index_ranking`,
   using `weight_table` (from the sheet's weight row) and
   `trait_correlation` (to account for trait relationships).
4. **Superiority over checks**: genotype performance relative to the
   selected/derived check set → `superiority_index`.
5. **Heritability & expected genetic gain** (`lpsi_selection_intensity`,
   default 0.10) → `heritability_gain`.
6. **Final decision**: ADVANCE / RETEST / DISCARD per genotype, from the
   index plus priority rules:
   - `advance_index_cutoff` / `retest_index_cutoff` — index-value thresholds
   - `priority_advance_cutoff_pct` / `priority_severe_weak_pct` — percentile-
     based overrides
   - `priority_weight_cutoff` — traits above this weight are treated as
     "priority traits" for override logic
   → `final_decision`.
7. Supplementary views: `si_lpsi_direct_selection()` (rank on one trait only,
   for comparison), `si_lpsi_method_comparison()` (index vs. direct
   selection agreement), `si_lpsi_pipeline_review()` (plain-language
   pipeline summary), `si_breeder_recommendation_table()` (shared format
   with MET).

Settings can change **after** the initial run without re-modeling: `app.R`'s
`lpsi_settings()` + `run_lpsi_with_settings()` re-run the pipeline cheaply
whenever checks/model/cutoffs change while `analysis_used() == "LPSI"`.

## 5. Multi-Environment Trial (MET) — statistical pipeline

Engine: `run_met_pipeline()` in `module_5_met.R`. Runs once per trait
(`run_met_all_traits()` loops it) into `met_by_trait[[trait]]`.

1. **Input prep** — `make_met_data()` reshapes to `(Genotype, Environment,
   Rep, Weight)` and reads that trait's `trait_direction` / `target_value` /
   `trait_weight`.
2. **Two-stage outlier cleaning**:
   - cell-level: per (Genotype × Environment) IQR fence
   - environment-level: per Environment IQR fence on the cell-cleaned data
   Both steps produce a before/after histogram and an `outlier_summary` row
   count table; if a step removes every row, the pipeline stops with a clear
   error rather than modeling empty data.
3. **Presence/control resolution** — builds a genotype × environment
   presence matrix (`presence`), resolves the check/control set from
   `check_varieties` or (fallback) genotypes present in *every* environment,
   and flags single-environment genotypes as high-uncertainty. Every
   fallback/gap is appended to a `notes` vector surfaced in the UI.
4. **Modeling** — user chooses one engine per run:
   - **Standard LMM** (`lme4`/`lmerTest`): Genotype, Environment, G×E as
     model terms (formula assembled by an internal `met_model_formula()`
     helper, replication/block included when present) → BLUPs.
     `met_safe_lrt()` runs likelihood-ratio tests for variance-component
     significance.
   - **Breeding LMM** (`sommer::mmes`): the same compound-symmetry genetic
     model (G + G×E) and design terms → BLUPs, PEV-based reliability, and
     prediction SEs that include covariance among fitted effects. The
     adapter requires `sommer >= 4.4.1` and normalizes its output to the same
     result fields used by the rest of the MET pipeline.
   - **ANOVA (RCBD)**, joint across environments → adjusted means (BLUEs).
   The pipeline picks correct labels throughout (`estimate_label`,
   `estimate_label_plural`, `estimate_matrix_source`) so tables/exports say
   "BLUP" vs. "adjusted mean (BLUE)" correctly for whichever model ran.
5. **Genotype & environment estimates** — `met_genotype_estimates_table()`,
   `met_location_estimates_table()` (overall and per-environment values),
   plus `build_met_qc_table()` for model diagnostics.
6. **Stability analysis** (only when `n_envs_total >= min_envs_for_biplot`,
   itself defaulted from `MET_MIN_ENVS_FOR_BIPLOT` / half the environment
   count):
   - **Finlay–Wickham regression** — genotype-specific regression on the
     environment mean (stability slope) → `fw_results`.
   - **AMMI** — `build_met_ammi_support()` → `ammi_genotype`, plus notes.
   - **GGE** — SVD-based biplot machinery (`met_gge_svd_scores()`,
     `met_gge_aec_axis()`, `met_gge_circle_data()`) feeding
     `build_met_gge_decision_views()`: mean-vs-stability ranking,
     which-won-where, environment grouping, PC variance explained.
7. **Per-trait selection ranking** — `build_met_selection_ranking()`
   combines mean performance with FW/ASV stability using breeder-set
   component weights (`normalize_met_component_weights()`, from the Mean/
   FW/ASV sliders), producing `met_selection` per trait.
8. **Cross-trait integration** — `build_met_integrated_ranking()` combines
   every trait's ranking (respecting each trait's own weight/direction/
   target) into one ranking (`met_integrated_ranking`) — the MGIDI-style
   multi-trait step, computed at the `app.R` level as
   `weighted_met_integrated()` since it needs *all* traits' results at once,
   not just one.
9. **Decision board** — `build_met_multitrait_methods()` +
   `build_met_decision_board()` turn the integrated ranking plus a
   breeder-set selection-intensity percentage (`met_selection_pct`) into a
   final ADVANCE / RETEST / DISCARD board, with the threshold rules shown
   alongside it (`met_decision_thresholds_table`).
10. Supplementary: `si_met_pipeline_review()` (plain-language summary),
    `si_met_trait_quality()` (per-trait QC/heritability-style quality
    check), `si_met_breeder_recommendation()` / `si_met_action()`.

Unlike LPSI, MET's cross-trait step (integrated ranking, decision board) is
**not** cached inside the per-trait result — it's recomputed reactively in
`app.R` (`weighted_met_integrated()`, `weighted_met_decision_board()`)
whenever the weight sliders or selection-intensity input change, without
re-running the (expensive) per-trait models.

## 6. State management (server-side reactives, `app.R`)

- `uploaded_data()` — raw data.frame from the current upload
- `diagnostic_data()` — Shapiro/residual-diagnostic view, independent of
  which analysis is selected
- `analysis_results()` / `analysis_used()` / `analysis_message()` — "last
  run" trio
- `saved_results` (reactiveValues: `MATING`, `BREEDING`, `DIVERSITY`,
  `LPSI`, `MET`) — one cached result per analysis type, reset on new upload
  (`si_reset_analysis_state()`)
- LPSI-specific: `lpsi_settings()`, `lpsi_result()`, `lpsi_direct_selection_r()`,
  `lpsi_method_comparison_r()` — cheap re-derivations from the cached result
- MET-specific: `met_result_for_table()` / `met_result_for_plot()` (per-trait
  lookup into `met_by_trait`), `met_component_weights()`,
  `weighted_met_selection_for_table/plot()`, `weighted_met_integrated()`,
  `met_selection_pct()`, `weighted_met_decision_board()`
- `breeder_recommendation()` — merges whatever LPSI/MET results exist (cached
  or active) into one shared recommendation table

## 7. Export system

`build_export_tables(analysis_type, results)` is the single source of truth
for "what goes in the downloadable workbook" — it's a big `if/else` per
analysis type, each branch calling a local `add_sheet(prefix, name, table)`
helper that also handles Excel's 31-character sheet-name limit and
de-duplication. `write_analysis_workbook()` wraps this + `writexl::write_xlsx()`
for the four single-analysis download buttons; `download_all` loops every
non-null `saved_results` entry (using the *live* result if that analysis is
currently active) plus a standalone breeder-recommendation sheet, and zips
everything with `zip::zipr()`.

**Rule for extenders:** if a value is computed for the Results tab, it
belongs in `build_export_tables()` too — the Results tab and the export
should never diverge.

## 8. Extension points

- **New trait-level MET stability metric:** add to `run_met_pipeline()`'s
  per-trait output list, then to `build_met_selection_ranking()`'s
  component-weight system if it should feed the ranking, then to
  `build_export_tables()`'s `MET` branch.
- **New LPSI decision rule:** add the cutoff to the constants block near the
  top of `app.R`, wire a matching input in the LPSI `conditionalPanel`, and
  extend `lpsi_settings()` / the LPSI engine.
- **New analysis type entirely:** follow the five-part pattern in
  `AGENTS.md` §3.
