# Architecture — Selection Index Pipeline

## 1. Shape of the app

### Missing traits in selection indices

The single-location mean-comparison selector uses SE for adjusted means plus
model-based SEs from emmeans under the selected CRD/RCBD/LMM fit, and SD for
observed arithmetic means plus within-genotype sample SDs. Statistical letters
remain based on model comparisons. The shared display/export builder applies
the same selection to downloads. Unavailable adjusted estimates (including
legacy nonparametric results) are labelled unavailable without a raw-mean fallback.

Single-trial and integrated MET indices require usable evidence for every
positive-weight trait for every entry. Missing individual plot measurements are
allowed when the fitted trait model remains estimable. A genotype with an entirely
unmeasured required trait stops index generation with the genotype/trait names;
the app does not impute values or redistribute weights across remaining traits.
Entirely missing required traits and missing required MET model results also stop
index generation. Resolve the records, explicitly remove an entry, or revise
weights consistently for all entries before rerunning. Non-estimable adjusted
means continue to stop the analysis. Zero-weight traits do not contribute to sums.
Constant observed traits standardize to zero; missing values remain missing.
The policy is included in existing weight/ranking tables and their exports.
The optional metan multivariate selection engine requires complete trait records
and now rejects partial records instead of silently deleting plots. MET stability
component availability remains separately reported in Component_notes.

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
significant. CV and per-pair LSD/Tukey thresholds are reported. Compact letters
start with `a` at the preferred adjusted mean, using the last DIRECTION metadata
row by default and the Results performance selector as an explicit override.
The selector displays only Higher better / Lower better, initialized from that
trait's metadata. Manual overrides are retained per trait and reset on upload.
The same direction applies to all effects, within-factor comparisons, charts,
and exports. Omnibus p-values remain in ANOVA rather than the mean table.

`mf_heritability()` estimates generalized broad-sense Cullis H² for Factor A,
assuming A contains genotypes. A separate REML model makes A and its interactions
random, keeps the remaining treatment effects fixed, and includes the design's
block/whole-plot/subplot random terms. The inverse joint mixed-model equations
provide the mean prediction-error variance of genotype differences, including
fixed-effect uncertainty; H² = 1 - mean PEV difference / (2 * genetic variance).
Zero genetic variance, convergence failures, and unidentifiable fits produce an
explicit unavailable value and a Summary note. Boundary nuisance components
are retained with a note. See the [genotype-difference heritability formulation](https://pmc.ncbi.nlm.nih.gov/articles/PMC6707473/).
The estimate is cached and exported in the Factor A mean-comparison footer;
assumptions and warnings are exported in Notes. B/C and interaction effects are
not labeled as genetic heritability.

Charts use bars for mean comparisons and interactions and a heatmap for superiority.
Mean-comparison chart selectors independently choose one factor (A/B/C), a
two-factor interaction (A:B/A:C/B:C), or the three-factor interaction (A:B:C),
using uploaded column names as labels. The preview and chart download use the
same selection, without changing the Results effect filter.
The Charts menu exposes Superiority and Mean comparison; interaction choices
live within Mean comparison rather than a separate Factor interaction entry.
Tables and charts follow the existing saved-result and workbook-export flow.

Multi-factor analysis accepts multiple response traits through `mf_run_traits()`.
The result retains the first trait's fields for compatibility and stores each
successful fit in `results_by_trait`; failed traits are named in Summary notes.
Results and Charts have independent trait selectors. Comparison effect,
direction, and method controls are visible only in the Mean comparison Results
view. Full workbook downloads combine all selected traits with a `Trait` column;
individual-view downloads use the currently displayed trait. The ANOVA table
includes computed mean squares and significance codes. Mixed-model residual
rows report REML variance without inventing classical residual SS or df.
Charts use a shared larger type scale and a fixed categorical palette for levels.
Results and Charts place the trait and contextual selectors in one compact,
full-width stack with consistent margins and no nested panel indentation.

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

The MET Multi Selection chart and Results table offer three methods:
Integrated Ranking (performance/FW/ASV blend), MGIDI (the existing workbook-weighted
distance in standardized trait space), and MGIDI with FA (`met_mgidi_fa()`, using
`metan::mgidi()` on direction-adjusted estimates from the selected MET model).
The FA method uses equal trait importance, 0–100 desirability rescaling, factors
with eigenvalues at least 1, and varimax rotation where multiple factors remain.
It does not use workbook weights or Mean/FW/ASV sliders. The selection percentage
controls the highlighted fraction, including ties at the cutoff.

`met_multitrait_methods$mgidi_fa` stores ranking, loadings/communalities, PCA
variance and retention, genotype scores, ideotype scores, factor contributions,
and notes. Results shows ranking, notes, loadings and variance; all seven tables
export through `build_export_tables()`. Constant/unusable traits are excluded
with notes. Incomplete genotypes remain visible with missing ranks; insufficient
or failed FA returns an explanatory unavailable result, never a substituted
distance method. FA notes also appear in pipeline review. FA selection does not
change the existing decision-board method-support counts.

Input workbook (first worksheet), long format, one row per plot/observation:

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
2. **Approved cleaning**: `modules/data_cleaning.R` provides trait-level screening
   and review. Original upload row numbers are retained. Only approved exclusions
   become missing trait measurements before analysis. MET IQR fences flag values
   for exploration and never delete them automatically. Historical outlier_summary
   fields remain for compatibility and report zero automatic removals. Applied
   exclusions and settings export separately in QC tables.
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

## Single-trial LMM CV and heritability

`si_lpsi_lmm_statistics()` in `advanced_analysis_extensions.R` separates two fits:
the selected comparison model `trait ~ ID + (1 | Rep)` supplies residual CV,
100 * residual SD / abs(observed trait mean), while a separate REML model
`trait ~ 1 + (1 | ID) + (1 | Rep)` supplies genetic variance and generalized
Cullis H2 = 1 - mean PEV of genotype differences / (2 * genetic variance).
The joint mixed-model equations include fixed-effect uncertainty and pairwise
covariance; the calculation supports unbalanced data. Mean comparisons retain
their original fixed-genotype model.

The mean-table footer labels LMM outputs Residual CV (%) and Broad-sense H2
(Cullis). The existing `heritability_gain` table, Results view and export carry
the estimates, separate CV-model residual variance, mean difference PEV, method,
model and diagnostic notes. Error_variance is from the H2 variance model; it
need not equal CV_model_residual_variance. ANOVA mean squares are not invented
for REML fits. Zero genetic variance with positive residual variation is reported
as H2=0 using a labelled zero-reliability limiting convention, rather than direct
division by zero in Cullis' formula. Nonconvergence and unidentifiable fits leave
H2 unavailable; a failed H2 fit does not discard an available CV.
Boundary nuisance variances are retained with notes. Near-zero means suppress CV.
Exact constant/additive score datasets with residual degrees of freedom are
recognized before optimization: residual CV is zero (unless the trait mean is
zero), while H2 is labelled Constant trait or Zero residual variation instead of
forcing a non-regular variance ratio. Regular models use standardized responses
with variances restored to original units and alternate optimizers on failure.
Per-trait CV_status and H2_status explain unavailable footer values and export
with the full diagnostic notes. This does not prove independent replication or
perfect heritability for score values copied across repetitions.
No conventional genetic-advance value is inferred from generalized H2.
See the [genotype-difference heritability formulation](https://pmc.ncbi.nlm.nih.gov/articles/PMC6707473/).

## Reviewed data cleaning

`modules/data_cleaning.R` supplies `si_clean_screen()` and `si_clean_apply()`.
The Analyze review panel accumulates per-trait decisions, requires exclusion reasons,
and preserves metadata rows and all other trait measurements. Residual screening
uses `si_trial_diagnostic_fit()` for both flags and diagnostic plots. Runtime
context comes from selected single-trial, MET or Multi-factor controls; it keeps
fixed/random and split-plot error strata consistent. MET includes environment,
GxE and nested replication/block terms. SOMMER preview diagnostics are explicitly
labelled preliminary lme4 CS approximations; final MET diagnostics use the actual
SOMMER fit. No spatial covariance model is introduced by this change.
The initial diagnostic chart always uses raw observations. The review uses the
existing trait selector, visible detection settings and flagged values by default.
Preview changes builds a draft dataset and raw/proposed comparison charts without
mutating approved data or clearing results. Accept reviewed data commits the draft
and invalidates saved analyses. Changed settings or selected decisions invalidate
the preview, so it must be regenerated before acceptance. Column mapping and
the full-observation view are inside Advanced settings.
Applying cleaning invalidates all saved results; the next run uses approved data.
QC sheets store review, settings, impact and, when a prior matching analysis exists,
before/after rank/action tables. Original and approved data are available in the
cleaning download. No automatic model refit-and-delete loop is performed.

## TYPE and UNIT input metadata

Below the observations, metadata labels belong in the configured genotype column
(`Variety` by default). Optional TYPE accepts `score` or `quantitative`;
`continuous` and `continues` are accepted quantitative aliases. Optional UNIT
accepts `score`, `gr`, `kg`, `number`, `cm`, `percent`, `brix`, or other unit text.
UNIT-only `score` declares a score; measurement units declare quantitative data.
TYPE overrides legacy value-based inference, and conflicting TYPE/UNIT declarations
or duplicate TYPE/UNIT rows produce errors. A missing declaration preserves the
legacy integer-1-to-5 heuristic; the metadata table labels this source as inferred.

Shared `si_trait_metadata()` returns Trait, Type, Unit and Metadata_source.
All WEIGHT/DIRECTION/TYPE/UNIT rows are excluded from observations. Entry/ID/plot
columns and common design fields are excluded from automatic trait detection.
Declared scores allow fractions within 1-5; out-of-range scores and percentages
are review flags, not automatic deletions. Units do not convert measurements or
introduce new model families. The Data and Results metadata tables and full
workbook exports expose declarations. A single comma or dot remains a decimal
separator; multiple/mixed separators fail numeric conversion rather than being
interpreted as thousands. Excel's stored numeric value is authoritative.

The app automatically reads the first worksheet through `si_read_excel_upload()`.
The raw preview preserves all rows/columns of the selected sheet, including input
metadata. Upload status displays filename, worksheet, ID columns and TYPE/UNIT
row counts, plus scoring and quantitative trait counts. The extra trait/unit table
and worksheet selector are not displayed. Uploading a new file resets cleaning
approvals and saved analyses.

## Analysis reliability: design preservation and estimation failures

`get_lpsi_normality_decision()` tests residuals of the selected CRD/RCBD/LMM model.
Shapiro-Wilk p-values and warnings are reported in `anova_full`; they never trigger
automatic unblocked Kruskal-Wallis. Existing rendering of historical KW/Dunn
results remains supported, but new single-trial runs retain the selected model.

Single-trial adjusted estimation fails explicitly on model/emmeans errors or
non-estimable means, before rankings are generated. There is no substitution of
raw means. Explicit descriptive raw-mean displays remain labelled raw. The shared
diversity BLUE/BLUP estimator also refuses to substitute raw means when replication
or the required estimation package is absent.

`residual_diagnostics` records observed/fitted values and raw/standardized residuals
from each final single-trial fit. Diagnostics are accessed in Analyze; duplicate
Model diagnostics entries are omitted from Results and Charts.
`build_export_tables()` retains the final diagnostics table for reproducibility.

### MET with two locations

Core MET modeling, integrated ranking, decisions, and workbook exports remain
available with two locations. AMMI/GGE require at least three locations under
the app policy. Optional GGE summaries retain typed empty columns before
grouping, so unavailable geometry does not abort otherwise valid MET results.
MET notifications strip terminal ANSI colors before displaying error messages.

MET replication selection never falls back to the first uploaded column. It
prefers a recognized replication name, rejects identity/environment choices,
and otherwise leaves an explicit unselected prompt. The pipeline validates
that replication and optional block columns are distinct from identity,
environment, and the response before fitting models.

Integrated MET ranking accepts an explicit trait scope from the analysis run.
Only selected traits are required; a failed selected positive-weight trait still
blocks the index. The result settings carry this scope into reactive ranking,
charts, and export recalculation. Unselected workbook traits are not treated as
failed models. Two-location regression checks render the non-AMMI/GGE charts
and exercise selections smaller than the uploaded trait set.
