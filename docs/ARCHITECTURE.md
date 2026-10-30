# Architecture — Selection Index Pipeline

### Selective GxE restoration — 30 October 2026

GxE/stability functions, ranking, FA-MGIDI, decision rules and plots are restored
from commit `68862510149ed9561e555a968c3175e1c1aa21a1` (26 October).
`met_run_gxe_20261026()` contains the historical single-trait engine.
`run_met_pipeline()` retains the current public mode contract and dispatches
GxE runs to that engine. The GxE Analyze controls again expose the historical
minimum-observed-location flag setting. Later joint prediction-SE calculations,
prediction-evidence gates, connectivity restrictions and two-location FW changes
do not apply to this restored path. Historical limitations are intentional.
The current SOMMER package API adapter is retained because installed releases
reject the October 26 optimizer arguments and changed variance-table columns;
the compound-symmetry model specification and downstream calculations are kept.

Mean Comparison remains current: `met_run_mean_comparison()` retains its
pre-restoration preparation and early-return workflow, with the same BLUP
comparison grid, covariance, tests, tables, boxplots, Environment slicer,
Observed/Adjusted slicer, CV-only footer and exports. Its connectivity helpers
remain available without changing GxE behavior. The two UI modes stay separate;
the legacy programmatic combined mode runs the historical GxE engine and then
the current comparison helper. No module is added and no result fields are
renamed. Source and numerical equivalence are checked by
`tests/test_met_selective_restore.R`. The pre-restoration source snapshot is
recorded in `tmp/met_gxe_restore_snapshot.txt`.

This restoration supersedes the later GxE-specific sections below. The current
Mean Comparison sections still apply.

### Single-stage mixed-model mean comparisons

Mean-comparison charts use boxplots of approved observed plot values, ordered
by adjusted performance, with adjusted-model comparison letters above the data.
Untested cells use diamond estimates and SE whiskers, never fabricated boxplots.
This supersedes earlier adjusted-mean bar-chart descriptions.
Interaction chart navigation offers a location slicer, defaulting to the first
available location. Both the preview and chart download use that single location;
changing location does not refit models or change Results tables. Programmatic
plot calls can omit `environment` to retain the full faceted overview.

Connected LMM mean-comparison runs now fit original plots with fixed environment,
random genotype and genotype-by-environment, and nested random design effects.
The SOMMER selection uses this explicitly labeled lme4 comparison engine; the
GxE SOMMER engine is unchanged. Fixed ANOVA retains its existing BLUE workflow.
Disconnected comparisons retain fixed-model estimability limitations and notes.

`met_comparison_blup_grid()` constructs a complete genotype/environment grid and
joint prediction-error covariance from the penalized mixed-model equations.
Unobserved interactions have zero predicted effect with their fitted interaction
variance added to prediction-error variance. Variance-component uncertainty is
not included. Predictions are never fitted again as observed data.

Overall genotype/GxE significance uses approximate random-effect LRTs; environment
uses a Satterthwaite fixed-effect test. Prediction contrasts use the joint grid
covariance with asymptotic normal/Tukey approximations. These are explicitly
labeled prediction comparisons, not classical fixed-genotype inference. Existing
means/pairwise/notes fields, Results, Charts and workbook exports are retained;
means additionally carry Estimate_type and Evidence_status. All earlier claims
that mixed-model mean comparison fixes the full interaction are superseded.

## 1. Shape of the app

### Complete MET prediction/evidence workflow

Both modes prepare original plot observations, approved exclusions and trait
metadata, then report connectivity from observed genotype/environment cells.
Selected checks are counted separately: all shared entries may link trials,
and selecting a check does not introduce artificial observations or a check-only
imputation model. `trial_connectivity` and `environment_links` are visible in
Results (QC and Mean Comparison) and exported per trait.

Check-superiority probabilities and rank-confidence simulations remain explicitly labelled exploratory marginal-normal approximations; they are not the separate Tukey/LSD inference outputs.

GxE estimates carry Confidence_status, Confidence_reason, Observed_locations,
Coverage_pct, Connected_group, Prediction_SE_ratio and Confidence_rule. Statuses
are Tested, Untested - predicted, Low confidence and Unsupported. Screening
flags are fewer than three observed locations, fewer than two replications in a
tested cell, missing prediction SE, or SE exceeding the observed trait SD.
Unconfirmed model convergence also flags supported cells. These are explicit
screening rules, not calibrated confidence probabilities. No unsupported
cross-group prediction is reported; original measurements remain unchanged.

For the current lme4 random-intercept model, prediction SE uses the inverse
penalized mixed-model precision in spherical coordinates. Fixed/random and
cross-effect covariance are included, including unseen interaction variance for
untested cells. Overall mean SE and random-genotype PEV for reliability are
separate quantities. Variance parameters are held at fitted values; parameter-
estimation uncertainty is not included. A failed uncertainty calculation returns
unavailable SEs and a note, never the old independent-component approximation.
Sommer keeps its covariance-aware PEV; fixed-model BLUEs use emmeans SEs.

Internal RCBD fallback values are not reportable and no longer enter stability
matrices. AMMI/GGE are withheld without a complete supported matrix. Both
heatmaps retain the same axes; predicted cells are faded, low-confidence cells
have dashed red borders, and unsupported cells are gray. FW regression points
are restricted to observed genotype/environment combinations; the environmental
index uses available model estimates. Two-location FW remains descriptive.

The decision board and passport expose tested/predicted/low-confidence/
unsupported cell counts. Low-confidence or unsupported evidence in positive-
weight traits gates candidates to RETEST. Missing random-genotype reliability
is not treated as not-applicable BLUE reliability. Unsupported integrated
rankings also withhold FA MGIDI rather than ranking disconnected trial groups.
The weighted non-FA option is labelled Weighted Ideotype Distance; MGIDI with FA
retains its factor-analysis meaning and equal-trait-importance behavior.

LMM predictions retain REML. Reduced LRT models are fitted by ML and one ML
full-model refit is reused across tests. Ordinary chi-square boundary p-values
remain approximate and are labelled in notes. Original observations, not filled
prediction tables, feed inference models. Mean Comparison keeps its independent
ANOVA/effect tests, Observe/Adjusted table and adjusted-only bars.

### Current MET supporting comparison outputs

Mean Comparison provides three views: Variety, Environment and Variety x
Environment. A single combined DT table aligns all completed/failed traits by
level; interaction rows also retain Environment. A top-right Observe / Adjusted
radio slicer switches raw mean plus observed SD versus adjusted mean plus SE,
without inference or refitting. Footer rows contain pooled comparison-model CV,
H2 only when supplied by an appropriate main MET model (otherwise Not estimated),
and a common or pair-specific critical difference. Non-significant/unavailable
tests have explicit Not applied footers and no letters.

Charts use adjusted-mean bars and SE whiskers only. Both Tukey and protected LSD
require a significant relevant test. Interaction displays compare varieties
within each environment, gated by the interaction test and each environment's
simple variety test. The former all-cell interaction contrast family is removed.
Non-significant means remain visible. Non-estimable values remain in tables.

One reference grid, overall joint tests, and within-environment joint tests are
reused. Family contrasts are formed only when gates pass. Letters are generated
from the same adjusted p-values as exported comparisons (coefficient indices,
not parsing genotype names). Critical differences reuse contrast confidence
limits, avoiding duplicate Tukey quantile computation. No model simplification,
imputation for inference, or change of degrees-of-freedom method is introduced.

Multi-location ANOVA and all-trait Variety/Environment/Interaction effect tests
are shown in Results and exported. Detailed workbook sheets retain effect_tests,
simple_tests, means, pairwise and notes; display downloads include Means, ANOVA
and Effects. These current rules supersede historical all-cell/four-view and
raw-chart behavior described below. GxE remains a separate analysis mode.

### MET analysis modes and mean comparisons

The Analyze sidebar places `met_analysis_mode` (GxE & Stability / Mean Comparison)
between Replication and Model. The default UI mode is `gxe`. The run records its
mode, comparison view and comparison method in settings and the result object.
Results/Charts navigation and workbook contents follow the completed mode,
not pending changes to Analyze controls.

`run_met_pipeline()` and `run_met_all_traits()` accept `analysis_mode`,
`comparison_view` and `comparison_method`. Their legacy `combined` default is
retained for existing programmatic callers; the UI always explicitly requests
`gxe` or `means`. GxE mode preserves the selected LMM/SOMMER/RCBD workflow and
both observed and model-estimated heatmaps, but never calls mean comparison.
Mean-only mode exits after common data validation, approved-data preparation,
coverage/connectivity checks and the requested comparison; it does not fit the
main prediction model, run stability extensions, or build integrated rankings.

Mean-only LMM and SOMMER selections use the existing lme4 companion inference
model: fixed Genotype * Environment with random replication/block terms nested
in environment. This actual inference engine is explicitly labeled. RCBD uses
fixed environment, replication/block, genotype and interaction terms. Missing
cells remain visible as non-estimable results. No main-model heritability is
reported for a mean-only run.

Mean Comparison calculates all four views with both Tukey and protected LSD in one run. Compare and Method selectors live in the Results and Charts sidebars and independently filter completed tables and charts. They are not shown in Analyze. Tukey adjustment is
within each displayed family. Protected LSD requires the corresponding
joint/simple-effect test at alpha 0.05; failed or unavailable gates suppress
letters, pairwise p-values and significance flags. It does not strongly control
familywise error for many means. Adjusted_mean_SE and Raw_mean_SD are available
in comparison exports. Results and Charts can switch between adjusted mean with
SE and raw mean with SD without rerunning models. Letters always refer to the
adjusted-model comparison. Changing view, method or raw/adjusted display does not refit a model. The Analysis mode stays in Analyze, so GxE continues to skip comparisons.
Non-significant families remain visible, and pair-specific thresholds are kept.

Direction/target metadata orders displayed means; letters are significance groups,
not breeding ranks. Singular fits, convergence issues, missing cells, and a
descriptive residual-SD ratio above two across environments are noted. Residual
homogeneity is an assumption, not automatically corrected by this feature.

Results and Charts have independent comparison/method selectors. Selected-table
downloads follow the Results selectors, chart downloads follow the Charts
selectors, and the full workbook exports both methods and all views per trait,
including companion ANOVA, simple-effect tests, and notes.

References: Lenth (2016), doi:10.18637/jss.v069.i01; Piepho (2004),
doi:10.1198/1061860043515; Hayter (1986), doi:10.1080/01621459.1986.10478364.

### Incomplete and two-location MET trials

Model choice remains explicit: LMM, SOMMER, or ANOVA_RCBD. The SOMMER adapter
supports both legacy mmes optimizer arguments and the newer mmes_control API.
Two-location trials retain core estimates and estimable comparisons, while FW,
AMMI and GGE stability outputs are unavailable. FW requires at least three
observed locations per genotype. Disconnected environment groups retain core
outputs for review, but AMMI/GGE and overall integrated ranking are withheld.
Missing adjusted performance under RCBD similarly withholds integrated ranking
without discarding trait models or supported within-environment comparisons.
Single-trait ranking uses common available stability components and reports
effective weights; non-estimable performance is never assigned an artificial rank.
The decision board routes unavailable stability evidence to RETEST. Comparison
letters do not enter decision scores. Existing model, selection, FW and integrated
weight tables carry these notes through Results and workbook exports.

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
MET means withhold overall ranking while preserving core trait results and
estimable comparisons. Zero-weight traits do not contribute to sums.
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

MET mean comparisons reuse the fitted main model for ANOVA_RCBD, preserving
fixed replication/block terms and sequential ANOVA. LMM and SOMMER use a
labelled companion model with fixed genotype, environment and GxE and random
design units; these adjusted means are distinct from breeding-model BLUPs.
The four views (Variety [overall], varieties across environments [separate
comparisons within each environment], environment, and variety x environment
combinations) retain equal-weight marginal means,
family-wise Tukey comparisons, protected LSD and non-estimable cells.
Mean comparison Results align all traits by Levels, with adjacent Adjusted_mean_SE
and Group columns under each trait header. Within-environment comparisons show
a separate wide table for each environment. The Results Trait selector is hidden
only for this view. Compare and Method apply across traits; selected-table and
full-workbook exports use the same wide layout (with an Environment column to
identify stacked families). Failed/unavailable traits retain their columns with
explicit status messages. Full
workbooks also include combined default-view tables and retain detailed
per-trait sheets for all comparison views and methods. Charts remain per trait.
The Results and Charts sidebars hold Compare and Method as visible radio
lists with all four comparison choices. The main comparison
panel shows ANOVA and the wide mean tables. Each comparison family ends with
CV (%), Heritability H2, and LSD (5%) or Tukey HSD (5%) rows. CV is 100 times
the comparison model residual sigma divided by the absolute observed grand mean;
zero means return unavailable. H2 comes from the main MET entry-mean model and
is shown only for overall genotype comparisons. Environment-specific CV and H2
are unavailable rather than copied from pooled results. A common critical
difference is shown only when all pairwise thresholds agree within numerical
tolerance; otherwise the footer says Pair-specific. Protected LSD reports when
its omnibus gate fails. Paging, ordering and searching are disabled to keep the
footers below their levels. Group letters still refer to adjusted model means.
Full workbook exports retain
N, model SE, pairwise tests and diagnostic notes for traceability.
Sparse mean-comparison designs retain aliased/missing genotype-environment
cells. ANOVA, omnibus, letter-group and pairwise failures are isolated so one
singular hypothesis does not discard other means or families. Failed tests are
reported in comparison notes, with unavailable p-values and blank letters where
appropriate; protected LSD never passes an unavailable omnibus test. Internal
emmeans singular-test console errors are replaced by readable notes when it
returns unavailable tests. The selected model and estimability rules are unchanged.
The combination view compares all genotype-environment cell means together;
protected LSD uses a joint test of equality of these cell means, rather than
the interaction-only test. Chart previews and downloads show observed mean
+/- observed SD with exact Group letters from the selected Results comparison,
only for families with a significant
overall/simple-effect test (p < 0.05). Other families are omitted; when no family
passes, the chart reports no significant difference. Tables retain all means
and comparison letters. The chart display selector has been removed. Chart
headings and axis titles use 10 pt; supporting text uses 9 pt and group letters
use 10.5 pt. Sidebar list options are indented 12 px beneath their headings.
Overall comparisons must be interpreted alongside GxE. No multi-factor
heritability footer is substituted for MET genetic statistics.

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

### MET Tukey-only update
The Tukey-only behavior described above supersedes older references below to
MET protected LSD or raw mean/SD comparison displays. Other analysis modules
retain their existing comparison methods.

### MET split-mode update
The analysis-mode description above supersedes earlier historical MET notes
about eager all-view/two-method computation or Tukey-only comparison displays.
No new module is introduced. Existing result field names and legacy combined
programmatic calls remain supported.
### Mean-comparison navigation update
Mean-only UI runs now request every comparison view and both methods. Settings record the full lists; detailed workbook sheets include every view/method. The combined workbook summary defaults to within-environment Tukey. Results and Chart downloads follow their respective display selectors. This supersedes earlier selected-family-only UI behavior; precomputation increases run time.

AMMI/GGE uses a fixed standard minimum of three locations, including the per-genotype coverage flag. The former half-of-locations default and adjustable UI threshold are removed. Legacy min_envs_for_biplot arguments remain accepted for API compatibility but do not override the standard. Connectivity and matrix/rank checks still apply; two-location MET retains core estimates and heatmaps.

FW now permits two observed locations per genotype in a connected design, provided the environmental index varies. Two-location slopes are labelled descriptive response rather than stable/responsive categories. AMMI and GGE retain their three-location minimum.
