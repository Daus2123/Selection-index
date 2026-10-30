# MET pipeline review — 30 October 2026

## Assessment

The current app is useful for exploratory trial analysis and breeder screening,
but should not yet be described as a validated automatic advancement/discard
system. Several reproducible issues can change decision outputs without a
corresponding change in valid biological evidence.

Scope: current local MET upload/preparation, approved cleaning, Mean Comparison,
restored October 26 GxE engine, uncertainty, FW/AMMI/GGE, integrated selection,
MGIDI/FA-MGIDI, decision logic, reactive controls, exports and selected tests.
The live shinyapps.io deployment and the user's paired trial exports were not
available for numerical comparison. This review does not resolve the previously
reported local/server discrepancy.

No application or analysis source was modified during this review. All runtime
checks used synthetic data. Review scripts, logs and synthetic results are in
`tmp/met_review_*` and `tmp/review_met_*`.

Restoration equivalence is not statistical validation: reproducing commit
6886251 also reproduces its limitations. Earlier passing restoration/export
checks should not be read as proof of calibrated confidence or correct final
breeding decisions.

## Findings requiring correction

### 1. P1 — Missing LMM uncertainty can pass advancement gates

Locations: `modules/module_5_met.R:1189`, `:1202`, `:1214`.

The reliability gate labels zero available reliability values as
`NOT_APPLICABLE_BLUE` without checking whether the model is actually ANOVA/BLUE.
Missing rank-confidence and check-probability values can also pass their gates.

Reproduction: starting from a valid LMM result, set SEs and reliabilities to
unavailable to simulate an uncertainty-calculation failure. With the remaining
conditions satisfied, G08 and G07 received ADVANCE despite missing SE and
reliability. The passport classified their reliability as NOT_APPLICABLE_BLUE.
Method membership was deliberately held fixed in this controlled guard test.

Required correction: distinguish unavailable random-genotype reliability from
genuinely inapplicable BLUE reliability. Required uncertainty unavailable should
route to RETEST with a specific reason. Validate required positive-weight traits
individually rather than accepting a partial set of available reliabilities.

### 2. P2 — The legacy recommendation helper contradicts the decision board

Locations: `modules/module_5_met.R:18`, `:166`, `:181`; `app.R:4971`;
`modules/shared_data_validation.R:597`.

The common/legacy breeder recommendation marks the top five integrated ranks as
ADVANCE candidates. It does not consume the decision board's precision,
coverage or check gates. The two outputs have different decision semantics.

Reproduction: force all entries below the decision board's precision requirement.
The board produces zero ADVANCE calls, while the common recommendation still
produces five. A top-five shortlist is acceptable as a shortlist; presenting it
as the same final action would be misleading. The app still defines a renderer
for this helper, but no corresponding visible `DTOutput` binding was found in
the current UI. This is an inconsistent legacy API rather than a reproduced
contradiction between two currently visible tables.

Required correction: one authoritative decision-policy output should supply
all action labels, reasons and exports. Keep rank-only lists explicitly labeled
as candidate shortlists, without final-action terminology.

### 3. P1 — Inadequate trial support does not stop pooled GxE conclusions

Locations: `modules/module_5_met.R:3947`, `:4543`, `:4546`, `:4607`, `:3017`.

The restored GxE path does not apply the current connectivity helper to gate
cross-group ranking or stability. FW uses every cell of the completed model
matrix, including untested cells, and does not require multiple observed
environments per genotype before calculating a slope.

Reproductions:

- A genotype observed at one environment received four regression points, a
  slope of 0.9842 and class Average.
- A trial network with two disconnected groups still produced eight GGE rows
  and an eight-genotype integrated ranking.

These outputs can be calculated under the model's distributional assumptions,
but do not establish observed cross-environment adaptation. A prediction is
not another independent trial.

Required correction: gate broad-adaptation claims using observed connectivity,
coverage, convergence and uncertainty. Retain supported within-group estimates.
If predicted-cell FW is retained for continuity, explicitly label it descriptive
model-implied response and prevent it from supplying observed stability evidence.

### 4. P1 — Internal ANOVA fallbacks enter FW, AMMI and GGE

Locations: `modules/module_5_met.R:4390`, `:4414`, `:4543`, `:4607`.

The ANOVA path correctly distinguishes reportable BLUEs from internal additive
fallbacks, but constructs stability matrices from BLUP_env rather than the
reportable estimate. Thus a value withheld in one table can still influence
stability and adaptation elsewhere.

Reproduction: three untested cells had zero finite reportable BLUEs but three
finite internal values. All three entered FW: 32 regression points versus 29
observed cells.

Required correction: use supported estimates consistently in downstream
calculations. Withhold methods requiring a complete matrix when that matrix
is not supported, or adopt an explicit validated incomplete-data method.

### 5. P1 — A Results browsing control changes the breeding decision trait

Location: `app.R:4968`.

The decision board receives `primary_trait = input$met_result_trait`.
The same selection used to browse trait results determines benchmark advantage,
rank confidence and adaptation evidence.

Reproduction: hold integrated ranks and method membership fixed. Switching the
primary trait between two opposing trait profiles changed G08 from ADVANCE to
RETEST, along with other candidates.

Required correction: add a deliberate decision-primary-trait setting stored
with the decision configuration. Ordinary Results navigation should not change
actions. Also separate the AMMI/GGE coverage-flag control from mandatory decision
coverage policy; the former currently influences the latter's default.

### 6. P1 — The restored GxE uncertainty is incomplete

Locations: `modules/module_5_met.R:4148`, `:4155`, `:4193`, `:1007`, `:313`.

The overall LMM estimate is grand mean plus genotype effect, but its reported
SE uses only the random-genotype conditional SD. Cell SEs sum component
variances without the required cross-component covariance. Check-difference
SEs omit cross-genotype covariance, and rank simulations assume independent
marginal normals. Target-distance comparisons need additional care because
absolute-distance transformations do not preserve normality.

Independent check: refit the same model and compare the old reported overall
SE with a joint mixed-model prediction-error calculation. The joint SE was
1.5749 times the reported SE in the balanced synthetic example. This is an
example, not a universal correction factor.

Required correction: retain the historical point-estimate model if desired,
but calculate joint uncertainty for the actual estimand. Use covariance-aware
contrasts and selection uncertainty. Until validated, label superiority and
rank probabilities exploratory and remove them as hard evidence gates.

The lme4 documentation distinguishes conditional random-effect covariance from
the full covariance structure, and notes that cross-group blocks need not be
zero: [ranef documentation](https://lme4.github.io/lme4/reference/ranef.html).

### 7. P2 — Zero-weight traits can veto a decision

Locations: `modules/module_5_met.R:420`, `:1032`, `:1046`, `:1189`.

Reliability and coverage summaries use all successful trait results, including
traits with zero selection weight. A diagnostic trait can therefore block
advancement although it contributes nothing to the index.

Reproduction: Yield had weight 1 and reliability 0.95; Height had weight 0 and
reliability 0.01. The decision reliability gate was FAIL.

Required correction: define explicit required selection traits, diagnostic
traits and mandatory biological constraints. Use the intended scope consistently
in ranking, precision gates and coverage gates.

### 8. P2 — Ranking agreement is a heuristic, not independent confirmation

Locations: `modules/module_5_met.R:853`, `:1212`, `:1029`.

The five method-support sets reuse the same fitted estimates and related
transformations. Two supporting methods are required for advancement, while
rank confidence concerns the primary trait's top five, not the integrated
selection percentage. Checks also count in the selection denominator.

Required correction: align the selection population, threshold and uncertainty
with the actual decision objective. Treat method agreement as sensitivity
analysis, not additional independent evidence. Add practical trait thresholds
or a breeder-defined meaningful check advantage; current defaults mainly use
relative rank and nonnegative estimated check advantage.

### 9. P2 — Mean Comparison is useful but its prediction inference is approximate

Locations: `modules/module_5_met.R:4813`, `:4855`, `:4872`, `:4990`.

Strengths: the model uses original observations; the prediction grid includes
joint fixed/random covariance and uncertainty for unseen interactions. It does
not refit predictions as observations.

Limitations: custom BLUP-grid contrasts use infinite degrees of freedom, fitted
variance parameters are treated as known, and variance-component significance
uses approximate boundary LRTs. Existing tests verify formulas, output consistency
and reference calculations; they do not establish Type I error or interval
coverage across sparse small trials. The notes acknowledge approximation, but
the current UI hides those notes and the chart caption was removed.

Required correction: retain the clean chart layout while exposing essential
model status in a compact status/diagnostic area. Validate small-sample and sparse
design behavior by simulation or bootstrap. Reserve conventional BLUE/Tukey
interpretations for the appropriate fixed-genotype inferential model.

See [lme4 inference guidance](https://lme4.github.io/lme4/reference/pvalues.html)
and [singular-fit concerns](https://lme4.github.io/lme4/reference/isSingular.html).
A singular fit is not automatically invalid, but needs an explicit review policy.

The two MET modes also have different models: GxE LMM treats environment as
random; the Mean Comparison LMM fixes environment and averages its cell
predictions. Same dropdown labels do not imply identical estimands. In one
balanced example, overall estimates differed by up to 0.08886. That difference
alone is not proof of a bug, but the distinction must be explained.

### 10. P2 — ANOVA effect tests and sparse-run failure handling need review

Locations: `modules/module_5_met.R:4282`, `:5335`, `:5418`.

The fixed ANOVA branch uses sequential sums of squares. With imbalance, effect
tests depend on the defined hypothesis and term order. In sparse fixed-interaction
models, overall means can be non-estimable; that limitation is legitimate.
However, the aggregate ranking runs outside the per-trait error boundary.

Reproduction: a sparse ANOVA run retained usable single-trait output, but
`run_met_all_traits()` stopped when G08's required adjusted mean was unavailable.
The Shiny run then loses otherwise usable core outputs.

Required correction: use explicitly defined estimable tests appropriate to the
model and preserve core models/comparisons when integrated ranking is unavailable.
Do not solve estimability by silently imputing values for inference.

### 11. P2 — The legacy combined wrapper falsely detects a real block

Locations: `modules/module_5_met.R:2712`, `:5357`, `:5358`.

When no block is selected, preparation creates a B1 placeholder. The combined
wrapper checks for any nonmissing Block value and therefore wrongly enables
block terms. Comparison failure is caught and reduced to a notes object.

Reproduction: the stored block-column setting was empty, yet has_block evaluated
TRUE. The comparison failed with `contrasts not defined for 0 degrees of freedom`.
`tests/test_met_wide_means.R` then failed because no means table existed.

This bug is in the compatibility wrapper added during selective restoration;
it does not affect the explicit Mean Comparison UI path, whose original block
detection was preserved. It nevertheless needs correction before a clean commit.

Required correction: derive design flags from the recorded selected columns,
not a placeholder factor. Test both explicit modes and legacy combined callers.

### 12. P2 — MGIDI labels and decision connections need explicit documentation

Locations: `modules/module_5_met.R:651`, `:783`, `:790`, `:853`;
`app.R:4951`, `:4957`, `:6254`, `:7948`.

As requested, the interface now calls the weighted Euclidean ideotype distance
MGIDI. Its calculation remains distinct from conventional FA-MGIDI.
MGIDI with FA uses the metan factor-analysis implementation with equal trait
importance; workbook weights and Mean/FW/ASV sliders do not apply to that branch.
Selecting either method changes the displayed ranking, but the decision board
continues using the integrated ranking. FA-MGIDI is not in the method-support
vote; the internal MGIDI vote is the weighted-distance result.

Required correction: retain the preferred short label if needed, but expose the
exact algorithm and active decision-ranking source in settings and exports.
Do not imply that selecting FA-MGIDI automatically makes it drive final actions.
Use the term MGIDI with FA when the conventional factor-space distance is intended.
Reference: [metan MGIDI definition](https://tiagoolivoto.github.io/metan/reference/mgidi.html).

## What passed and is useful

- Approved-cleaning and metadata regression checks passed. Extremes are flagged
  for review rather than automatically removed.
- Sparse Mean Comparison regression checks passed across LMM, SOMMER companion
  and ANOVA choices, including covariance properties, charts, disconnected
  safeguards and exports.
- The FA-MGIDI test passed against metan, including directions/targets, missing
  records, plots and all-model exports.
- In an independent balanced-data check with the same random-intercept model,
  LMM and SOMMER grand means agreed; maximum genotype-estimate difference was
  8.96e-8 and SOMMER reported convergence.
- Required positive-weight trait coverage is checked before integrated ranking.
  This is a useful safeguard; core results should be retained when it fails.
- Observed boxplots and untested predictions use different marks, and Environment
  slicing does not refit the model.

These are evidence of useful components, not blanket certification for every
trial design, package version or missing-data pattern.

## Validation status

Executed during this review:

| Check | Status |
|---|---|
| `tests/test_met_blup_comparison.R` | Passed |
| `tests/test_met_mgidi_fa.R` | Passed |
| `tests/test_data_cleaning.R` | Passed |
| `tests/test_input_metadata.R` | Passed |
| `tests/test_met_mean_ui.R` | Passed |
| `tests/test_met_complete_workflow.R` | Failed: expects prediction-confidence fields removed by restoration |
| `tests/test_met_wide_means.R` | Failed: combined-mode means absent; stale footer expectations also remain |
| Targeted reliability and independent-engine probes | Completed; findings reproduced above |

The test suite is therefore not green. Reconcile tests with intentional
restoration/UI changes while retaining tests that expose real decision defects.
Do not delete failing reliability tests merely to obtain a passing suite.

## Recommended release sequence

1. Fix missing-evidence gates, unify recommendation/action output, and decouple
   primary decision settings from browsing controls.
2. Enforce observed evidence scope for pooled ranking/stability; keep predictions
   visibly model-dependent. Remove internal ANOVA fallbacks from reportable analyses.
3. Restore valid joint uncertainty without unnecessarily changing historical
   point-estimate models. Validate the approximate comparison path with simulations.
4. Fix the combined-wrapper block flag and isolate integrated-ranking failures.
5. Clarify the active ranking method, required traits, practical thresholds,
   target-environment weights and method-agreement interpretation.
6. Reconcile and run the release regression suite. Validate the same synthetic
   reference dataset locally and on shinyapps.io before claiming reproducibility.
7. Record source/version ID, package versions, actual formula, workbook/worksheet
   identity or hash, exclusions and completed decision settings in each export.

A checkpoint commit of the current restoration can be useful if labeled as such.
It should not be labeled a statistically validated decision-support release yet.
Before staging, inspect untracked `tmp/` and `outputs/`: neither directory is
currently ignored, and existing workspace artifacts should not be included via
an unreviewed `git add .`, especially if they contain real breeding-trial data.
