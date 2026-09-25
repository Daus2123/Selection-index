# PRD — Selection Index Pipeline (Single-Environment & Multi-Environment)

## 1. Problem

Plant breeders run trials with many genotypes, replications, and traits, and
need to decide which genotypes to **advance, retest, or discard** — for one
trial location (single environment) and, once a genotype has been through
several sites/seasons, across all of them (multi-environment). Doing this
correctly requires ANOVA/mixed-model statistics, multi-trait index
construction, and stability analysis (AMMI/GGE/Finlay–Wickham) that most
breeders don't want to run by hand in R for every trial. Mistakes here
(picking the wrong model, ignoring trait direction, letting one loud trait
dominate the index) lead to advancing the wrong material.

## 2. Users

- **Primary:** plant breeders / trial coordinators who have clean trial data
  in Excel and want a defensible ranking + decision, not a stats lecture.
- **Secondary:** breeding-program analysts/statisticians who want to inspect
  the ANOVA, variance components, heritability, and stability tables behind
  the ranking, and who set trait weights/directions.

## 3. Goals

1. Turn one uploaded trial workbook into a **weighted, multi-trait selection
   index** and an **advance/retest/discard recommendation** per genotype,
   for a **single environment**.
2. When the same genotypes have data from **multiple environments**, produce
   an **integrated, stability-aware ranking** and decision board, so a
   genotype that's great in one place but wildly unstable elsewhere isn't
   blindly advanced.
3. Make every number traceable: every ranking table sits next to the ANOVA /
   variance-component / heritability table that justifies it.
4. Let the breeder tune the decision without re-uploading data: check-variety
   selection, trait weights, advance/retest/discard cutoffs, and
   selection-intensity are all live controls.
5. Every result table is downloadable (per-analysis `.xlsx`, or everything as
   one `.zip`).

## 4. Non-goals

- Not a general-purpose statistics tool — only the analyses listed in
  `analysis_method` (Pre-Breeding, Genetic Diversity, Mating, LPSI, MET).
- Not a trial-design tool (randomization, plot layout) — input is assumed to
  already be a valid RCBD (or similar) design.
- Not multi-tenant / multi-user data storage — one workbook per session,
  nothing persisted server-side beyond the session.
- Not responsible for raw-data QA beyond the built-in upload validation and
  the automatic cell-/environment-level outlier filtering.

## 5. Feature: Single-Location Trial (LPSI)

**User story:** *As a breeder with one trial location, I want to rank
genotypes on a weighted combination of my traits and get an
advance/retest/discard call, so I know what to keep without doing the stats
myself.*

Requirements:
- Compute per-trait ANOVA (RCBD) and, when requested, LSD or Tukey mean
  comparisons with compact letter display.
- Read trait **weight** and **direction** (higher/lower-is-better) from the
  uploaded sheet, so the index reflects breeding priorities, not raw scale.
- Standardize traits before combining them (so a trait measured in the
  thousands doesn't dominate one measured in single digits).
- Produce: trait summary, ANOVA, mean comparison, superiority-over-check
  index, the combined selection index ranking, standardized scores,
  per-trait weighted contributions, trait correlation matrix, heritability
  and expected genetic gain, and a final ADVANCE/RETEST/DISCARD decision
  per genotype, using configurable cutoffs
  (`advance_index_cutoff`, `retest_index_cutoff`, `priority_advance_cutoff_pct`,
  `priority_severe_weak_pct`, `priority_weight_cutoff`).
- Support re-running the decision instantly when the breeder changes checks,
  model, or cutoffs — no full re-analysis / re-upload needed.
- Offer a single-trait "direct selection" view and a method-comparison view
  (index selection vs. direct selection on one trait) so a breeder can sanity
  check the index against simple intuition.

## 6. Feature: Multi-Environment Trial (MET)

**User story:** *As a breeder with the same genotypes tested across several
locations/seasons, I want to know which ones perform well **and** reliably
across environments, so I don't advance something that only got lucky once.*

Requirements:
- Accept data with an `Environment` column; fit either a **linear mixed
  model (LMM)**, giving BLUPs, or a **joint ANOVA (RCBD)**, giving adjusted
  means (BLUEs) — user-selectable per run.
- Two-stage outlier cleaning (cell-level, then environment-level, via
  IQR fencing) with before/after diagnostic plots and a row-count summary.
- Handle genotypes with missing environments; flag genotypes present in only
  one environment as high-uncertainty rather than silently ranking them
  alongside fully-tested ones.
- Let the breeder pick check/control genotypes explicitly, or fall back to
  genotypes present in *every* environment.
- Compute, per trait: genotype summary + QC table, model summary and
  variance components, genotype estimates (overall and per-environment),
  Finlay–Wickham regression (stability), AMMI, and GGE (mean-vs-stability,
  which-won-where, environment grouping, PC variance) — but only when there
  are enough environments (configurable minimum) for a meaningful biplot.
- Combine mean performance with stability metrics (Finlay–Wickham slope,
  AMMI Stability Value) into one **weighted per-trait selection ranking**,
  with breeder-adjustable weights (mean / FW / ASV sliders).
- Combine **across traits** into one **integrated ranking** (MGIDI-style
  multi-trait index) respecting each trait's direction/target and weight.
- Produce a final **decision board** (ADVANCE / RETEST / DISCARD) driven by
  a breeder-adjustable selection-intensity percentage, with the decision
  rules themselves shown alongside the board.
- Share a common "breeder recommendation" view/table with the LPSI module so
  a breeder moving from single-site to multi-site data sees a consistent
  recommendation format.

## 7. Cross-cutting requirements

- One Excel upload (≤ 50 MB) per session; first-sheet data only.
- Every result table must be exportable (xlsx, and CSV for single tables).
- Every chart must be downloadable as PNG (300 DPI) or PDF, at a
  user-chosen size.
- Changing a decision setting (checks, cutoffs, weights) must update results
  without a full re-run of the underlying models where avoidable (LPSI
  decision settings re-run cheaply; full MET remodeling only on explicit
  "Run analysis").
- Errors (bad columns, all rows filtered as outliers, invalid model choice)
  must produce a clear, actionable message, not a crash.

## 8. Success metrics (qualitative, since this is an internal analysis tool)

- A breeder can go from "workbook in hand" to "ranked list + decision" in
  under a few minutes without reading R code.
- Every advance/retest/discard call is explainable by pointing at a specific
  table (ANOVA, heritability, stability) already in the Results tab.
- Adding a new trait or environment to the workbook requires no app changes.

## 9. Open questions / future scope

- Should LPSI and MET recommendations be reconciled automatically (e.g. a
  genotype ADVANCE in LPSI but RETEST in MET) rather than shown separately?
- Should decision cutoffs be trait-family-specific (e.g. yield vs. disease
  score) rather than a single global set?
- Persisting settings/results across sessions (currently session-only).
