# Using reviewed data cleaning

1. Upload the workbook and choose the analysis and trial settings.
2. Optionally run the original-data analysis first. Applying cleaning will retain
   that run as a baseline for available rank/action comparisons, then clear results.
3. In Analyze, expand **Data cleaning review** below the diagnostic preview.
4. Inspect the raw chart, then set the detection method, flag threshold and optional
   trait minimum/maximum. Click **Check selected trait**. The default is residual
   screening at threshold 3, with detected identifier/design columns.
5. Only unusual numeric values, values outside your limits, or already excluded
   values for the selected trait appear by default. Missing values and repeated
   identifier/design combinations can be included under Advanced settings.
   Row is the original uploaded data-row number, not the Excel worksheet row number.
6. If needed, expand **Advanced settings** to check identifier/design columns
   or show all observations. Include
   relevant block/treatment columns; do not select response traits as design columns.
7. Select table rows. Choose Keep or Exclude, enter a reason, and click
   **Preview changes**. Exclusions require reasons. A flag is not proof of error.
8. Compare the raw and proposed charts. Previewing does not change analysis data.
   Adjust settings/choices and preview again if needed. When satisfied, click
   **Accept reviewed data**. To retain all values, preview with no rows selected
   and accept. Repeat for other traits; earlier accepted decisions are retained.
9. Previous results are cleared only after acceptance; click **Run analysis** to
   refit approved data. Impact details are inside Advanced settings.
10. The first chart always shows raw data. The comparison below shows raw versus
    proposed data. QQ plots show residuals; choose Observed trait histogram to
    inspect measurements. Only accepted exclusions feed the final analysis.
11. Results contains an Applied cleaning record. Available rank/action comparisons
    show the run before the latest cleaning change versus the new run. The baseline
    is only comparable if analysis/model/trait and decision settings are held fixed.
12. Download the cleaning workbook for original data, approved data, review, settings
    and impact. Full analysis workbooks also contain QC sheets. Records are session
    scoped: download before closing or uploading a new workbook.

## Interpretation and limits

- Exclude affects only the chosen trait measurement; other traits and metadata stay.
  For a whole-plot error, review/exclude that row in every affected trait.
- Missing identifiers or duplicates require investigating the source workbook.
  Excluding a trait does not repair identifiers. Correct source errors and reupload.
- Cleaning flags and diagnostic charts use one shared fit service driven by the
  selected analysis settings. CRD/RCBD/LMM retain their selected fixed/random terms;
  MET includes environment, GxE and nested replication/block terms. Multi-factor
  screening uses its analysis engine, including split-plot strata. SOMMER screening
  uses a labelled preliminary lme4 compound-symmetry approximation; final MET
  diagnostics remain based on the actual SOMMER fit. Too little replication can
  make screening unavailable; notes say so.
- After a single-trial run, use **Model diagnostics** in Results or Charts to inspect
  residuals from the final approved-data model. These observations export with the
  analysis workbook. Shapiro-Wilk results are advisory and never automatically
  replace the selected design model with Kruskal-Wallis.
- IQR screening pools observations for the chosen trait. Use it for exploration,
  especially cautiously when environments differ; final MET charts also show
  environment-specific IQR flags without automatic removal.
- Keep plausible genetic extremes. Investigate field notes and measurement errors.
  Do not delete values simply to improve QQ appearance or Shapiro p-values.
- Accepted data remains ready for analysis when browsing traits or changing controls.
  To change exclusions, explicitly start a new check, preview and accept it. An
  unfinished new review does not block analysis of a previously accepted snapshot;
  the app says that unaccepted draft changes are not included. Before the first
  acceptance, an unfinished review must be accepted before analysis. Merely editing
  settings does not revoke previously accepted exclusions. If no cleaning review is used,
  analyses retain valid available observations.
- Approved exclusions are shared across analyses in the session. Applying cleaning
  clears all cached analyses so they must be rerun on the new approved data.
