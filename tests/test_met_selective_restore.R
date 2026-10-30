suppressMessages(invisible(source('app.R')))
snapshot <- if (file.exists('tmp/met_gxe_restore_snapshot.txt'))
  readLines('tmp/met_gxe_restore_snapshot.txt', warn = FALSE)[1] else NULL
have_snapshot <- length(snapshot) == 1 && file.exists(file.path(snapshot, 'modules/module_5_met.R'))
before <- if (have_snapshot) new.env(parent = .GlobalEnv) else NULL
historical <- new.env(parent = .GlobalEnv)
if (have_snapshot) sys.source(file.path(snapshot, 'modules/module_5_met.R'), envir = before)
historical_source <- system2('git', c('show', '68862510149ed9561e555a968c3175e1c1aa21a1:modules/module_5_met.R'), stdout = TRUE)
stopifnot(length(historical_source) > 4000)
historical_file <- tempfile(fileext = '.R')
writeLines(historical_source, historical_file)
sys.source(historical_file, envir = historical)
unlink(historical_file)
preserved <- c('met_trial_evidence', 'met_design_connected', 'met_mean_comparison_views',
  'met_comparison_blup_grid', 'met_build_mean_comparison', 'met_comparison_table',
  'met_combined_comparison_tables', 'met_wide_mean_table', 'met_mean_display_table',
  'met_plot_mean_comparison')
if (have_snapshot) for (name in preserved) {
  stopifnot(identical(body(get(name)), body(before[[name]])),
    identical(formals(get(name)), formals(before[[name]])))
}
for (name in ls(historical)) {
  original <- historical[[name]]
  if (!is.function(original) || name %in% c('run_met_pipeline', 'run_met_all_traits', 'run_met_sommer_cs')) next
  stopifnot(identical(body(get(name)), body(original)), identical(formals(get(name)), formals(original)))
}
stopifnot(identical(body(met_run_gxe_20261026), body(historical$run_met_pipeline)),
  identical(formals(met_run_gxe_20261026), formals(historical$run_met_pipeline)))
if (have_snapshot) stopifnot(identical(body(run_met_sommer_cs), body(before$run_met_sommer_cs)))
# The historical call uses removed sommer arguments. Only its package API adapter
# is retained; use the same adapter in the reference environment for parity.
historical$run_met_sommer_cs <- run_met_sommer_cs
cat('PASS historical GxE functions and preserved Mean Comparison functions match source\n')

set.seed(926)
trial <- expand.grid(Variety = LETTERS[1:6], Environment = paste0('E', 1:3), Rep = 1:3)
g <- match(trial$Variety, LETTERS); e <- match(trial$Environment, paste0('E', 1:3))
trial$Yield <- 25 + 2*g + 3*e + .5*g*e + rnorm(nrow(trial), sd = .7)
sparse <- subset(trial, !(Variety == 'F' & Environment == 'E3'))
pdf(NULL)
for (engine in c('LMM', 'ANOVA_RCBD', 'SOMMER')) {
  # Mean Comparison must retain its complete pre-restoration output.
  old_mean <- if (have_snapshot) before$run_met_pipeline(sparse, 'Yield', replication_col = 'Rep',
    model_type = engine, analysis_mode = 'means') else NULL
  new_mean <- run_met_pipeline(sparse, 'Yield', replication_col = 'Rep',
    model_type = engine, analysis_mode = 'means')
  if (have_snapshot) for (field in c('anova', 'effect_tests', 'simple_tests', 'means', 'pairwise', 'notes'))
    stopifnot(isTRUE(all.equal(new_mean$mean_comparison[[field]], old_mean$mean_comparison[[field]])))
  for (location in c('E1', 'E3')) {
    chart <- met_plot_mean_comparison(new_mean, 'interaction', 'tukey', environment = location)
    stopifnot(all(chart$data$Family == location), inherits(chart$layers[[1]]$geom, 'GeomBoxplot'))
    ggplot2::ggplot_build(chart)
  }
  means <- list(analysis_mode = 'means', settings = data.frame(Model = engine), met_trait_names = 'Yield',
    met_by_trait = list(Yield = new_mean), met_failed_traits = data.frame(Trait = character(), Error = character()))
  wide <- met_combined_comparison_tables(means)$Means
  stopifnot(identical(tail(wide$Levels, 1), 'CV (%)'),
    !any(grepl('Heritability|LSD|HSD', wide$Levels)))
  stopifnot(any(grepl('03c_mc_means', names(build_export_tables('MET', means)))))
  cat('PASS preserved sparse Mean Comparison, environment charts and export:', engine, '\n')

  # Compare restored results against the original commit on the same observations.
  set.seed(19)
  invisible(capture.output(expected <- historical$run_met_pipeline(trial, 'Yield',
    replication_col = 'Rep', model_type = engine)))
  set.seed(19)
  invisible(capture.output(actual <- run_met_pipeline(trial, 'Yield',
    replication_col = 'Rep', model_type = engine, analysis_mode = 'gxe')))
  for (field in c('blups_main', 'blups_environment', 'variance_components', 'lrt_table',
      'fw_results', 'ammi_genotype', 'gge_genotype', 'met_selection', 'model_summary'))
    stopifnot(isTRUE(all.equal(actual[[field]], expected[[field]])))
  stopifnot(actual$analysis_mode == 'gxe', is.null(actual$mean_comparison$means))
  for (name in grep('^p_', names(actual), value = TRUE))
    if (inherits(actual[[name]], 'ggplot')) ggplot2::ggplot_build(actual[[name]])
  gx <- list(analysis_mode = 'gxe', settings = data.frame(Model = engine), met_trait_names = 'Yield',
    met_by_trait = list(Yield = actual), met_failed_traits = data.frame(Trait = character(), Error = character()))
  workbook <- build_export_tables('MET', gx)
  stopifnot(length(workbook) > 5, !any(grepl('03c_mc_|01c_connectivity', names(workbook))))
  path <- tempfile(fileext = '.xlsx'); writexl::write_xlsx(workbook, path)
  stopifnot(file.info(path)$size > 0); unlink(path)
  cat('PASS historical GxE equivalence, charts and export:', engine, '\n')
}
dev.off()
cat('Selective MET restoration passed.\n')
