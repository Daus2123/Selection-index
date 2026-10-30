suppressMessages(invisible(source('app.R')))
set.seed(701)
d <- expand.grid(Variety = LETTERS[1:6], Environment = paste0('E', 1:3), Rep = 1:3)
g <- as.integer(factor(d$Variety)); e <- as.integer(factor(d$Environment))
d$Yield <- 20 + 3*g + 2*e + .6*g*e + rnorm(nrow(d), sd = 1)
d <- subset(d, !(Variety == 'F' & Environment == 'E3'))
for (engine in c('LMM', 'SOMMER', 'ANOVA_RCBD')) {
  r <- run_met_all_traits(d, trait_cols = 'Yield', replication_col = 'Rep',
    model_type = engine, analysis_mode = 'means')
  m <- r$met_by_trait$Yield$mean_comparison
  stopifnot(nrow(r$met_by_trait$Yield$met_cleaned_data) == nrow(d))
  if (engine != 'ANOVA_RCBD') {
    stopifnot(all(is.finite(m$means$Adjusted_mean)), all(is.finite(m$means$Adjusted_SE)),
      any(m$means$Evidence_status == 'Untested - predicted'),
      all(m$means$Estimate_type == 'BLUP-based prediction'),
      all(is.finite(m$effect_tests$p.value)),
      all(c('Variety', 'Environment', 'Variety x Environment') %in% m$effect_tests$Effect))
    missing <- subset(m$means, Family == 'E3' & Level == 'F')
    stopifnot(nrow(missing) > 0, all(is.na(missing$Raw_mean)), all(is.na(missing$N)))
    pair <- subset(m$pairwise, is.finite(estimate))
    stopifnot(nrow(pair) > 0, all(is.finite(pair$SE)), all(pair$SE > 0))
  } else stopifnot(any(!m$means$Estimable))
  wb <- build_export_tables('MET', r)
  stopifnot(any(grepl('03c_mc_means', names(wb))))
  chart <- met_plot_mean_comparison(r$met_by_trait$Yield, 'interaction', 'tukey')
  ggplot2::ggplot_build(chart)
  stopifnot(inherits(chart$layers[[1]]$geom, 'GeomBoxplot'),
    nrow(chart$layers[[1]]$data) == nrow(d),
    nrow(chart$layers[[2]]$data) == 1)
  if (engine == 'LMM') ggplot2::ggsave('tmp/met_boxplot_preview.png', chart, width = 12, height = 6)
  cat('PASS sparse means, exports and charts:', engine, '\n')
}
# Prediction covariance must agree with direct mixed-model equations, including
# the extra interaction variance of an untested cell.
dat <- make_met_data(d, 'Yield', replication_col = 'Rep')$data
b <- met_comparison_blup_grid(dat, TRUE, FALSE)
v <- b$grid@V
stopifnot(max(abs(v - t(v))) < 1e-8, min(eigen(v, symmetric = TRUE)$values) > -1e-8)
grid <- b$grid@grid
stopifnot(stats::nobs(b$model) == nrow(dat))
missing_id <- which(grid$Genotype == 'F' & grid$Environment == 'E3')
fixed_design <- model.matrix(~ Environment, grid,
  contrasts.arg = attr(lme4::getME(b$model, 'X'), 'contrasts'))
expected <- as.numeric(fixed_design[missing_id, , drop = FALSE] %*% lme4::fixef(b$model)) +
  lme4::ranef(b$model)$Genotype['F', 1]
# Missing GxE uses the fitted environment and genotype effects, never a fake plot.
stopifnot(abs(b$grid@bhat[missing_id] - expected) < 1e-8)
raw <- dat %>% dplyr::group_by(Genotype, Environment) %>%
  dplyr::summarise(n = dplyr::n(), .groups = 'drop')
stopifnot(nrow(grid) == 18, nrow(raw) == 17)
# Disconnected trials must not produce a pooled BLUP table.
dis <- subset(d, (Variety %in% c('A','B')) == (Environment == 'E1'))
rr <- run_met_pipeline(dis, 'Yield', replication_col = 'Rep', analysis_mode = 'means')
stopifnot(!any(rr$mean_comparison$means$Estimate_type == 'BLUP-based prediction'),
  any(grepl('Disconnected', rr$mean_comparison$notes$Note)))
cat('PASS prediction covariance and disconnected safeguards\n')
full <- met_build_mean_comparison(dat, TRUE, FALSE, methods = c('tukey', 'lsd'))
for (method in c('tukey', 'lsd')) {
  one <- met_build_mean_comparison(dat, TRUE, FALSE, views = 'interaction', methods = method)
  expected <- subset(full$pairwise, View == 'interaction' & Method == method)
  stopifnot(isTRUE(all.equal(one$pairwise, expected, check.attributes = FALSE)))
}
cat('PASS selected prediction contrasts equal workbook slices\n')
