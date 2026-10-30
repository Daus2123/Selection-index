suppressMessages(invisible(source('app.R')))
set.seed(810)
d <- expand.grid(Variety = LETTERS[1:5], Environment = paste0('E', 1:3), Rep = 1:3)
d$Yield <- 10 + as.integer(factor(d$Variety)) + as.integer(factor(d$Environment)) + rnorm(nrow(d), sd=.5)
# Prove mean-only mode never calls GxE model fitting, extensions or ranking.
original_lmer <- lmer
lmer <- function(...) stop('Unexpected GxE fit')
original_extensions <- run_met_single_trait_extensions
run_met_single_trait_extensions <- function(...) stop('Unexpected extensions')
original_integrated <- build_met_integrated_ranking
build_met_integrated_ranking <- function(...) stop('Unexpected ranking')
for (engine in c('LMM','ANOVA_RCBD','SOMMER')) {
 for (sparse in c(FALSE, TRUE)) {
  trial <- if (sparse) subset(d, !(Variety == 'E' & Environment == 'E3')) else d
  for (method in c('tukey','lsd')) {
   start <- proc.time()[['elapsed']]
   r <- run_met_all_traits(trial, trait_cols='Yield', replication_col='Rep', model_type=engine,
     analysis_mode='means', comparison_view='interaction', comparison_method=method)
   m <- r$met_by_trait$Yield$mean_comparison
   stopifnot(r$analysis_mode=='means', setequal(m$means$Method, method),
     setequal(m$means$View,'interaction'), nrow(m$means)==15,
     nrow(r$met_integrated_ranking)==0, is.null(r$met_by_trait$Yield$blups_main),
     all(c('Raw_mean','Observed_SD','Adjusted_mean','Adjusted_SE') %in% names(m$means)))
   if (sparse) stopifnot(if(engine == 'ANOVA_RCBD') any(!m$means$Estimable) else all(m$means$Estimable))
   for (display in c('adjusted_se','observed_sd')) {
    p <- met_plot_mean_comparison(r$met_by_trait$Yield,'interaction',method,display)
    ggplot2::ggplot_build(p)
    stopifnot(isTRUE(all.equal(p$data$Plot_mean, p$data$Adjusted_mean)))
    tab <- met_combined_comparison_tables(r,'interaction',method,display)$Means
    stopifnot(any(grepl(if(display=='observed_sd') 'Raw_mean_SD' else 'Adjusted_mean_SE', names(tab))))
   }
   wb <- build_export_tables('MET',r)
   stopifnot(sum(grepl('03c_mc_',names(wb)))==6, !any(grepl('05_fw|06_ammi|09_overall',names(wb))))
   path <- tempfile(fileext='.xlsx'); writexl::write_xlsx(wb,path); unlink(path)
   cat('MEANS',engine,'sparse=',sparse,method,round(proc.time()[['elapsed']]-start,2),'seconds PASS\n')
  }
 }
}
lmer <- original_lmer
run_met_single_trait_extensions <- original_extensions
build_met_integrated_ranking <- original_integrated
# Prove GxE mode never invokes the comparison model.
original_comparison <- met_build_mean_comparison
met_build_mean_comparison <- function(...) stop('Unexpected comparison')
for (engine in c('LMM','ANOVA_RCBD','SOMMER')) {
 for (sparse in c(FALSE,TRUE)) {
  trial <- if(sparse) subset(d, !(Variety=='E' & Environment=='E3')) else d
  r <- run_met_pipeline(trial,'Yield',replication_col='Rep',model_type=engine,analysis_mode='gxe')
  stopifnot(is.null(r$mean_comparison$means), inherits(r$p_perf_heatmap,'ggplot'),
    inherits(r$p_perf_heatmap_observed,'ggplot'), nrow(r$blups_main)==5,
    !any(grepl('Unexpected comparison',r$model_summary$Notes)))
  cat('GXE',engine,'sparse=',sparse,'PASS\n')
 }
}
met_build_mean_comparison <- original_comparison
# Selected family/method must equal the same slice of a complete comparison run.
dat <- make_met_data(d,'Yield',replication_col='Rep')$data
full <- met_build_mean_comparison(dat,TRUE,FALSE,methods=c('tukey','lsd'))
for(method in c('tukey','lsd')) {
 one <- met_build_mean_comparison(dat,TRUE,FALSE,views='interaction',methods=method)
 expected <- subset(full$means, View=='interaction' & Method==method)
 stopifnot(isTRUE(all.equal(one$means,expected,check.attributes=FALSE)))
 expected <- subset(full$pairwise, View=='interaction' & Method==method)
 stopifnot(isTRUE(all.equal(one$pairwise,expected,check.attributes=FALSE)))
}
cat('Split MET modes, selected inference, heatmaps, all engines and exports passed.\n')
