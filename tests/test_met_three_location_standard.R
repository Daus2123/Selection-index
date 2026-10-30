suppressMessages(invisible(source('app.R')))
set.seed(3)
d <- expand.grid(Variety=LETTERS[1:5],Environment=paste0('E',1:3),Rep=1:3)
d$Y <- 10+as.integer(factor(d$Variety))+rnorm(nrow(d))
for(engine in c('LMM','ANOVA_RCBD','SOMMER')) {
 for(nenv in c(2,3)) {
  trial <- subset(d,Environment %in% paste0('E',seq_len(nenv)))
  r <- run_met_pipeline(trial,'Y',replication_col='Rep',model_type=engine,analysis_mode='gxe',min_envs_for_biplot=5)
  stopifnot(r$model_summary$AMMI_GGE_min_observed_locations==3,
    r$ammi_notes$Min_locations_required==3)
  if(nenv==2) stopifnot(nrow(r$ammi_genotype)==0,nrow(r$gge_genotype)==0) else
    stopifnot(nrow(r$ammi_genotype)>0,nrow(r$gge_genotype)>0)
 }
}
cat('Fixed three-location standard passed for LMM, RCBD and SOMMER; two locations skip AMMI/GGE.\n')
