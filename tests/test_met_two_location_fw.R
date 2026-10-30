suppressMessages(invisible(source('app.R')))
set.seed(3)
d <- expand.grid(Variety=LETTERS[1:5],Environment=paste0('E',1:3),Rep=1:3)
d$Y <- 10+as.integer(factor(d$Variety))+2*as.integer(factor(d$Environment))+rnorm(nrow(d))
for(engine in c('LMM','ANOVA_RCBD','SOMMER')) {
 r <- run_met_pipeline(subset(d,Environment!='E3'),'Y',replication_col='Rep',model_type=engine,analysis_mode='gxe')
 stopifnot(all(is.finite(r$fw_results$Sens)),all(r$fw_results$b_interp=='Two-location response'),
   nrow(r$ammi_genotype)==0,nrow(r$gge_genotype)==0)
 ggplot2::ggplot_build(r$p_fw_regression); ggplot2::ggplot_build(r$p_fw_mean_sens)
 wb <- build_export_tables('MET',list(analysis_mode='gxe',met_trait_names='Y',met_by_trait=list(Y=r)))
 stopifnot(any(grepl('05_fw',names(wb))))
 cat(engine,'two-location FW passed\n')
}
cat('FW minimum two, AMMI/GGE minimum three passed.\n')
