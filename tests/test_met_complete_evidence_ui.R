suppressMessages(invisible(source('app.R')))
set.seed(92)
d <- expand.grid(Variety=LETTERS[1:4],Environment=paste0('E',1:3),Rep=1:3)
d$Y <- 20+as.integer(factor(d$Variety))+rnorm(nrow(d))
d <- subset(d,!(Variety=='D' & Environment=='E3'))
r <- run_met_all_traits(d,trait_cols='Y',replication_col='Rep',analysis_mode='gxe',check_varieties='A')
shiny::testServer(server,{
 analysis_results(r);analysis_used('MET');saved_results$MET<-r
 session$setInputs(result_module='met',result_view='met_qc',met_result_trait='Y',chart_module='met',plot_view='met_performance_heatmap',met_plot_trait='Y',met_performance_heatmap_view='estimated',chart_width=12,chart_height=7)
 stopifnot(nrow(selected_result_tables()$tables$Trial_connectivity)==3,
   nrow(selected_result_tables()$tables$Environment_links)==3)
 invisible(output$met_trial_connectivity);invisible(output$met_environment_links)
 session$setInputs(result_view='met_env_blup')
 stopifnot('Confidence_status' %in% names(selected_result_tables()$tables$Genotype_by_location))
 ggplot2::ggsave('tmp/met_complete_evidence_heatmap.png',selected_chart()$plot,width=10,height=6)
})
cat('Complete MET evidence UI, selected downloads and heatmap passed.\n')
