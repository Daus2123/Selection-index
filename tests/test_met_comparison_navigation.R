suppressMessages(invisible(source('app.R')))
set.seed(81)
d <- expand.grid(Variety=LETTERS[1:4],Environment=paste0('E',1:3),Rep=1:3)
d$Yield <- 10 + as.integer(factor(d$Variety)) + rnorm(nrow(d))
for(engine in c('LMM','SOMMER','ANOVA_RCBD')) {
 for(sparse in c(FALSE,TRUE)) {
  trial <- if(sparse) subset(d,!(Variety=='D' & Environment=='E3')) else d
  r <- run_met_all_traits(trial,trait_cols='Yield',replication_col='Rep',model_type=engine,analysis_mode='means')
  m <- r$met_by_trait$Yield$mean_comparison
  stopifnot(setequal(m$means$View,unname(met_mean_comparison_views())),
    setequal(m$means$Method,c('tukey','lsd')),nrow(m$means)==38,
    nrow(r$settings)==1,nrow(r$met_integrated_ranking)==0)
  if(sparse) stopifnot(if(engine == 'ANOVA_RCBD') any(!m$means$Estimable) else all(m$means$Estimable))
  wb <- build_export_tables('MET',r)
  stopifnot(sum(grepl('03c_mc_',names(wb)))==6)
 }
}
# Reading tables/charts must never refit a model.
met_build_mean_comparison <- function(...) stop('Unexpected refit')
shiny::testServer(server, {
 analysis_results(r); analysis_used('MET'); saved_results$MET <- r
 session$setInputs(result_module='met',result_view='met_means',met_result_trait='Yield',
   chart_module='met',plot_view='met_mean_plot',met_plot_trait='Yield',chart_width=12,chart_height=7)
 stopifnot(grepl('Protected LSD',output$met_mean_controls$html,fixed=TRUE),
   grepl('met_mean_chart_method',output$chart_met_detail$html,fixed=TRUE))
 for(view in unname(met_mean_comparison_views())) {
  for(method in c('tukey','lsd')) {
   session$setInputs(met_mean_view=view,met_mean_method=method,
     met_mean_chart_view=view,met_mean_chart_method=method)
   stopifnot(nrow(selected_result_tables()$tables$Means)>0,
     inherits(selected_chart()$plot,'ggplot'),
     identical(met_completed_method(TRUE),method))
   for(display in c('adjusted_se','observed_sd')) {
    session$setInputs(met_mean_display=display,met_mean_chart_display=display)
    plot <- selected_chart()$plot
    if(is.data.frame(plot$data)) stopifnot(all(plot$data$Method==method),all(plot$data$View==view))
    stopifnot(any(grepl(if(display=='observed_sd') 'Raw_mean_SD' else 'Adjusted_mean_SE',names(selected_result_tables()$tables$Means))))
    invisible(output$met_mean_table); invisible(output$met_mean_combined)
   }
  }
 }
 # Results and Charts selectors are independent.
 session$setInputs(met_mean_view='genotypes',met_mean_method='lsd',
   met_mean_chart_view='interaction',met_mean_chart_method='tukey')
 stopifnot(met_completed_view()=='genotypes',met_completed_view(TRUE)=='interaction',
   met_completed_method()=='lsd',met_completed_method(TRUE)=='tukey')
 path <- tempfile(fileext='.xlsx'); writexl::write_xlsx(d,path)
 session$setInputs(excel_file=list(name='synthetic.xlsx',size=file.info(path)$size,type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',datapath=path),analysis_method='MET',met_analysis_mode='means')
 controls <- output$check_variety_inputs$html
 stopifnot(grepl('met_analysis_mode',controls,fixed=TRUE),
   !grepl('met_run_comparison_view',controls,fixed=TRUE),!grepl('met_run_comparison_method',controls,fixed=TRUE))
 unlink(path)
})
cat('All-view/method mean-only execution, settings, exports and display-only selectors passed.\n')
