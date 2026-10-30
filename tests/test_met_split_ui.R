suppressMessages(invisible(source('app.R')))
set.seed(345)
d <- expand.grid(Variety=LETTERS[1:4], Environment=c('E1','E2'), Rep=1:3)
d$Yield <- 10 + as.integer(factor(d$Variety)) + rnorm(nrow(d))
r <- run_met_all_traits(d,trait_cols='Yield',replication_col='Rep',analysis_mode='means')
shiny::testServer(server, {
 analysis_results(r); analysis_used('MET'); saved_results$MET <- r
 session$setInputs(result_module='met',result_view='met_means',met_result_trait='Yield',
   chart_module='met',plot_view='met_mean_plot',met_plot_trait='Yield',chart_width=12,chart_height=7,met_mean_method='lsd',met_mean_chart_method='lsd')
 stopifnot(grepl('Compare',output$met_mean_controls$html,fixed=TRUE),
   grepl('Observe',output$met_mean_table$html,fixed=TRUE),
   !grepl('Decision Board',output$result_met_detail$html,fixed=TRUE),
   !grepl('AMMI',output$chart_met_detail$html,fixed=TRUE),
   grepl('Protected LSD',selected_chart()$plot$labels$subtitle,fixed=TRUE))
 for(display in c('adjusted_se','observed_sd')) {
  session$setInputs(met_mean_display=display,met_mean_chart_display=display)
  p <- selected_chart()$plot
  stopifnot(p$labels$y=='Adjusted mean \u00b1 SE',
    any(grepl(if(display=='observed_sd') 'Raw_mean_SD' else 'Adjusted_mean_SE',names(selected_result_tables()$tables$Means))))
  invisible(output$met_mean_table); invisible(output$met_mean_combined)
  stopifnot(length(output$met_mean_plot$src)==1)
 }
 # Pending Analyze changes cannot relabel an already completed result.
 session$setInputs(met_run_comparison_method='tukey',met_run_comparison_view='interaction')
 stopifnot(grepl('Protected LSD',selected_chart()$plot$labels$subtitle,fixed=TRUE),
   identical(met_completed_view(),'genotypes'))
 # Read the actual Analyze controls using a synthetic uploaded workbook.
 path <- tempfile(fileext='.xlsx'); writexl::write_xlsx(d,path)
 session$setInputs(excel_file=list(name='synthetic.xlsx',size=file.info(path)$size,type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',datapath=path),
   analysis_method='MET',met_analysis_mode='means')
 controls <- output$check_variety_inputs$html
 stopifnot(regexpr('met_rep_col',controls)[1] < regexpr('met_analysis_mode',controls)[1],
   regexpr('met_analysis_mode',controls)[1] < regexpr('met_model_type',controls)[1],
   !grepl('met_run_comparison_method',controls,fixed=TRUE))
 g <- r; g$analysis_mode <- 'gxe'; analysis_results(g)
 session$setInputs(result_view='met_summary',plot_view='met_env_cor')
 stopifnot(!grepl('met_means',output$result_met_detail$html,fixed=TRUE),
   !grepl('met_mean_plot',output$chart_met_detail$html,fixed=TRUE))
 unlink(path)
})
cat('Split MET Shiny controls, mode navigation, selected results and displays passed.\n')
