suppressMessages(invisible(source('app.R')))
set.seed(615)
d <- expand.grid(Variety=LETTERS[1:4],Environment=paste0('E',1:3),Rep=1:3)
g <- match(d$Variety,LETTERS); e <- match(d$Environment,paste0('E',1:3))
d$Yield <- 20+g+e+g*e+rnorm(nrow(d),sd=.2)
d$Height <- 40+2*g+e+2*g*e+rnorm(nrow(d),sd=.2)
r <- run_met_all_traits(d,trait_cols=c('Yield','Height'),replication_col='Rep',analysis_mode='means')
met_build_mean_comparison <- function(...) stop('Unexpected model refit')
shiny::testServer(server, {
 analysis_results(r); analysis_used('MET'); saved_results$MET <- r
 session$setInputs(result_module='met',result_view='met_means',met_result_trait='Yield',chart_module='met',
   plot_view='met_mean_plot',met_plot_trait='Yield',chart_width=12,chart_height=7)
 stopifnot(identical(unname(met_mean_comparison_views()),c('genotypes','environments','interaction')),
   !grepl('met_mean_chart_display',output$chart_met_detail$html,fixed=TRUE))
 corner <- output$met_mean_table$html
 stopifnot(grepl('Observe',corner,fixed=TRUE),grepl('Adjusted',corner,fixed=TRUE),
   grepl('met_mean_combined',corner,fixed=TRUE),!grepl('met_mean_family_',corner,fixed=TRUE),
   !grepl('met_mean_display',output$met_mean_controls$html,fixed=TRUE))
 for(view in unname(met_mean_comparison_views())) for(method in c('tukey','lsd')) {
  session$setInputs(met_mean_view=view,met_mean_method=method,met_mean_chart_view=view,met_mean_chart_method=method)
  for(display in c('observed_sd','adjusted_se')) {
   session$setInputs(met_mean_display=display,met_mean_chart_display='observed_sd')
   tables <- selected_result_tables()$tables
   stopifnot(setequal(names(tables),c('Means','ANOVA','Effects')),
     nrow(tables$Effects)==6, all(c('Yield','Height') %in% tables$ANOVA$Trait),
     sum(grepl(if(display=='observed_sd') 'Raw_mean_SD' else 'Adjusted_mean_SE',names(tables$Means)))==2)
   invisible(output$met_mean_combined); invisible(output$met_mean_anova); invisible(output$met_mean_tests)
   p <- selected_chart()$plot
   stopifnot(inherits(p$layers[[1]]$geom,'GeomCol'),identical(p$data$Plot_mean,p$data$Adjusted_mean),
     identical(p$data$Plot_error,p$data$Adjusted_SE))
  }
 }
 session$setInputs(met_mean_chart_view='interaction')
 ggplot2::ggsave('tmp/met_support_bar_preview.png',selected_chart()$plot,width=10,height=5)
})
cat('Support MET single-table corner slicer, effect tests, all-trait downloads and adjusted-only bars passed.\n')
