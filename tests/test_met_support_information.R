suppressMessages(invisible(source('app.R')))
set.seed(615)
d <- expand.grid(Variety=LETTERS[1:5],Environment=paste0('E',1:3),Rep=1:4)
g <- match(d$Variety,LETTERS); e <- match(d$Environment,paste0('E',1:3))
d$Yield <- 20 + g + e + g*e + rnorm(nrow(d),sd=.2)
d$Height <- 40 + 2*g + e + 2*g*e + rnorm(nrow(d),sd=.3)
for(engine in c('LMM','ANOVA_RCBD','SOMMER')) {
 for(sparse in c(FALSE,TRUE)) {
  trial <- if(sparse) subset(d,!(Variety=='E' & Environment=='E3')) else d
  r <- run_met_all_traits(trial,trait_cols=c('Yield','Height'),replication_col='Rep',model_type=engine,analysis_mode='means')
  stopifnot(length(r$met_by_trait)==2, nrow(r$met_integrated_ranking)==0)
  for(trait in r$met_trait_names) {
   m <- r$met_by_trait[[trait]]$mean_comparison
   stopifnot(nrow(m$effect_tests)==3,nrow(m$means)==46,
     setequal(m$means$View,c('genotypes','environments','interaction')),
     setequal(m$means$Method,c('lsd','tukey')))
   if(sparse) stopifnot(any(!m$means$Estimable))
   blocked <- subset(m$means,Posthoc_status!='Applied')
   stopifnot(all(blocked$Group==''))
   if(!sparse) stopifnot(all(m$effect_tests$p.value<.05),any(nzchar(m$means$Group)))
   for(view in unname(met_mean_comparison_views())) {
    for(method in c('tukey','lsd')) {
     p <- met_plot_mean_comparison(r$met_by_trait[[trait]],view,method,display='observed_sd')
     if(is.data.frame(p$data)) stopifnot(inherits(p$layers[[1]]$geom,'GeomCol'),
       identical(p$data$Plot_mean,p$data$Adjusted_mean),identical(p$data$Plot_error,p$data$Adjusted_SE))
     ggplot2::ggplot_build(p)
    }
   }
  }
  for(view in unname(met_mean_comparison_views())) {
   for(display in c('observed_sd','adjusted_se')) {
    tab <- met_combined_comparison_tables(r,view,'tukey',display)
    stopifnot(nrow(tab$Effects)==6,all(c('CV (%)','Heritability H\u00b2') %in% tab$Means$Levels),
      any(grepl('Yield',names(tab$Means))),any(grepl('Height',names(tab$Means))))
    if(view=='interaction') stopifnot('Environment' %in% names(tab$Means))
   }
  }
  wb <- build_export_tables('MET',r)
  stopifnot(sum(grepl('03c_mc_',names(wb)))==12, any(grepl('effect_tests',names(wb))))
  path <- tempfile(fileext='.xlsx'); writexl::write_xlsx(wb,path); unlink(path)
  cat('PASS',engine,'sparse=',sparse,'\n')
 }
}
# Exactly balanced zero treatment effects: charts must remain but letters/tests omitted.
null <- expand.grid(Genotype=factor(LETTERS[1:4]),Environment=factor(paste0('E',1:3)),Rep=factor(1:4))
noise <- c(-.3,.1,.4,-.2)
null$Weight <- 10 + as.integer(null$Rep) + noise[as.integer(null$Genotype)] * c(-1,1,-1,1)[as.integer(null$Rep)]
m <- met_build_mean_comparison(null,TRUE,FALSE,methods=c('lsd','tukey'))
stopifnot(all(m$means$Group==''),all(is.na(m$pairwise$p.value)),all(is.na(m$pairwise$Critical_difference)))
for(view in unname(met_mean_comparison_views())) {
 p <- met_plot_mean_comparison(list(mean_comparison=m),view)
 stopifnot(is.data.frame(p$data),any(is.finite(p$data$Plot_mean)),all(p$data$Group==''))
}
# Validate letters against the same p-values, without depending on alphabet choices.
m <- r$met_by_trait$Yield$mean_comparison
for(view in unique(m$means$View)) for(method in c('tukey','lsd')) {
 for(family in unique(subset(m$means,View==view & Method==method)$Family)) {
  a <- subset(m$means,View==view & Method==method & Family==family)
  b <- subset(m$pairwise,View==view & Method==method & Family==family)
  for(k in which(is.finite(b$p.value))) {
   ids <- strsplit(as.character(b$contrast[k]),' - ',fixed=TRUE)[[1]]
   if(length(ids)!=2) next
   letters <- a$Group[match(ids,a$Level)]
   if(any(is.na(letters)) || any(!nzchar(letters))) next
   shared <- length(intersect(strsplit(letters[1],'')[[1]],strsplit(letters[2],'')[[1]]))>0
   stopifnot(shared == (b$p.value[k]>=.05))
  }
 }
}
cat('Support MET ANOVA/effects, significance gates, bar charts, combined tables, all engines and exports passed.\n')
