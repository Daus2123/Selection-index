suppressMessages(invisible(source('app.R')))
set.seed(882)
d <- expand.grid(Variety=LETTERS[1:5],Environment=paste0('E',1:3),Rep=1:3)
g <- match(d$Variety,LETTERS); e <- match(d$Environment,paste0('E',1:3))
d$Yield <- 20+g+e+.3*g*e+rnorm(nrow(d),sd=.5)
d$Height <- 30+2*g+e+.2*g*e+rnorm(nrow(d),sd=.5)
for(engine in c('LMM','SOMMER','ANOVA_RCBD')) {
 for(kind in c('balanced','sparse','disconnected')) {
  trial <- switch(kind,balanced=d,sparse=subset(d,!(Variety=='E' & Environment!='E1')),
    disconnected=subset(d,(Variety %in% c('A','B'))==(Environment=='E1')))
  r <- run_met_all_traits(trial,trait_cols=c('Yield','Height'),replication_col='Rep',model_type=engine,
    analysis_mode='gxe',check_varieties=c('A','MissingCheck'))
  for(trait in r$met_trait_names) {
   result <- r$met_by_trait[[trait]]
   cells <- result$blups_environment
   stopifnot(nrow(cells)==15,nrow(result$trial_connectivity)==3,nrow(result$environment_links)==3,
     all(c('Confidence_status','Confidence_reason','Uncertainty_Method','Prediction_SE') %in% names(cells)),
     nrow(result$met_cleaned_data)==nrow(trial),is.null(result$mean_comparison$means))
   if(kind=='sparse') stopifnot(any(cells$Confidence_status=='Low confidence'),
     sum(cells$Cell_Status=='UNTESTED')==2)
   if(kind=='disconnected') stopifnot(length(unique(result$trial_connectivity$Connected_group))==2,
     any(cells$Confidence_status=='Unsupported'),all(is.na(cells$Reportable_Estimate[cells$Confidence_status=='Unsupported'])))
   if(kind=='sparse' && engine=='ANOVA_RCBD') stopifnot(nrow(result$ammi_genotype)==0,
     all(is.na(cells$BLUP_env[cells$Cell_Status=='UNTESTED'])))
   ggplot2::ggplot_build(result$p_perf_heatmap)
   stopifnot(all(c('Confidence_status','Confidence_rule') %in% names(met_location_estimates_table(result))))
  }
  if(kind=='disconnected' || (kind=='sparse' && engine=='ANOVA_RCBD')) {
   stopifnot(nrow(r$met_integrated_ranking)==0,
     !any(is.finite(r$met_multitrait_methods$mgidi_fa$ranking$MGIDI_distance)))
  } else {
   stopifnot(nrow(r$met_integrated_ranking)==5)
   board <- build_met_decision_board(trial,r$met_by_trait,r$met_integrated_ranking,
     method_membership=r$met_multitrait_methods$method_membership,check_genotypes='A')
   if(kind=='sparse') stopifnot(board$board$Prediction_evidence_gate[board$board$Genotype=='E']=='RETEST',
     board$board$Final_action[board$board$Genotype=='E']=='RETEST')
  }
  wb <- build_export_tables('MET',r)
  stopifnot(sum(grepl('01c_connectivity',names(wb)))==2,sum(grepl('01d_environment_links',names(wb)))==2)
  path <- tempfile(fileext='.xlsx');writexl::write_xlsx(wb,path);unlink(path)
  cat('PASS',engine,kind,'\n')
 }
 # Separate inference mode must never receive artificial prediction rows.
 trial <- subset(d,!(Variety=='E' & Environment!='E1'))
 m <- run_met_all_traits(trial,trait_cols='Yield',replication_col='Rep',model_type=engine,analysis_mode='means',check_varieties='A')
 stopifnot(nrow(m$met_by_trait$Yield$met_cleaned_data)==nrow(trial),
   sum(m$met_by_trait$Yield$mean_comparison$means$N,na.rm=TRUE)>0,
   is.null(m$met_by_trait$Yield$blups_environment),nrow(m$met_integrated_ranking)==0)
}
cat('Complete MET workflow: evidence, all models, modes, missing cells, connectivity, decisions and exports passed.\n')
