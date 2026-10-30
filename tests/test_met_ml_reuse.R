suppressMessages(invisible(source('app.R')))
set.seed(25)
d <- expand.grid(Variety=LETTERS[1:5],Environment=paste0('E',1:3),Rep=1:3)
d$Y <- 10+as.integer(factor(d$Variety))+as.integer(factor(d$Environment))+rnorm(nrow(d))
r <- run_met_pipeline(d,'Y',replication_col='Rep',analysis_mode='gxe',model_type='LMM')
x <- r$met_cleaned_data
full <- lme4::lmer(Weight~(1|Environment)+(1|Environment:Rep)+(1|Genotype)+(1|Genotype:Environment),data=x,control=lme4::lmerControl(optimizer='bobyqa'))
reduced <- lme4::lmer(Weight~(1|Environment)+(1|Environment:Rep)+(1|Genotype),data=x,control=lme4::lmerControl(optimizer='bobyqa'))
old <- as.data.frame(anova(full,reduced))
new <- subset(r$lrt_table,Test=='GxE')
stopifnot(abs(max(old$Chisq,na.rm=TRUE)-max(new$Chisq,na.rm=TRUE))<1e-4,
  all(is.finite(r$blups_main$SE_G)))
cat('ML reuse matches REML-to-ML reference LRT; REML prediction uncertainty retained.\n')
