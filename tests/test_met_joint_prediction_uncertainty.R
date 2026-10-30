suppressMessages(invisible(source('app.R')))
set.seed(741)
d <- expand.grid(Genotype=factor(LETTERS[1:4]),Environment=factor(paste0('E',1:3)),Rep=factor(1:4))
d <- subset(d,!(Genotype=='D' & Environment=='E3'))
d$Weight <- 20 + rep(rnorm(4,sd=2),length.out=nrow(d)) + 2*as.integer(d$Environment) + rnorm(nrow(d),sd=1)
m <- lme4::lmer(Weight~1+(1|Genotype)+(1|Environment)+(1|Genotype:Environment)+(1|Environment:Rep),data=d)
a <- met_lmm_prediction_uncertainty(m,d,as.numeric(lme4::VarCorr(m)$`Genotype:Environment`))
# Independent observation-space universal-kriging calculation.
v <- diag(sigma(m)^2,nrow(d)); target <- a$cells; cross <- matrix(0,nrow(target),nrow(d)); prior <- rep(0,nrow(target))
for(term in names(lme4::VarCorr(m))) {
 cols <- strsplit(term,':',fixed=TRUE)[[1]]
 key <- function(data) do.call(paste,c(lapply(data[cols],as.character),sep='|'))
 obs <- key(d); tau <- as.numeric(lme4::VarCorr(m)[[term]])
 v <- v + tau*outer(obs,obs,'==')
 if(all(cols %in% c('Genotype','Environment'))) {
  cross <- cross + tau*outer(key(target),obs,'=='); prior <- prior+tau
 }
}
vi <- solve(v); x <- matrix(1,nrow(d),1); bc <- solve(t(x)%*%vi%*%x)
rem <- matrix(1,nrow(target),1)-cross%*%vi%*%x
reference <- prior-diag(cross%*%vi%*%t(cross))+diag(rem%*%bc%*%t(rem))
stopifnot(max(abs(a$cells$Prediction_SE^2-reference))<1e-7)
cat('Joint mixed-model prediction SE agrees with independent universal-kriging calculation.\n')
