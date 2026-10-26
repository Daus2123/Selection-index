suppressMessages(invisible(source("app.R")))
set.seed(274)
d <- expand.grid(ID = as.character(1:12), Rep = as.character(1:6), stringsAsFactors = FALSE)
genotype <- rnorm(12, sd = 3)
replication <- rnorm(6, sd = 1.5)
d$Yield <- 30 + genotype[as.integer(d$ID)] + replication[as.integer(d$Rep)] + rnorm(nrow(d), sd = 1)
d$ID <- factor(d$ID); d$Rep <- factor(d$Rep)
stats <- si_lpsi_lmm_statistics(d, "Yield", .1)
comparison <- fit_single_location_model(d, "Yield", "LMM")
genetic <- lme4::lmer(Yield ~ 1 + (1 | ID) + (1 | Rep), d, REML = TRUE,
  control = lme4::lmerControl(optimizer = "bobyqa"))
vg <- as.numeric(lme4::VarCorr(genetic)$ID)
ve <- sigma(genetic)^2
# For a balanced RCBD, generalized H2 equals the standard entry-mean expression.
stopifnot(stats$CV_pct == round(100 * sigma(comparison) / abs(mean(d$Yield)), 2),
  abs(stats$Broad_sense_H2 - vg / (vg + ve / 6)) < 1e-4,
  is.finite(stats$Mean_PEV_difference), is.na(stats$Genetic_advance))
# Independent dense MME calculation for an unbalanced trial.
u <- d[-c(1, 3, 7, 20, 32, 56), ]
us <- si_lpsi_lmm_statistics(u, "Yield", .1)
fit <- lme4::lmer(Yield ~ 1 + (1 | ID) + (1 | Rep), u, REML = TRUE,
  control = lme4::lmerControl(optimizer = "bobyqa"))
x <- model.matrix(~ 1, u)
zg <- model.matrix(~ 0 + ID, u); zr <- model.matrix(~ 0 + Rep, u)
z <- cbind(zg, zr)
vg <- as.numeric(lme4::VarCorr(fit)$ID); vr <- as.numeric(lme4::VarCorr(fit)$Rep)
ve <- sigma(fit)^2
design <- cbind(x, z)
precision <- crossprod(design) / ve + diag(c(0, rep(1 / vg, ncol(zg)), rep(1 / vr, ncol(zr))))
pev <- solve(precision)[1 + seq_len(ncol(zg)), 1 + seq_len(ncol(zg))]
pair_variance <- outer(diag(pev), diag(pev), "+") - 2 * pev
expected <- 1 - mean(pair_variance[upper.tri(pair_variance)]) / (2 * vg)
stopifnot(abs(us$Broad_sense_H2 - expected) < 1e-4)
# Units and genotype order do not change CV or H2.
scaled <- u; scaled$Yield <- scaled$Yield * 1000
ss <- si_lpsi_lmm_statistics(scaled, "Yield", .1)
reordered <- u[rev(seq_len(nrow(u))), ]; reordered$ID <- factor(reordered$ID, levels = rev(levels(u$ID)))
rs <- si_lpsi_lmm_statistics(reordered, "Yield", .1)
stopifnot(abs(ss$Broad_sense_H2 - us$Broad_sense_H2) < 1e-4,
  abs(ss$CV_pct - us$CV_pct) < .02, abs(rs$Broad_sense_H2 - us$Broad_sense_H2) < 1e-4)
# A zero genetic-variance boundary must not remove the separately estimable CV.
zero <- d
noise <- matrix(rnorm(nrow(zero)), nrow = 12)
noise <- sweep(noise, 1, rowMeans(noise), "-")
noise <- sweep(noise, 2, colMeans(noise), "-")
zero$Yield <- 30 + as.vector(noise)
zs <- si_lpsi_lmm_statistics(zero, "Yield", .1)
stopifnot(is.finite(zs$CV_pct), zs$Broad_sense_H2 == 0,
  zs$H2_status == "Zero genetic variance", grepl("zero boundary", zs$Note))
# Rounded/constant scores must report legitimate zero CV without forcing H2.
constant <- d; constant$Yield <- 5
cs <- si_lpsi_lmm_statistics(constant, "Yield", .1)
stopifnot(cs$CV_pct == 0, is.na(cs$Broad_sense_H2), cs$H2_status == "Constant trait")
repeated <- d; repeated$Yield <- rep(seq(2, 4.75, length.out = 12), 6)
ps <- si_lpsi_lmm_statistics(repeated, "Yield", .1)
stopifnot(ps$CV_pct == 0, is.na(ps$Broad_sense_H2), ps$H2_status == "Zero residual variation")
zero_mean <- d; zero_mean$Yield <- 0
ms <- si_lpsi_lmm_statistics(zero_mean, "Yield", .1)
stopifnot(is.na(ms$CV_pct), ms$CV_status == "Zero mean", is.na(ms$Broad_sense_H2))
trial <- d; names(trial)[names(trial) == "ID"] <- id_col
trial$Constant_score <- 5
trial$Repeated_score <- repeated$Yield
result <- run_selection_pipeline(trial, check_varieties = "1", model_type = "LMM")
view <- si_lpsi_mean_comparison_view(result)
yield_stats <- result$heritability_gain[result$heritability_gain$Trait == "Yield", ]
stopifnot(is.finite(yield_stats$CV_pct), is.finite(yield_stats$Broad_sense_H2),
  grepl("%", view$summary$Yield[1]), view$summary$Yield[2] != "Not estimated",
  grepl("Cullis", view$summary$Original_ID[2]),
  view$summary$Constant_score[1] == "0.00%", view$summary$Constant_score[2] == "Constant trait",
  view$summary$Repeated_score[1] == "0.00%", view$summary$Repeated_score[2] == "Zero residual variation")
tables <- build_export_tables("LPSI", result)
gain <- tables[[which(grepl("heritability_gain", names(tables)))[1]]]
stopifnot(all(is.finite(gain$CV_pct)), is.finite(gain$Broad_sense_H2[gain$Trait == "Yield"]),
  "Heritability_method" %in% names(gain), "H2_status" %in% names(gain))
cat("LMM residual CV, balanced/unbalanced Cullis H2, boundary handling, footer and export checks passed.\n")
