suppressMessages(invisible(source("app.R")))
set.seed(531)
d <- expand.grid(Genotype = LETTERS[1:4], Environment = paste0("E", 1:3), Rep = 1:4)
d$Genotype <- factor(d$Genotype); d$Environment <- factor(d$Environment); d$Rep <- factor(d$Rep)
rep_effect <- rnorm(12, sd = 0.6)
d$Weight <- 10 + as.integer(d$Genotype) * 1.5 + as.integer(d$Environment) +
  rep_effect[as.integer(interaction(d$Environment, d$Rep))] + rnorm(nrow(d), sd = 0.4)
result <- met_build_mean_comparison(d, TRUE, FALSE)
stopifnot(nrow(result$means) == 62,
  setequal(result$means$View, unname(met_mean_comparison_views())),
  all(result$means$Estimable), all(is.finite(result$means$Adjusted_SE)))
# Balanced design: fixed cell means reproduce observed cell means.
within <- subset(result$means, View == "within_environment" & Method == "tukey")
for (i in seq_len(nrow(within))) {
  raw <- subset(d, Genotype == within$Level[i] & Environment == within$Family[i])$Weight
  stopifnot(abs(within$Adjusted_mean[i] - mean(raw)) < 1e-6,
    abs(within$Raw_mean[i] - mean(raw)) < 1e-10,
    abs(within$Observed_SD[i] - sd(raw)) < 1e-10)
}
pair <- subset(result$pairwise, Method == "lsd" & Protected_test_passed)
stopifnot(all(abs(pair$Critical_difference - qt(0.975, pair$df) * pair$SE) < 1e-10))
# Remove treatment effects exactly: protected LSD must suppress significance/letters.
null <- d
noise <- c(-0.3, 0.1, 0.4, -0.2)
null$Weight <- 10 + as.integer(null$Rep) + noise[as.integer(null$Genotype)] *
  c(-1, 1, -1, 1)[as.integer(null$Rep)]
n <- met_build_mean_comparison(null, TRUE, FALSE)
blocked <- subset(n$pairwise, Method == "lsd" & !Protected_test_passed)
stopifnot(nrow(blocked) > 0, all(is.na(blocked$p.value)), all(is.na(blocked$Significant_0.05)))
stopifnot(all(subset(n$means, Method == "lsd" & Omnibus_p >= 0.05)$Group == ""))
# Missing cells remain visible and are explicitly noted.
sparse <- subset(d, !(Genotype == "D" & Environment == "E3"))
s <- met_build_mean_comparison(sparse, TRUE, FALSE)
stopifnot(any(!s$means$Estimable), any(grepl("missing", s$notes$Note)))
for (view in unname(met_mean_comparison_views())) {
  for (display in c("adjusted_se", "observed_sd")) {
    plot <- met_plot_mean_comparison(list(mean_comparison = result), view, display = display)
    built <- ggplot2::ggplot_build(plot)
    if (!is.data.frame(plot$data)) next
    stopifnot(nrow(built$data[[1]]) > 0,
      isTRUE(all.equal(plot$data$Plot_mean, plot$data$Raw_mean)),
      isTRUE(all.equal(plot$data$Plot_error, plot$data$Observed_SD)))
  }
}
# Full pipeline and workbook contract under each existing engine.
trial <- data.frame(Variety = d$Genotype, Environment = d$Environment, Rep = d$Rep, Yield = d$Weight)
for (engine in c("LMM", "ANOVA_RCBD", "SOMMER")) {
  r <- run_met_pipeline(trial, "Yield", replication_col = "Rep", model_type = engine)
  stopifnot(nrow(r$mean_comparison$means) == 62)
  if (engine == "ANOVA_RCBD") {
    stopifnot(any(grepl("fitted MET ANOVA model", r$mean_comparison$notes$Note)))
    cleaned <- r$met_cleaned_data
    main <- lm(Weight ~ Environment + Environment:Rep + Genotype + Genotype:Environment, data = cleaned)
    expected <- as.data.frame(emmeans::emmeans(main, ~ Genotype | Environment))
    actual <- subset(r$mean_comparison$means, View == "within_environment" & Method == "tukey")
    key <- paste(expected$Genotype, expected$Environment)
    idx <- match(paste(actual$Level, actual$Family), key)
    stopifnot(max(abs(actual$Adjusted_mean - expected$emmean[idx])) < 1e-8,
      max(abs(actual$Adjusted_SE - expected$SE[idx])) < 1e-8)
  }
  workbook <- build_export_tables("MET", list(met_trait_names = "Yield", met_by_trait = list(Yield = r)))
  stopifnot(sum(grepl("03c_mc_", names(workbook))) == 5)
  path <- tempfile(fileext = ".xlsx")
  writexl::write_xlsx(workbook, path)
  stopifnot(file.info(path)$size > 0)
  unlink(path)
}
cat("MET mean comparisons, protection, missing cells, charts and all-engine exports passed.\n")
