suppressMessages(invisible(source("app.R")))
set.seed(531)
d <- expand.grid(Genotype = LETTERS[1:4], Environment = paste0("E", 1:3), Rep = 1:4)
d$Genotype <- factor(d$Genotype); d$Environment <- factor(d$Environment); d$Rep <- factor(d$Rep)
rep_effect <- rnorm(12, sd = 0.6)
d$Weight <- 10 + as.integer(d$Genotype) * 1.5 + as.integer(d$Environment) +
  rep_effect[as.integer(interaction(d$Environment, d$Rep))] + rnorm(nrow(d), sd = 0.4)
fixed_model <- lm(Weight ~ Environment + Environment:Rep + Genotype * Environment, data = d)
result <- met_build_mean_comparison(d, TRUE, FALSE, fitted_model = fixed_model)
stopifnot(nrow(result$means) == 19,
  setequal(result$means$View, unname(met_mean_comparison_views())),
  all(result$means$Estimable), all(is.finite(result$means$Adjusted_SE)))
# Balanced design: fixed cell means reproduce observed cell means.
within <- subset(result$means, View == "interaction" & Method == "tukey")
for (i in seq_len(nrow(within))) {
  raw <- subset(d, Genotype == within$Level[i] & Environment == within$Family[i])$Weight
  stopifnot(abs(within$Adjusted_mean[i] - mean(raw)) < 1e-6)
}
stopifnot(all(result$means$Method == "tukey"), all(result$pairwise$Method == "tukey"),
  all(c("Raw_mean", "Observed_SD", "Raw_mean_SD") %in% names(result$means)))
pair <- result$pairwise
valid <- is.finite(pair$Critical_difference)
stopifnot(all(abs(pair$Critical_difference[valid] - (pair$upper.CL[valid] - pair$lower.CL[valid]) / 2) < 1e-10))
# Missing cells remain visible and are explicitly noted.
sparse <- subset(d, !(Genotype == "D" & Environment == "E3"))
s <- met_build_mean_comparison(sparse, TRUE, FALSE)
stopifnot(all(s$means$Estimable), any(grepl("missing", s$notes$Note)),
  any(s$means$Evidence_status == "Untested - predicted"))
for (view in unname(met_mean_comparison_views())) {
  for (display in "adjusted_se") {
    plot <- met_plot_mean_comparison(list(mean_comparison = result), view, display = display)
    built <- ggplot2::ggplot_build(plot)
    if (!is.data.frame(plot$data)) next
    stopifnot(nrow(built$data[[1]]) > 0,
      isTRUE(all.equal(plot$data$Plot_mean, plot$data$Adjusted_mean)),
      isTRUE(all.equal(plot$data$Plot_error, plot$data$Adjusted_SE)))
  }
}
# Full pipeline and workbook contract under each existing engine.
trial <- data.frame(Variety = d$Genotype, Environment = d$Environment, Rep = d$Rep, Yield = d$Weight)
for (engine in c("LMM", "ANOVA_RCBD", "SOMMER")) {
  r <- run_met_pipeline(trial, "Yield", replication_col = "Rep", model_type = engine)
  stopifnot(nrow(r$mean_comparison$means) == 19)
  if (engine == "ANOVA_RCBD") {
    stopifnot(any(grepl("fitted MET ANOVA model", r$mean_comparison$notes$Note)))
    cleaned <- r$met_cleaned_data
    main <- lm(Weight ~ Environment + Environment:Rep + Genotype + Genotype:Environment, data = cleaned)
    expected <- as.data.frame(emmeans::emmeans(main, ~ Genotype | Environment))
    actual <- subset(r$mean_comparison$means, View == "interaction" & Method == "tukey")
    key <- paste(expected$Genotype, expected$Environment)
    idx <- match(paste(actual$Level, actual$Family), key)
    stopifnot(max(abs(actual$Adjusted_mean - expected$emmean[idx])) < 1e-8,
      max(abs(actual$Adjusted_SE - expected$SE[idx])) < 1e-8)
  }
  workbook <- build_export_tables("MET", list(met_trait_names = "Yield", met_by_trait = list(Yield = r)))
  stopifnot(sum(grepl("03c_mc_", names(workbook))) == 6)
  path <- tempfile(fileext = ".xlsx")
  writexl::write_xlsx(workbook, path)
  stopifnot(file.info(path)$size > 0)
  unlink(path)
}
cat("MET mean comparisons, protection, missing cells, charts and all-engine exports passed.\n")
