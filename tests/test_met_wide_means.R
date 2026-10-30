suppressMessages(invisible(source("app.R")))
set.seed(814)
trial <- expand.grid(Variety = LETTERS[1:4], Environment = paste0("E", 1:3), Rep = factor(1:3))
trial$Yield <- 10 + as.integer(factor(trial$Variety)) * 2 + rnorm(nrow(trial), sd = 0.4)
for (engine in c("LMM", "ANOVA_RCBD")) {
  result <- run_met_pipeline(trial, "Yield", replication_col = "Rep", model_type = engine,
    min_envs_for_biplot = 4)
  other <- result
  other$mean_comparison$means <- other$mean_comparison$means[nrow(other$mean_comparison$means):1, ]
  results <- list(met_trait_names = c("Yield", "Other"), met_by_trait = list(Yield = result, Other = other))
  for (method in "tukey") {
    wide <- met_combined_comparison_tables(results, "genotypes", method)$Means
    stopifnot(identical(names(wide), c("Levels", "Yield | Adjusted_mean_SE", "Yield | Group", "Other | Adjusted_mean_SE", "Other | Group")),
      identical(wide[[2]], wide[[4]]), identical(wide[[3]], wide[[5]]),
      identical(tail(wide$Levels, 3), c("CV (%)", "Heritability H\u00b2", "Tukey HSD (5%)")),
      is.finite(as.numeric(wide[[2]][wide$Levels == "CV (%)"])),
      identical(wide[[2]][wide$Levels == "Heritability H\u00b2"],
        if (is.finite(result$model_summary$Broad_sense_H2[1]))
          formatC(result$model_summary$Broad_sense_H2[1], digits = 3, format = "f") else "Not available"),
      (tail(wide[[2]], 1) == "Pair-specific" || is.finite(as.numeric(tail(wide[[2]], 1)))))
    within <- met_combined_comparison_tables(results, "interaction", method)$Means
    stopifnot(setequal(within$Environment, paste0("E", 1:3)),
      all(is.finite(as.numeric(within[[3]][within$Levels == "CV (%)"]))))
    altered <- results
    altered$met_by_trait$Yield$mean_comparison$pairwise$Critical_difference <-
      seq_len(nrow(result$mean_comparison$pairwise))
    pair_specific <- met_combined_comparison_tables(altered, "genotypes", method)$Means
    stopifnot(tail(pair_specific[[2]], 1) == "Pair-specific")
    altered$met_by_trait$Yield$mean_comparison$pairwise$Critical_difference <- 1.25
    common <- met_combined_comparison_tables(altered, "genotypes", method)$Means
    stopifnot(tail(common[[2]], 1) == "1.25")
  }
  workbook <- build_export_tables("MET", results)
  exported <- workbook[[grep("combined_mean_comparison", names(workbook))[1]]]
  stopifnot(identical(exported, met_combined_comparison_tables(results)$Means))
  path <- tempfile(fileext = ".xlsx")
  writexl::write_xlsx(workbook, path)
  stopifnot(file.info(path)$size > 0)
  unlink(path)
}
cat("Wide MET alignment, footer statistics, thresholds and LMM/ANOVA exports passed.\n")
