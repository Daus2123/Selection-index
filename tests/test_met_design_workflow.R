suppressMessages(invisible(source("app.R")))
set.seed(925)
d <- expand.grid(Variety = paste0("G", 1:6), Environment = paste0("E", 1:4), Rep = 1:3)
d$Yield <- 20 + match(d$Variety, paste0("G", 1:6)) +
  match(d$Environment, paste0("E", 1:4)) + rnorm(nrow(d))
designs <- list(balanced = d,
  unequal = d[-c(1, 10, 20), ],
  missing_cells = subset(d, !(Variety == "G6" & Environment != "E4")),
  disconnected = subset(d, (Variety %in% paste0("G", 1:3)) == (Environment %in% c("E1", "E2"))),
  two_locations = subset(d, Environment %in% c("E1", "E2")))
for (engine in c("LMM", "ANOVA_RCBD", "SOMMER")) {
  for (design in names(designs)) {
    cat("RUN", engine, design, "\n")
    result <- run_met_pipeline(designs[[design]], "Yield", replication_col = "Rep", model_type = engine)
    stopifnot(nrow(result$blups_main) == 6, nrow(result$met_selection) == 6)
    if (all(is.finite(result$blups_main$BLUP_G)) && design != "disconnected") {
      stopifnot(all(is.finite(result$met_selection$Combined_score)))
    } else stopifnot(all(is.na(result$met_selection$Final_rank)))
    for (display in "adjusted_se") {
      ggplot2::ggplot_build(met_plot_mean_comparison(result, display = display))
    }
    if (design %in% c("two_locations", "disconnected")) stopifnot(nrow(result$ammi_genotype) == 0,
      nrow(result$gge_genotype) == 0)
    if (design == "disconnected") stopifnot(all(is.na(result$fw_results$Sens)))
    if (design == "two_locations") stopifnot(all(is.finite(result$fw_results$Sens)))
    integrated <- build_met_integrated_ranking(designs[[design]], list(Yield = result), trait_cols = "Yield")
    if (design == "disconnected" || any(!is.finite(result$blups_main$BLUP_G))) {
      stopifnot(nrow(integrated$ranking) == 0, nrow(integrated$trait_weights) > 0)
    } else stopifnot(nrow(integrated$ranking) == 6)
    if (design == "two_locations") {
      board <- build_met_decision_board(designs[[design]], list(Yield = result), integrated$ranking)
      stopifnot(all(board$board$Final_action == "RETEST"))
    }
    stopifnot(length(build_export_tables("MET", list(met_by_trait = list(Yield = result),
      met_trait_names = "Yield"))) > 0)
  }
}
# Application orchestration must keep trait results when overall BLUEs cannot
# be estimated, and an optional comparison failure must never erase core MET.
all_traits <- run_met_all_traits(designs$missing_cells, trait_cols = "Yield",
  replication_col = "Rep", model_type = "ANOVA_RCBD")
stopifnot(length(all_traits$met_by_trait) == 1, nrow(all_traits$met_integrated_ranking) == 0)
original_comparison <- met_build_mean_comparison
met_build_mean_comparison <- function(...) stop("Injected comparison failure")
protected <- run_met_pipeline(designs$two_locations, "Yield", replication_col = "Rep", model_type = "LMM")
met_build_mean_comparison <- original_comparison
stopifnot(nrow(protected$blups_main) == 6,
  any(grepl("Injected comparison failure", protected$mean_comparison$notes$Note)))
cat("MET design workflow passed.\n")
