suppressMessages(invisible(source("app.R")))

set.seed(7)
trial <- expand.grid(
  Variety = c("A", "B", "C", "D"),
  Rep = as.character(1:3),
  Environment = c("E1", "E2", "E3")
)
trial$Yield <- 10 +
  match(trial$Variety, c("A", "B", "C", "D")) * 1.2 +
  match(trial$Environment, c("E1", "E2", "E3")) * 0.8 +
  rnorm(nrow(trial), sd = 0.3)

for (model_type in c("LMM", "ANOVA_RCBD", "SOMMER")) {
  trait_result <- run_met_pipeline(
    trial,
    trait_used = "Yield",
    check_varieties = "A",
    replication_col = "Rep",
    min_envs_for_biplot = 4,
    model_type = model_type
  )
  result <- list(
    settings = data.frame(Model = model_type),
    met_by_trait = list(Yield = trait_result),
    met_trait_names = "Yield",
    met_failed_traits = data.frame(Trait = character(), Error = character()),
    met_integrated_ranking = NULL,
    met_decision_board = NULL
  )
  stopifnot(
    length(result$met_trait_names) == 1,
    length(build_export_tables("MET", result)) > 0
  )
  # Multi-trait aggregation must retain valid estimates and decision rows.
  # The dedicated two-location test covers unavailable GGE results.
  multi_trial <- trial
  multi_trial$Height <- multi_trial$Yield
  multi_results <- list(Yield = trait_result, Height = trait_result)
  integrated <- build_met_integrated_ranking(multi_trial, multi_results)
  methods <- build_met_multitrait_methods(multi_trial, multi_results, integrated$ranking)
  stopifnot(is.list(methods))
  decision <- build_met_decision_board(multi_trial, multi_results, integrated$ranking,
                                      check_genotypes = "A")
  stopifnot(nrow(decision$board) == 3, nrow(decision$passport) == 3,
            setequal(decision$board$Genotype, setdiff(trial$Variety, "A")))
  chart <- trait_result$p_fw_regression
  stopifnot(identical(chart$labels$x, "Environment / environmental index (mean of other genotypes)"),
            setequal(unique(as.character(chart$data$Environment)), unique(trial$Environment)))
  axis <- chart$scales$get_scales("x")
  stopifnot(length(axis$breaks) == 3,
            setequal(axis$labels, unique(as.character(trial$Environment))),
            is.null(chart$labels$caption))
  built <- ggplot2::ggplot_build(chart)
  stopifnot(nrow(built$data[[2]]) == 12,
            all(is.finite(built$data[[2]]$x)),
            all(is.finite(built$data[[2]]$y)))
  for (point_group in unique(built$data[[2]]$group)) {
    points <- built$data[[2]][built$data[[2]]$group == point_group, ]
    line <- built$data[[1]][built$data[[1]]$group == point_group, ]
    genotype_levels <- sort(unique(as.character(chart$data$Genotype)))
    genotype_colors <- built$plot$scales$get_scales("colour")$map(genotype_levels)
    genotype <- genotype_levels[match(points$colour[1], genotype_colors)]
    expected_slope <- trait_result$fw_results$Sens[trait_result$fw_results$Genotype == genotype]
    plotted_slope <- (tail(line$y, 1) - line$y[1]) / (tail(line$x, 1) - line$x[1])
    stopifnot(isTRUE(all.equal(unname(plotted_slope), unname(expected_slope), tolerance = 1e-7)))
  }
}

cat("MET LMM, ANOVA, and sommer export checks passed.\n")
