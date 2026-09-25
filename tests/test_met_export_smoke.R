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
}

cat("MET LMM, ANOVA, and sommer export checks passed.\n")
