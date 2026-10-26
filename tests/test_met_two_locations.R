suppressMessages(invisible(source("app.R")))
set.seed(208)
d <- expand.grid(Variety = LETTERS[1:6], Rep = 1:3,
                 Location = c("Rain", "Winter"), stringsAsFactors = FALSE)
d$Yield <- 10 + match(d$Variety, LETTERS) + rnorm(nrow(d))
d$Height <- 30 + match(d$Variety, LETTERS) + rnorm(nrow(d))
# An unselected workbook trait must not prevent the chosen traits from running.
d$Unselected <- NA_real_
pdf(file = NULL)
for (model in c("LMM", "ANOVA_RCBD", "SOMMER")) {
  result <- run_met_all_traits(d, trait_cols = c("Yield", "Height"), replication_col = "Rep", model_type = model)
  stopifnot(length(result$met_trait_names) == 2, nrow(result$met_failed_traits) == 0,
            nrow(result$met_integrated_ranking) == 6,
            length(build_export_tables("MET", result)) > 0)
  stopifnot(setequal(result$met_integrated_trait_weights$Trait, c("Yield", "Height")))
  for (trait_result in result$met_by_trait) {
    stopifnot(all(trait_result$ammi_notes$Status == "NOT AVAILABLE"),
              nrow(trait_result$ammi_genotype) == 0,
              nrow(trait_result$gge_genotype) == 0)
    plot_names <- grep("^p_", names(trait_result), value = TRUE)
    plot_names <- plot_names[!grepl("ammi|gge", plot_names)]
    for (plot_name in plot_names) print(trait_result[[plot_name]])
  }
  methods <- result$met_multitrait_methods
  for (plot in list(result$p_met_integrated_ranking, methods$plot_gt,
                    methods$plot_gyt, methods$plot_method_agreement)) print(plot)
  reranked <- build_met_integrated_ranking(d, result$met_by_trait,
                                          trait_cols = result$settings$Trait)
  stopifnot(nrow(reranked$ranking) == 6)
  missing_model <- tryCatch({
    build_met_integrated_ranking(d, result$met_by_trait["Yield"], trait_cols = c("Yield", "Height"))
    "NO ERROR"
  }, error = conditionMessage)
  stopifnot(grepl("required MET trait models are missing: Height", missing_model, fixed = TRUE))
  board <- build_met_decision_board(d, result$met_by_trait, result$met_integrated_ranking)
  stopifnot(nrow(board$board) == 6, nrow(board$passport) == 6)
}
dev.off()
stopifnot(identical(si_plain_error_message(simpleError("\033[31mColumn missing\033[39m")), "Column missing"))
cat("Two-location MET models, decisions, exports, and plain errors passed.\n")

prepared <- prepare_met_trait_settings(d)
stopifnot(identical(met_default_replication_col(prepared, "Genotype"), "Rep"),
          !any(c("Genotype", "Environment") %in% met_design_column_choices(prepared)))
no_rep <- d[, setdiff(names(d), "Rep")]
stopifnot(identical(met_default_replication_col(prepare_met_trait_settings(no_rep)), ""))
for (bad_rep in c("Genotype", id_col, "Environment", "Location", "Yield")) {
  message <- tryCatch({make_met_data(d, "Yield", replication_col = bad_rep); "NO ERROR"},
                      error = conditionMessage)
  stopifnot(grepl("Invalid MET replication column", message, fixed = TRUE))
}
cat("MET replication selector and invalid-column guards passed.\n")
