suppressMessages(invisible(source("app.R")))
set.seed(812)
stopifnot(identical(si_clean_value_flag(c("Repeated identifier/design combination: verify plot identity",
  "Missing or non-numeric/non-finite value", "Standardized residual flag", "IQR flag; Repeated identifier",
  "Above trait maximum", "")), c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE)))
trial <- expand.grid(Variety = c("A", "B", "C"), Rep = as.character(1:8))
trial$Yield <- 10 + match(trial$Variety, c("A", "B", "C")) + rnorm(nrow(trial), sd = .2)
trial$Height <- 50 + rnorm(nrow(trial))
trial$Yield[1] <- 100
screen <- si_clean_screen(trial, "Yield", id_col, c(id_col, rep_col), "iqr", 1.5)
stopifnot(nrow(screen$review) == nrow(trial),
          grepl("IQR", screen$review$Flag[1]), all(screen$review$Action == "Keep"))
screen$review$Action[1] <- "Exclude"
stopifnot(inherits(try(si_clean_apply(trial, screen$review), silent = TRUE), "try-error"))
screen$review$Reason[1] <- "Confirmed weighing error"
approved <- si_clean_apply(trial, screen$review)
stopifnot(is.na(approved$Yield[1]), identical(approved$Height, trial$Height), trial$Yield[1] == 100)
res <- si_clean_screen(trial, "Yield", id_col, c(id_col, rep_col), "residual", 3)
stopifnot(grepl("residual", res$review$Flag[1]))
metadata <- trial[1, ]; metadata$Variety <- "WEIGHT"
with_metadata <- rbind(metadata, trial)
stopifnot(!1L %in% si_clean_screen(with_metadata, "Yield", id_col, character(), "iqr", 1.5)$review$Row)
# Model diagnostics consume supplied data without a second hidden IQR pass.
stopifnot(si_diagnostic_preview(make_diagnostic_data(trial), "Yield", "trimmed")$n_used == nrow(trial),
          si_diagnostic_preview(make_diagnostic_data(approved), "Yield", "trimmed")$n_used == nrow(trial) - 1L)
# End-to-end reactive approval, raw/trim source, stale-result invalidation, and export.
path <- tempfile(fileext = ".xlsx")
writexl::write_xlsx(trial, path)
shiny::testServer(server, {
  session$setInputs(excel_file = list(name = "synthetic.xlsx", datapath = path, size = file.info(path)$size),
    analysis_method = "LPSI", clean_identity = id_col, clean_trait = "Yield",
    clean_predictors = c(id_col, rep_col), clean_method = "iqr", clean_threshold = 1.5,
    clean_min = "", clean_max = "", diagnostic_plot_type = "Raw trait histogram")
  session$setInputs(clean_check = 1)
  stopifnot(nrow(cleaning$review) == nrow(trial), cleaning$pending)
  stopifnot(identical(clean_visible_rows(), which(si_clean_value_flag(cleaning$review$Flag))))
  session$setInputs(clean_review_rows_selected = 1L, clean_action = "Exclude", clean_reason = "Weighing error")
  analysis_results(list(index_ranking = data.frame(ID = c("A", "B"), Rank = c(1, 2))))
  analysis_used("LPSI")
  session$setInputs(clean_apply = 1)
  stopifnot(cleaning$pending, analysis_data()$Yield[1] == 100)
  session$setInputs(clean_preview = 1)
  stopifnot(!is.null(cleaning$preview), analysis_data()$Yield[1] == 100,
    is.na(si_clean_apply(uploaded_data(), cleaning$preview)$Yield[1]))
  # Render the proposed chart before acceptance.
  stopifnot(!is.null(output$clean_preview_plot))
  session$setInputs(clean_apply = 2)
  stopifnot(!cleaning$pending, is.na(analysis_data()$Yield[1]),
    cleaning_impact()$Excluded == 1L, is.null(analysis_results()))
  stopifnot(diagnostic_data()$data$Yield[1] == 100)
  result <- attach_cleaning(list(settings = data.frame(Test = TRUE),
    index_ranking = data.frame(ID = c("A", "B"), Rank = c(2, 1))))
  stopifnot(nrow(result$cleaning_rank_changes) == 2,
    result$cleaning_rank_changes$Rank_before[1] == 1,
    result$cleaning_rank_changes$Rank_after[1] == 2)
  sheets <- build_export_tables("BREEDING", result)
  stopifnot(any(grepl("cleaning_review", names(sheets))), any(grepl("cleaning_impact", names(sheets))))
  session$setInputs(clean_threshold = 2)
  stopifnot(!cleaning$pending, cleaning$needs_check)
  session$setInputs(clean_apply = 3)
  stopifnot(!cleaning$pending, cleaning$approved$Reason[1] == "Weighing error")
  # Browsing other traits, changing analysis controls, or table redraws must not
  # revoke approval or prevent a run on the already accepted dataset.
  session$setInputs(eval_trait = "Height", lpsi_trial_model = "RCBD",
    met_model_type = "LMM", clean_review_rows_selected = integer(), clean_reason = "")
  stopifnot(!cleaning$pending, is.na(analysis_data()$Yield[1]))
  session$setInputs(run_analysis = 1)
  stopifnot(!is.null(analysis_results()), is.na(analysis_results()$cleaned_data$Yield[1]),
    sum(analysis_results()$cleaning_review$Action == "Exclude") == 1)
  # An explicitly started new review still needs preview and acceptance.
  session$setInputs(clean_check = 2)
  stopifnot(cleaning$pending)
  # Further trait checks cannot block a run on the last accepted snapshot.
  session$setInputs(run_analysis = 2)
  stopifnot(!is.null(analysis_results()), cleaning$pending, cleaning$has_approved,
    is.na(analysis_results()$cleaned_data$Yield[1]),
    identical(analysis_results()$cleaning_review, cleaning$approved))
})
unlink(path)
cat("Cleaning helper, reactive approval, preview, and export checks passed.\n")
