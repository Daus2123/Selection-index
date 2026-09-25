source(file.path("modules", "module_multi_factor.R"))

set.seed(23)
two <- expand.grid(A = c("check", "new", "new2"), B = c("low", "high"), Rep = 1:5)
two$Y <- 10 + 2 * (two$A == "new") + 3 * (two$A == "new2") +
  1.5 * (two$B == "high") + stats::rnorm(nrow(two), sd = 0.7)
three <- expand.grid(A = c("check", "new"), B = c("low", "high"),
                     C = c("early", "late"), Rep = 1:5)
three$Y <- 10 + 2 * (three$A == "new") + (three$B == "high") +
  (three$C == "late") + stats::rnorm(nrow(three), sd = 0.7)

cases <- list(
  list(two, "factorial", "CRD", "ANOVA"),
  list(two, "factorial", "RCBD", "ANOVA"),
  list(two, "factorial", "RCBD", "LMM"),
  list(two, "split_plot", "CRD", "LMM"),
  list(two, "split_plot", "RCBD", "LMM"),
  list(three, "factorial", "CRD", "ANOVA"),
  list(three, "factorial", "RCBD", "ANOVA"),
  list(three, "factorial", "RCBD", "LMM"),
  list(three, "split_plot", "CRD", "LMM"),
  list(three, "split_plot", "RCBD", "LMM")
)
for (case in cases) {
  dat <- case[[1]]
  fit <- run_multifactor_pipeline(
    dat, "Y", if ("C" %in% names(dat)) c("A", "B", "C") else c("A", "B"),
    design = case[[2]], randomization = case[[3]], model_type = case[[4]],
    replication_col = if (case[[2]] == "split_plot" || case[[3]] == "RCBD") "Rep" else NULL,
    comparison_factor = "A", checks = "check"
  )
  stopifnot(nrow(fit$anova) > 0, nrow(fit$mean_comparison) > 0,
            nrow(fit$pairwise) > 0, nrow(fit$superiority) > 0,
            all(is.finite(fit$superiority$Check_baseline)),
            all(nzchar(fit$mean_comparison$Group)),
            all(c("Residual CV (%)", "Heritability (h²)", "Tukey HSD (0.05)") %in%
                  fit$effect_comparison$Level))
  plots <- list(mf_plot_mean_comparison(fit), mf_plot_superiority(fit),
                mf_plot_interaction(fit))
  stopifnot(inherits(plots[[1]]$layers[[1]]$geom, "GeomCol"),
            inherits(plots[[2]]$layers[[1]]$geom, "GeomTile"),
            inherits(plots[[3]]$layers[[1]]$geom, "GeomCol"))
  for (plot in plots) invisible(ggplot2::ggplot_build(plot))
  print(c(case[[2]], case[[3]], case[[4]], nrow(fit$anova), nrow(fit$superiority)))
}
lsd_fit <- run_multifactor_pipeline(two, "Y", c("A", "B"), "factorial", "RCBD",
                                    "ANOVA", "Rep", "A", "check",
                                    comparison_method = "lsd")
stopifnot(all(is.finite(lsd_fit$pairwise$Critical_difference)))
interaction_trial <- two
interaction_trial$Y <- 10 + 8 * (interaction_trial$A == "new" & interaction_trial$B == "high") +
  stats::rnorm(nrow(interaction_trial), sd = 0.15)
interaction_fit <- run_multifactor_pipeline(interaction_trial, "Y", c("A", "B"),
                                             "factorial", "RCBD", "ANOVA", "Rep", "A", "check")
stopifnot(all(c("A", "B", "A:B") %in% interaction_fit$effect_comparison$Effect),
          interaction_fit$anova$p_value[interaction_fit$anova$Source == "A:B"] < 0.05,
          all(nzchar(interaction_fit$effect_comparison$Group[
            interaction_fit$effect_comparison$Effect == "A:B" &
              interaction_fit$effect_comparison$Row_type == "Level"])),
          all(interaction_fit$effect_comparison$Raw_mean_SD[
            interaction_fit$effect_comparison$Row_type == "Level"] != ""))
missing_cell <- two[!(two$A == "new" & two$B == "high"), ]
stopifnot(inherits(try(run_multifactor_pipeline(missing_cell, "Y", c("A", "B"),
                                                 "factorial", "CRD", "ANOVA",
                                                 comparison_factor = "A", checks = "check"),
                       silent = TRUE), "try-error"))
app_env <- new.env(parent = globalenv())
suppressMessages(invisible(source("app.R", local = app_env)))
export <- app_env$build_export_tables("MULTIFACTOR", fit)
stopifnot(length(export) == 7L,
          all(c("00_settings", "01_summary", "02_anova", "03_mean_comparison",
                "04_pairwise", "05_superiority", "06_notes") %in% names(export)))
upload_path <- tempfile(fileext = ".xlsx")
writexl::write_xlsx(two, upload_path)
shiny::testServer(app_env$server, {
  session$setInputs(excel_file = list(name = "synthetic.xlsx", datapath = upload_path,
                                      size = file.info(upload_path)$size,
                                      type = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"),
                    analysis_method = "MULTIFACTOR")
  stopifnot(any(grepl("Number of factors", unlist(output$multifactor_controls), fixed = TRUE)))
  session$setInputs(mf_factor_count = "2", mf_design = "factorial", mf_factor_a = "A",
                    mf_factor_b = "B", mf_response = "Y", mf_randomization = "RCBD",
                    mf_model = "ANOVA", mf_rep = "Rep", mf_checks = "check",
                    run_analysis = 1)
  stopifnot(!is.null(saved_results$MULTIFACTOR),
            nrow(saved_results$MULTIFACTOR$anova) > 0,
            any(grepl("Multi-factor", unlist(output$result_multifactor_detail), fixed = TRUE)))
  session$setInputs(result_module = "multifactor", result_view = "mf_means",
                    chart_module = "multifactor", plot_view = "mf_superiority_plot")
  stopifnot(length(selected_result_tables()$tables) == 2L,
            inherits(selected_chart()$plot, "ggplot"),
            any(grepl("mf_superiority_plot", unlist(output$chart_preview_ui), fixed = TRUE)))
  session$setInputs(mf_result_effect = "B", mf_result_method = "lsd",
                    mf_result_direction = "Lower better")
  selected <- selected_result_tables()$tables$Mean_comparison
  stopifnot(all(selected$Effect == "B"), "LSD (0.05)" %in% selected$Level,
            all(multifactor_result()$superiority$Superiority_pct < 0))
})
user_three <- three
names(user_three)[match(c("A", "B", "C", "Y"), names(user_three))] <-
  c("Variety", "Season", "Location", "Percent_severity")
user_three$Percent_severity <- gsub("\\.", ",", sprintf("%.2f", user_three$Percent_severity))
three_path <- tempfile(fileext = ".xlsx")
writexl::write_xlsx(user_three, three_path)
shiny::testServer(app_env$server, {
  session$setInputs(excel_file = list(name = "three_factors.xlsx", datapath = three_path,
                                      size = file.info(three_path)$size,
                                      type = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"),
                    analysis_method = "MULTIFACTOR")
  stopifnot(any(grepl("Percent_severity", unlist(output$multifactor_controls), fixed = TRUE)))
  session$setInputs(mf_factor_count = "3", mf_design = "factorial",
                    mf_factor_a = "Variety", mf_factor_b = "Season",
                    mf_factor_c = "Location", mf_response = "Percent_severity",
                    mf_randomization = "RCBD", mf_model = "ANOVA", mf_rep = "Rep",
                    mf_checks = "check", run_analysis = 1)
  stopifnot(!is.null(saved_results$MULTIFACTOR),
            all(c("A", "B", "C", "A:B", "A:C", "B:C", "A:B:C") %in%
                  saved_results$MULTIFACTOR$effect_comparison$Effect))
})
cat("Multi-factor models, charts, and export checks passed.\n")
