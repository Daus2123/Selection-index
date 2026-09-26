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
            is.finite(fit$heritability$value), fit$heritability$value >= 0,
            fit$heritability$value <= 1,
            !any(c("Omnibus_p", "Row_type") %in% names(fit$effect_comparison)),
            nrow(fit$pairwise) > 0, nrow(fit$superiority) > 0,
            all(is.finite(fit$superiority$Check_baseline)),
            all(nzchar(fit$mean_comparison$Group)),
            all(c("Residual CV (%)", "Heritability (h\u00b2)", "Tukey HSD (0.05)") %in%
                  fit$effect_comparison$Level))
  plots <- list(mf_plot_mean_comparison(fit), mf_plot_superiority(fit),
                mf_plot_interaction(fit))
  stopifnot(inherits(plots[[1]]$layers[[1]]$geom, "GeomCol"),
            inherits(plots[[2]]$layers[[1]]$geom, "GeomTile"),
            inherits(plots[[3]]$layers[[1]]$geom, "GeomCol"))
  for (plot in plots) invisible(ggplot2::ggplot_build(plot))
  residual <- fit$anova[fit$anova$Source == "Residuals", ]
  stopifnot(nrow(residual) == 1, is.finite(residual$Mean_Sq),
            "Significance" %in% names(fit$anova))
  if (case[[4]] == "ANOVA") {
    stopifnot(all(is.finite(fit$anova$Mean_Sq)),
      isTRUE(all.equal(fit$anova$Mean_Sq, fit$anova$Sum_Sq / fit$anova$Num_df)),
      residual$Num_df == stats::df.residual(fit$model))
  } else stopifnot(is.na(residual$Sum_Sq), is.na(residual$Num_df),
    isTRUE(all.equal(residual$Mean_Sq, stats::sigma(fit$model)^2)))
  stopifnot(length(unique(ggplot2::ggplot_build(plots[[1]])$data[[1]]$fill)) > 1,
            plots[[1]]$theme$axis.text$size >= 12)
  print(c(case[[2]], case[[3]], case[[4]], nrow(fit$anova), nrow(fit$superiority)))
}
lsd_fit <- run_multifactor_pipeline(two, "Y", c("A", "B"), "factorial", "RCBD",
                                    "ANOVA", "Rep", "A", "check",
                                    comparison_method = "lsd")
stopifnot(all(is.finite(lsd_fit$pairwise$Critical_difference)))
for (method in c("lsd", "tukey")) {
  for (direction in c("Higher better", "Lower better")) {
    refreshed <- mf_refresh_result(lsd_fit, method, direction)
    for (effect in unique(refreshed$effect_comparison$Effect)) {
      rows <- subset(refreshed$effect_comparison, Effect == effect & !is.na(N))
      if (any(nzchar(rows$Group))) {
        best <- if (direction == "Higher better") which.max(rows$Adjusted_mean) else which.min(rows$Adjusted_mean)
        stopifnot(grepl("a", rows$Group[best], fixed = TRUE))
      }
    }
    for (rows in split(refreshed$mean_comparison, refreshed$mean_comparison$B)) {
      best <- if (direction == "Higher better") which.max(rows$Adjusted_mean) else which.min(rows$Adjusted_mean)
      stopifnot(grepl("a", rows$Group[best], fixed = TRUE))
    }
  }
}
# In a balanced two-factor trial, Cullis H2 agrees with the variance-ratio
# entry-mean formula, including genotype-by-factor variance.
h2_model <- lme4::lmer(Y ~ B + (1 | A) + (1 | A:B) + (1 | Rep), data = two)
h2_v <- lme4::VarCorr(h2_model)
h2_expected <- as.numeric(h2_v$A) / (as.numeric(h2_v$A) +
  as.numeric(h2_v[["A:B"]]) / nlevels(two$B) + stats::sigma(h2_model)^2 / 10)
stopifnot(isTRUE(all.equal(lsd_fit$heritability$value, h2_expected, tolerance = 1e-6)))
interaction_trial <- two
interaction_trial$Y <- 10 + 8 * (interaction_trial$A == "new" & interaction_trial$B == "high") +
  stats::rnorm(nrow(interaction_trial), sd = 0.15)
interaction_fit <- run_multifactor_pipeline(interaction_trial, "Y", c("A", "B"),
                                             "factorial", "RCBD", "ANOVA", "Rep", "A", "check")
stopifnot(all(c("A", "B", "A:B") %in% interaction_fit$effect_comparison$Effect),
          interaction_fit$anova$p_value[interaction_fit$anova$Source == "A:B"] < 0.05,
          all(nzchar(interaction_fit$effect_comparison$Group[
            interaction_fit$effect_comparison$Effect == "A:B" &
              !is.na(interaction_fit$effect_comparison$N)])),
          all(interaction_fit$effect_comparison$Raw_mean_SD[
            !is.na(interaction_fit$effect_comparison$N)] != ""))
missing_cell <- two[!(two$A == "new" & two$B == "high"), ]
stopifnot(inherits(try(run_multifactor_pipeline(missing_cell, "Y", c("A", "B"),
                                                 "factorial", "CRD", "ANOVA",
                                                 comparison_factor = "A", checks = "check"),
                       silent = TRUE), "try-error"))
app_env <- new.env(parent = globalenv())
suppressMessages(invisible(source("app.R", local = app_env)))
export <- app_env$build_export_tables("MULTIFACTOR", fit)
stopifnot(length(export) == 7L,
          !any(c("Omnibus_p", "Row_type") %in% names(export$`03_mean_comparison`)),
          any(grepl("%", export$`03_mean_comparison`$Raw_mean_SD[
            export$`03_mean_comparison`$Effect == "A" &
              grepl("Heritability", export$`03_mean_comparison`$Level)])),
          all(c("00_settings", "01_summary", "02_anova", "03_mean_comparison",
                "04_pairwise", "05_superiority", "06_notes") %in% names(export)))
upload_path <- tempfile(fileext = ".xlsx")
upload <- two
upload$Y2 <- 100 - 2 * upload$Y
writexl::write_xlsx(upload, upload_path)
shiny::testServer(app_env$server, {
  session$setInputs(excel_file = list(name = "synthetic.xlsx", datapath = upload_path,
                                      size = file.info(upload_path)$size,
                                      type = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"),
                    analysis_method = "MULTIFACTOR")
  stopifnot(any(grepl("Number of factors", unlist(output$multifactor_controls), fixed = TRUE)))
  session$setInputs(mf_factor_count = "2", mf_design = "factorial", mf_factor_a = "A",
                    mf_factor_b = "B", mf_response = c("Y", "Y2"), mf_randomization = "RCBD",
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
  session$setInputs(plot_view = "mf_mean_plot", mf_chart_order = "2", mf_chart_effect = "A:B")
  stopifnot(all(selected_chart()$plot$data$Effect == "A:B"))
  session$setInputs(mf_chart_order = "1", mf_chart_effect = "A")
  stopifnot(all(selected_chart()$plot$data$Effect == "A"),
            all(selected_result_tables()$tables$Mean_comparison$Effect == "B"))
  session$setInputs(result_view = "mf_anova", mf_result_trait = "Y2")
  stopifnot(multifactor_result()$response == "Y2",
    isTRUE(all.equal(selected_result_tables()$tables$ANOVA$Mean_Sq,
      saved_results$MULTIFACTOR$results_by_trait$Y2$anova$Mean_Sq)))
  session$setInputs(result_view = "mf_superiority", mf_result_direction = "From data")
  stopifnot(all(selected_result_tables()$tables$Superiority$Superiority_pct < 0))
  session$setInputs(mf_result_trait = "Y", mf_chart_trait = "Y2")
  stopifnot(multifactor_result()$response == "Y", multifactor_chart_result()$response == "Y2",
            selected_chart()$plot$labels$y == "Y2",
            mf_resolved_direction("Y") == "Lower better",
            mf_resolved_direction("Y2") == "Higher better")
  all_traits_export <- app_env$build_export_tables("MULTIFACTOR", multifactor_export_result())
  stopifnot(all(vapply(all_traits_export, function(table) {
    setequal(unique(table$Trait), c("Y", "Y2"))
  }, logical(1))))
  controls <- paste(unlist(output$result_multifactor_detail), collapse = "")
  stopifnot(grepl("input.result_view == &#39;mf_means&#39;", controls, fixed = TRUE) ||
              grepl("input.result_view == 'mf_means'", controls, fixed = TRUE))
  stopifnot(!grepl("From DIRECTION row", controls, fixed = TRUE),
            !grepl("Factor interaction", paste(unlist(output$chart_multifactor_detail), collapse = ""), fixed = TRUE))
})
user_three <- three
names(user_three)[match(c("A", "B", "C", "Y"), names(user_three))] <-
  c("Variety", "Season", "Location", "Percent_severity")
user_three$Percent_severity <- gsub("\\.", ",", sprintf("%.2f", user_three$Percent_severity))
metadata <- user_three[1:2, ]
metadata[,] <- NA
metadata$Variety <- "DIRECTION"
metadata$Percent_severity <- c("high", "low")
user_three <- rbind(user_three, metadata)
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
  stopifnot(mf_resolved_direction("Percent_severity") == "Lower better",
            all(multifactor_result()$superiority$Superiority_pct < 0))
  session$setInputs(mf_result_direction = "Higher better", plot_view = "mf_mean_plot",
                    chart_module = "multifactor", mf_chart_order = "3", mf_chart_effect = "A:B:C")
  stopifnot(all(multifactor_result()$superiority$Superiority_pct > 0),
            all(selected_chart()$plot$data$Effect == "A:B:C"))
  for (effect in c("A:B", "A:C", "B:C")) {
    session$setInputs(mf_chart_order = "2", mf_chart_effect = effect)
    stopifnot(all(selected_chart()$plot$data$Effect == effect))
  }
  for (effect in c("A", "B", "C")) {
    session$setInputs(mf_chart_order = "1", mf_chart_effect = effect)
    stopifnot(all(selected_chart()$plot$data$Effect == effect))
  }
})
cat("Multi-factor models, charts, and export checks passed.\n")
