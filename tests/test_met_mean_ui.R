suppressMessages(invisible(source("app.R")))
set.seed(345)
d <- expand.grid(Genotype = LETTERS[1:4], Environment = c("E1", "E2"), Rep = factor(1:3))
d$Weight <- 10 + as.integer(d$Genotype) + rnorm(nrow(d))
r <- list(mean_comparison = met_build_mean_comparison(d, TRUE, FALSE),
          model_summary = data.frame(Model = "LMM"))
stopifnot(identical(unname(met_mean_comparison_views()),
  c("genotypes", "within_environment", "environments", "interaction")))
gated <- r
gated$mean_comparison$means$Omnibus_p <- 0.5
empty <- met_plot_mean_comparison(gated)
stopifnot(grepl("No significant difference", ggplot2::ggplot_build(empty)$data[[1]]$label))
gated$mean_comparison$means$Omnibus_p <- 0.01
significant <- met_plot_mean_comparison(gated)
stopifnot(all(significant$data$Plot_mean == significant$data$Raw_mean),
  any(vapply(significant$layers, function(x) inherits(x$geom, "GeomText"), logical(1))))
for (method in c("tukey", "lsd")) {
  plot <- met_plot_mean_comparison(gated, method = method)
  expected <- met_comparison_table(gated, "means", method = method)
  stopifnot(identical(as.character(ggplot2::ggplot_build(plot)$data[[3]]$label), expected$Group),
    identical(plot$data$Plot_error, expected$Observed_SD),
    plot$theme$plot.title$size == 10, plot$theme$plot.subtitle$size == 9)
}
shiny::testServer(server, {
  analysis_results(list(met_by_trait = list(Yield = r, Height = r), met_trait_names = c("Yield", "Height")))
  analysis_used("MET")
  session$setInputs(result_module = "met", result_view = "met_means", met_result_trait = "Yield",
    chart_module = "met", plot_view = "met_mean_plot", met_plot_trait = "Yield",
    chart_width = 12, chart_height = 7)
  controls <- output$met_mean_controls$html
  stopifnot(grepl("Variety", controls, fixed = TRUE),
    grepl('type="radio"', controls, fixed = TRUE),
    !grepl("<select", controls, fixed = TRUE))
  for (method in c("tukey", "lsd")) {
    for (view in unname(met_mean_comparison_views())) {
      session$setInputs(met_mean_method = method, met_mean_view = view,
        met_mean_chart_method = method, met_mean_chart_view = view)
      stopifnot(length(output$met_mean_plot$src) == 1,
        grepl("met_mean_plot", output$chart_preview_ui$html, fixed = TRUE),
        nrow(selected_result_tables()$tables$Means) > 0,
        inherits(selected_chart()$plot, "ggplot"))
      stopifnot(identical(names(selected_result_tables()$tables$Means),
        c("Trait", "Level", "Raw_mean_SD", "Group", "Adjusted_mean")),
        setequal(names(selected_result_tables()$tables), c("Means", "ANOVA")))
      stopifnot(setequal(selected_result_tables()$tables$Means$Trait, c("Yield", "Height")),
        setequal(selected_result_tables()$tables$ANOVA$Trait, c("Yield", "Height")))
      invisible(output$met_mean_table)
      invisible(output$met_mean_pairwise)
      invisible(output$met_mean_anova)
      invisible(output$met_mean_tests)
      for (display in c("adjusted_se", "observed_sd")) {
        session$setInputs(met_mean_chart_display = display)
        plot <- selected_chart()$plot
        stopifnot(identical(plot$labels$y, "Observed mean \u00b1 SD"),
          length(output$met_mean_plot$src) == 1)
      }
    }
  }
  ggplot2::ggsave("tmp/met_mean_comparison_preview.png", selected_chart()$plot, width = 9, height = 5)
})
# Synthetic previews for all four comparison families, using the same chart
# function as Shiny previews and chart downloads.
for (view in unname(met_mean_comparison_views())) {
  for (display in c("adjusted_se", "observed_sd")) {
    plot <- met_plot_mean_comparison(r, view, "tukey", display)
    ggplot2::ggsave(paste0("tmp/met_mean_", view, "_", display, ".png"),
      plot, width = 11, height = 7)
  }
}
cat("Shiny MET comparison controls, tables, previews and downloads passed.\n")
