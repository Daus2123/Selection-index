suppressMessages(invisible(source("app.R")))

set.seed(42)
trial <- expand.grid(ID = c("1", "2", "3"), Rep = as.character(1:4))
trial$Original_ID <- paste0("Variety ", trial$ID)
trial$Yield <- c(10, 13, 16)[match(trial$ID, c("1", "2", "3"))] +
  c(-0.3, 0.2, 0.1, -0.1)[match(trial$Rep, as.character(1:4))] +
  rnorm(nrow(trial), sd = 0.2)
trial$Score <- rep(c(2, 3, 4), 4)

make_result <- function(data, model = "RCBD") {
  list(
    cleaned_data = data,
    trait_info = data.frame(
      Trait = c("Yield", "Score"),
      Direction = c("Higher better", "Lower better")
    ),
    anova_full = data.frame(
      Trait = c("Yield", "Score"),
      Source = "ID",
      Test = c("ANOVA", "Kruskal-Wallis"),
      p_value = c(0.001, 0.01)
    ),
    heritability_gain = data.frame(
      Trait = c("Yield", "Score"),
      CV_pct = c(7.9, 10.2),
      Broad_sense_H2 = if (model == "LMM") c(NA_real_, NA_real_) else c(0.68, 0.55)
    ),
    decision_settings = data.frame(Model = model)
  )
}

rcbd <- si_lpsi_mean_comparison_view(make_result(trial), "lsd")
stopifnot(
  nrow(rcbd$means) == 3,
  nrow(rcbd$summary) == 3,
  nrow(rcbd$combined) == 6,
  identical(names(rcbd$means), names(rcbd$summary)),
  identical(names(rcbd$means), c("ID", "Original_ID", "Yield", "Yield Group", "Score", "Score Group")),
  grepl(" \u00b1 ", rcbd$means$Yield[1]),
  grepl("\u00b1", rcbd$means$Score[1]),
  rcbd$summary$Yield[1] == "7.90%",
  rcbd$summary$Yield[2] == "0.680",
  is.finite(as.numeric(rcbd$summary$Yield[3])),
  grepl("^p=0.01 \\(KW; Dunn: [1-9]", rcbd$summary$Score[3]),
  any(nzchar(rcbd$means$`Score Group`)),
  rcbd$means$Yield[1] == sprintf("%.2f \u00b1 %.2f", mean(trial$Yield[trial$ID == "1"]), stats::sd(trial$Yield[trial$ID == "1"]) / 2),
  all(rcbd$summary$`Yield Group` == ""),
  all(rcbd$summary$`Score Group` == "")
)

tukey <- si_lpsi_mean_comparison_view(make_result(trial), "tukey")
stopifnot(
  tukey$summary$Original_ID[3] == "Tukey HSD (0.05)",
  is.finite(as.numeric(tukey$summary$Yield[3]))
)

crd <- si_lpsi_mean_comparison_view(make_result(trial, "CRD"), "lsd")
stopifnot(grepl(" \u00b1 ", crd$means$Yield[1]))

sd_view <- si_lpsi_mean_comparison_view(make_result(trial), "lsd", "sd")
stopifnot(grepl(sprintf("%.2f", stats::sd(trial$Yield[trial$ID == "1"])), sd_view$means$Yield[1], fixed = TRUE))

unbalanced <- trial[!(trial$ID == "3" & trial$Rep == "4"), , drop = FALSE]
unequal <- si_lpsi_mean_comparison_view(make_result(unbalanced), "lsd")
stopifnot(unequal$summary$Yield[3] == "Varies by pair")

lmm <- si_lpsi_mean_comparison_view(make_result(trial, "LMM"), "lsd")
stopifnot(lmm$summary$Yield[2] == "Not estimated")

nonsignificant <- make_result(trial)
nonsignificant$anova_full$p_value[2] <- 0.50
blank_letters <- si_lpsi_mean_comparison_view(nonsignificant, "lsd")
stopifnot(
  blank_letters$summary$Score[3] == "",
  !nzchar(blank_letters$means$`Score Group`[1])
)

nonsignificant$anova_full$p_value[1] <- 0.50
blank_lsd <- si_lpsi_mean_comparison_view(nonsignificant, "lsd")
blank_tukey <- si_lpsi_mean_comparison_view(nonsignificant, "tukey")
stopifnot(
  blank_lsd$summary$Yield[3] == "",
  blank_tukey$summary$Yield[3] == "",
  blank_tukey$summary$Score[3] == "",
  all(blank_lsd$means$`Yield Group` == "")
)

flat <- trial
flat$Score <- rep(1:4, each = 3)
no_pairwise_difference <- si_lpsi_mean_comparison_view(make_result(flat), "lsd")
stopifnot(
  grepl("^p=0.01 \\(KW; Dunn: no significant pairs\\)$", no_pairwise_difference$summary$Score[3]),
  !any(nzchar(no_pairwise_difference$means$`Score Group`))
)

footer <- htmltools::tags$tfoot(lapply(seq_len(nrow(rcbd$summary)), function(i) {
  htmltools::tags$tr(lapply(rcbd$summary[i, , drop = TRUE], htmltools::tags$th))
}))
container <- htmltools::tags$table(
  htmltools::tags$thead(htmltools::tags$tr(lapply(names(rcbd$means), htmltools::tags$th))),
  footer
)
widget <- DT::datatable(rcbd$means, container = container, rownames = FALSE)
stopifnot(inherits(widget, "datatables"))

export_path <- tempfile(fileext = ".xlsx")
writexl::write_xlsx(list(Mean_comparison = rcbd$combined), export_path)
exported <- readxl::read_xlsx(export_path, col_types = "text")
stopifnot(
  nrow(exported) == 6,
  identical(exported$Original_ID[4:6], rcbd$summary$Original_ID),
  grepl("^p=0.01 \\(KW; Dunn:", exported$Score[6]),
  grepl("\u00b1", exported$Yield[1]),
  identical(exported$`Score Group`[1:3], rcbd$means$`Score Group`)
)

plot_result <- make_result(trial)
plot_result$lsd_long <- data.frame(
  Trait = "Yield", ID = c("1", "2", "3"),
  Original_ID = paste("Variety", 1:3), emmean = c(10, 13, 16),
  SE = c(0.1, 0.2, 0.3),
  LSD_group = c("c", "b", "a")
)
plot_result$mean_comparison_method <- "lsd"
mean_plot <- plot_lpsi_mean_comparison(plot_result, "Yield")
stopifnot(inherits(mean_plot$layers[[1]]$geom, "GeomBoxplot"))
invisible(ggplot2::ggplot_build(mean_plot))
annotation_data <- mean_plot$layers[[2]]$data
stopifnot(all(c("Top", "Group") %in% names(annotation_data)),
          !any(c("SE.x", "SE.y") %in% names(annotation_data)),
          all(vapply(seq_len(nrow(annotation_data)), function(i) {
            annotation_data$Top[i] == max(mean_plot$data$Value[
              mean_plot$data$Label == annotation_data$Label[i]])
          }, logical(1))))

# Exercise the real Charts reactive data (including adjusted SE), preview,
# and download plot with the punctuation used in uploaded trait names.
shiny::testServer(server, {
  session$setInputs(chart_module = "selection_index", chart_lpsi_mode = "single",
                    plot_view = "lpsi_mean_comparison_plot",
                    lpsi_chart_trait = "No.ofFruit/plant", lpsi_chart_mean_comparison_method = "lsd")
  for (model in c("CRD", "RCBD", "LMM")) {
    result <- make_result(trial, model)
    names(result$cleaned_data)[names(result$cleaned_data) == "Yield"] <- "No.ofFruit/plant"
    for (field in c("trait_info", "anova_full", "heritability_gain")) {
      result[[field]]$Trait[result[[field]]$Trait == "Yield"] <- "No.ofFruit/plant"
    }
    saved_results$LPSI <- result
    for (method in c("lsd", "tukey")) {
      session$setInputs(lpsi_chart_trait = "No.ofFruit/plant", lpsi_chart_mean_comparison_method = method)
      comparison <- lpsi_selected_mean_comparison()
      stopifnot("SE" %in% names(comparison$plot_data),
                length(output$lpsi_mean_comparison_plot$src) == 1)
      chart <- selected_chart()$plot
      built <- ggplot2::ggplot_build(chart)
      stopifnot(nrow(built$data[[2]]) == 3, all(is.finite(built$data[[2]]$y)))
      png_path <- tempfile(fileext = ".png")
      ggplot2::ggsave(png_path, chart, width = 8, height = 5)
      stopifnot(file.info(png_path)$size > 0)
      unlink(png_path)
    }
    session$setInputs(lpsi_chart_trait = "Score")
    stopifnot(lpsi_selected_mean_comparison()$method == "Dunn_Holm",
              length(output$lpsi_mean_comparison_plot$src) == 1)
    invisible(ggplot2::ggplot_build(selected_chart()$plot))
    result$anova_full$p_value <- 0.5
    saved_results$LPSI <- result
    session$setInputs(lpsi_chart_trait = "No.ofFruit/plant")
    stopifnot(nrow(ggplot2::ggplot_build(selected_chart()$plot)$data[[2]]) == 0,
              length(output$lpsi_mean_comparison_plot$src) == 1)
  }
})

cat("LPSI mean comparison checks passed.\n")
