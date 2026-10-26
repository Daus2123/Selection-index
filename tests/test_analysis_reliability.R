suppressMessages(invisible(source("app.R")))
set.seed(819)
d <- expand.grid(Variety = LETTERS[1:8], Rep = as.character(1:6), stringsAsFactors = FALSE)
d$Yield <- 20 + match(d$Variety, LETTERS[1:8]) + 2 * as.numeric(d$Rep) + rexp(nrow(d), rate = 1)
internal <- d; names(internal)[names(internal) == id_col] <- "ID"
internal$ID <- factor(internal$ID); internal$Rep <- factor(internal$Rep)
for (model in c("CRD", "RCBD", "LMM")) {
  context <- list(kind = "single", identity = id_col, rep = rep_col, model = model)
  fit <- si_trial_diagnostic_fit(d, "Yield", context)
  expected <- fit_single_location_model(internal, "Yield", model)
  stopifnot(max(abs(fit$observations$Residual - residuals(expected))) < 1e-5)
  review <- si_clean_screen(d, "Yield", id_col, c(id_col, rep_col), "residual", 2, context = context)
  flagged <- grepl("Standardized residual flag", review$review$Flag)
  stopifnot(identical(flagged, is.finite(fit$observations$Standardized) & abs(fit$observations$Standardized) > 2))
  diagnostic <- make_diagnostic_data(d, context)
  plot <- make_diagnostic_plot(diagnostic, "Yield", "Residual histogram")
  stopifnot(max(abs(plot$data$Residual - fit$observations$Residual)) < 1e-5)
  decision <- get_lpsi_normality_decision(internal, "Yield", model)
  stopifnot(decision$use_anova, grepl("selected design model retained", decision$note))
}
# Environment names and nested replication are consistently used in MET.
met <- d[rep(seq_len(nrow(d)), 3), ]
met$Location <- rep(c("E1", "E2", "E3"), each = nrow(d))
met$Yield <- met$Yield + rep(c(-10, 0, 12), each = nrow(d)) + rnorm(nrow(met))
for (model in c("ANOVA_RCBD", "LMM", "SOMMER")) {
  context <- list(kind = "met", identity = id_col, environment = "Location", rep = rep_col, block = "", model = model)
  fit <- si_trial_diagnostic_fit(met, "Yield", context)
  diagnostic <- make_diagnostic_data(met, context)
  plot <- make_diagnostic_plot(diagnostic, "Yield", "QQ plot")
  stopifnot(max(abs(plot$data$Residual - fit$observations$Residual)) < 1e-4,
    grepl("Environment", fit$label, fixed = TRUE), grepl("Replication", fit$label, fixed = TRUE))
  if (model == "SOMMER") stopifnot(grepl("approximation", fit$label))
}
# Split-plot diagnostics use the same error strata as the Multi-factor engine.
mf <- expand.grid(Variety = c("A1", "A2"), Treatment = c("B1", "B2"), Rep = as.character(1:5), stringsAsFactors = FALSE)
mf$Yield <- 10 + 2 * (mf$Variety == "A2") + (mf$Treatment == "B2") + rnorm(nrow(mf))
context <- list(kind = "multifactor", identity = id_col, factors = c(id_col, "Treatment"),
  rep = rep_col, design = "split_plot", randomization = "RCBD", model = "LMM")
fit <- si_trial_diagnostic_fit(mf, "Yield", context)
stopifnot(grepl("WholePlot", fit$label), nrow(fit$observations) == nrow(mf))
# A skewed trait must retain its RCBD model, not switch to unblocked KW.
skew <- d
skew$Yield <- 30 + 2 * as.numeric(skew$Rep) + match(skew$Variety, LETTERS[1:8]) + c(rep(0, nrow(d) - 1), 40)
result <- run_selection_pipeline(skew, check_varieties = "A", model_type = "RCBD")
stopifnot(any(result$anova_full$Normality_p < .05, na.rm = TRUE),
  all(result$anova_full$Test == "ANOVA"), !any(result$anova_full$Test == "Kruskal-Wallis"),
  all(grepl("Rep", result$residual_diagnostics$Model)),
  inherits(si_lpsi_diagnostic_plot(result, "Yield"), "patchwork"),
  any(grepl("model_diagnostics", names(build_export_tables("LPSI", result)))))
# A failed adjusted model must never become a raw-mean ranking.
bad <- d; bad$Rep <- "1"
failure <- tryCatch(run_selection_pipeline(bad, check_varieties = "A", model_type = "LMM"), error = identity)
stopifnot(inherits(failure, "error"), grepl("No raw-mean fallback", conditionMessage(failure)))
blue_failure <- tryCatch(si_ext_estimate_trait_values(d, "Yield", id_col, NULL, "blue"), error = identity)
stopifnot(inherits(blue_failure, "error"))
cat("Design-consistent screening/plots, advisory Shapiro, explicit estimation failures and final diagnostics checks passed.\n")
