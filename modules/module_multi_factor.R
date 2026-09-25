# Multi-factor factorial and split-plot trial analysis.
# A, B, and C are treatment factors; Rep is a block for RCBD and a
# whole-plot replicate identifier for CRD split-plot designs.

mf_required_column <- function(df, column, label) {
  if (length(column) != 1 || is.na(column) || !nzchar(column) || !column %in% names(df)) {
    stop(label, " must be a column in the uploaded file.", call. = FALSE)
  }
  column
}

mf_prepare_data <- function(df, response_col, factor_cols, replication_col,
                            design, randomization, metadata_col = NULL,
                            metadata_labels = character(0)) {
  if (!is.data.frame(df) || nrow(df) == 0) stop("Upload a nonempty trial table.", call. = FALSE)
  if (anyDuplicated(names(df))) stop("Column names must be unique.", call. = FALSE)
  response_col <- mf_required_column(df, response_col, "Response trait")
  factor_cols <- as.character(factor_cols)
  if (!length(factor_cols) %in% c(2, 3) || anyDuplicated(factor_cols)) {
    stop("Select two or three different factor columns.", call. = FALSE)
  }
  invisible(lapply(seq_along(factor_cols), function(i) {
    mf_required_column(df, factor_cols[i], paste("Factor", LETTERS[i]))
  }))
  if (response_col %in% factor_cols) stop("The response trait cannot be a factor column.", call. = FALSE)
  needs_rep <- randomization == "RCBD" || design == "split_plot"
  if (needs_rep) {
    replication_col <- mf_required_column(df, replication_col, "Block / whole-plot replicate")
    if (replication_col %in% c(response_col, factor_cols)) {
      stop("Block / whole-plot replicate must differ from the response and factors.", call. = FALSE)
    }
  }
  notes <- character(0)
  rows <- rep(TRUE, nrow(df))
  if (!is.null(metadata_col) && metadata_col %in% names(df) && length(metadata_labels) > 0) {
    metadata <- toupper(trimws(as.character(df[[metadata_col]]))) %in% metadata_labels
    metadata[is.na(metadata)] <- FALSE
    if (any(metadata)) {
      rows <- rows & !metadata
      notes <- c(notes, paste(sum(metadata), "trait-metadata row(s) excluded."))
    }
  }
  data <- df[rows, , drop = FALSE]
  value <- suppressWarnings(as.numeric(gsub(",", ".", as.character(data[[response_col]]))))
  if (!any(is.finite(value))) stop("The selected response has no numeric observations.", call. = FALSE)
  factor_values <- lapply(factor_cols, function(column) trimws(as.character(data[[column]])))
  valid <- is.finite(value)
  for (x in factor_values) valid <- valid & !is.na(x) & nzchar(x)
  if (needs_rep) {
    rep_value <- trimws(as.character(data[[replication_col]]))
    valid <- valid & !is.na(rep_value) & nzchar(rep_value)
  }
  if (any(!valid)) {
    notes <- c(notes, paste(sum(!valid), "row(s) excluded for missing response, factor, or replication."))
  }
  data <- data[valid, , drop = FALSE]
  if (nrow(data) == 0) stop("No complete observations remain for this analysis.", call. = FALSE)
  out <- data.frame(Y = value[valid])
  for (i in seq_along(factor_cols)) {
    values <- factor_values[[i]][valid]
    out[[LETTERS[i]]] <- factor(values, levels = unique(values))
    if (nlevels(out[[LETTERS[i]]]) < 2) {
      stop("Factor ", LETTERS[i], " needs at least two observed levels.", call. = FALSE)
    }
  }
  if (needs_rep) {
    out$Rep <- factor(rep_value[valid])
    if (nlevels(out$Rep) < 2) {
      stop("At least two block or whole-plot replicate levels are required.", call. = FALSE)
    }
  }
  factors <- LETTERS[seq_along(factor_cols)]
  observed <- unique(out[, factors, drop = FALSE])
  expected <- prod(vapply(out[factors], nlevels, integer(1)))
  if (nrow(observed) < expected) {
    stop("At least one factor combination is missing. Include every A",
         if (length(factors) == 3) " x B x C" else " x B",
         " combination before running Multi-factor analysis.", call. = FALSE)
  }
  cell_counts <- stats::aggregate(out$Y, out[factors], length)$x
  if (length(unique(cell_counts)) > 1) {
    notes <- c(notes, "Factor combinations have unequal observation counts; estimated means and pairwise thresholds may differ in precision.")
  }
  if (design == "split_plot") {
    unit_counts <- stats::aggregate(rep(1, nrow(out)), out[c("Rep", factors)], length)$x
    if (any(unit_counts > 1)) {
      stop("Split-plot data contain multiple rows for a smallest experimental unit. Aggregate subsamples to one plot value first.", call. = FALSE)
    }
    out$WholePlot <- interaction(out$Rep, out$A, drop = TRUE)
    whole_plot_counts <- stats::aggregate(out$B, list(WholePlot = out$WholePlot),
                                          function(x) length(unique(x)))$x
    if (any(whole_plot_counts < nlevels(out$B))) {
      stop("Every whole plot must contain every Factor B subplot level.", call. = FALSE)
    }
    if (length(factors) == 3) {
      out$SubPlot <- interaction(out$Rep, out$A, out$B, drop = TRUE)
      sub_plot_counts <- stats::aggregate(out$C, list(SubPlot = out$SubPlot),
                                          function(x) length(unique(x)))$x
      if (any(sub_plot_counts < nlevels(out$C))) {
        stop("Every subplot must contain every Factor C sub-subplot level.", call. = FALSE)
      }
    }
  } else if (randomization == "RCBD") {
    per_block <- stats::aggregate(rep(1, nrow(out)), out[c("Rep", factors)], length)
    expected_block_cells <- nlevels(out$Rep) * expected
    if (nrow(per_block) < expected_block_cells) {
      notes <- c(notes, "Some block-by-treatment combinations are missing; inspect the model estimates and standard errors.")
    }
  }
  list(data = out, factor_names = factor_cols, factors = factors,
       response = response_col, replication = if (needs_rep) replication_col else "",
       notes = notes)
}

mf_fit_model <- function(data, factors, design, randomization, model_type) {
  fixed <- paste(factors, collapse = "*")
  if (design == "factorial" && model_type == "ANOVA") {
    formula <- stats::as.formula(paste("Y ~", if (randomization == "RCBD") "Rep +" else "", fixed))
    contrasts <- stats::setNames(lapply(factors, function(x) stats::contr.sum(nlevels(data[[x]]))), factors)
    if (randomization == "RCBD") contrasts$Rep <- stats::contr.sum(nlevels(data$Rep))
    model <- stats::lm(formula, data = data, contrasts = contrasts)
    if (anyNA(stats::coef(model))) stop("The factorial model is rank deficient; check factor combinations and replication.", call. = FALSE)
    test <- as.data.frame(car::Anova(model, type = 3))
    test$Source <- rownames(test)
    return(list(model = model, anova = test, method = "Fixed-effects ANOVA (Type III)"))
  }
  if (design == "factorial" && randomization == "CRD") {
    stop("Factorial LMM needs a block / replication column. Choose ANOVA for CRD.", call. = FALSE)
  }
  random_terms <- if (design == "factorial") {
    "(1 | Rep)"
  } else {
    paste(c(if (randomization == "RCBD") "(1 | Rep)",
            "(1 | WholePlot)",
            if (length(factors) == 3) "(1 | SubPlot)"), collapse = " + ")
  }
  formula <- stats::as.formula(paste("Y ~", fixed, "+", random_terms))
  contrasts <- stats::setNames(lapply(factors, function(x) stats::contr.sum(nlevels(data[[x]]))), factors)
  model <- lmerTest::lmer(formula, data = data, REML = TRUE, contrasts = contrasts)
  if (anyNA(lme4::fixef(model, add.dropped = TRUE))) {
    stop("The mixed model is rank deficient; check factor combinations and replication.", call. = FALSE)
  }
  test <- as.data.frame(stats::anova(model, type = 3, ddf = "Satterthwaite"))
  test$Source <- rownames(test)
  list(model = model, anova = test, method = "LMM Type III (Satterthwaite)")
}

mf_anova_table <- function(test) {
  get_column <- function(names_to_try) {
    column <- intersect(names_to_try, names(test))[1]
    if (is.na(column)) rep(NA_real_, nrow(test)) else suppressWarnings(as.numeric(test[[column]]))
  }
  data.frame(
    Source = as.character(test$Source),
    Num_df = get_column(c("NumDF", "Df")),
    Den_df = get_column(c("DenDF")),
    Sum_Sq = get_column(c("Sum Sq", "Sum_Sq")),
    Mean_Sq = get_column(c("Mean Sq", "Mean_Sq")),
    F_value = get_column(c("F value", "F.value")),
    p_value = get_column(c("Pr(>F)", "Pr(>Chisq)")),
    check.names = FALSE
  ) |> dplyr::filter(Source != "(Intercept)")
}

mf_comparison_tables <- function(model, factors, comparison_factor, method) {
  other <- setdiff(factors, comparison_factor)
  formula <- stats::as.formula(paste("~", comparison_factor, "|", paste(other, collapse = "*")))
  emm <- emmeans::emmeans(model, specs = formula)
  means <- as.data.frame(emm)
  group_cols <- c(comparison_factor, other)
  letter_table <- suppressMessages(as.data.frame(multcomp::cld(
    emm, by = other, Letters = c(base::letters, LETTERS),
    adjust = if (method == "tukey") "tukey" else "none",
    alpha = 0.05, sort = FALSE
  )))
  means$Group <- ""
  key <- function(x) do.call(paste, c(lapply(x[group_cols], as.character), sep = "\r"))
  means$Group <- trimws(as.character(letter_table$.group)[match(key(means), key(letter_table))])
  if (anyNA(means$Group)) stop("Could not align compact letters with the mean estimates.", call. = FALSE)
  names(means)[names(means) == "emmean"] <- "Adjusted_mean"
  names(means)[names(means) == "lower.CL"] <- "Lower_95_CI"
  names(means)[names(means) == "upper.CL"] <- "Upper_95_CI"
  means$Adjusted_mean_SE <- sprintf("%.3f ± %.3f", means$Adjusted_mean, means$SE)
  first <- c(group_cols, "Adjusted_mean_SE", "Group")
  means <- means[, c(first, setdiff(names(means), first)), drop = FALSE]
  pairwise <- as.data.frame(summary(pairs(
    emm, by = other, adjust = if (method == "tukey") "tukey" else "none"
  )))
  names(pairwise)[names(pairwise) == "p.value"] <- "Adjusted_p_value"
  pairwise$Critical_difference <- if (method == "lsd") {
    stats::qt(0.975, pairwise$df) * pairwise$SE
  } else {
    stats::qtukey(0.95, nlevels(means[[comparison_factor]]), pairwise$df) * pairwise$SE / sqrt(2)
  }
  pairwise$Significant_0.05 <- is.finite(pairwise$Adjusted_p_value) & pairwise$Adjusted_p_value < 0.05
  list(means = means, pairwise = pairwise)
}

mf_effect_terms <- function(factors) {
  unlist(lapply(seq_along(factors), function(n) {
    vapply(utils::combn(factors, n, simplify = FALSE), paste, character(1), collapse = ":")
  }), use.names = FALSE)
}

mf_effect_comparisons <- function(model, data, factors, anova, method) {
  terms <- mf_effect_terms(factors)
  tables <- lapply(terms, function(effect) {
    columns <- strsplit(effect, ":", fixed = TRUE)[[1]]
    emm <- suppressMessages(emmeans::emmeans(model, specs = stats::as.formula(paste("~", effect))))
    adjusted <- as.data.frame(emm)
    raw <- data |>
      dplyr::group_by(dplyr::across(dplyr::all_of(columns))) |>
      dplyr::summarise(N = dplyr::n(), Raw_mean = mean(Y), Raw_SD = stats::sd(Y), .groups = "drop")
    adjusted <- dplyr::left_join(adjusted, raw, by = columns)
    p <- anova$p_value[match(effect, anova$Source)]
    if (length(p) == 0 || !is.finite(p)) p <- NA_real_
    adjusted$Group <- ""
    if (is.finite(p) && p < 0.05) {
      letters <- suppressMessages(as.data.frame(multcomp::cld(
        emm, Letters = c(base::letters, LETTERS),
        adjust = if (method == "tukey") "tukey" else "none",
        alpha = 0.05, sort = FALSE
      )))
      key <- function(x) do.call(paste, c(lapply(x[columns], as.character), sep = "\r"))
      adjusted$Group <- trimws(as.character(letters$.group)[match(key(adjusted), key(letters))])
      if (anyNA(adjusted$Group)) stop("Could not align effect comparison letters.", call. = FALSE)
    }
    level <- apply(adjusted[columns], 1, function(x) paste(x, collapse = " × "))
    means <- data.frame(
      Effect = effect, Level = level,
      Raw_mean_SD = sprintf("%.3f ± %.3f", adjusted$Raw_mean, adjusted$Raw_SD),
      Group = adjusted$Group, N = adjusted$N,
      Raw_mean = adjusted$Raw_mean, Raw_SD = adjusted$Raw_SD,
      Adjusted_mean = adjusted$emmean, Adjusted_SE = adjusted$SE,
      Omnibus_p = p, Row_type = "Level", stringsAsFactors = FALSE
    )
    contrast <- as.data.frame(summary(pairs(
      emm, adjust = if (method == "tukey") "tukey" else "none"
    )))
    names(contrast)[names(contrast) == "p.value"] <- "Adjusted_p_value"
    contrast$Effect <- effect
    contrast$Critical_difference <- if (method == "lsd") {
      stats::qt(0.975, contrast$df) * contrast$SE
    } else {
      stats::qtukey(0.95, nrow(adjusted), contrast$df) * contrast$SE / sqrt(2)
    }
    contrast$Significant_0.05 <- is.finite(contrast$Adjusted_p_value) & contrast$Adjusted_p_value < 0.05
    critical <- contrast$Critical_difference[is.finite(contrast$Critical_difference)]
    critical_label <- if (!length(critical)) "Not estimable" else if (diff(range(critical)) < 1e-6) {
      sprintf("%.3f", critical[1])
    } else "Varies by pair; see pairwise table"
    cv <- 100 * stats::sigma(model) / abs(mean(data$Y))
    footer <- data.frame(
      Effect = effect,
      Level = c("Residual CV (%)", "Heritability (h²)",
                if (method == "lsd") "LSD (0.05)" else "Tukey HSD (0.05)"),
      Raw_mean_SD = c(if (is.finite(cv)) sprintf("%.2f%%", cv) else "Not estimable",
                      "Not estimated", critical_label),
      Group = "", N = NA_integer_, Raw_mean = NA_real_, Raw_SD = NA_real_,
      Adjusted_mean = NA_real_, Adjusted_SE = NA_real_, Omnibus_p = p,
      Row_type = "Statistic", stringsAsFactors = FALSE
    )
    list(means = dplyr::bind_rows(means, footer), pairwise = contrast)
  })
  list(means = dplyr::bind_rows(lapply(tables, `[[`, "means")),
       pairwise = dplyr::bind_rows(lapply(tables, `[[`, "pairwise")))
}

mf_superiority_table <- function(means, factors, comparison_factor, checks, direction) {
  other <- setdiff(factors, comparison_factor)
  strata <- interaction(means[other], drop = TRUE, lex.order = TRUE)
  parts <- split(means, strata, drop = TRUE)
  rows <- lapply(parts, function(part) {
    baseline_rows <- part[as.character(part[[comparison_factor]]) %in% checks &
                            is.finite(part$Adjusted_mean), , drop = FALSE]
    candidates <- part[!as.character(part[[comparison_factor]]) %in% checks, , drop = FALSE]
    if (nrow(candidates) == 0) return(NULL)
    baseline <- if (nrow(baseline_rows) > 0) mean(baseline_rows$Adjusted_mean) else NA_real_
    difference <- if (direction == "Higher better") {
      candidates$Adjusted_mean - baseline
    } else {
      baseline - candidates$Adjusted_mean
    }
    out <- candidates[, factors, drop = FALSE]
    out$Candidate_mean <- candidates$Adjusted_mean
    out$Check_baseline <- baseline
    out$Advantage <- difference
    out$Superiority_pct <- if (is.finite(baseline) && abs(baseline) > 1e-8) {
      100 * difference / abs(baseline)
    } else {
      rep(NA_real_, nrow(candidates))
    }
    out$Check_levels <- paste(as.character(baseline_rows[[comparison_factor]]), collapse = ", ")
    out$Note <- if (nrow(baseline_rows) == 0) "No check estimate in this factor combination" else if (abs(baseline) <= 1e-8) "Check baseline is zero; percent unavailable" else ""
    out
  })
  dplyr::bind_rows(rows)
}

run_multifactor_pipeline <- function(df, response_col, factor_cols,
                                     design = c("factorial", "split_plot"),
                                     randomization = c("CRD", "RCBD"),
                                     model_type = c("ANOVA", "LMM"),
                                     replication_col = NULL,
                                     comparison_factor = "A", checks = NULL,
                                     direction = c("Higher better", "Lower better"),
                                     comparison_method = c("tukey", "lsd"),
                                     metadata_col = NULL,
                                     metadata_labels = character(0)) {
  design <- match.arg(design)
  randomization <- match.arg(randomization)
  model_type <- match.arg(model_type)
  direction <- match.arg(direction)
  comparison_method <- match.arg(comparison_method)
  if (design == "split_plot" && model_type != "LMM") {
    stop("Split-plot and split-split-plot designs use the LMM error strata.", call. = FALSE)
  }
  if (design == "factorial" && randomization == "CRD" && model_type == "LMM") {
    stop("Factorial CRD uses ANOVA; select RCBD for a block-random LMM.", call. = FALSE)
  }
  prepared <- mf_prepare_data(df, response_col, factor_cols, replication_col,
                              design, randomization, metadata_col, metadata_labels)
  if (!comparison_factor %in% prepared$factors) {
    stop("Choose a comparison factor among A, B, and C.", call. = FALSE)
  }
  checks <- unique(trimws(as.character(checks)))
  checks <- checks[!is.na(checks) & nzchar(checks)]
  available <- levels(prepared$data[[comparison_factor]])
  if (length(checks) == 0 || !all(checks %in% available)) {
    stop("Choose at least one check level from the comparison factor.", call. = FALSE)
  }
  if (length(setdiff(available, checks)) == 0) {
    stop("At least one non-check candidate level is required.", call. = FALSE)
  }
  fitted <- mf_fit_model(prepared$data, prepared$factors, design, randomization, model_type)
  summary <- prepared$data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(prepared$factors))) |>
    dplyr::summarise(N = dplyr::n(), Raw_mean = mean(Y), Raw_SD = stats::sd(Y),
                     Raw_SE = Raw_SD / sqrt(N), .groups = "drop")
  comparison <- mf_comparison_tables(fitted$model, prepared$factors,
                                     comparison_factor, comparison_method)
  anova <- mf_anova_table(fitted$anova)
  effects <- mf_effect_comparisons(fitted$model, prepared$data, prepared$factors,
                                   anova, comparison_method)
  superiority <- mf_superiority_table(comparison$means, prepared$factors,
                                      comparison_factor, checks, direction)
  notes <- prepared$notes
  notes <- c(notes, "Main-effect comparisons average across the other factors; interpret them alongside interaction tests.")
  if (inherits(fitted$model, "merMod") && lme4::isSingular(fitted$model)) {
    notes <- c(notes, "The mixed model has a singular random-effects fit; inspect variance estimates and interpret tests cautiously.")
  }
  if (anyNA(superiority$Check_baseline)) {
    notes <- c(notes, "Some candidate combinations have no estimable check baseline; their superiority is unavailable.")
  }
  settings <- data.frame(
    Setting = c("Response", "Factor A", "Factor B", if (length(factor_cols) == 3) "Factor C",
                "Design", "Base randomization", "Model", "Block / replicate",
                "Comparison factor", "Checks", "Trait direction", "Mean comparison"),
    Value = c(response_col, factor_cols, design, randomization, model_type,
              prepared$replication, comparison_factor, paste(checks, collapse = ", "),
              direction, comparison_method),
    stringsAsFactors = FALSE
  )
  list(settings = settings, summary = as.data.frame(summary),
       anova = anova, effect_comparison = effects$means,
       effect_pairwise = effects$pairwise,
       mean_comparison = comparison$means, pairwise = comparison$pairwise,
       superiority = superiority,
       notes = data.frame(Note = if (length(notes)) notes else "No data exclusions or design warnings."),
       factor_names = factor_cols, factors = prepared$factors,
       comparison_factor = comparison_factor, response = response_col,
       method_label = fitted$method, model = fitted$model,
       analysis_data = prepared$data, checks = checks)
}

mf_refresh_result <- function(result, method = "tukey", direction = "Higher better",
                              effect = "All effects") {
  comparison <- mf_comparison_tables(result$model, result$factors, "A", method)
  effects <- mf_effect_comparisons(result$model, result$analysis_data, result$factors,
                                   result$anova, method)
  result$mean_comparison <- comparison$means
  result$pairwise <- comparison$pairwise
  result$effect_comparison <- effects$means
  result$effect_pairwise <- effects$pairwise
  result$superiority <- mf_superiority_table(comparison$means, result$factors,
                                              "A", result$checks, direction)
  result$comparison_factor <- "A"
  result$selected_effect <- effect
  result$settings$Value[result$settings$Setting == "Mean comparison"] <- method
  result$settings$Value[result$settings$Setting == "Trait direction"] <- direction
  result
}

mf_selected_effect_table <- function(result, field) {
  table <- result[[field]]
  effect <- if (is.null(result$selected_effect)) "All effects" else result$selected_effect
  if (!identical(effect, "All effects")) table <- table[table$Effect == effect, , drop = FALSE]
  table
}

mf_plot_mean_comparison <- function(result) {
  d <- mf_selected_effect_table(result, "effect_comparison")
  d <- d[d$Row_type == "Level", , drop = FALSE]
  gap <- max(diff(range(c(d$Raw_mean - d$Raw_SD, d$Raw_mean + d$Raw_SD), na.rm = TRUE)) * 0.04, 0.03)
  ggplot2::ggplot(d, ggplot2::aes(x = Level, y = Raw_mean, fill = Effect)) +
    ggplot2::geom_col(width = 0.72, show.legend = FALSE) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Raw_mean - Raw_SD, ymax = Raw_mean + Raw_SD),
                           width = 0.16, na.rm = TRUE) +
    ggplot2::geom_text(ggplot2::aes(y = Raw_mean + Raw_SD + gap, label = Group),
                       fontface = "bold", na.rm = TRUE) +
    ggplot2::facet_wrap(~ Effect, scales = "free_x") +
    ggplot2::scale_fill_manual(values = rep(c("#3498DB", "#2ECC71", "#E74C3C"), length.out = length(unique(d$Effect)))) +
    ggplot2::labs(title = "Multi-factor mean comparison",
                  subtitle = "Bars and error bars: raw mean ± SD; letters: adjusted model comparisons",
                  x = "Factor level or combination", y = result$response) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

mf_plot_superiority <- function(result) {
  d <- result$superiority
  other <- setdiff(result$factors, "A")
  d$Context <- apply(d[other], 1, function(x) paste(x, collapse = " × "))
  d$Label <- ifelse(is.finite(d$Superiority_pct), sprintf("%+.1f%%", d$Superiority_pct), "")
  limit <- max(abs(d$Superiority_pct), na.rm = TRUE)
  if (!is.finite(limit) || limit < 1) limit <- 1
  ggplot2::ggplot(d, ggplot2::aes(x = Context, y = A, fill = Superiority_pct)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5, na.rm = TRUE) +
    ggplot2::geom_text(ggplot2::aes(label = Label), color = "gray15", size = 3) +
    ggplot2::scale_fill_gradient2(low = "#E74C3C", mid = "white", high = "#2ECC71",
                                  midpoint = 0, limits = c(-limit, limit), na.value = "gray90") +
    ggplot2::labs(title = "Superiority versus checks", x = paste(other, collapse = " × "),
                  y = result$factor_names[1], fill = "Advantage (%)") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

mf_plot_interaction <- function(result) {
  d <- result$summary
  effect <- paste(result$factors, collapse = ":")
  labels <- result$effect_comparison
  labels <- labels[labels$Effect == effect & labels$Row_type == "Level", c("Level", "Group"), drop = FALSE]
  d$Level <- apply(d[result$factors], 1, function(x) paste(x, collapse = " × "))
  d <- dplyr::left_join(d, labels, by = "Level")
  gap <- max(diff(range(c(d$Raw_mean - d$Raw_SE, d$Raw_mean + d$Raw_SE), na.rm = TRUE)) * 0.04, 0.03)
  dodge <- ggplot2::position_dodge(width = 0.8)
  ggplot2::ggplot(d, ggplot2::aes(x = A, y = Raw_mean, fill = B)) +
    ggplot2::geom_col(position = dodge, width = 0.72) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Raw_mean - Raw_SE,
                                       ymax = Raw_mean + Raw_SE),
                            position = dodge, width = 0.14, na.rm = TRUE) +
    ggplot2::geom_text(ggplot2::aes(y = Raw_mean + Raw_SE + gap, label = Group),
                       position = dodge, fontface = "bold", na.rm = TRUE) +
    {if ("C" %in% result$factors) ggplot2::facet_wrap(~ C) else NULL} +
    ggplot2::scale_fill_manual(values = rep(c("#3498DB", "#2ECC71", "#E74C3C"), length.out = length(unique(d$B)))) +
    ggplot2::labs(title = "Factor interaction", x = result$factor_names[1],
                  y = paste(result$response, "raw mean ± SE"), fill = result$factor_names[2],
                  subtitle = "Letters from adjusted interaction comparisons when the interaction is significant") +
    ggplot2::theme_bw()
}
