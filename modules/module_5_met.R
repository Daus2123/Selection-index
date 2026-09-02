# Module 5 - Multi-Environment Trial (MET)
#
# This is the single production source for MET preparation, modelling,
# Finlay-Wilkinson stability, AMMI, GGE, and multi-trait selection support.
# Test code remains outside this module so sourcing it cannot execute tests.
MODULE_5_NAME <- "Multi-environment trial (MET) analysis"

# MET package adapter and reporting helpers -------------------------------

si_met_action <- function(recommendation) {
  recommendation <- as.character(recommendation)
  ifelse(grepl("^Advance", recommendation, ignore.case = TRUE), "ADVANCE",
         ifelse(grepl("^Retest|review", recommendation, ignore.case = TRUE), "RETEST", "DISCARD"))
}

si_met_breeder_recommendation <- function(met_results, top_n = 20) {
  ranking <- si_safe_table(if (!is.null(met_results)) met_results$met_integrated_ranking else NULL)
  rec <- si_met_recommendation_table(ranking, top_n = top_n)
  if (nrow(rec) == 0) return(data.frame())

  out <- data.frame(
    Genotype = as.character(rec$Genotype),
    Source = "MET",
    Final_action = si_met_action(rec$MET_recommendation),
    Rank = suppressWarnings(as.numeric(rec$MET_rank)),
    Index = suppressWarnings(as.numeric(rec$MET_index)),
    Index_advantage = NA_real_,
    Weakness_trait = "",
    Reason = as.character(rec$MET_recommendation),
    stringsAsFactors = FALSE
  )
  out[order(match(out$Final_action, c("ADVANCE", "RETEST", "DISCARD")), out$Rank, na.last = TRUE), , drop = FALSE]
}
si_ext_run_metan_met <- function(data, env, gen, rep, trait) {
  si_ext_require("metan")
  d <- data
  d$ENV <- factor(si_ext_chr(d[[env]]))
  d$GEN <- factor(si_ext_chr(d[[gen]]))
  d$REP <- factor(si_ext_chr(d[[rep]]))
  d$TRAIT <- si_ext_num(d[[trait]])
  d <- d[stats::complete.cases(d[, c("ENV", "GEN", "REP", "TRAIT")]), , drop = FALSE]
  if (nrow(d) == 0) stop("No complete ENV/GEN/REP/trait rows are available.", call. = FALSE)

  list(
    analysis = "Expanded MET stability via metan",
    trait = trait,
    descriptive = si_ext_safe(metan::desc_stat(d, stats = "all")),
    means_genotype = si_ext_safe(metan::means_by(d, GEN)),
    means_environment = si_ext_safe(metan::means_by(d, ENV)),
    anova_individual = si_ext_safe(metan::anova_ind(d, ENV, GEN, REP, resp = TRAIT)),
    anova_joint = si_ext_safe(metan::anova_joint(d, ENV, GEN, REP, TRAIT)),
    bartlett = si_ext_safe(stats::bartlett.test(d$TRAIT ~ d$ENV)),
    annicchiarico = si_ext_safe(metan::Annicchiarico(d, ENV, GEN, REP, TRAIT)),
    ecovalence = si_ext_safe(metan::ecovalence(d, ENV, GEN, REP, TRAIT)),
    shukla = si_ext_safe(metan::Shukla(d, ENV, GEN, REP, TRAIT)),
    regression_stability = si_ext_safe(metan::ge_reg(d, ENV, GEN, REP, TRAIT)),
    superiority = si_ext_safe(metan::superiority(d, ENV, GEN, TRAIT)),
    fox = si_ext_safe(metan::Fox(d, ENV, GEN, TRAIT)),
    factor_analysis = si_ext_safe(metan::ge_factanal(d, ENV, GEN, REP, TRAIT)),
    stability_wrap = si_ext_safe(metan::ge_stats(d, ENV, GEN, REP, TRAIT)),
    ammi = si_ext_safe(metan::performs_ammi(d, ENV, GEN, REP, TRAIT)),
    waas = si_ext_safe(metan::waas(d, ENV, GEN, REP, TRAIT)),
    gge_environment = si_ext_safe(metan::gge(d, ENV, GEN, TRAIT, svp = "environment")),
    gge_genotype = si_ext_safe(metan::gge(d, ENV, GEN, TRAIT, svp = "genotype")),
    gge_symmetrical = si_ext_safe(metan::gge(d, ENV, GEN, TRAIT, svp = "symmetrical"))
  )
}

si_met_trait_quality <- function(results) {
  met_by_trait <- results$met_by_trait
  if (is.null(met_by_trait) || length(met_by_trait) == 0) {
    return(data.frame())
  }

  rows <- lapply(names(met_by_trait), function(trait) {
    result <- met_by_trait[[trait]]
    ms <- si_safe_table(result$model_summary)
    top <- si_safe_table(result$met_selection)
    top <- if (nrow(top) > 0) top[1, , drop = FALSE] else top
    genotype_col <- si_first_existing(top, c("Genotype", "GEN", "Original_ID", "ID"))
    score_col <- si_first_existing(top, c("Combined_score", "MET_selection_score", "Selection_score", "Integrated_MET_Index"))

    data.frame(
      Trait = trait,
      N_genotypes = if ("N_genotypes" %in% names(ms)) ms$N_genotypes[1] else NA,
      N_environments = if ("N_environments" %in% names(ms)) ms$N_environments[1] else NA,
      Broad_sense_H2 = if ("Broad_sense_H2" %in% names(ms)) ms$Broad_sense_H2[1] else NA,
      Stability_ratio = if ("Stability_ratio_G_over_G_plus_GxE" %in% names(ms)) ms$Stability_ratio_G_over_G_plus_GxE[1] else NA,
      Controls_used = if ("Controls_used" %in% names(ms)) ms$Controls_used[1] else NA,
      Top_genotype = if (!is.na(genotype_col) && nrow(top) > 0) as.character(top[[genotype_col]][1]) else NA,
      Top_score = if (!is.na(score_col) && nrow(top) > 0) suppressWarnings(as.numeric(top[[score_col]][1])) else NA_real_,
      Notes = if ("Notes" %in% names(ms)) ms$Notes[1] else "",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)

  failures <- si_safe_table(results$met_failed_traits)
  if (nrow(failures) > 0) {
    fail_rows <- data.frame(
      Trait = as.character(failures$Trait),
      N_genotypes = NA,
      N_environments = NA,
      Broad_sense_H2 = NA,
      Stability_ratio = NA,
      Controls_used = NA,
      Top_genotype = NA,
      Top_score = NA,
      Notes = paste("Failed:", failures$Error),
      stringsAsFactors = FALSE
    )
    out <- rbind(out, fail_rows)
  }
  out
}

si_met_pipeline_review <- function(results) {
  quality <- si_met_trait_quality(results)
  if (nrow(quality) == 0) {
    return(data.frame(Module = "MET", Check = "Trait coverage", Result = "No MET traits found",
                      Why_this_matters = "MET needs successful per-trait models before stability ranking.",
                      stringsAsFactors = FALSE))
  }
  successful <- quality[!grepl("^Failed:", quality$Notes %||% "", ignore.case = TRUE), , drop = FALSE]
  failed_n <- nrow(quality) - nrow(successful)
  h2 <- suppressWarnings(as.numeric(successful$Broad_sense_H2))
  stability <- suppressWarnings(as.numeric(successful$Stability_ratio))
  low_conf <- sum(grepl("Single-environment genotypes", successful$Notes, ignore.case = TRUE), na.rm = TRUE)

  data.frame(
    Module = "MET",
    Check = c(
      "Trait model coverage",
      "Median heritability",
      "Median stability ratio",
      "Low-confidence genotype warning",
      "Failed trait count"
    ),
    Result = c(
      paste0(nrow(successful), " successful trait model(s)"),
      if (all(is.na(h2))) "H2 not available" else round(stats::median(h2, na.rm = TRUE), 4),
      if (all(is.na(stability))) "Stability ratio not available" else round(stats::median(stability, na.rm = TRUE), 4),
      paste0(low_conf, " trait(s) flagged single-environment genotypes"),
      failed_n
    ),
    Why_this_matters = c(
      "MET conclusions are only as broad as the traits that fitted successfully.",
      "Higher H2 supports more reliable genotype ranking.",
      "Higher G/(G+GxE) means genotype performance is less dominated by interaction.",
      "Single-environment genotypes can look promising but have weak stability evidence.",
      "Failed traits should be reviewed before using the integrated MET ranking."
    ),
    stringsAsFactors = FALSE
  )
}

si_met_recommendation_table <- function(integrated_ranking, top_n = 20) {
  x <- si_safe_table(integrated_ranking)
  if (nrow(x) == 0) return(data.frame())
  genotype_col <- si_first_existing(x, c("Genotype", "GEN", "Original_ID", "ID"))
  rank_col <- si_first_existing(x, c("Integrated_rank", "Rank", "MET_rank"))
  index_col <- si_first_existing(x, c("Integrated_MET_Index", "MET_Index", "Selection_score"))
  if (is.na(genotype_col)) return(x)
  out <- data.frame(
    Genotype = as.character(x[[genotype_col]]),
    MET_rank = if (!is.na(rank_col)) suppressWarnings(as.numeric(x[[rank_col]])) else seq_len(nrow(x)),
    MET_index = if (!is.na(index_col)) suppressWarnings(as.numeric(x[[index_col]])) else NA_real_,
    stringsAsFactors = FALSE
  )
  out <- out[order(out$MET_rank, -out$MET_index, na.last = TRUE), , drop = FALSE]
  out$MET_recommendation <- ifelse(
    out$MET_rank <= max(1, min(5, nrow(out))),
    "Advance candidate across environments",
    ifelse(out$MET_rank <= top_n, "Retest / review candidate", "Lower priority")
  )
  out
}

# Decision-support helpers -------------------------------------------------
# MET decision-support helpers.
# These functions keep breeder-facing decision logic separate from the Shiny UI.

met_ds_null <- function(a, b) if (!is.null(a)) a else b

met_ds_empty_plot <- function(message_text) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = message_text, size = 4) +
    ggplot2::theme_void()
}

met_ds_standardize <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0) return(numeric(0))
  if (all(is.na(x))) return(rep(0, length(x)))
  sx <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sx) || sx == 0) return(ifelse(is.na(x), NA_real_, 0))
  as.numeric((x - mean(x, na.rm = TRUE)) / sx)
}

met_ds_round_numeric <- function(x, digits = 4) {
  if (!is.data.frame(x) || nrow(x) == 0) return(x)
  num_cols <- vapply(x, is.numeric, logical(1))
  x[num_cols] <- lapply(x[num_cols], round, digits = digits)
  x
}

met_result_uses_anova <- function(result) {
  model_summary <- tryCatch(as.data.frame(result$model_summary), error = function(e) data.frame())
  model_name <- if (nrow(model_summary) > 0 && "Model" %in% names(model_summary)) {
    as.character(model_summary$Model[1])
  } else {
    "LMM"
  }
  isTRUE(grepl("^ANOVA", model_name, ignore.case = TRUE))
}

met_genotype_estimates_table <- function(result) {
  estimates <- tryCatch(as.data.frame(result$blups_main), error = function(e) data.frame())
  if (ncol(estimates) == 0) return(estimates)

  if (met_result_uses_anova(result)) {
    rename_map <- c(
      BLUP_G = "BLUE_G",
      SE_G = "SE_BLUE",
      CI_lower = "BLUE_CI_Lower",
      CI_upper = "BLUE_CI_Upper",
      Rank_BLUP = "Rank_BLUE"
    )
    estimates$Reliability <- NULL
    preferred_order <- c(
      "Genotype", "Raw_Mean", "SE_Raw_Mean", "BLUE_G", "SE_BLUE",
      "BLUE_CI_Lower", "BLUE_CI_Upper", "Trait_direction", "Rank_BLUE"
    )
  } else {
    rename_map <- c(
      SE_G = "SE_BLUP",
      CI_lower = "BLUP_CI_Lower",
      CI_upper = "BLUP_CI_Upper"
    )
    preferred_order <- c(
      "Genotype", "Raw_Mean", "SE_Raw_Mean", "BLUP_G", "SE_BLUP",
      "Reliability", "BLUP_CI_Lower", "BLUP_CI_Upper", "Trait_direction", "Rank_BLUP"
    )
  }

  for (old_name in names(rename_map)) {
    if (old_name %in% names(estimates)) {
      names(estimates)[names(estimates) == old_name] <- unname(rename_map[[old_name]])
    }
  }
  estimates[, c(intersect(preferred_order, names(estimates)), setdiff(names(estimates), preferred_order)), drop = FALSE]
}

met_location_estimates_table <- function(result) {
  estimates <- tryCatch(as.data.frame(result$blups_environment), error = function(e) data.frame())
  if (nrow(estimates) == 0) return(estimates)

  trait_name <- tryCatch(as.character(result$model_summary$Trait_used[1]), error = function(e) NA_character_)
  estimate_method <- if (met_result_uses_anova(result)) "BLUE" else "BLUP"
  estimated_value_column <- if ("Reportable_Estimate" %in% names(estimates)) "Reportable_Estimate" else "BLUP_env"
  value_or_na <- function(column, mode = c("numeric", "character", "integer")) {
    mode <- match.arg(mode)
    if (!column %in% names(estimates)) {
      return(switch(mode,
        numeric = rep(NA_real_, nrow(estimates)),
        integer = rep(NA_integer_, nrow(estimates)),
        character = rep(NA_character_, nrow(estimates))
      ))
    }
    switch(mode,
      numeric = suppressWarnings(as.numeric(estimates[[column]])),
      integer = suppressWarnings(as.integer(estimates[[column]])),
      character = as.character(estimates[[column]])
    )
  }

  data.frame(
    Trait = trait_name,
    Genotype = value_or_na("Genotype", "character"),
    Location = value_or_na("Environment", "character"),
    Estimate_Method = estimate_method,
    Cell_Status = value_or_na("Cell_Status", "character"),
    Observed_Mean = value_or_na("Observed_Mean"),
    N_Replications = value_or_na("N_Replications", "integer"),
    Estimated_Value = value_or_na(estimated_value_column),
    Prediction_SE = value_or_na("Prediction_SE"),
    Prediction_CI_Lower = value_or_na("Prediction_CI_Lower"),
    Prediction_CI_Upper = value_or_na("Prediction_CI_Upper"),
    Uncertainty_Method = value_or_na("Uncertainty_Method", "character"),
    Benchmark_Check = value_or_na("Benchmark_Check", "character"),
    Benchmark_Estimate = value_or_na("Benchmark_Estimate"),
    Check_Advantage = value_or_na("Check_Advantage"),
    Check_Advantage_CI_Lower = value_or_na("Check_Advantage_CI_Lower"),
    Check_Advantage_CI_Upper = value_or_na("Check_Advantage_CI_Upper"),
    Probability_Superior = value_or_na("Probability_Superior"),
    Evidence_Flag = value_or_na("Evidence_Flag", "character"),
    Trait_Direction = value_or_na("Trait_Direction", "character"),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(Location, dplyr::desc(Estimated_Value), Genotype) %>%
    met_ds_round_numeric()
}

met_ds_rank_confidence <- function(estimates, direction = "Higher better", target_value = NA_real_,
                                   top_k = 5L, simulations = 2000L, seed = 24051L) {
  estimates <- as.data.frame(estimates)
  required <- c("Genotype", "Genotype_estimate", "Estimate_SE")
  if (nrow(estimates) == 0 || !all(required %in% names(estimates))) {
    return(data.frame(Genotype = character(), Rank_confidence = numeric(), stringsAsFactors = FALSE))
  }
  means <- suppressWarnings(as.numeric(estimates$Genotype_estimate))
  ses <- suppressWarnings(as.numeric(estimates$Estimate_SE))
  valid <- is.finite(means) & is.finite(ses) & ses >= 0
  confidence <- rep(NA_real_, nrow(estimates))
  if (!any(valid)) {
    return(data.frame(Genotype = as.character(estimates$Genotype), Rank_confidence = confidence, stringsAsFactors = FALSE))
  }

  top_k <- min(max(1L, suppressWarnings(as.integer(top_k)[1])), sum(valid))
  simulations <- max(250L, suppressWarnings(as.integer(simulations)[1]))
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)

  draws <- vapply(which(valid), function(i) {
    if (ses[i] == 0) rep(means[i], simulations) else stats::rnorm(simulations, means[i], ses[i])
  }, numeric(simulations))
  if (is.null(dim(draws))) draws <- matrix(draws, ncol = 1)
  scores <- met_ds_direction_score(draws, direction, target_value)
  dim(scores) <- dim(draws)
  simulated_ranks <- t(apply(scores, 1, function(row) rank(-row, ties.method = "average")))
  confidence[valid] <- colMeans(simulated_ranks <= top_k)
  data.frame(
    Genotype = as.character(estimates$Genotype),
    Rank_confidence = confidence,
    stringsAsFactors = FALSE
  )
}

met_ds_direction_score <- function(value, direction = "Higher better", target_value = NA_real_) {
  value <- suppressWarnings(as.numeric(value))
  direction <- as.character(met_ds_null(direction, "Higher better"))[1]
  target_value <- suppressWarnings(as.numeric(target_value)[1])
  if (identical(direction, "Lower better")) {
    return(-value)
  }
  if (identical(direction, "Target trait") && is.finite(target_value)) {
    return(-abs(value - target_value))
  }
  value
}

met_ds_safe_table <- function(x, label = "object") {
  if (is.null(x)) {
    return(data.frame(Status = "Not available", Note = paste(label, "was not generated"), stringsAsFactors = FALSE))
  }
  if (inherits(x, "si_ext_error")) {
    note <- as.character(met_ds_null(x$error, paste(label, "failed")))
    return(data.frame(Status = "Not available", Note = note, stringsAsFactors = FALSE))
  }
  out <- tryCatch(as.data.frame(x), error = function(e) NULL)
  if (!is.null(out) && ncol(out) > 0) {
    return(out)
  }
  if (is.list(x)) {
    component_names <- names(x)
    if (is.null(component_names)) component_names <- paste0("component_", seq_along(x))
    return(data.frame(
      Component = component_names,
      Class = vapply(x, function(v) paste(class(v), collapse = "/"), character(1)),
      Rows = vapply(x, function(v) if (is.data.frame(v)) nrow(v) else NA_integer_, integer(1)),
      Columns = vapply(x, function(v) if (is.data.frame(v)) ncol(v) else NA_integer_, integer(1)),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(Value = as.character(x), stringsAsFactors = FALSE)
}

met_ds_threshold_defaults <- function(n_envs = NA_integer_, min_envs = NA_integer_) {
  n_envs <- suppressWarnings(as.integer(n_envs)[1])
  min_envs <- suppressWarnings(as.integer(min_envs)[1])
  min_coverage <- if (is.finite(n_envs) && n_envs > 0 && is.finite(min_envs)) {
    round(100 * min_envs / n_envs, 1)
  } else {
    50
  }
  list(
    min_reliability = 0.50,
    min_coverage_pct = min(max(min_coverage, 0), 100),
    min_check_advantage = 0,
    selection_pct = 20
  )
}

met_ds_clean_checks <- function(checks, genotypes = character(0)) {
  checks <- unique(trimws(as.character(met_ds_null(checks, character(0)))))
  checks <- checks[!is.na(checks) & checks != ""]
  if (length(genotypes) > 0) {
    checks <- intersect(checks, as.character(genotypes))
  }
  checks
}

met_ds_trait_matrix <- function(df_raw, met_results) {
  prepared <- prepare_met_trait_settings(df_raw)
  traits <- intersect(as.character(prepared$trait_cols), names(met_results))
  traits <- traits[vapply(traits, function(tr) {
    !is.null(met_results[[tr]]$blups_main) && nrow(as.data.frame(met_results[[tr]]$blups_main)) > 0
  }, logical(1))]
  if (length(traits) == 0) {
    return(list(
      traits = character(0),
      trait_info = data.frame(),
      long = data.frame(),
      raw_wide = data.frame(),
      adjusted_wide = data.frame(),
      standardized_wide = data.frame(),
      weights = numeric(0)
    ))
  }

  long <- dplyr::bind_rows(lapply(traits, function(tr) {
    result <- met_results[[tr]]
    blup <- as.data.frame(result$blups_main)
    direction <- as.character(met_ds_null(result$trait_direction, met_ds_null(prepared$trait_direction[[tr]], "Higher better")))
    target <- suppressWarnings(as.numeric(met_ds_null(result$target_value, met_ds_null(prepared$target_traits[[tr]], NA_real_))))
    raw_weight <- suppressWarnings(as.numeric(met_ds_null(result$trait_weight, met_ds_null(prepared$weights_raw_used[[tr]], 1))))
    data.frame(
      Genotype = as.character(blup$Genotype),
      Trait = tr,
      Raw_BLUP = suppressWarnings(as.numeric(blup$BLUP_G)),
      Adjusted_score = met_ds_direction_score(blup$BLUP_G, direction, target),
      Reliability = suppressWarnings(as.numeric(blup$Reliability)),
      Trait_direction = direction,
      Target_value = ifelse(identical(direction, "Target trait"), target, NA_real_),
      Raw_weight = raw_weight,
      stringsAsFactors = FALSE
    )
  }))
  long <- long %>%
    dplyr::group_by(Trait) %>%
    dplyr::mutate(
      Standardized_score = met_ds_standardize(Adjusted_score),
      Trait_rank = rank(-Adjusted_score, ties.method = "average", na.last = "keep")
    ) %>%
    dplyr::ungroup()

  weights_raw <- prepared$weights_raw_used[traits]
  weights_raw[is.na(weights_raw)] <- 1
  weights_raw[weights_raw < 0] <- 0
  if (sum(weights_raw, na.rm = TRUE) <= 0) weights_raw <- stats::setNames(rep(1, length(traits)), traits)
  weights <- weights_raw / sum(weights_raw, na.rm = TRUE)

  raw_wide <- long %>%
    dplyr::select(Genotype, Trait, Raw_BLUP) %>%
    tidyr::pivot_wider(names_from = Trait, values_from = Raw_BLUP)
  adjusted_wide <- long %>%
    dplyr::select(Genotype, Trait, Adjusted_score) %>%
    tidyr::pivot_wider(names_from = Trait, values_from = Adjusted_score)
  standardized_wide <- long %>%
    dplyr::select(Genotype, Trait, Standardized_score) %>%
    tidyr::pivot_wider(names_from = Trait, values_from = Standardized_score)

  list(
    traits = traits,
    trait_info = prepared$trait_info[match(traits, prepared$trait_info$Trait), , drop = FALSE],
    long = long,
    raw_wide = raw_wide,
    adjusted_wide = adjusted_wide,
    standardized_wide = standardized_wide,
    weights = weights
  )
}

met_ds_biplot <- function(score_wide, title, subtitle = "", arrow_label = "Trait") {
  if (is.null(score_wide) || nrow(score_wide) < 2) {
    return(list(plot = met_ds_empty_plot(paste(title, "needs at least two genotypes.")), genotype_scores = data.frame(), trait_scores = data.frame()))
  }
  trait_cols <- setdiff(names(score_wide), "Genotype")
  if (length(trait_cols) < 2) {
    return(list(plot = met_ds_empty_plot(paste(title, "needs at least two traits.")), genotype_scores = data.frame(), trait_scores = data.frame()))
  }
  mat <- as.matrix(score_wide[, trait_cols, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- as.character(score_wide$Genotype)
  mat[!is.finite(mat)] <- 0
  if (sum(abs(mat), na.rm = TRUE) <= sqrt(.Machine$double.eps)) {
    return(list(plot = met_ds_empty_plot(paste(title, "has no usable variation.")), genotype_scores = data.frame(), trait_scores = data.frame()))
  }
  fit <- stats::prcomp(mat, center = FALSE, scale. = FALSE)
  var_pct <- fit$sdev^2 / sum(fit$sdev^2)
  pcs <- min(2, ncol(fit$x), ncol(fit$rotation))
  geno <- as.data.frame(fit$x[, seq_len(pcs), drop = FALSE])
  names(geno) <- paste0("PC", seq_len(pcs))
  if (pcs == 1) geno$PC2 <- 0
  geno$Genotype <- rownames(mat)
  trait <- as.data.frame(fit$rotation[, seq_len(pcs), drop = FALSE])
  names(trait) <- paste0("PC", seq_len(pcs))
  if (pcs == 1) trait$PC2 <- 0
  trait[[arrow_label]] <- rownames(fit$rotation)
  scale_factor <- max(abs(c(geno$PC1, geno$PC2)), na.rm = TRUE)
  if (!is.finite(scale_factor) || scale_factor == 0) scale_factor <- 1
  trait$PC1 <- trait$PC1 * scale_factor * 0.85
  trait$PC2 <- trait$PC2 * scale_factor * 0.85
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = trait,
      ggplot2::aes(x = 0, y = 0, xend = PC1, yend = PC2),
      arrow = ggplot2::arrow(length = grid::unit(0.22, "cm"), type = "closed"),
      color = "#2E7D32",
      linewidth = 0.7
    ) +
    ggplot2::geom_text(
      data = trait,
      ggplot2::aes(x = PC1 * 1.12, y = PC2 * 1.12, label = .data[[arrow_label]]),
      color = "#2E7D32",
      size = 3.3,
      fontface = "bold"
    ) +
    ggplot2::geom_point(data = geno, ggplot2::aes(x = PC1, y = PC2), color = "#244F9E", size = 2.6) +
    ggplot2::geom_text(data = geno, ggplot2::aes(x = PC1, y = PC2, label = Genotype), color = "#244F9E", vjust = -0.8, size = 2.8) +
    ggplot2::geom_hline(yintercept = 0, color = "gray55", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 0, color = "gray55", linewidth = 0.4) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = paste0("PC1 (", round(var_pct[1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(ifelse(length(var_pct) >= 2, var_pct[2], 0) * 100, 1), "%)")
    ) +
    ggplot2::theme_bw()
  list(
    plot = p,
    genotype_scores = met_ds_round_numeric(geno[, c("Genotype", "PC1", "PC2"), drop = FALSE]),
    trait_scores = met_ds_round_numeric(trait)
  )
}


build_met_ammi_support <- function(gxe_matrix, ammi_genotype = NULL, metan_single = NULL,
                                   trait_name = "Trait") {
  pc_table <- data.frame()
  waas_fallback <- data.frame()
  if (!is.null(gxe_matrix) && nrow(as.data.frame(gxe_matrix)) >= 2 && ncol(as.data.frame(gxe_matrix)) >= 2) {
    x <- as.matrix(gxe_matrix)
    storage.mode(x) <- "double"
    row_means <- rowMeans(x, na.rm = TRUE)
    col_means <- colMeans(x, na.rm = TRUE)
    grand <- mean(x, na.rm = TRUE)
    centered <- x - outer(row_means, rep(1, ncol(x))) - outer(rep(1, nrow(x)), col_means) + grand
    fit <- svd(centered)
    k <- min(nrow(centered) - 1, ncol(centered) - 1, length(fit$d))
    if (k >= 1) {
      pc_ss <- fit$d[seq_len(k)]^2
      total_ss <- sum(centered^2, na.rm = TRUE)
      pc_table <- data.frame(
        Source = paste0("IPCA", seq_len(k)),
        Sum_sq = pc_ss,
        Percent_GxE_SS = if (total_ss > 0) 100 * pc_ss / total_ss else NA_real_,
        Cumulative_percent = if (total_ss > 0) 100 * cumsum(pc_ss) / total_ss else NA_real_,
        Interpretation = ifelse(seq_len(k) <= 2, "Shown in AMMI biplot", "Higher-order interaction not shown in 2D chart"),
        stringsAsFactors = FALSE
      )
    }
  }
  ammi_genotype <- as.data.frame(met_ds_null(ammi_genotype, data.frame()))
  pc_cols <- grep("^PC[0-9]+$", names(ammi_genotype), value = TRUE)
  if (length(pc_cols) > 0 && nrow(pc_table) > 0) {
    weights <- pc_table$Percent_GxE_SS[seq_along(pc_cols)]
    weights <- weights / sum(weights, na.rm = TRUE)
    waas_score <- rowSums(abs(as.matrix(ammi_genotype[, pc_cols, drop = FALSE])) * rep(weights, each = nrow(ammi_genotype)), na.rm = TRUE)
    waas_fallback <- data.frame(
      Genotype = as.character(ammi_genotype$Genotype),
      WAAS_like = waas_score,
      WAAS_like_rank = rank(waas_score, ties.method = "average", na.last = "keep"),
      Source = "Calculated from AMMI IPCA scores when metan WAAS is unavailable",
      stringsAsFactors = FALSE
    )
  }
  list(
    ammi_pc_anova = met_ds_round_numeric(pc_table),
    ammi_anova = met_ds_safe_table(if (!is.null(metan_single)) metan_single$anova_joint else NULL, "metan::anova_joint"),
    ammi_metan = met_ds_safe_table(if (!is.null(metan_single)) metan_single$ammi else NULL, "metan::performs_ammi"),
    waas = met_ds_safe_table(if (!is.null(metan_single)) metan_single$waas else NULL, "metan::waas"),
    waas_fallback = met_ds_round_numeric(waas_fallback)
  )
}

run_met_single_trait_extensions <- function(dat_clean) {
  if (!requireNamespace("metan", quietly = TRUE) || !exists("si_ext_run_metan_met")) {
    return(structure(list(error = "Package metan is not available."), class = "si_ext_error"))
  }
  tryCatch(
    si_ext_run_metan_met(dat_clean, env = "Environment", gen = "Genotype", rep = "Rep", trait = "Weight"),
    error = function(e) structure(list(error = conditionMessage(e)), class = "si_ext_error")
  )
}

run_met_multitrait_extensions <- function(df_raw, trait_cols, replication_col, selection_intensity = 20) {
  if (!requireNamespace("metan", quietly = TRUE) || !exists("si_ext_run_multitrait_selection")) {
    return(structure(list(error = "Package metan is not available."), class = "si_ext_error"))
  }
  prepared <- prepare_met_trait_settings(df_raw)
  traits <- intersect(as.character(trait_cols), as.character(prepared$trait_cols))
  traits <- traits[traits %in% names(prepared$data)]
  if (length(traits) < 2) {
    return(structure(list(error = "Multi-trait methods need at least two successful traits."), class = "si_ext_error"))
  }
  if (is.null(replication_col) || !replication_col %in% names(prepared$data)) {
    return(structure(list(error = "Multi-trait methods need a replication column."), class = "si_ext_error"))
  }
  goals <- ifelse(as.character(prepared$trait_direction[traits]) == "Lower better", "l", "h")
  names(goals) <- traits
  weights <- prepared$weights_raw_used[traits]
  high_traits <- traits[as.character(prepared$trait_direction[traits]) != "Lower better"]
  candidate_traits <- if (length(high_traits) > 0) high_traits else traits
  yield_trait <- candidate_traits[which.max(weights[candidate_traits])]
  tryCatch(
    si_ext_run_multitrait_selection(
      prepared$data,
      env = "Environment",
      gen = "Genotype",
      rep = replication_col,
      traits = traits,
      goals = goals,
      selection_intensity = selection_intensity,
      yield_trait = yield_trait,
      economic_weights = weights
    ),
    error = function(e) structure(list(error = conditionMessage(e)), class = "si_ext_error")
  )
}

build_met_multitrait_methods <- function(df_raw, met_results, integrated_ranking,
                                         metan_multitrait = NULL, selection_pct = 20) {
  matrix_info <- met_ds_trait_matrix(df_raw, met_results)
  traits <- matrix_info$traits
  if (length(traits) < 2 || nrow(matrix_info$standardized_wide) < 2) {
    empty <- data.frame(Status = "Not available", Note = "Multi-trait methods need at least two successful MET traits.", stringsAsFactors = FALSE)
    return(list(
      gt_table = empty,
      gt_trait_scores = data.frame(),
      gyt_table = empty,
      gyt_scores = data.frame(),
      mgidi = empty,
      waasb = empty,
      method_membership = empty,
      method_agreement = empty,
      metan_mgidi = met_ds_safe_table(if (!is.null(metan_multitrait)) metan_multitrait$mgidi else NULL, "metan::mgidi"),
      metan_waasb = met_ds_safe_table(if (!is.null(metan_multitrait)) metan_multitrait$waasb else NULL, "metan::waasb"),
      plot_gt = met_ds_empty_plot("GT biplot needs at least two MET traits."),
      plot_gyt = met_ds_empty_plot("GYT biplot needs at least two MET traits."),
      plot_method_agreement = met_ds_empty_plot("Method agreement needs at least two MET traits.")
    ))
  }

  weights <- matrix_info$weights
  std_wide <- matrix_info$standardized_wide
  trait_info <- matrix_info$trait_info
  gt <- met_ds_biplot(
    std_wide,
    "GT Biplot - direction-adjusted MET trait scores",
    "Every trait is transformed so larger values are favorable before plotting.",
    "Trait"
  )
  gt_table <- matrix_info$long %>%
    dplyr::mutate(
      Normalized_weight = as.numeric(weights[Trait]),
      Direction_adjusted_rank = rank(-Adjusted_score, ties.method = "average", na.last = "keep")
    ) %>%
    dplyr::arrange(Trait, Direction_adjusted_rank)

  high_traits <- trait_info$Trait[as.character(trait_info$Direction) != "Lower better"]
  candidate_traits <- if (length(high_traits) > 0) high_traits else traits
  anchor <- candidate_traits[which.max(weights[candidate_traits])]
  if (is.na(anchor) || !anchor %in% traits) anchor <- traits[which.max(weights)]
  combo_traits <- setdiff(traits, anchor)
  gyt_wide <- data.frame(Genotype = std_wide$Genotype, stringsAsFactors = FALSE)
  for (tr in combo_traits) {
    combo_name <- paste(anchor, tr, sep = "_x_")
    gyt_wide[[combo_name]] <- rowMeans(cbind(std_wide[[anchor]], std_wide[[tr]]), na.rm = TRUE)
  }
  if (length(combo_traits) == 0) {
    gyt_wide[[paste0(anchor, "_single_trait")]] <- std_wide[[anchor]]
  }
  gyt <- met_ds_biplot(
    gyt_wide,
    "GYT Biplot - direction-aware trait combinations",
    paste0("Anchor trait: ", anchor, ". Larger combination scores are favorable."),
    "GYT"
  )
  gyt_table <- gyt_wide %>%
    tidyr::pivot_longer(-Genotype, names_to = "GYT_combination", values_to = "GYT_score") %>%
    dplyr::group_by(GYT_combination) %>%
    dplyr::mutate(GYT_rank = rank(-GYT_score, ties.method = "average", na.last = "keep")) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(GYT_combination, GYT_rank)

  mat <- as.matrix(std_wide[, traits, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- std_wide$Genotype
  mat[!is.finite(mat)] <- NA_real_
  ideal <- apply(mat, 2, max, na.rm = TRUE)
  ideal[!is.finite(ideal)] <- 0
  weighted_distance <- vapply(seq_len(nrow(mat)), function(i) {
    available <- is.finite(mat[i, ]) & is.finite(ideal) & is.finite(weights)
    if (!any(available)) return(NA_real_)
    sqrt(sum(weights[available] * (ideal[available] - mat[i, available])^2, na.rm = TRUE) / sum(weights[available], na.rm = TRUE))
  }, numeric(1))
  mgidi <- data.frame(
    Genotype = rownames(mat),
    MGIDI_distance = weighted_distance,
    MGIDI_rank = rank(weighted_distance, ties.method = "average", na.last = "keep"),
    MGIDI_source = "Direction-aware weighted distance to the favorable ideotype",
    stringsAsFactors = FALSE
  ) %>% dplyr::arrange(MGIDI_rank)

  waas_rows <- dplyr::bind_rows(lapply(traits, function(tr) {
    result <- met_results[[tr]]
    waas_source <- as.data.frame(met_ds_null(result$ammi_support$waas_fallback, data.frame()))
    if (nrow(waas_source) == 0 && !is.null(result$ammi_genotype) && nrow(as.data.frame(result$ammi_genotype)) > 0) {
      waas_source <- as.data.frame(result$ammi_genotype) %>%
        dplyr::transmute(Genotype = as.character(Genotype), WAAS_like = suppressWarnings(as.numeric(ASV)))
    }
    if (nrow(waas_source) == 0) return(data.frame())
    waas_source %>%
      dplyr::transmute(
        Genotype = as.character(Genotype),
        Trait = tr,
        Stability_proxy_raw = suppressWarnings(as.numeric(WAAS_like)),
        Stability_proxy_score = met_ds_standardize(-suppressWarnings(as.numeric(WAAS_like)))
      )
  }))
  waasb <- if (nrow(waas_rows) > 0) {
    waas_rows %>%
      dplyr::mutate(Normalized_weight = as.numeric(weights[Trait])) %>%
      dplyr::group_by(Genotype) %>%
      dplyr::summarise(
        Stability_proxy_index = sum(Stability_proxy_score * Normalized_weight, na.rm = TRUE) /
          sum(Normalized_weight[!is.na(Stability_proxy_score)], na.rm = TRUE),
        N_traits_with_stability = sum(!is.na(Stability_proxy_score)),
        .groups = "drop"
      ) %>%
      dplyr::mutate(Stability_proxy_rank = rank(-Stability_proxy_index, ties.method = "average", na.last = "keep")) %>%
      dplyr::arrange(Stability_proxy_rank)
  } else {
    data.frame(Status = "Not available", Note = "WAASB-like stability could not be calculated.", stringsAsFactors = FALSE)
  }

  integrated <- as.data.frame(met_ds_null(integrated_ranking, data.frame()))
  n_gen <- nrow(std_wide)
  top_n <- max(1, ceiling(n_gen * suppressWarnings(as.numeric(selection_pct)[1]) / 100))
  if (!is.finite(top_n)) top_n <- max(1, ceiling(n_gen * 0.20))
  gyt_index <- gyt_wide %>%
    dplyr::mutate(GYT_index = rowMeans(dplyr::across(-Genotype), na.rm = TRUE)) %>%
    dplyr::arrange(dplyr::desc(GYT_index))
  gge_score <- dplyr::bind_rows(lapply(traits, function(tr) {
    tab <- as.data.frame(met_ds_null(met_results[[tr]]$gge_mean_stability_decision, data.frame()))
    if (!all(c("Genotype", "Direction_adjusted_mean_rank", "Stability_rank") %in% names(tab))) return(data.frame())
    tab %>%
      dplyr::transmute(
        Genotype = as.character(Genotype),
        Trait = tr,
        GGE_score = -rowMeans(cbind(Direction_adjusted_mean_rank, Stability_rank), na.rm = TRUE),
        Weight = as.numeric(weights[tr])
      )
  })) %>%
    dplyr::group_by(Genotype) %>%
    dplyr::summarise(GGE_broad_score = sum(GGE_score * Weight, na.rm = TRUE) / sum(Weight, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(GGE_broad_score))

  method_sets <- list(
    Integrated = if ("Integrated_rank" %in% names(integrated)) as.character(integrated$Genotype[order(integrated$Integrated_rank)][seq_len(min(top_n, nrow(integrated)))]) else character(0),
    MGIDI = as.character(mgidi$Genotype[seq_len(min(top_n, nrow(mgidi)))]),
    GYT = as.character(gyt_index$Genotype[seq_len(min(top_n, nrow(gyt_index)))]),
    Stability_proxy = if (all(c("Genotype", "Stability_proxy_rank") %in% names(waasb))) as.character(waasb$Genotype[order(waasb$Stability_proxy_rank)][seq_len(min(top_n, nrow(waasb)))]) else character(0),
    GGE_Broad = if (nrow(gge_score) > 0) as.character(gge_score$Genotype[seq_len(min(top_n, nrow(gge_score)))]) else character(0)
  )
  all_genotypes <- sort(unique(c(std_wide$Genotype, unlist(method_sets, use.names = FALSE))))
  membership <- data.frame(Genotype = all_genotypes, stringsAsFactors = FALSE)
  for (nm in names(method_sets)) {
    membership[[nm]] <- membership$Genotype %in% method_sets[[nm]]
  }
  method_cols <- names(method_sets)
  membership$Method_support_count <- rowSums(membership[, method_cols, drop = FALSE])
  membership$Agreement_pattern <- apply(membership[, method_cols, drop = FALSE], 1, function(x) {
    hit <- names(x)[as.logical(x)]
    if (length(hit) == 0) "No top-method support" else paste(hit, collapse = " + ")
  })
  membership$Method_recommendation <- dplyr::case_when(
    membership$Method_support_count >= 3 ~ "Strong multi-method support",
    membership$Method_support_count == 2 ~ "Moderate multi-method support",
    membership$Method_support_count == 1 ~ "Single-method support; review",
    TRUE ~ "Low method support"
  )
  membership <- membership[order(-membership$Method_support_count, membership$Genotype), , drop = FALSE]

  agreement <- expand.grid(Method1 = method_cols, Method2 = method_cols, stringsAsFactors = FALSE)
  agreement$Coincidence_pct <- vapply(seq_len(nrow(agreement)), function(i) {
    a <- method_sets[[agreement$Method1[i]]]
    b <- method_sets[[agreement$Method2[i]]]
    if (top_n == 0) return(NA_real_)
    100 * length(intersect(a, b)) / top_n
  }, numeric(1))
  p_agreement <- ggplot2::ggplot(agreement, ggplot2::aes(x = Method1, y = Method2, fill = Coincidence_pct)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = paste0(round(Coincidence_pct), "%")), size = 3.5) +
    ggplot2::scale_fill_gradient(low = "#FDEBD0", high = "#1D9E75", limits = c(0, 100), na.value = "gray85") +
    ggplot2::labs(title = "MET method agreement", subtitle = paste0("Top ", top_n, " genotype(s) per method"), x = NULL, y = NULL, fill = "Coincidence") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

  list(
    gt_table = met_ds_round_numeric(gt_table),
    gt_trait_scores = met_ds_round_numeric(gt$trait_scores),
    gyt_table = met_ds_round_numeric(gyt_table),
    gyt_scores = met_ds_round_numeric(gyt$trait_scores),
    mgidi = met_ds_round_numeric(mgidi),
    waasb = met_ds_round_numeric(waasb),
    method_membership = met_ds_round_numeric(membership),
    method_agreement = met_ds_round_numeric(agreement),
    metan_mgidi = met_ds_safe_table(if (!is.null(metan_multitrait)) metan_multitrait$mgidi else NULL, "metan::mgidi"),
    metan_waasb = met_ds_safe_table(if (!is.null(metan_multitrait)) metan_multitrait$waasb else NULL, "metan::waasb"),
    plot_gt = gt$plot,
    plot_gyt = gyt$plot,
    plot_method_agreement = p_agreement
  )
}

build_met_decision_board <- function(df_raw, met_results, integrated_ranking, method_membership = NULL,
                                     check_genotypes = NULL, thresholds = NULL, primary_trait = NULL) {
  integrated <- as.data.frame(met_ds_null(integrated_ranking, data.frame()))
  if (nrow(integrated) == 0 || !"Genotype" %in% names(integrated)) {
    empty <- data.frame(Status = "Not available", Note = "Decision board needs an integrated MET ranking.", stringsAsFactors = FALSE)
    return(list(board = empty, passport = empty, thresholds = data.frame()))
  }
  matrix_info <- met_ds_trait_matrix(df_raw, met_results)
  traits <- matrix_info$traits
  n_envs <- suppressWarnings(max(vapply(met_results, function(r) {
    if (!is.null(r$presence)) length(setdiff(names(as.data.frame(r$presence)), c("Genotype", "n_envs"))) else NA_integer_
  }, integer(1)), na.rm = TRUE))
  min_envs <- suppressWarnings(max(vapply(met_results, function(r) {
    as.integer(met_ds_null(r$model_summary$AMMI_GGE_min_observed_locations[1], NA_integer_))
  }, integer(1)), na.rm = TRUE))
  defaults <- met_ds_threshold_defaults(n_envs, min_envs)
  thresholds <- modifyList(defaults, met_ds_null(thresholds, list()))

  genotypes <- as.character(integrated$Genotype)
  checks <- met_ds_clean_checks(check_genotypes, genotypes)
  if (is.null(primary_trait) || length(primary_trait) == 0 || !primary_trait[1] %in% traits) {
    primary_trait <- if (length(traits) > 0) traits[1] else NA_character_
  } else {
    primary_trait <- as.character(primary_trait[1])
  }
  check_index <- if (length(checks) > 0 && "Integrated_MET_Index" %in% names(integrated)) {
    max(suppressWarnings(as.numeric(integrated$Integrated_MET_Index[integrated$Genotype %in% checks])), na.rm = TRUE)
  } else {
    NA_real_
  }
  if (!is.finite(check_index)) check_index <- NA_real_

  primary_result <- if (!is.na(primary_trait) && primary_trait %in% names(met_results)) met_results[[primary_trait]] else NULL
  primary_blups <- if (!is.null(primary_result)) {
    as.data.frame(met_ds_null(primary_result$blups_main, data.frame()))
  } else {
    data.frame()
  }
  primary_direction <- if (!is.null(primary_result)) as.character(met_ds_null(primary_result$trait_direction, "Higher better"))[1] else "Higher better"
  primary_target <- if (!is.null(primary_result)) suppressWarnings(as.numeric(met_ds_null(primary_result$target_value, NA_real_))[1]) else NA_real_
  primary_method <- if (!is.null(primary_result) && met_result_uses_anova(primary_result)) "BLUE" else "BLUP"
  primary_estimates <- if (nrow(primary_blups) > 0 && all(c("Genotype", "BLUP_G") %in% names(primary_blups))) {
    primary_numeric <- function(column) {
      if (column %in% names(primary_blups)) suppressWarnings(as.numeric(primary_blups[[column]])) else rep(NA_real_, nrow(primary_blups))
    }
    data.frame(
      Genotype = as.character(primary_blups$Genotype),
      Primary_trait = primary_trait,
      Estimate_method = primary_method,
      Raw_Mean = primary_numeric("Raw_Mean"),
      SE_Raw_Mean = primary_numeric("SE_Raw_Mean"),
      Genotype_estimate = primary_numeric("BLUP_G"),
      Estimate_SE = primary_numeric("SE_G"),
      Estimate_reliability = primary_numeric("Reliability"),
      Estimate_CI_Lower = primary_numeric("CI_lower"),
      Estimate_CI_Upper = primary_numeric("CI_upper"),
      Trait_direction = primary_direction,
      stringsAsFactors = FALSE
    ) %>%
      dplyr::mutate(
        .Decision_score = met_ds_direction_score(Genotype_estimate, primary_direction, primary_target)
      )
  } else {
    data.frame(
      Genotype = character(), Primary_trait = character(), Estimate_method = character(),
      Raw_Mean = numeric(), SE_Raw_Mean = numeric(), Genotype_estimate = numeric(),
      Estimate_SE = numeric(), Estimate_reliability = numeric(), Estimate_CI_Lower = numeric(),
      Estimate_CI_Upper = numeric(), Trait_direction = character(), Best_check = character(),
      Best_check_estimate = numeric(), Favorable_check_advantage = numeric(),
      Check_advantage_SE = numeric(), Check_advantage_CI_lower = numeric(),
      Check_advantage_CI_upper = numeric(), Probability_superior = numeric(), stringsAsFactors = FALSE
    )
  }
  best_check_row <- if (nrow(primary_estimates) > 0 && length(checks) > 0) {
    primary_estimates %>%
      dplyr::filter(Genotype %in% checks, is.finite(.Decision_score)) %>%
      dplyr::arrange(dplyr::desc(.Decision_score)) %>%
      head(1)
  } else {
    data.frame()
  }
  best_check <- if (nrow(best_check_row) > 0) as.character(best_check_row$Genotype[1]) else NA_character_
  best_check_estimate <- if (nrow(best_check_row) > 0) as.numeric(best_check_row$Genotype_estimate[1]) else NA_real_
  best_check_score <- if (nrow(best_check_row) > 0) as.numeric(best_check_row$.Decision_score[1]) else NA_real_
  best_check_se <- if (nrow(best_check_row) > 0) as.numeric(best_check_row$Estimate_SE[1]) else NA_real_
  if (nrow(primary_estimates) > 0) {
    primary_estimates <- primary_estimates %>%
      dplyr::mutate(
        Best_check = best_check,
        Best_check_estimate = best_check_estimate,
        Favorable_check_advantage = if (is.finite(best_check_score)) .Decision_score - best_check_score else NA_real_,
        Check_advantage_SE = ifelse(
          Genotype == best_check,
          0,
          sqrt(Estimate_SE^2 + best_check_se^2)
        ),
        Check_advantage_CI_lower = ifelse(
          is.na(Check_advantage_SE), NA_real_, Favorable_check_advantage - 1.96 * Check_advantage_SE
        ),
        Check_advantage_CI_upper = ifelse(
          is.na(Check_advantage_SE), NA_real_, Favorable_check_advantage + 1.96 * Check_advantage_SE
        ),
        Probability_superior = dplyr::case_when(
          Genotype == best_check ~ 0.5,
          is.na(Check_advantage_SE) ~ NA_real_,
          Check_advantage_SE == 0 ~ as.numeric(Favorable_check_advantage > 0),
          TRUE ~ stats::pnorm(Favorable_check_advantage / Check_advantage_SE)
        )
      ) %>%
      dplyr::select(-.Decision_score)
  }

  rank_confidence <- met_ds_rank_confidence(
    primary_estimates,
    direction = primary_direction,
    target_value = primary_target,
    top_k = min(5L, max(1L, nrow(primary_estimates)))
  )

  reliability <- matrix_info$long %>%
    dplyr::group_by(Genotype) %>%
    dplyr::summarise(
      Mean_reliability = if (any(!is.na(Reliability))) mean(Reliability, na.rm = TRUE) else NA_real_,
      Min_reliability = if (any(!is.na(Reliability))) min(Reliability, na.rm = TRUE) else NA_real_,
      N_traits_with_reliability = sum(!is.na(Reliability)),
      .groups = "drop"
    )
  reliability$Mean_reliability[!is.finite(reliability$Mean_reliability)] <- NA_real_
  reliability$Min_reliability[!is.finite(reliability$Min_reliability)] <- NA_real_

  coverage <- dplyr::bind_rows(lapply(traits, function(tr) {
    pres <- as.data.frame(met_ds_null(met_results[[tr]]$presence, data.frame()))
    if (!all(c("Genotype", "n_envs") %in% names(pres))) return(data.frame())
    env_cols <- setdiff(names(pres), c("Genotype", "n_envs"))
    data.frame(
      Genotype = as.character(pres$Genotype),
      Trait = tr,
      Observed_locations = suppressWarnings(as.numeric(pres$n_envs)),
      Total_locations = length(env_cols),
      Coverage_pct = if (length(env_cols) > 0) 100 * suppressWarnings(as.numeric(pres$n_envs)) / length(env_cols) else NA_real_,
      stringsAsFactors = FALSE
    )
  })) %>%
    dplyr::group_by(Genotype) %>%
    dplyr::summarise(
      Mean_coverage_pct = mean(Coverage_pct, na.rm = TRUE),
      Min_coverage_pct = min(Coverage_pct, na.rm = TRUE),
      Min_observed_locations = min(Observed_locations, na.rm = TRUE),
      .groups = "drop"
    )
  coverage$Mean_coverage_pct[!is.finite(coverage$Mean_coverage_pct)] <- NA_real_
  coverage$Min_coverage_pct[!is.finite(coverage$Min_coverage_pct)] <- NA_real_
  coverage$Min_observed_locations[!is.finite(coverage$Min_observed_locations)] <- NA_real_

  primary_presence <- if (!is.null(primary_result)) {
    as.data.frame(met_ds_null(primary_result$presence, data.frame()))
  } else {
    data.frame()
  }
  primary_coverage <- if (nrow(primary_presence) > 0 && all(c("Genotype", "n_envs") %in% names(primary_presence))) {
    primary_envs <- setdiff(names(primary_presence), c("Genotype", "n_envs"))
    observed_environment_text <- apply(primary_presence[, primary_envs, drop = FALSE], 1, function(row) {
      paste(primary_envs[!is.na(row) & suppressWarnings(as.numeric(row)) > 0], collapse = ", ")
    })
    data.frame(
      Genotype = as.character(primary_presence$Genotype),
      Environments_tested = paste0(primary_presence$n_envs, "/", length(primary_envs)),
      Coverage_pct = if (length(primary_envs) > 0) 100 * suppressWarnings(as.numeric(primary_presence$n_envs)) / length(primary_envs) else NA_real_,
      Tested_environment_names = observed_environment_text,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      Genotype = character(), Environments_tested = character(), Coverage_pct = numeric(),
      Tested_environment_names = character(), stringsAsFactors = FALSE
    )
  }

  gge_summary <- dplyr::bind_rows(lapply(traits, function(tr) {
    mean_tab <- as.data.frame(met_ds_null(met_results[[tr]]$gge_mean_stability_decision, data.frame()))
    win_tab <- as.data.frame(met_ds_null(met_results[[tr]]$gge_winner_decision, data.frame()))
    mean_part <- if (nrow(mean_tab) > 0 && all(c("Genotype", "Decision_class") %in% names(mean_tab))) {
      mean_tab %>%
        dplyr::transmute(
          Genotype = as.character(Genotype),
          Trait = tr,
          Broad_stable = grepl("broad-adaptation|stable", Decision_class, ignore.case = TRUE),
          Broad_high = grepl("Strong broad-adaptation|High performance", Decision_class, ignore.case = TRUE),
          Broad_class = as.character(Decision_class)
        )
    } else data.frame()
    winner_part <- if (nrow(win_tab) > 0 && "Winner" %in% names(win_tab)) {
      win_tab %>%
        dplyr::count(Genotype = as.character(Winner), name = "Winning_locations")
    } else data.frame(Genotype = character(), Winning_locations = integer())
    if (nrow(mean_part) == 0) return(data.frame())
    mean_part %>% dplyr::left_join(winner_part, by = "Genotype") %>%
      dplyr::mutate(Winning_locations = tidyr::replace_na(Winning_locations, 0L))
  })) %>%
    dplyr::group_by(Genotype) %>%
    dplyr::summarise(
      GGE_broad_stable_traits = sum(Broad_stable, na.rm = TRUE),
      GGE_high_traits = sum(Broad_high, na.rm = TRUE),
      GGE_winning_locations = sum(Winning_locations, na.rm = TRUE),
      GGE_class_summary = paste(unique(Broad_class), collapse = " | "),
      .groups = "drop"
    )

  primary_mean_stability <- if (!is.null(primary_result)) {
    as.data.frame(met_ds_null(primary_result$gge_mean_stability_decision, data.frame()))
  } else {
    data.frame()
  }
  primary_stability <- if (nrow(primary_mean_stability) > 0 && all(c("Genotype", "Decision_class") %in% names(primary_mean_stability))) {
    primary_mean_stability %>%
      dplyr::transmute(
        Genotype = as.character(Genotype),
        GGE_stability_evidence = as.character(Decision_class),
        Stability_class = dplyr::case_when(
          grepl("Strong broad-adaptation|Stable;", Decision_class, ignore.case = TRUE) ~ "STABLE",
          grepl("environment-sensitive|confirm stability", Decision_class, ignore.case = TRUE) ~ "RESPONSIVE",
          grepl("Lower broad-adaptation priority|stability unavailable", Decision_class, ignore.case = TRUE) ~ "UNSTABLE",
          TRUE ~ "RESPONSIVE"
        )
      )
  } else {
    data.frame(Genotype = character(), GGE_stability_evidence = character(), Stability_class = character(), stringsAsFactors = FALSE)
  }
  primary_winners <- if (!is.null(primary_result)) {
    as.data.frame(met_ds_null(primary_result$gge_winner_decision, data.frame()))
  } else {
    data.frame()
  }
  primary_targets <- if (nrow(primary_winners) > 0 && all(c("Environment", "Winner") %in% names(primary_winners))) {
    primary_winners %>%
      dplyr::transmute(Genotype = as.character(Winner), Environment = as.character(Environment)) %>%
      dplyr::filter(!is.na(Genotype), Genotype != "") %>%
      dplyr::group_by(Genotype) %>%
      dplyr::summarise(Target_environments = paste(sort(unique(Environment)), collapse = ", "), .groups = "drop")
  } else {
    data.frame(Genotype = character(), Target_environments = character(), stringsAsFactors = FALSE)
  }
  primary_adaptation <- dplyr::full_join(primary_stability, primary_targets, by = "Genotype") %>%
    dplyr::mutate(
      Target_environments = tidyr::replace_na(Target_environments, ""),
      Adaptation_type = dplyr::case_when(
        grepl("Strong broad-adaptation", GGE_stability_evidence, ignore.case = TRUE) ~ "BROAD",
        Target_environments != "" ~ "SPECIFIC",
        Stability_class == "STABLE" ~ "BROAD",
        Stability_class %in% c("RESPONSIVE", "UNSTABLE") ~ "SPECIFIC",
        TRUE ~ "UNCLASSIFIED"
      )
    )

  method_membership <- as.data.frame(met_ds_null(method_membership, data.frame()))
  method_small <- if (nrow(method_membership) > 0 && "Genotype" %in% names(method_membership)) {
    method_membership %>%
      dplyr::select(Genotype, Method_support_count, Agreement_pattern, Method_recommendation)
  } else {
    data.frame(Genotype = character(), Method_support_count = integer(), Agreement_pattern = character(), Method_recommendation = character())
  }

  board <- integrated %>%
    dplyr::mutate(Genotype = as.character(Genotype)) %>%
    dplyr::left_join(primary_estimates, by = "Genotype") %>%
    dplyr::left_join(rank_confidence, by = "Genotype") %>%
    dplyr::left_join(reliability, by = "Genotype") %>%
    dplyr::left_join(coverage, by = "Genotype") %>%
    dplyr::left_join(primary_coverage, by = "Genotype") %>%
    dplyr::left_join(gge_summary, by = "Genotype") %>%
    dplyr::left_join(primary_adaptation, by = "Genotype") %>%
    dplyr::left_join(method_small, by = "Genotype") %>%
    dplyr::mutate(
      Is_check = Genotype %in% checks,
      Best_check_integrated_index = check_index,
      Integrated_check_advantage = if ("Integrated_MET_Index" %in% names(.)) suppressWarnings(as.numeric(Integrated_MET_Index)) - check_index else NA_real_,
      Reliability_gate = dplyr::case_when(
        tidyr::replace_na(N_traits_with_reliability, 0) == 0 ~ "NOT_APPLICABLE_BLUE",
        !is.na(Min_reliability) & Min_reliability >= thresholds$min_reliability ~ "PASS",
        TRUE ~ "FAIL"
      ),
      Pass_reliability = Reliability_gate != "FAIL",
      Pass_coverage = !is.na(Min_coverage_pct) & Min_coverage_pct >= thresholds$min_coverage_pct,
      Pass_check = ifelse(
        length(checks) == 0,
        TRUE,
        !is.na(Favorable_check_advantage) &
          Favorable_check_advantage >= thresholds$min_check_advantage &
          (is.na(Probability_superior) | Probability_superior >= 0.50)
      ),
      Mandatory_eligibility = dplyr::case_when(
        Is_check ~ "CHECK",
        !Pass_reliability ~ "FAIL: reliability",
        !Pass_coverage ~ "FAIL: coverage",
        !Pass_check ~ "FAIL: below check",
        TRUE ~ "PASS"
      ),
      Final_action = dplyr::case_when(
        Is_check ~ NA_character_,
        Mandatory_eligibility == "PASS" & Integrated_rank <= max(1, ceiling(dplyr::n() * thresholds$selection_pct / 100)) &
          tidyr::replace_na(Method_support_count, 0) >= 2 &
          (is.na(Rank_confidence) | Rank_confidence >= 0.50) ~ "ADVANCE",
        !Pass_reliability | !Pass_coverage ~ "RETEST",
        tidyr::replace_na(Method_support_count, 0) >= 1 |
          tidyr::replace_na(GGE_winning_locations, 0) > 0 |
          tidyr::replace_na(GGE_high_traits, 0) > 0 |
          (!is.na(Probability_superior) & Probability_superior > 0.20) ~ "RETEST",
        TRUE ~ "DISCARD"
      ),
      Decision_reason = dplyr::case_when(
        Final_action == "ADVANCE" ~ paste0(
          "Passes precision, coverage, check, and ranking gates; adaptation = ",
          tidyr::replace_na(Adaptation_type, "unclassified"), "."
        ),
        Final_action == "RETEST" & grepl("^FAIL", Mandatory_eligibility) ~ paste(
          "Promising evidence remains, but more testing is required:", Mandatory_eligibility
        ),
        Final_action == "RETEST" ~ paste0(
          "Some favorable evidence is present, but confidence is insufficient for advancement; adaptation = ",
          tidyr::replace_na(Adaptation_type, "unclassified"), "."
        ),
        Final_action == "DISCARD" ~ "Adequate evidence is available, but current performance and selection support are insufficient.",
        TRUE ~ "Benchmark check; excluded from candidate decisions."
      )
    ) %>%
    dplyr::arrange(match(Final_action, c("ADVANCE", "RETEST", "DISCARD")), Integrated_rank)

  format_decision_number <- function(value, digits = 3L) {
    value <- suppressWarnings(as.numeric(value)[1])
    if (!is.finite(value)) "not available" else format(round(value, digits), trim = TRUE)
  }
  board$Decision_reason <- vapply(seq_len(nrow(board)), function(i) {
    if (isTRUE(board$Is_check[i])) return("Benchmark check; excluded from candidate decisions.")
    evidence <- paste0(
      "estimate ", format_decision_number(board$Genotype_estimate[i]),
      "; advantage over ", ifelse(is.na(board$Best_check[i]), "selected check", board$Best_check[i]),
      " = ", format_decision_number(board$Favorable_check_advantage[i]),
      "; probability superior = ", format_decision_number(board$Probability_superior[i], 2),
      "; reliability = ", format_decision_number(board$Estimate_reliability[i], 2),
      "; coverage = ", format_decision_number(board$Coverage_pct[i], 1), "%",
      "; rank confidence = ", format_decision_number(board$Rank_confidence[i], 2),
      "; stability = ", tolower(ifelse(is.na(board$Stability_class[i]), "unclassified", board$Stability_class[i])), "."
    )
    switch(
      board$Final_action[i],
      ADVANCE = paste("Advance:", evidence),
      RETEST = paste("Retest: evidence is promising or incomplete;", evidence),
      DISCARD = paste("Discard: adequate evidence does not support advancement;", evidence),
      evidence
    )
  }, character(1))

  trait_strengths <- matrix_info$long %>%
    dplyr::group_by(Genotype) %>%
    dplyr::summarise(
      Strength_traits = paste(head(Trait[order(-Standardized_score)], 3), collapse = ", "),
      Weakness_traits = paste(head(Trait[order(Standardized_score)], 3), collapse = ", "),
      .groups = "drop"
    )
  winner_envs <- dplyr::bind_rows(lapply(traits, function(tr) {
    win <- as.data.frame(met_ds_null(met_results[[tr]]$gge_winner_decision, data.frame()))
    if (!all(c("Environment", "Winner") %in% names(win))) return(data.frame())
    data.frame(Genotype = as.character(win$Winner), Win = paste(tr, win$Environment, sep = "@"), stringsAsFactors = FALSE)
  })) %>%
    dplyr::group_by(Genotype) %>%
    dplyr::summarise(Location_wins = paste(unique(Win), collapse = ", "), .groups = "drop")
  passport <- board %>%
    dplyr::filter(!Is_check) %>%
    dplyr::select(
      Genotype, Final_action, Mandatory_eligibility, Primary_trait, Estimate_method,
      Genotype_estimate, Estimate_SE, Estimate_reliability, Best_check,
      Favorable_check_advantage, Check_advantage_CI_lower, Check_advantage_CI_upper,
      Probability_superior, Rank_confidence, Environments_tested, Coverage_pct,
      Tested_environment_names, Integrated_rank, Integrated_MET_Index,
      Mean_reliability, Min_reliability, Reliability_gate, Mean_coverage_pct, Min_coverage_pct,
      Stability_class, Adaptation_type, Target_environments, GGE_broad_stable_traits,
      GGE_high_traits, GGE_winning_locations, Method_support_count, Agreement_pattern, Decision_reason
    ) %>%
    dplyr::left_join(trait_strengths, by = "Genotype") %>%
    dplyr::left_join(winner_envs, by = "Genotype") %>%
    dplyr::mutate(
      Location_wins = tidyr::replace_na(Location_wins, ""),
      Breeder_readout = paste(
        Final_action,
        "-", Primary_trait, Estimate_method, "=", round(Genotype_estimate, 3),
        "| favorable advantage over", tidyr::replace_na(Best_check, "no selected check"), "=", round(Favorable_check_advantage, 3),
        "| adaptation:", tidyr::replace_na(Adaptation_type, "unclassified"),
        "- strengths:", Strength_traits,
        "| weaknesses:", Weakness_traits,
        "|", Decision_reason
      )
    )
  board_display <- board %>%
    dplyr::filter(!Is_check) %>%
    dplyr::mutate(
      Check_advantage_CI = dplyr::if_else(
        is.na(Check_advantage_CI_lower) | is.na(Check_advantage_CI_upper),
        NA_character_,
        paste0("[", round(Check_advantage_CI_lower, 3), ", ", round(Check_advantage_CI_upper, 3), "]")
      )
    ) %>%
    dplyr::transmute(
      Genotype,
      Decision_Trait = Primary_trait,
      Estimate_Method = Estimate_method,
      Genotype_Estimate = Genotype_estimate,
      SE = Estimate_SE,
      Reliability = Estimate_reliability,
      Benchmark_Check = Best_check,
      Check_Advantage = Favorable_check_advantage,
      Check_Advantage_CI = Check_advantage_CI,
      Probability_Superior = Probability_superior,
      Environments_Tested = Environments_tested,
      Coverage_Pct = Coverage_pct,
      Stability_Class = Stability_class,
      Adaptation_Type = Adaptation_type,
      Target_Environments = Target_environments,
      Rank_Confidence = Rank_confidence,
      Final_Action = Final_action,
      Reason = Decision_reason
    )
  threshold_table <- data.frame(
    Gate = c(
      "Minimum reliability", "Minimum coverage pct", "Minimum check advantage",
      "Minimum probability for an uncertain candidate", "Minimum rank confidence for advance",
      "Selection pct for method support"
    ),
    Value = c(
      thresholds$min_reliability, thresholds$min_coverage_pct, thresholds$min_check_advantage,
      0.20, 0.50, thresholds$selection_pct
    ),
    Use = c(
      "Blocks advance when LMM genotype reliability is too low; not applied to ANOVA BLUEs.",
      "Blocks advance when genotype has weak location coverage.",
      "Requires candidate to equal or exceed the best selected check on the decision trait.",
      "Candidates above this probability remain in RETEST rather than DISCARD.",
      "ADVANCE requires at least a 50% simulated probability of ranking in the top five.",
      "Defines top-method sets for method agreement."
    ),
    stringsAsFactors = FALSE
  )
  list(
    board = met_ds_round_numeric(board_display),
    passport = met_ds_round_numeric(passport),
    thresholds = met_ds_round_numeric(threshold_table)
  )
}

# Direction-aware GGE views -----------------------------------------------
# Direction-aware GGE decision views for the MET pipeline.
#
# The functions in this file deliberately keep the statistical calculation
# separate from the Shiny wiring so the geometry can be tested independently.

met_gge_validate_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (length(dim(x)) != 2 || nrow(x) < 2 || ncol(x) < 2) {
    stop("GGE decision views require at least 2 genotypes and 2 environments.")
  }
  if (is.null(rownames(x)) || anyNA(rownames(x)) || any(rownames(x) == "")) {
    rownames(x) <- paste0("G", seq_len(nrow(x)))
  }
  if (is.null(colnames(x)) || anyNA(colnames(x)) || any(colnames(x) == "")) {
    colnames(x) <- paste0("E", seq_len(ncol(x)))
  }
  if (any(!is.finite(x))) {
    stop("The GGE decision matrix contains missing or non-finite values.")
  }
  x
}

met_gge_direction_matrix <- function(x, direction = "Higher better", target_value = NA_real_) {
  x <- met_gge_validate_matrix(x)
  if (is.null(direction) || length(direction) == 0 || is.na(direction[1]) || direction[1] == "") {
    direction <- "Higher better"
  }
  direction <- as.character(direction)[1]
  if (identical(direction, "Lower better")) {
    return(-x)
  }
  if (identical(direction, "Target trait")) {
    target_value <- suppressWarnings(as.numeric(target_value)[1])
    if (!is.finite(target_value)) {
      stop("A finite target value is required for a target-direction GGE analysis.")
    }
    return(-abs(x - target_value))
  }
  x
}

met_gge_direction_text <- function(direction, target_value = NA_real_) {
  if (identical(direction, "Lower better")) {
    return("Lower is favorable; values were sign-reversed for decision geometry")
  }
  if (identical(direction, "Target trait")) {
    return(paste0(
      "Closest to target ", format(target_value, trim = TRUE),
      " is favorable; geometry uses negative absolute target distance"
    ))
  }
  "Higher is favorable"
}

# Singular-value partitioning follows M = U D V':
# genotype-focused: U D and V; environment-focused: U and V D;
# symmetrical: U sqrt(D) and V sqrt(D).
met_gge_svd_scores <- function(x, svp = c("symmetrical", "genotype", "environment")) {
  svp <- match.arg(svp)
  x <- met_gge_validate_matrix(x)
  centered <- sweep(x, 2, colMeans(x), FUN = "-")
  fit <- svd(centered)
  d <- fit$d
  total_ss <- sum(d^2)
  variance_pct <- if (is.finite(total_ss) && total_ss > .Machine$double.eps) {
    100 * d^2 / total_ss
  } else {
    rep(0, length(d))
  }
  tolerance <- if (length(d) == 0 || max(d) == 0) {
    .Machine$double.eps
  } else {
    max(dim(centered)) * max(d) * .Machine$double.eps
  }
  effective_rank <- sum(d > tolerance)
  alpha <- switch(svp, genotype = 1, environment = 0, symmetrical = 0.5)
  k <- length(d)
  genotype_all <- sweep(fit$u[, seq_len(k), drop = FALSE], 2, d^alpha, FUN = "*")
  environment_all <- sweep(fit$v[, seq_len(k), drop = FALSE], 2, d^(1 - alpha), FUN = "*")
  rownames(genotype_all) <- rownames(x)
  rownames(environment_all) <- colnames(x)

  take_two <- function(scores) {
    result <- matrix(0, nrow = nrow(scores), ncol = 2)
    rownames(result) <- rownames(scores)
    colnames(result) <- c("PC1", "PC2")
    available <- min(2, ncol(scores))
    if (available > 0) {
      result[, seq_len(available)] <- scores[, seq_len(available), drop = FALSE]
    }
    result
  }

  list(
    centered = centered,
    genotype = take_two(genotype_all),
    environment = take_two(environment_all),
    genotype_all = genotype_all,
    environment_all = environment_all,
    singular_values = d,
    variance_pct = variance_pct,
    pc1_pct = if (length(variance_pct) >= 1) variance_pct[1] else 0,
    pc2_pct = if (length(variance_pct) >= 2) variance_pct[2] else 0,
    pc12_pct = sum(head(variance_pct, 2)),
    effective_rank = effective_rank,
    svp = svp
  )
}

met_gge_aec_axis <- function(genotype_scores, environment_scores, favorable_means) {
  genotype_scores <- as.matrix(genotype_scores)
  environment_scores <- as.matrix(environment_scores)
  if (ncol(genotype_scores) != ncol(environment_scores) || ncol(genotype_scores) < 1) {
    stop("Genotype and environment GGE scores must have the same non-zero number of components.")
  }
  axis_vector <- colMeans(environment_scores)
  axis_norm <- sqrt(sum(axis_vector^2))
  axis_defined <- is.finite(axis_norm) && axis_norm > sqrt(.Machine$double.eps)
  if (!axis_defined) {
    axis_vector <- rep(NA_real_, ncol(environment_scores))
    projected <- rep(NA_real_, nrow(genotype_scores))
    stability_distance <- rep(NA_real_, nrow(genotype_scores))
  } else {
    axis_vector <- axis_vector / axis_norm
    projected <- as.numeric(genotype_scores %*% axis_vector)
    orientation <- suppressWarnings(stats::cor(projected, favorable_means, use = "complete.obs"))
    if (is.finite(orientation) && orientation < 0) {
      axis_vector <- -axis_vector
      projected <- -projected
    }
    projected_coordinates <- outer(projected, axis_vector)
    stability_distance <- sqrt(rowSums((genotype_scores - projected_coordinates)^2))
  }
  perpendicular <- if (ncol(genotype_scores) == 2 && axis_defined) {
    c(-axis_vector[2], axis_vector[1])
  } else {
    rep(NA_real_, ncol(genotype_scores))
  }
  list(
    axis = axis_vector,
    perpendicular = perpendicular,
    defined = axis_defined,
    genotype_projection = projected,
    genotype_stability = stability_distance,
    signed_2d_stability = if (ncol(genotype_scores) == 2 && axis_defined) {
      as.numeric(genotype_scores %*% perpendicular)
    } else {
      rep(NA_real_, nrow(genotype_scores))
    }
  )
}

met_gge_circle_data <- function(radii, n = 240) {
  radii <- radii[is.finite(radii) & radii > 0]
  if (length(radii) == 0) {
    return(data.frame(x = numeric(), y = numeric(), radius = numeric()))
  }
  angles <- seq(0, 2 * pi, length.out = n)
  do.call(rbind, lapply(radii, function(radius) {
    data.frame(
      x = radius * cos(angles),
      y = radius * sin(angles),
      radius = radius
    )
  }))
}

met_gge_empty_plot <- function(message) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = message, size = 5, color = "gray35") +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::theme_void()
}

met_gge_plot_theme <- function() {
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(color = "gray35", size = 10),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(8, 14, 8, 8)
    )
}

met_gge_display_group <- function(genotypes, confidence_flag, check_genotypes = character(0), winners = character(0)) {
  ifelse(
    confidence_flag == "Low confidence", "Low confidence",
    ifelse(
      genotypes %in% check_genotypes, "Check",
      ifelse(genotypes %in% winners, "Location winner", "Candidate")
    )
  )
}

met_gge_confidence_frame <- function(genotypes, confidence = NULL, total_environments) {
  defaults <- data.frame(
    Genotype = genotypes,
    Observed_locations = total_environments,
    Total_locations = total_environments,
    Coverage_pct = 100,
    Confidence_flag = "OK",
    stringsAsFactors = FALSE
  )
  if (is.null(confidence) || nrow(as.data.frame(confidence)) == 0 ||
      !"Genotype" %in% names(confidence)) {
    return(defaults)
  }
  confidence <- as.data.frame(confidence, stringsAsFactors = FALSE)
  keep <- intersect(
    c("Genotype", "Observed_locations", "Total_locations", "Coverage_pct", "Confidence_flag"),
    names(confidence)
  )
  confidence <- confidence[, keep, drop = FALSE]
  idx <- match(genotypes, as.character(confidence$Genotype))
  for (column in setdiff(names(defaults), "Genotype")) {
    if (column %in% names(confidence)) {
      replacement <- confidence[[column]][idx]
      defaults[[column]][!is.na(replacement)] <- replacement[!is.na(replacement)]
    }
  }
  defaults
}

met_gge_label_subset <- function(table, id_col, required, limit = 25) {
  ids <- as.character(table[[id_col]])
  if (length(ids) <= limit) {
    return(ids)
  }
  unique(as.character(required[!is.na(required) & required != ""]))
}

build_met_gge_decision_views <- function(
    raw_matrix,
    direction = "Higher better",
    target_value = NA_real_,
    trait_name = "Trait",
    trait_weight = 1,
    confidence = NULL,
    cell_source = NULL,
    check_genotypes = character(0),
    estimate_matrix_source = "LMM predicted BLUP matrix") {
  raw_matrix <- met_gge_validate_matrix(raw_matrix)
  if (is.null(direction) || length(direction) == 0 || is.na(direction[1]) || direction[1] == "") {
    direction <- "Higher better"
  }
  direction <- as.character(direction)[1]
  decision_matrix <- met_gge_direction_matrix(raw_matrix, direction, target_value)
  genotype_names <- rownames(raw_matrix)
  environment_names <- colnames(raw_matrix)
  n_genotypes <- nrow(raw_matrix)
  n_environments <- ncol(raw_matrix)
  check_genotypes <- intersect(as.character(check_genotypes), genotype_names)
  trait_weight <- suppressWarnings(as.numeric(trait_weight)[1])
  if (!is.finite(trait_weight) || trait_weight < 0) trait_weight <- 1
  confidence_frame <- met_gge_confidence_frame(
    genotype_names,
    confidence = confidence,
    total_environments = n_environments
  )

  raw_means <- rowMeans(raw_matrix)
  favorable_means <- rowMeans(decision_matrix)
  target_distance <- if (identical(direction, "Target trait")) {
    rowMeans(abs(raw_matrix - target_value))
  } else {
    rep(NA_real_, n_genotypes)
  }
  direction_note <- met_gge_direction_text(direction, target_value)

  genotype_fit <- met_gge_svd_scores(decision_matrix, "genotype")
  symmetrical_fit <- met_gge_svd_scores(decision_matrix, "symmetrical")
  environment_fit <- met_gge_svd_scores(decision_matrix, "environment")
  pc1_pct <- symmetrical_fit$pc1_pct
  pc2_pct <- symmetrical_fit$pc2_pct
  pc12_pct <- symmetrical_fit$pc12_pct
  fidelity <- if (pc12_pct >= 70) {
    "Good two-PC summary"
  } else if (pc12_pct >= 50) {
    "Moderate two-PC summary; confirm with decision tables"
  } else {
    "Low two-PC summary; use decision tables as primary evidence"
  }
  pc_subtitle <- paste0(
    direction_note, " | PC1 + PC2 = ", round(pc12_pct, 1), "% (", fidelity, ")"
  )

  # 1. Mean-versus-stability (genotype-focused SVP = 1).
  mean_axis_full <- met_gge_aec_axis(
    genotype_fit$genotype_all,
    genotype_fit$environment_all,
    favorable_means
  )
  mean_axis <- met_gge_aec_axis(
    genotype_fit$genotype,
    genotype_fit$environment,
    favorable_means
  )
  full_performance_rank <- if (mean_axis_full$defined) {
    rank(-mean_axis_full$genotype_projection, ties.method = "min")
  } else {
    rep(NA_integer_, n_genotypes)
  }
  full_stability_rank <- if (mean_axis_full$defined) {
    rank(mean_axis_full$genotype_stability, ties.method = "min")
  } else {
    rep(NA_integer_, n_genotypes)
  }
  mean_genotype <- data.frame(
    Genotype = genotype_names,
    PC1 = genotype_fit$genotype[, 1],
    PC2 = genotype_fit$genotype[, 2],
    Mean_BLUP = raw_means,
    Mean_distance_to_target = target_distance,
    Direction_adjusted_mean = favorable_means,
    Direction_adjusted_mean_rank = rank(-favorable_means, ties.method = "min"),
    AEC_performance = mean_axis_full$genotype_projection,
    AEC_performance_rank = full_performance_rank,
    Stability_distance = mean_axis_full$genotype_stability,
    Stability_rank = full_stability_rank,
    PC1_PC2_AEC_performance = mean_axis$genotype_projection,
    PC1_PC2_stability_distance = mean_axis$genotype_stability,
    stringsAsFactors = FALSE
  )
  mean_genotype <- merge(mean_genotype, confidence_frame, by = "Genotype", all.x = TRUE, sort = FALSE)
  mean_genotype <- mean_genotype[match(genotype_names, mean_genotype$Genotype), , drop = FALSE]
  top_cutoff <- max(1, ceiling(n_genotypes * 0.25))
  stable_cutoff <- if (mean_axis_full$defined) {
    stats::median(mean_genotype$Stability_distance, na.rm = TRUE)
  } else {
    NA_real_
  }
  mean_genotype$Decision_class <- if (!mean_axis_full$defined) {
    ifelse(
      mean_genotype$Direction_adjusted_mean_rank <= top_cutoff,
      "High performance; confirm stability with FW/AMMI",
      "AEC stability unavailable; use FW/AMMI"
    )
  } else {
    ifelse(
      mean_genotype$Direction_adjusted_mean_rank <= top_cutoff & mean_genotype$Stability_distance <= stable_cutoff,
      "Strong broad-adaptation candidate",
      ifelse(
        mean_genotype$Direction_adjusted_mean_rank <= top_cutoff,
        "High performance; environment-sensitive",
        ifelse(
          mean_genotype$Stability_distance <= stable_cutoff,
          "Stable; performance not in top quartile",
          "Lower broad-adaptation priority"
        )
      )
    )
  }
  mean_genotype$Is_check <- mean_genotype$Genotype %in% check_genotypes
  mean_genotype$Trait_direction <- direction
  mean_genotype$Target_value <- if (identical(direction, "Target trait")) target_value else NA_real_
  mean_genotype$Integrated_trait_weight <- trait_weight
  mean_genotype$Display_group <- met_gge_display_group(
    mean_genotype$Genotype,
    mean_genotype$Confidence_flag,
    check_genotypes
  )
  mean_genotype$Projection_x <- mean_genotype$PC1_PC2_AEC_performance * mean_axis$axis[1]
  mean_genotype$Projection_y <- mean_genotype$PC1_PC2_AEC_performance * mean_axis$axis[2]
  mean_genotype <- mean_genotype[order(mean_genotype$Direction_adjusted_mean_rank, mean_genotype$Stability_rank), ]

  mean_environment <- data.frame(
    Environment = environment_names,
    PC1 = genotype_fit$environment[, 1],
    PC2 = genotype_fit$environment[, 2],
    stringsAsFactors = FALSE
  )
  mean_limit <- max(abs(c(
    mean_genotype$PC1, mean_genotype$PC2,
    mean_environment$PC1, mean_environment$PC2
  )), na.rm = TRUE)
  if (!is.finite(mean_limit) || mean_limit <= 0) mean_limit <- 1
  axis_segment <- data.frame(
    x = -1.1 * mean_limit * mean_axis$axis[1],
    y = -1.1 * mean_limit * mean_axis$axis[2],
    xend = 1.1 * mean_limit * mean_axis$axis[1],
    yend = 1.1 * mean_limit * mean_axis$axis[2]
  )
  mean_labels <- met_gge_label_subset(
    mean_genotype,
    "Genotype",
    c(
      head(mean_genotype$Genotype, 10),
      mean_genotype$Genotype[mean_genotype$Stability_rank <= 5],
      mean_genotype$Genotype[mean_genotype$Confidence_flag == "Low confidence"],
      check_genotypes
    )
  )
  p_mean_stability <- if (genotype_fit$effective_rank < 1 || !mean_axis$defined) {
    met_gge_empty_plot("The average-environment axis is undefined; use the numeric mean plus FW/AMMI stability evidence.")
  } else {
    ggplot2::ggplot() +
      ggplot2::geom_hline(yintercept = 0, color = "gray75", linewidth = 0.35) +
      ggplot2::geom_vline(xintercept = 0, color = "gray75", linewidth = 0.35) +
      ggplot2::geom_segment(
        data = axis_segment,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        color = "#34495E", linewidth = 0.8,
        arrow = grid::arrow(length = grid::unit(0.2, "cm"), type = "closed")
      ) +
      ggplot2::geom_segment(
        data = mean_genotype,
        ggplot2::aes(x = PC1, y = PC2, xend = Projection_x, yend = Projection_y),
        color = "#6C8EBF", linetype = "dotted", linewidth = 0.5, alpha = 0.8
      ) +
      ggplot2::geom_segment(
        data = mean_environment,
        ggplot2::aes(x = 0, y = 0, xend = PC1, yend = PC2),
        color = "#2E8B57", linewidth = 0.6,
        arrow = grid::arrow(length = grid::unit(0.13, "cm"), type = "closed")
      ) +
      ggplot2::geom_text(
        data = mean_environment,
        ggplot2::aes(x = PC1 * 1.08, y = PC2 * 1.08, label = Environment),
        color = "#2E8B57", fontface = "bold", size = 3.2
      ) +
      ggplot2::geom_point(
        data = mean_genotype,
        ggplot2::aes(x = PC1, y = PC2, color = Display_group, shape = Display_group),
        size = 2.8, stroke = 0.9
      ) +
      ggplot2::geom_text(
        data = mean_genotype[mean_genotype$Genotype %in% mean_labels, , drop = FALSE],
        ggplot2::aes(x = PC1, y = PC2, label = Genotype, color = Display_group),
        vjust = -0.85, size = 3, show.legend = FALSE
      ) +
      ggplot2::annotate(
        "text",
        x = axis_segment$xend,
        y = axis_segment$yend,
        label = "More favorable mean",
        color = "#34495E", fontface = "bold", size = 3, vjust = -0.7
      ) +
      ggplot2::scale_color_manual(values = c(
        "Candidate" = "#244F9E", "Check" = "#111111",
        "Location winner" = "#C0392B", "Low confidence" = "#D97706"
      )) +
      ggplot2::scale_shape_manual(values = c(
        "Candidate" = 16, "Check" = 18,
        "Location winner" = 17, "Low confidence" = 1
      )) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(
        title = paste0("GGE Mean vs Stability - ", trait_name),
        subtitle = pc_subtitle,
        x = paste0("PC1 (", round(genotype_fit$pc1_pct, 1), "%)"),
        y = paste0("PC2 (", round(genotype_fit$pc2_pct, 1), "%)"),
        color = NULL, shape = NULL
      ) +
      met_gge_plot_theme()
  }
  genotype_ideal_projection <- if (mean_axis$defined && any(is.finite(mean_genotype$PC1_PC2_AEC_performance))) {
    max(mean_genotype$PC1_PC2_AEC_performance, na.rm = TRUE)
  } else {
    NA_real_
  }
  genotype_ideal <- data.frame(
    x = genotype_ideal_projection * mean_axis$axis[1],
    y = genotype_ideal_projection * mean_axis$axis[2]
  )
  genotype_rank_dist <- if (is.finite(genotype_ideal$x) && is.finite(genotype_ideal$y)) {
    sqrt((mean_genotype$PC1 - genotype_ideal$x)^2 + (mean_genotype$PC2 - genotype_ideal$y)^2)
  } else {
    rep(NA_real_, nrow(mean_genotype))
  }
  mean_genotype$Genotype_ranking_distance <- genotype_rank_dist
  rank_circle_max <- max(genotype_rank_dist, na.rm = TRUE)
  if (!is.finite(rank_circle_max) || rank_circle_max <= 0) rank_circle_max <- mean_limit
  genotype_rank_circles <- met_gge_circle_data(seq(rank_circle_max / 5, rank_circle_max, length.out = 5))
  if (nrow(genotype_rank_circles) > 0 && is.finite(genotype_ideal$x) && is.finite(genotype_ideal$y)) {
    genotype_rank_circles$x <- genotype_rank_circles$x + genotype_ideal$x
    genotype_rank_circles$y <- genotype_rank_circles$y + genotype_ideal$y
  }
  genotype_rank_labels <- met_gge_label_subset(
    mean_genotype[order(mean_genotype$Genotype_ranking_distance), , drop = FALSE],
    "Genotype",
    c(head(mean_genotype$Genotype[order(mean_genotype$Genotype_ranking_distance)], 12), check_genotypes)
  )
  p_genotype_ranking <- if (genotype_fit$effective_rank < 1 || !mean_axis$defined || !is.finite(genotype_ideal$x)) {
    met_gge_empty_plot("GGE genotype ranking needs a defined average-environment axis.")
  } else {
    ggplot2::ggplot() +
      ggplot2::geom_path(
        data = genotype_rank_circles,
        ggplot2::aes(x = x, y = y, group = radius),
        color = "gray70", linewidth = 0.35, alpha = 0.65
      ) +
      ggplot2::geom_hline(yintercept = 0, color = "gray75", linewidth = 0.35) +
      ggplot2::geom_vline(xintercept = 0, color = "gray75", linewidth = 0.35) +
      ggplot2::geom_segment(
        data = axis_segment,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        color = "#34495E", linewidth = 0.8,
        arrow = grid::arrow(length = grid::unit(0.2, "cm"), type = "closed")
      ) +
      ggplot2::geom_point(
        data = genotype_ideal,
        ggplot2::aes(x = x, y = y),
        shape = 21, size = 4.5, stroke = 1, fill = NA, color = "#34495E"
      ) +
      ggplot2::annotate(
        "text", x = genotype_ideal$x, y = genotype_ideal$y,
        label = "Ideal", vjust = -1.2, color = "#34495E", fontface = "bold", size = 3
      ) +
      ggplot2::geom_point(
        data = mean_genotype,
        ggplot2::aes(x = PC1, y = PC2, color = Display_group, shape = Display_group),
        size = 2.8, stroke = 0.9
      ) +
      ggplot2::geom_text(
        data = mean_genotype[mean_genotype$Genotype %in% genotype_rank_labels, , drop = FALSE],
        ggplot2::aes(x = PC1, y = PC2, label = Genotype, color = Display_group),
        vjust = -0.85, size = 3, show.legend = FALSE
      ) +
      ggplot2::scale_color_manual(values = c(
        "Candidate" = "#244F9E", "Check" = "#111111",
        "Location winner" = "#C0392B", "Low confidence" = "#D97706"
      )) +
      ggplot2::scale_shape_manual(values = c(
        "Candidate" = 16, "Check" = 18,
        "Location winner" = 17, "Low confidence" = 1
      )) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(
        title = paste0("GGE Genotype Ranking - ", trait_name),
        subtitle = paste0(pc_subtitle, " | Smaller distance to the ideal point is better"),
        x = paste0("PC1 (", round(genotype_fit$pc1_pct, 1), "%)"),
        y = paste0("PC2 (", round(genotype_fit$pc2_pct, 1), "%)"),
        color = NULL, shape = NULL
      ) +
      met_gge_plot_theme()
  }

  # 2. Which-won-where (symmetrical SVP = 3).
  which_genotype <- data.frame(
    Genotype = genotype_names,
    PC1 = symmetrical_fit$genotype[, 1],
    PC2 = symmetrical_fit$genotype[, 2],
    stringsAsFactors = FALSE
  )
  which_environment <- data.frame(
    Environment = environment_names,
    PC1 = symmetrical_fit$environment[, 1],
    PC2 = symmetrical_fit$environment[, 2],
    stringsAsFactors = FALSE
  )
  winner_rows <- lapply(seq_len(n_environments), function(j) {
    values <- decision_matrix[, j]
    ord <- order(values, decreasing = TRUE, na.last = NA)
    winner_index <- ord[1]
    runner_index <- if (length(ord) >= 2) ord[2] else ord[1]
    favorable_margin <- values[winner_index] - values[runner_index]
    tie_tolerance <- max(1, max(abs(values), na.rm = TRUE)) * 1e-8
    check_index <- which(genotype_names %in% check_genotypes)
    best_check_index <- if (length(check_index) > 0) {
      check_index[which.max(values[check_index])]
    } else {
      NA_integer_
    }
    approximation <- as.numeric(symmetrical_fit$genotype %*% symmetrical_fit$environment[j, ])
    biplot_index <- which.max(approximation)
    data.frame(
      Environment = environment_names[j],
      Winner = genotype_names[winner_index],
      Winner_BLUP = raw_matrix[winner_index, j],
      Winner_favorable_score = values[winner_index],
      Runner_up = genotype_names[runner_index],
      Runner_up_BLUP = raw_matrix[runner_index, j],
      Favorable_margin_trait_units = favorable_margin,
      Near_tie = ifelse(favorable_margin <= tie_tolerance, "Yes - treat as co-winners", "No"),
      Winner_target_distance = if (identical(direction, "Target trait")) abs(raw_matrix[winner_index, j] - target_value) else NA_real_,
      Runner_up_target_distance = if (identical(direction, "Target trait")) abs(raw_matrix[runner_index, j] - target_value) else NA_real_,
      Best_check = if (is.finite(best_check_index)) genotype_names[best_check_index] else NA_character_,
      Best_check_BLUP = if (is.finite(best_check_index)) raw_matrix[best_check_index, j] else NA_real_,
      Favorable_advantage_over_best_check = if (is.finite(best_check_index)) values[winner_index] - values[best_check_index] else NA_real_,
      PC1_PC2_winner = genotype_names[biplot_index],
      Full_vs_biplot = ifelse(winner_index == biplot_index, "Agree", "Different - trust table"),
      stringsAsFactors = FALSE
    )
  })
  winner_table <- do.call(rbind, winner_rows)
  if (!is.null(cell_source) && nrow(as.data.frame(cell_source)) > 0 &&
      all(c("Genotype", "Environment", "Source") %in% names(cell_source))) {
    source_key <- paste(as.character(cell_source$Genotype), as.character(cell_source$Environment), sep = "\r")
    winner_key <- paste(winner_table$Winner, winner_table$Environment, sep = "\r")
    winner_table$Winner_cell_source <- as.character(cell_source$Source[match(winner_key, source_key)])
  } else {
    winner_table$Winner_cell_source <- NA_character_
  }
  winner_table$Environment_observed_cell_pct <- vapply(winner_table$Environment, function(environment) {
    if (is.null(cell_source) || nrow(as.data.frame(cell_source)) == 0 ||
        !all(c("Environment", "Source") %in% names(cell_source))) {
      return(NA_real_)
    }
    source_values <- as.character(cell_source$Source[as.character(cell_source$Environment) == environment])
    if (length(source_values) == 0) return(NA_real_)
    100 * mean(source_values == "Observed", na.rm = TRUE)
  }, numeric(1))
  winner_confidence_index <- match(winner_table$Winner, confidence_frame$Genotype)
  winner_table$Winner_coverage_pct <- confidence_frame$Coverage_pct[winner_confidence_index]
  winner_table$Winner_confidence <- confidence_frame$Confidence_flag[winner_confidence_index]
  winner_table$Winner_is_check <- winner_table$Winner %in% check_genotypes
  winner_table$Direction <- direction
  winner_table$Integrated_trait_weight <- trait_weight
  winner_table$PC1_PC2_explained_pct <- round(pc12_pct, 2)
  winner_table <- winner_table[order(winner_table$Environment), , drop = FALSE]
  location_winners <- unique(winner_table$Winner)

  hull_index <- integer(0)
  hull_path <- data.frame()
  sector_lines <- data.frame()
  if (symmetrical_fit$effective_rank >= 2 && n_genotypes >= 3) {
    hull_index <- unique(grDevices::chull(which_genotype$PC1, which_genotype$PC2))
    if (length(hull_index) >= 3) {
      closed_index <- c(hull_index, hull_index[1])
      hull_path <- which_genotype[closed_index, , drop = FALSE]
      plot_radius <- max(abs(c(
        which_genotype$PC1, which_genotype$PC2,
        which_environment$PC1, which_environment$PC2
      )), na.rm = TRUE)
      if (!is.finite(plot_radius) || plot_radius <= 0) plot_radius <- 1
      sector_parts <- lapply(seq_along(hull_index), function(k) {
        first <- which_genotype[hull_index[k], c("PC1", "PC2")]
        second <- which_genotype[hull_index[if (k == length(hull_index)) 1 else k + 1], c("PC1", "PC2")]
        edge <- c(second$PC1 - first$PC1, second$PC2 - first$PC2)
        normal <- c(-edge[2], edge[1])
        normal_norm <- sqrt(sum(normal^2))
        if (!is.finite(normal_norm) || normal_norm <= sqrt(.Machine$double.eps)) return(NULL)
        normal <- normal / normal_norm
        data.frame(
          x = -1.35 * plot_radius * normal[1],
          y = -1.35 * plot_radius * normal[2],
          xend = 1.35 * plot_radius * normal[1],
          yend = 1.35 * plot_radius * normal[2]
        )
      })
      sector_lines <- do.call(rbind, sector_parts)
      if (is.null(sector_lines)) sector_lines <- data.frame()
    }
  }
  which_genotype <- merge(which_genotype, confidence_frame, by = "Genotype", all.x = TRUE, sort = FALSE)
  which_genotype <- which_genotype[match(genotype_names, which_genotype$Genotype), , drop = FALSE]
  which_genotype$Display_group <- met_gge_display_group(
    which_genotype$Genotype,
    which_genotype$Confidence_flag,
    check_genotypes,
    location_winners
  )
  which_labels <- met_gge_label_subset(
    which_genotype,
    "Genotype",
    c(which_genotype$Genotype[hull_index], location_winners, check_genotypes,
      which_genotype$Genotype[which_genotype$Confidence_flag == "Low confidence"])
  )
  p_which_won <- if (symmetrical_fit$effective_rank < 2 || nrow(hull_path) == 0) {
    met_gge_empty_plot("Which-won-where needs two non-zero PCs and at least three non-collinear genotypes.")
  } else {
    ggplot2::ggplot() +
      ggplot2::geom_hline(yintercept = 0, color = "gray75", linewidth = 0.35) +
      ggplot2::geom_vline(xintercept = 0, color = "gray75", linewidth = 0.35) +
      ggplot2::geom_segment(
        data = sector_lines,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        color = "#7A7DB8", linetype = "dotted", linewidth = 0.6
      ) +
      ggplot2::geom_path(
        data = hull_path,
        ggplot2::aes(x = PC1, y = PC2),
        color = "#244F9E", linewidth = 0.9
      ) +
      ggplot2::geom_segment(
        data = which_environment,
        ggplot2::aes(x = 0, y = 0, xend = PC1, yend = PC2),
        color = "#2E8B57", linewidth = 0.55,
        arrow = grid::arrow(length = grid::unit(0.12, "cm"), type = "closed")
      ) +
      ggplot2::geom_text(
        data = which_environment,
        ggplot2::aes(x = PC1 * 1.08, y = PC2 * 1.08, label = Environment),
        color = "#2E8B57", fontface = "bold", size = 3.2
      ) +
      ggplot2::geom_point(
        data = which_genotype,
        ggplot2::aes(x = PC1, y = PC2, color = Display_group, shape = Display_group),
        size = 2.8, stroke = 0.9
      ) +
      ggplot2::geom_text(
        data = which_genotype[which_genotype$Genotype %in% which_labels, , drop = FALSE],
        ggplot2::aes(x = PC1, y = PC2, label = Genotype, color = Display_group),
        vjust = -0.85, size = 3, show.legend = FALSE
      ) +
      ggplot2::scale_color_manual(values = c(
        "Candidate" = "#244F9E", "Check" = "#111111",
        "Location winner" = "#C0392B", "Low confidence" = "#D97706"
      )) +
      ggplot2::scale_shape_manual(values = c(
        "Candidate" = 16, "Check" = 18,
        "Location winner" = 17, "Low confidence" = 1
      )) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(
        title = paste0("GGE Which-Won-Where - ", trait_name),
        subtitle = paste0(pc_subtitle, " | Red = full-matrix location winner"),
        x = paste0("PC1 (", round(pc1_pct, 1), "%)"),
        y = paste0("PC2 (", round(pc2_pct, 1), "%)"),
        color = NULL, shape = NULL
      ) +
      met_gge_plot_theme()
  }

  # 3. Discriminativeness-versus-representativeness. The paper-like styling is
  # retained, but environment-focused SVP = 2 is used so vector length has its
  # strongest environment-metric interpretation. Full-space values drive the
  # decision table; the two-PC values explain the chart.
  environment_axis_full <- met_gge_aec_axis(
    environment_fit$genotype_all,
    environment_fit$environment_all,
    favorable_means
  )
  environment_axis <- met_gge_aec_axis(
    environment_fit$genotype,
    environment_fit$environment,
    favorable_means
  )
  env_scores <- environment_fit$environment
  env_scores_full <- environment_fit$environment_all
  env_length_2d <- sqrt(rowSums(env_scores^2))
  env_length_full <- sqrt(rowSums(env_scores_full^2))
  env_projection_full <- if (environment_axis_full$defined) {
    as.numeric(env_scores_full %*% environment_axis_full$axis)
  } else {
    rep(NA_real_, n_environments)
  }
  env_cosine_full <- ifelse(
    environment_axis_full$defined & env_length_full > sqrt(.Machine$double.eps),
    env_projection_full / env_length_full,
    NA_real_
  )
  env_cosine_full <- pmax(-1, pmin(1, env_cosine_full))
  env_angle_full <- acos(env_cosine_full) * 180 / pi
  env_projection_2d <- if (environment_axis$defined) {
    as.numeric(env_scores %*% environment_axis$axis)
  } else {
    rep(NA_real_, n_environments)
  }
  env_cosine_2d <- ifelse(
    environment_axis$defined & env_length_2d > sqrt(.Machine$double.eps),
    env_projection_2d / env_length_2d,
    NA_real_
  )
  env_cosine_2d <- pmax(-1, pmin(1, env_cosine_2d))
  env_angle_2d <- acos(env_cosine_2d) * 180 / pi
  signal_captured <- ifelse(
    env_length_full > sqrt(.Machine$double.eps),
    100 * env_length_2d^2 / env_length_full^2,
    NA_real_
  )
  max_length_2d <- max(env_length_2d, na.rm = TRUE)
  if (!is.finite(max_length_2d) || max_length_2d <= 0) max_length_2d <- 1
  max_length_full <- max(env_length_full, na.rm = TRUE)
  if (!is.finite(max_length_full) || max_length_full <= 0) max_length_full <- 1
  ideal_point <- if (environment_axis$defined) {
    max_length_2d * environment_axis$axis
  } else {
    c(NA_real_, NA_real_)
  }
  ideal_distance_2d <- if (environment_axis$defined) {
    sqrt((env_scores[, 1] - ideal_point[1])^2 + (env_scores[, 2] - ideal_point[2])^2)
  } else {
    rep(NA_real_, n_environments)
  }
  cor_matrix <- suppressWarnings(stats::cor(decision_matrix, use = "pairwise.complete.obs"))
  diag(cor_matrix) <- NA_real_
  closest_environment <- rep(NA_character_, n_environments)
  closest_correlation <- rep(NA_real_, n_environments)
  for (j in seq_len(n_environments)) {
    candidates <- cor_matrix[j, ]
    if (any(is.finite(candidates))) {
      best <- which.max(replace(candidates, !is.finite(candidates), -Inf))
      closest_environment[j] <- environment_names[best]
      closest_correlation[j] <- candidates[best]
    }
  }
  discrimination_median <- stats::median(env_length_full, na.rm = TRUE)
  angle_median <- stats::median(env_angle_full, na.rm = TRUE)
  if (!is.finite(discrimination_median)) discrimination_median <- 0
  if (!is.finite(angle_median)) angle_median <- 90
  core_score_raw <- pmax(env_projection_full, 0)
  core_score <- if (any(is.finite(core_score_raw)) && max(core_score_raw, na.rm = TRUE) > 0) {
    100 * core_score_raw / max(core_score_raw, na.rm = TRUE)
  } else {
    rep(0, n_environments)
  }
  environment_role <- if (!environment_axis_full$defined) {
    rep("Average-environment axis undefined", n_environments)
  } else {
    ifelse(
      env_length_full >= discrimination_median & env_angle_full <= angle_median,
      "Core testing-environment candidate",
      ifelse(
        env_length_full >= discrimination_median,
        "Discriminating specialist environment",
        ifelse(
          env_angle_full <= angle_median,
          "Representative but lower discrimination",
          "Limited broad-network information"
        )
      )
    )
  }
  recommended_use <- ifelse(
    environment_role == "Core testing-environment candidate",
    "Prioritize for broad-adaptation screening",
    ifelse(
      environment_role == "Discriminating specialist environment",
      "Retain when its target population is important",
      ifelse(
        environment_role == "Representative but lower discrimination",
        "Use as confirmatory evidence",
        ifelse(
          environment_role == "Average-environment axis undefined",
          "Use discrimination and correlation only; AEC representativeness is unavailable",
          "Review cost and unique breeding value before retention"
        )
      )
    )
  )
  redundancy_flag <- is.finite(closest_correlation) & closest_correlation >= 0.90
  recommended_use[redundancy_flag] <- paste0(
    recommended_use[redundancy_flag],
    "; compare with ", closest_environment[redundancy_flag], " for redundancy"
  )
  observed_genotypes <- rep(NA_integer_, n_environments)
  environment_coverage_pct <- rep(NA_real_, n_environments)
  if (!is.null(cell_source) && nrow(as.data.frame(cell_source)) > 0 &&
      all(c("Genotype", "Environment", "Source") %in% names(cell_source))) {
    for (j in seq_len(n_environments)) {
      source_rows <- as.data.frame(cell_source)[
        as.character(cell_source$Environment) == environment_names[j], , drop = FALSE
      ]
      if (nrow(source_rows) > 0) {
        observed_genotypes[j] <- length(unique(as.character(
          source_rows$Genotype[as.character(source_rows$Source) == "Observed"]
        )))
        environment_coverage_pct[j] <- 100 * observed_genotypes[j] / n_genotypes
      }
    }
  }
  core_rank_order <- order(-core_score, ideal_distance_2d, environment_names, na.last = TRUE)
  core_rank <- rep(NA_integer_, n_environments)
  core_rank[core_rank_order] <- seq_along(core_rank_order)
  environment_table <- data.frame(
    Environment = environment_names,
    PC1 = env_scores[, 1],
    PC2 = env_scores[, 2],
    Full_discriminating_power = env_length_full,
    PC1_PC2_vector_length = env_length_2d,
    Environment_signal_captured_pct = signal_captured,
    Discrimination_pct_of_max = 100 * env_length_full / max_length_full,
    Discrimination_rank = rank(-env_length_full, ties.method = "min"),
    Full_representativeness_angle_deg = env_angle_full,
    PC1_PC2_representativeness_angle_deg = env_angle_2d,
    Full_AEC_cosine = env_cosine_full,
    Ideal_environment_distance_2D = ideal_distance_2d,
    Core_value_score_0_100 = core_score,
    Core_environment_rank = core_rank,
    Closest_environment = closest_environment,
    Max_response_correlation = closest_correlation,
    Potential_redundancy = ifelse(redundancy_flag, "Review pair", "No strong duplicate signal"),
    Observed_genotypes = observed_genotypes,
    Total_genotypes = n_genotypes,
    Observed_cell_coverage_pct = environment_coverage_pct,
    Imputed_cell_pct = 100 - environment_coverage_pct,
    Trait_direction = direction,
    Integrated_trait_weight = trait_weight,
    Environment_role = environment_role,
    Recommended_use = recommended_use,
    stringsAsFactors = FALSE
  )
  environment_table <- environment_table[order(environment_table$Core_environment_rank, environment_table$Discrimination_rank), ]
  circles <- met_gge_circle_data(seq(max_length_2d / 5, max_length_2d, length.out = 5))
  environment_axis_segment <- data.frame(
    x = -1.15 * max_length_2d * environment_axis$axis[1],
    y = -1.15 * max_length_2d * environment_axis$axis[2],
    xend = 1.15 * max_length_2d * environment_axis$axis[1],
    yend = 1.15 * max_length_2d * environment_axis$axis[2]
  )
  environment_genotype <- data.frame(
    Genotype = genotype_names,
    PC1 = environment_fit$genotype[, 1],
    PC2 = environment_fit$genotype[, 2],
    stringsAsFactors = FALSE
  )
  p_environment <- if (environment_fit$effective_rank < 1 || !environment_axis$defined) {
    met_gge_empty_plot("The average-environment axis is undefined; use discrimination and correlation in the table.")
  } else {
    ggplot2::ggplot() +
      ggplot2::geom_path(
        data = circles,
        ggplot2::aes(x = x, y = y, group = radius),
        color = "gray70", linewidth = 0.35, alpha = 0.65
      ) +
      ggplot2::geom_hline(yintercept = 0, color = "gray75", linewidth = 0.35) +
      ggplot2::geom_vline(xintercept = 0, color = "gray75", linewidth = 0.35) +
      ggplot2::geom_segment(
        data = environment_axis_segment,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        color = "#34495E", linewidth = 0.75,
        arrow = grid::arrow(length = grid::unit(0.18, "cm"), type = "closed")
      ) +
      ggplot2::geom_point(
        data = environment_genotype,
        ggplot2::aes(x = PC1, y = PC2),
        color = "#244F9E", alpha = 0.38, size = 1.7
      ) +
      ggplot2::geom_segment(
        data = environment_table,
        ggplot2::aes(x = 0, y = 0, xend = PC1, yend = PC2),
        color = "#2E8B57", linewidth = 0.75,
        arrow = grid::arrow(length = grid::unit(0.14, "cm"), type = "closed")
      ) +
      ggplot2::geom_point(
        data = environment_table,
        ggplot2::aes(x = PC1, y = PC2, size = Core_value_score_0_100),
        shape = 23, fill = "#74B95A", color = "#1D6B3A", stroke = 0.6
      ) +
      ggplot2::geom_text(
        data = environment_table,
        ggplot2::aes(x = PC1 * 1.08, y = PC2 * 1.08, label = Environment),
        color = "#1D6B3A", fontface = "bold", size = 3.2
      ) +
      ggplot2::annotate(
        "point", x = ideal_point[1], y = ideal_point[2],
        shape = 21, size = 4.5, stroke = 1, fill = NA, color = "#34495E"
      ) +
      ggplot2::annotate(
        "text", x = ideal_point[1], y = ideal_point[2],
        label = "Ideal", vjust = -1.2, color = "#34495E", fontface = "bold", size = 3
      ) +
      ggplot2::scale_size_continuous(range = c(2.3, 4.2), guide = "none") +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(
        title = paste0("GGE Environment Value - ", trait_name),
        subtitle = paste0(pc_subtitle, " | SVP = 2; long vector = discriminating; small AEC angle = representative"),
        x = paste0("PC1 (", round(pc1_pct, 1), "%)"),
        y = paste0("PC2 (", round(pc2_pct, 1), "%)")
      ) +
      met_gge_plot_theme()
  }
  environment_rank_labels <- met_gge_label_subset(
    environment_table[order(environment_table$Ideal_environment_distance_2D), , drop = FALSE],
    "Environment",
    head(environment_table$Environment[order(environment_table$Ideal_environment_distance_2D)], 12)
  )
  p_environment_ranking <- if (environment_fit$effective_rank < 1 || !environment_axis$defined) {
    met_gge_empty_plot("GGE environment ranking needs a defined average-environment axis.")
  } else {
    ggplot2::ggplot() +
      ggplot2::geom_path(
        data = circles,
        ggplot2::aes(x = x + ideal_point[1], y = y + ideal_point[2], group = radius),
        color = "gray70", linewidth = 0.35, alpha = 0.65
      ) +
      ggplot2::geom_hline(yintercept = 0, color = "gray75", linewidth = 0.35) +
      ggplot2::geom_vline(xintercept = 0, color = "gray75", linewidth = 0.35) +
      ggplot2::geom_segment(
        data = environment_axis_segment,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        color = "#34495E", linewidth = 0.75,
        arrow = grid::arrow(length = grid::unit(0.18, "cm"), type = "closed")
      ) +
      ggplot2::geom_segment(
        data = environment_table,
        ggplot2::aes(x = 0, y = 0, xend = PC1, yend = PC2),
        color = "#2E8B57", linewidth = 0.75,
        arrow = grid::arrow(length = grid::unit(0.14, "cm"), type = "closed")
      ) +
      ggplot2::geom_point(
        data = environment_table,
        ggplot2::aes(x = PC1, y = PC2, size = Core_value_score_0_100),
        shape = 23, fill = "#74B95A", color = "#1D6B3A", stroke = 0.6
      ) +
      ggplot2::geom_text(
        data = environment_table[environment_table$Environment %in% environment_rank_labels, , drop = FALSE],
        ggplot2::aes(x = PC1 * 1.08, y = PC2 * 1.08, label = Environment),
        color = "#1D6B3A", fontface = "bold", size = 3.2
      ) +
      ggplot2::annotate(
        "point", x = ideal_point[1], y = ideal_point[2],
        shape = 21, size = 4.5, stroke = 1, fill = NA, color = "#34495E"
      ) +
      ggplot2::annotate(
        "text", x = ideal_point[1], y = ideal_point[2],
        label = "Ideal", vjust = -1.2, color = "#34495E", fontface = "bold", size = 3
      ) +
      ggplot2::scale_size_continuous(range = c(2.3, 4.2), guide = "none") +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(
        title = paste0("GGE Environment Ranking - ", trait_name),
        subtitle = paste0(pc_subtitle, " | Ranked by distance to ideal test environment"),
        x = paste0("PC1 (", round(pc1_pct, 1), "%)"),
        y = paste0("PC2 (", round(pc2_pct, 1), "%)")
      ) +
      met_gge_plot_theme()
  }

  # Technical overview uses environment-focused SVP = 2 without arbitrary
  # post-hoc vector scaling. It remains secondary to the three decision views.
  overview_genotype <- data.frame(
    Genotype = genotype_names,
    PC1 = environment_fit$genotype[, 1],
    PC2 = environment_fit$genotype[, 2],
    stringsAsFactors = FALSE
  )
  overview_genotype <- merge(overview_genotype, confidence_frame, by = "Genotype", all.x = TRUE, sort = FALSE)
  overview_genotype <- overview_genotype[match(genotype_names, overview_genotype$Genotype), , drop = FALSE]
  overview_genotype$Display_group <- met_gge_display_group(
    overview_genotype$Genotype,
    overview_genotype$Confidence_flag,
    check_genotypes,
    location_winners
  )
  overview_environment <- data.frame(
    Environment = environment_names,
    PC1 = environment_fit$environment[, 1],
    PC2 = environment_fit$environment[, 2],
    stringsAsFactors = FALSE
  )
  overview_labels <- met_gge_label_subset(
    overview_genotype,
    "Genotype",
    c(location_winners, check_genotypes,
      overview_genotype$Genotype[overview_genotype$Confidence_flag == "Low confidence"])
  )
  p_overview <- if (environment_fit$effective_rank < 1) {
    met_gge_empty_plot("No G+GxE variation is available for a GGE overview.")
  } else {
    ggplot2::ggplot() +
      ggplot2::geom_hline(yintercept = 0, color = "gray75", linewidth = 0.35) +
      ggplot2::geom_vline(xintercept = 0, color = "gray75", linewidth = 0.35) +
      ggplot2::geom_segment(
        data = overview_environment,
        ggplot2::aes(x = 0, y = 0, xend = PC1, yend = PC2),
        color = "#2E8B57", linewidth = 0.7,
        arrow = grid::arrow(length = grid::unit(0.14, "cm"), type = "closed")
      ) +
      ggplot2::geom_text(
        data = overview_environment,
        ggplot2::aes(x = PC1 * 1.08, y = PC2 * 1.08, label = Environment),
        color = "#2E8B57", fontface = "bold", size = 3.2
      ) +
      ggplot2::geom_point(
        data = overview_genotype,
        ggplot2::aes(x = PC1, y = PC2, color = Display_group, shape = Display_group),
        size = 2.8, stroke = 0.9
      ) +
      ggplot2::geom_text(
        data = overview_genotype[overview_genotype$Genotype %in% overview_labels, , drop = FALSE],
        ggplot2::aes(x = PC1, y = PC2, label = Genotype, color = Display_group),
        vjust = -0.85, size = 3, show.legend = FALSE
      ) +
      ggplot2::scale_color_manual(values = c(
        "Candidate" = "#244F9E", "Check" = "#111111",
        "Location winner" = "#C0392B", "Low confidence" = "#D97706"
      )) +
      ggplot2::scale_shape_manual(values = c(
        "Candidate" = 16, "Check" = 18,
        "Location winner" = 17, "Low confidence" = 1
      )) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(
        title = paste0("GGE Technical Overview - ", trait_name),
        subtitle = paste0(pc_subtitle, " | Environment-focused SVP = 2"),
        x = paste0("PC1 (", round(pc1_pct, 1), "%)"),
        y = paste0("PC2 (", round(pc2_pct, 1), "%)"),
        color = NULL, shape = NULL
      ) +
      met_gge_plot_theme()
  }

  winner_agreement_n <- sum(winner_table$Full_vs_biplot == "Agree", na.rm = TRUE)
  predicted_winners_n <- sum(grepl("Imputed", winner_table$Winner_cell_source), na.rm = TRUE)
  notes <- data.frame(
    Decision_item = c(
      "Primary chart",
      "Direction treatment",
      "Integrated trait weight",
      "Data source",
      "PC1 explained",
      "PC2 explained",
      "PC1 + PC2 fidelity",
      "Full matrix vs two-PC winners",
      "Winning cells predicted by the model",
      "Reading order"
    ),
    Value = c(
      "GGE Mean vs Stability",
      direction_note,
      format(trait_weight, trim = TRUE),
      paste("Environment-centered", estimate_matrix_source),
      paste0(round(pc1_pct, 1), "%"),
      paste0(round(pc2_pct, 1), "%"),
      paste0(round(pc12_pct, 1), "% - ", fidelity),
      paste0(winner_agreement_n, " of ", n_environments, " environments agree"),
      paste0(predicted_winners_n, " of ", n_environments),
      "1) broad adaptation, 2) location winners, 3) environment value, then confirm uncertainty"
    ),
    Breeder_interpretation = c(
      "Start here for the clearest combined performance-and-stability decision.",
      "All chart directions use a common rule: larger decision score is more favorable.",
      if (trait_weight > 0) {
        "Weight affects the cross-trait integrated recommendation, not this single-trait geometry."
      } else {
        "This trait is diagnostic only and does not contribute to the cross-trait integrated recommendation."
      },
      "Unobserved cells are model predictions and remain flagged in the winner table.",
      "Share of G+GxE pattern shown on the horizontal axis.",
      "Additional share shown on the vertical axis.",
      "When fidelity is limited, the numeric decision tables outrank visual sectors.",
      "Differences occur because the chart displays only two components; trust the full-matrix winner column.",
      "Predicted winners need confirmation before advancement.",
      "This sequence separates broad choice, specific adaptation, and testing-network design."
    ),
    stringsAsFactors = FALSE
  )

  round_columns <- function(data, digits = 4) {
    numeric_columns <- vapply(data, is.numeric, logical(1))
    data[numeric_columns] <- lapply(data[numeric_columns], round, digits = digits)
    data
  }
  mean_output <- mean_genotype[, c(
    "Genotype", "Mean_BLUP", "Mean_distance_to_target", "Direction_adjusted_mean",
    "Direction_adjusted_mean_rank", "AEC_performance", "AEC_performance_rank",
    "Stability_distance", "Stability_rank", "PC1_PC2_AEC_performance",
    "PC1_PC2_stability_distance", "Observed_locations", "Total_locations",
    "Coverage_pct", "Confidence_flag", "Is_check", "Trait_direction", "Target_value",
    "Integrated_trait_weight", "Decision_class"
  ), drop = FALSE]

  list(
    direction = direction,
    target_value = target_value,
    decision_matrix = decision_matrix,
    notes = notes,
    mean_stability = round_columns(mean_output),
    winners = round_columns(winner_table),
    environments = round_columns(environment_table),
    overview_genotype = round_columns(overview_genotype[, setdiff(names(overview_genotype), "Display_group"), drop = FALSE]),
    overview_environment = round_columns(overview_environment),
    pc_variance = data.frame(
      Component = paste0("PC", seq_along(symmetrical_fit$variance_pct)),
      Explained_pct = round(symmetrical_fit$variance_pct, 4),
      Cumulative_pct = round(cumsum(symmetrical_fit$variance_pct), 4),
      stringsAsFactors = FALSE
    ),
    p_mean_stability = p_mean_stability,
    p_which_won = p_which_won,
    p_environment = p_environment,
    p_genotype_ranking = p_genotype_ranking,
    p_environment_ranking = p_environment_ranking,
    p_overview = p_overview,
    fits = list(
      genotype = genotype_fit,
      symmetrical = symmetrical_fit,
      environment = environment_fit
    )
  )
}

# MET pipeline -------------------------------------------------------------
MET_W_YIELD <- 0.50
MET_W_FW <- 0.25
MET_W_ASV <- 0.25
MET_MIN_ENVS_FOR_BIPLOT <- 3
met_rep_col_candidates <- c("Rep", "REP", "Replication", "Replicate", "Rep_No", "RepNo", "Replication_No")
met_block_col_candidates <- c("Block", "BLOCK", "Blk", "Block_No", "BlockNo", "Incomplete_Block")
prepare_met_input_frame <- function(df_raw) {
  df <- as.data.frame(df_raw)
  names(df) <- make.unique(trimws(names(df)), sep = "_")
  if (!"Environment" %in% names(df) && "Location" %in% names(df)) {
    names(df)[names(df) == "Location"] <- "Environment"
  }
  if (!"Genotype" %in% names(df) && id_col %in% names(df)) {
    names(df)[names(df) == id_col] <- "Genotype"
  }
  if (!"Genotype" %in% names(df)) {
    stop("MET requires a Genotype column. If your file uses Variety, it will be converted automatically.")
  }
  if (!"Environment" %in% names(df)) {
    stop("MET requires an Environment column. If your file uses Location, it will be converted automatically.")
  }
  df
}
pick_met_numeric_traits <- function(df, protected_cols) {
  candidate_cols <- setdiff(names(df), protected_cols)
  numeric_count <- sapply(candidate_cols, function(tr) sum(!is.na(to_number(df[[tr]]))))
  trait_cols <- candidate_cols[numeric_count > 0]
  if (length(trait_cols) == 0) {
    stop("MET requires at least one numeric trait column.")
  }
  trait_cols
}
get_met_trait_cols <- function(df_raw) {
  prepare_met_trait_settings(df_raw)$trait_cols
}
prepare_met_trait_settings <- function(df_raw) {
  df <- prepare_met_input_frame(df_raw)
  genotype_values_upper <- toupper(clean_text(df$Genotype))
  weight_rows <- which(genotype_values_upper %in% weight_row_labels)
  direction_rows <- which(genotype_values_upper %in% direction_row_labels)
  weight_row <- if (length(weight_rows) > 0) tail(weight_rows, 1) else integer(0)
  direction_row <- if (length(direction_rows) > 0) tail(direction_rows, 1) else integer(0)
  metadata_rows <- sort(unique(c(weight_rows, direction_rows)))
  df_data <- df
  if (length(metadata_rows) > 0) {
    df_data <- df_data[-metadata_rows, , drop = FALSE]
  }
  trait_cols <- pick_met_numeric_traits(
    df_data,
    c("Genotype", "Environment", met_rep_col_candidates, met_block_col_candidates)
  )
  weights_raw_used <- setNames(rep(1, length(trait_cols)), trait_cols)
  if (length(weight_row) > 0) {
    weight_values <- to_number(unlist(
      df[weight_row, trait_cols, drop = FALSE],
      use.names = FALSE
    ))
    names(weight_values) <- trait_cols
    weights_raw_used <- weight_values
    weights_raw_used[is.na(weights_raw_used)] <- 1
    weights_raw_used[weights_raw_used < 0] <- 0
    if (sum(weights_raw_used, na.rm = TRUE) <= 0) {
      weights_raw_used <- setNames(rep(1, length(trait_cols)), trait_cols)
    }
  }
  trait_direction <- setNames(rep("Higher better", length(trait_cols)), trait_cols)
  target_traits <- numeric(0)
  if (length(direction_row) > 0) {
    direction_values <- unlist(
      df[direction_row, trait_cols, drop = FALSE],
      use.names = FALSE
    )
    names(direction_values) <- trait_cols
    for (tr in trait_cols) {
      parsed <- parse_trait_direction(direction_values[[tr]])
      trait_direction[[tr]] <- parsed$direction
      if (parsed$direction == "Target trait") {
        target_traits[[tr]] <- parsed$target_value
      }
    }
  }
  trait_info <- data.frame(
    Trait = trait_cols,
    Direction = as.character(trait_direction[trait_cols]),
    Target_value = ifelse(
      trait_cols %in% names(target_traits),
      as.numeric(target_traits[trait_cols]),
      NA_real_
    ),
    Raw_weight = as.numeric(weights_raw_used[trait_cols]),
    stringsAsFactors = FALSE
  )
  list(
    data = df_data,
    trait_cols = trait_cols,
    weights_raw_used = weights_raw_used,
    trait_direction = trait_direction,
    target_traits = target_traits,
    trait_info = trait_info
  )
}
make_met_data <- function(df_raw, trait_used = NULL, replication_col = NULL, block_col = NULL) {
  prepared <- prepare_met_trait_settings(df_raw)
  df <- prepared$data
  trait_cols <- prepared$trait_cols
  if (is.null(trait_used)) {
    trait_used <- if ("Weight" %in% trait_cols) "Weight" else trait_cols[1]
  }
  if (!trait_used %in% trait_cols) {
    stop("Selected MET trait was not found or is not numeric: ", trait_used)
  }
  df$Weight <- to_number(df[[trait_used]])
  use_replication <- !is.null(replication_col) && replication_col != "" && replication_col %in% names(df)
  if (!use_replication) {
    stop("MET requires a replication column. Select a valid Rep/Replication column before running MET.")
  }
  use_block <- !is.null(block_col) && block_col != "" && block_col %in% names(df)
  dat <- df %>%
    mutate(
      Genotype = trimws(as.character(Genotype)),
      Environment = trimws(as.character(Environment)),
      Rep = if (use_replication) trimws(as.character(.data[[replication_col]])) else NA_character_,
      Block = if (use_block) trimws(as.character(.data[[block_col]])) else NA_character_
    ) %>%
    filter(
      !is.na(Genotype), Genotype != "",
      !is.na(Environment), Environment != "",
      !is.na(Weight)
    )
  dat <- dat %>% filter(!is.na(Rep), Rep != "")
  if (use_block) {
    dat <- dat %>% filter(!is.na(Block), Block != "")
  }
  dat <- dat %>%
    mutate(
      Genotype = factor(Genotype),
      Environment = factor(Environment),
      Rep = factor(Rep),
      Block = if (use_block) factor(Block) else factor("B1")
    )
  if (nrow(dat) == 0) {
    stop("No valid MET rows remained after cleaning Genotype, Environment, and Weight.")
  }
  if (n_distinct(dat$Genotype) < 2) {
    stop("MET requires at least 2 genotypes.")
  }
  if (n_distinct(dat$Environment) < 2) {
    stop("MET requires at least 2 environments.")
  }
  list(
    data = dat,
    trait_used = trait_used,
    trait_direction = as.character(prepared$trait_direction[[trait_used]] %||% "Higher better"),
    target_value = if (trait_used %in% names(prepared$target_traits)) {
      as.numeric(prepared$target_traits[[trait_used]])
    } else {
      NA_real_
    },
    trait_weight = as.numeric(prepared$weights_raw_used[[trait_used]] %||% 1),
    replication_col = replication_col,
    block_col = if (use_block) block_col else NULL
  )
}
met_safe_lrt <- function(model_a, model_b, test_name) {
  if (inherits(model_a, "met_model_error") || inherits(model_b, "met_model_error")) {
    note_parts <- character(0)
    if (inherits(model_a, "met_model_error")) note_parts <- c(note_parts, model_a$error)
    if (inherits(model_b, "met_model_error")) note_parts <- c(note_parts, model_b$error)
    note <- paste(note_parts, collapse = " | ")
    return(data.frame(Test = test_name, Model = "LRT model failed", npar = NA, AIC = NA, BIC = NA, logLik = NA, deviance = NA, Chisq = NA, Df = NA, `Pr(>Chisq)` = NA, Note = note, check.names = FALSE))
  }
  tryCatch({
    as.data.frame(anova(model_a, model_b)) %>% rownames_to_column("Model") %>% mutate(Test = test_name, .before = 1)
  }, error = function(e) {
    data.frame(Test = test_name, Model = "LRT failed", npar = NA, AIC = NA, BIC = NA, logLik = NA, deviance = NA, Chisq = NA, Df = NA, `Pr(>Chisq)` = NA, Note = e$message, check.names = FALSE)
  })
}
build_met_qc_table <- function(result) {
  if (is.null(result)) {
    return(data.frame())
  }
  rows <- list()
  add_qc <- function(area, metric, value, status = "OK", note = "") {
    rows[[length(rows) + 1]] <<- data.frame(
      QC_area = area,
      Metric = metric,
      Value = as.character(value),
      Status = status,
      Note = as.character(note),
      stringsAsFactors = FALSE
    )
  }
  presence <- as.data.frame(result$presence)
  genotype_summary <- as.data.frame(result$genotype_summary)
  outliers <- as.data.frame(result$outlier_summary)
  model_summary <- as.data.frame(result$model_summary)
  ammi_notes <- as.data.frame(result$ammi_notes)
  model_is_anova <- met_result_uses_anova(result)
  estimate_name <- if (model_is_anova) "adjusted-mean (BLUE)" else "BLUP"
  env_cols <- setdiff(names(presence), c("Genotype", "n_envs"))
  n_genotypes <- nrow(presence)
  n_envs <- length(env_cols)
  if (n_genotypes > 0 && n_envs > 0) {
    presence_matrix <- as.matrix(presence[, env_cols, drop = FALSE])
    observed_cells <- sum(suppressWarnings(as.numeric(presence_matrix)) == 1, na.rm = TRUE)
    expected_cells <- n_genotypes * n_envs
    missing_cells <- expected_cells - observed_cells
    missing_pct <- if (expected_cells > 0) round(100 * missing_cells / expected_cells, 1) else NA_real_
    coverage_status <- if (missing_cells == 0) "OK" else if (missing_pct <= 20) "Review" else "High risk"
    add_qc(
      "Coverage",
      "Observed genotype-location cells",
      paste0(observed_cells, " of ", expected_cells, " (", round(100 - missing_pct, 1), "%)"),
      coverage_status,
      paste0("Missing cells are model-estimated in the ", estimate_name, " matrix used by AMMI/GGE.")
    )
    min_envs_observed <- min(presence$n_envs, na.rm = TRUE)
    max_envs_observed <- max(presence$n_envs, na.rm = TRUE)
    add_qc(
      "Coverage",
      "Locations observed per genotype",
      paste0(min_envs_observed, " to ", max_envs_observed, " of ", n_envs),
      if (min_envs_observed < n_envs) "Review" else "OK",
      "Large differences mean some genotypes have weaker across-location evidence."
    )
    if (missing_cells > 0) {
      missing_examples <- presence %>%
        dplyr::select(Genotype, dplyr::all_of(env_cols)) %>%
        tidyr::pivot_longer(-Genotype, names_to = "Environment", values_to = "Present") %>%
        mutate(Present = suppressWarnings(as.numeric(Present))) %>%
        filter(is.na(Present) | Present == 0) %>%
        transmute(Cell = paste0(Genotype, "@", Environment)) %>%
        pull(Cell)
      add_qc(
        "Coverage",
        "Missing cells",
        missing_cells,
        coverage_status,
        if (length(missing_examples) == 0) "" else paste(head(missing_examples, 12), collapse = " | ")
      )
    }
  }
  min_required <- suppressWarnings(as.integer(ammi_notes$Min_locations_required[1] %||% NA_integer_))
  if (nrow(presence) > 0 && is.finite(min_required)) {
    low_conf <- presence %>%
      filter(n_envs < min_required) %>%
      pull(Genotype) %>%
      as.character()
    add_qc(
      "AMMI/GGE confidence",
      "Minimum observed locations required",
      min_required,
      "OK",
      "Genotypes below this threshold are shown as low confidence in AMMI/GGE."
    )
    add_qc(
      "AMMI/GGE confidence",
      "Low-confidence genotypes",
      length(low_conf),
      if (length(low_conf) == 0) "OK" else "Review",
      if (length(low_conf) == 0) "" else paste(low_conf, collapse = ", ")
    )
  }
  if (nrow(ammi_notes) > 0 && "Data_source" %in% names(ammi_notes)) {
    add_qc(
      "AMMI/GGE confidence",
      "Biplot data source",
      ammi_notes$Data_source[1],
      "Review",
      ammi_notes$Imputed_cell_warning[1] %||% ""
    )
  }
  if (nrow(outliers) > 0) {
    removed <- suppressWarnings(as.numeric(outliers$Removed_total[1] %||% NA_real_))
    original <- suppressWarnings(as.numeric(outliers$Rows_original[1] %||% NA_real_))
    removed_pct <- if (is.finite(original) && original > 0) round(100 * removed / original, 1) else NA_real_
    outlier_status <- if (!is.finite(removed) || removed == 0) "OK" else if (removed_pct <= 5) "Review" else "High risk"
    add_qc(
      "Outliers",
      "Rows removed by automatic filtering",
      paste0(removed, " of ", original, " (", removed_pct, "%)"),
      outlier_status,
      paste0("High removal changes may affect ", estimate_name, " and ranking decisions.")
    )
  }
  if (nrow(genotype_summary) > 0 && "N_observations" %in% names(genotype_summary)) {
    n_obs <- suppressWarnings(as.numeric(genotype_summary$N_observations))
    add_qc(
      "Replication balance",
      "Observations per genotype",
      paste0(min(n_obs, na.rm = TRUE), " to ", max(n_obs, na.rm = TRUE)),
      if (length(unique(stats::na.omit(n_obs))) > 1) "Review" else "OK",
      if (model_is_anova) {
        "Unequal observations can increase differences between raw and adjusted genotype means."
      } else {
        "Unequal observations can increase shrinkage differences among genotype BLUPs."
      }
    )
  }
  if (nrow(model_summary) > 0) {
    add_qc(
      "Controls",
      "Controls used",
      model_summary$Controls_used[1] %||% "",
      "Review",
      model_summary$Notes[1] %||% ""
    )
    add_qc(
      "Model",
      "Model formula",
      model_summary$Model_formula[1] %||% "",
      "OK",
      paste0("Broad-sense H2 = ", model_summary$Broad_sense_H2[1] %||% "")
    )
  }
  if (length(rows) == 0) {
    return(data.frame())
  }
  dplyr::bind_rows(rows)
}
met_fill_rank <- function(x) {
  if (length(x) == 0) return(x)
  if (all(is.na(x))) return(rep((length(x) + 1) / 2, length(x)))
  replace(x, is.na(x), max(x, na.rm = TRUE) + 1)
}
normalize_met_component_weights <- function(mean_weight, fw_weight, asv_weight) {
  scalar_weight <- function(value) {
    value <- suppressWarnings(as.numeric(value))
    if (length(value) == 0 || !is.finite(value[1]) || value[1] < 0) {
      return(0)
    }
    value[1]
  }
  weights <- c(
    mean = scalar_weight(mean_weight),
    fw = scalar_weight(fw_weight),
    asv = scalar_weight(asv_weight)
  )
  if (sum(weights) <= 0) {
    weights <- c(mean = 50, fw = 25, asv = 25)
  }
  weights / sum(weights)
}
build_met_selection_ranking <- function(result, component_weights = c(mean = 0.5, fw = 0.25, asv = 0.25)) {
  asv_table <- if (nrow(result$ammi_genotype) > 0) {
    result$ammi_genotype[, c("Genotype", "ASV"), drop = FALSE]
  } else {
    data.frame(Genotype = character(), ASV = numeric())
  }
  selection_base <- result$blups_main %>%
    left_join(result$fw_results[, c("Genotype", "Sens", "b_interp")], by = "Genotype") %>%
    left_join(asv_table, by = "Genotype")
  trait_direction <- as.character(result$trait_direction %||% "Higher better")
  target_value <- suppressWarnings(as.numeric(result$target_value %||% NA_real_))
  favorable_blup <- if (identical(trait_direction, "Lower better")) {
    -selection_base$BLUP_G
  } else if (identical(trait_direction, "Target trait") && is.finite(target_value)) {
    -abs(selection_base$BLUP_G - target_value)
  } else {
    selection_base$BLUP_G
  }
  rank_blup <- rank(-favorable_blup, ties.method = "average")
  rank_stability <- met_fill_rank(rank(abs(selection_base$Sens - 1), ties.method = "average", na.last = "keep"))
  rank_asv <- met_fill_rank(rank(selection_base$ASV, ties.method = "average", na.last = "keep"))
  selection_base %>%
    mutate(
      Rank_BLUP = rank_blup,
      Favorable_BLUP_score = favorable_blup,
      Trait_direction = trait_direction,
      Target_value = ifelse(identical(trait_direction, "Target trait"), target_value, NA_real_),
      Rank_stability = rank_stability,
      Rank_ASV = rank_asv,
      Mean_weight = round(component_weights[["mean"]] * 100, 1),
      FW_weight = round(component_weights[["fw"]] * 100, 1),
      ASV_weight = round(component_weights[["asv"]] * 100, 1),
      Combined_score = component_weights[["mean"]] * Rank_BLUP +
        component_weights[["fw"]] * Rank_stability +
        component_weights[["asv"]] * Rank_ASV,
      Final_rank = rank(Combined_score, ties.method = "first")
    ) %>%
    arrange(Final_rank) %>%
    mutate(CI_overlap_flag = CI_upper > lead(CI_lower))
}
plot_met_selection_ranking <- function(selection, trait_used, estimate_label = "BLUP") {
  ggplot(selection, aes(x = reorder(Genotype, -Final_rank), y = BLUP_G, fill = b_interp)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.3, color = "#2C3E50", linewidth = 0.6) +
    geom_text(aes(y = BLUP_G / 2, label = paste0("#", Final_rank)), vjust = 0.5, size = 3.5, fontface = "bold", color = "white") +
    scale_fill_manual(values = c("Responsive" = "#E74C3C", "Average" = "#F39C12", "Stable" = "#2ECC71"), na.value = "gray70") +
    labs(
      title = paste0("Hybrid selection - ", trait_used, " performance & stability"),
      subtitle = paste0(
        "Direction: ", unique(selection$Trait_direction)[1], " | ",
        "Weights: Mean=", round(unique(selection$Mean_weight)[1], 1),
        "%, FW=", round(unique(selection$FW_weight)[1], 1),
        "%, ASV=", round(unique(selection$ASV_weight)[1], 1), "%"
      ),
      x = "Genotype",
      y = paste0(estimate_label, " for ", trait_used),
      fill = "Stability (FW)"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8), plot.margin = margin(t = 10, r = 24, b = 24, l = 12))
}
build_met_integrated_ranking <- function(df_raw, met_results, component_weights = c(mean = 1, fw = 0, asv = 0)) {
  prepared <- prepare_met_trait_settings(df_raw)
  successful_traits <- intersect(prepared$trait_cols, names(met_results))
  if (length(successful_traits) == 0) {
    return(list(
      ranking = data.frame(),
      trait_weights = data.frame(),
      adjusted_performance = data.frame(),
      standardized_scores = data.frame(),
      plot = empty_plot("Integrated ranking needs at least one successful MET trait.")
    ))
  }
  weights_raw <- prepared$weights_raw_used[successful_traits]
  weights_raw[is.na(weights_raw)] <- 1
  weights_raw[weights_raw < 0] <- 0
  if (sum(weights_raw, na.rm = TRUE) <= 0) {
    weights_raw <- setNames(rep(1, length(successful_traits)), successful_traits)
  }
  weights <- weights_raw / sum(weights_raw)
  trait_blups <- purrr::map_dfr(successful_traits, function(tr) {
    fw_scores <- met_results[[tr]]$fw_results %>%
      transmute(
        Genotype = as.character(Genotype),
        FW_stability_raw = -abs(as.numeric(Sens) - 1)
      )
    asv_scores <- if (nrow(met_results[[tr]]$ammi_genotype) > 0) {
      met_results[[tr]]$ammi_genotype %>%
        transmute(
          Genotype = as.character(Genotype),
          ASV_stability_raw = -as.numeric(ASV)
        )
    } else {
      data.frame(Genotype = character(), ASV_stability_raw = numeric())
    }
    met_results[[tr]]$blups_main %>%
      transmute(
        Genotype = as.character(Genotype),
        Trait = tr,
        Raw_BLUP = as.numeric(BLUP_G)
      ) %>%
      left_join(fw_scores, by = "Genotype") %>%
      left_join(asv_scores, by = "Genotype")
  })
  adjust_value <- function(value, trait_name) {
    if (trait_name %in% names(prepared$target_traits)) {
      return(-abs(value - prepared$target_traits[[trait_name]]))
    }
    if (prepared$trait_direction[[trait_name]] == "Lower better") {
      return(-value)
    }
    value
  }
  trait_blups$Adjusted_performance <- mapply(
    adjust_value,
    trait_blups$Raw_BLUP,
    trait_blups$Trait
  )
  trait_blups <- trait_blups %>%
    group_by(Trait) %>%
    mutate(
      Mean_component = standardize_trait(Adjusted_performance),
      FW_component = standardize_trait(FW_stability_raw),
      ASV_component = standardize_trait(ASV_stability_raw),
      Component_weight_coverage =
        ifelse(!is.na(Mean_component), component_weights[["mean"]], 0) +
        ifelse(!is.na(FW_component), component_weights[["fw"]], 0) +
        ifelse(!is.na(ASV_component), component_weights[["asv"]], 0),
      Standardized_score = ifelse(
        Component_weight_coverage > 0,
        (
          replace_na(Mean_component, 0) * ifelse(!is.na(Mean_component), component_weights[["mean"]], 0) +
            replace_na(FW_component, 0) * ifelse(!is.na(FW_component), component_weights[["fw"]], 0) +
            replace_na(ASV_component, 0) * ifelse(!is.na(ASV_component), component_weights[["asv"]], 0)
        ) / Component_weight_coverage,
        NA_real_
      )
    ) %>%
    ungroup()
  adjusted_wide <- trait_blups %>%
    dplyr::select(Genotype, Trait, Adjusted_performance) %>%
    pivot_wider(
      names_from = Trait,
      values_from = Adjusted_performance
    )
  standardized_wide <- trait_blups %>%
    dplyr::select(Genotype, Trait, Standardized_score) %>%
    pivot_wider(
      names_from = Trait,
      values_from = Standardized_score
    )
  score_long <- trait_blups %>%
    dplyr::select(Genotype, Trait, Standardized_score) %>%
    mutate(Normalized_weight = as.numeric(weights[Trait]))
  integrated <- score_long %>%
    group_by(Genotype) %>%
    summarise(
      N_traits_used = sum(!is.na(Standardized_score)),
      Weight_coverage = sum(Normalized_weight[!is.na(Standardized_score)], na.rm = TRUE),
      Integrated_MET_Index = ifelse(
        Weight_coverage > 0,
        sum(Standardized_score * Normalized_weight, na.rm = TRUE) / Weight_coverage,
        NA_real_
      ),
      .groups = "drop"
    ) %>%
    arrange(desc(Integrated_MET_Index)) %>%
    mutate(Integrated_rank = row_number(), .before = 1)
  adjusted_output <- adjusted_wide %>%
    rename_with(~ paste0("Perf_", .x), all_of(successful_traits))
  standardized_output <- standardized_wide %>%
    rename_with(~ paste0("Std_", .x), all_of(successful_traits))
  ranking <- integrated %>%
    left_join(adjusted_output, by = "Genotype") %>%
    left_join(standardized_output, by = "Genotype") %>%
    mutate(
      Mean_weight = round(component_weights[["mean"]] * 100, 1),
      FW_weight = round(component_weights[["fw"]] * 100, 1),
      ASV_weight = round(component_weights[["asv"]] * 100, 1),
      Integrated_MET_Index = round(Integrated_MET_Index, 4),
      Weight_coverage = round(Weight_coverage, 4)
    )
  trait_weights <- prepared$trait_info %>%
    filter(Trait %in% successful_traits) %>%
    mutate(
      Normalized_weight = round(as.numeric(weights[Trait]), 4),
      Used_in_integrated_ranking = "YES"
    )
  p_integrated <- if (nrow(ranking) == 0) {
    empty_plot("Integrated ranking needs at least one genotype.")
  } else {
    plot_dat <- ranking %>%
      mutate(
        Rank_group = case_when(
          Integrated_rank <= ceiling(n() * 0.20) ~ "Top 20%",
          Integrated_MET_Index >= 0 ~ "Above average",
          TRUE ~ "Below average"
        )
      )
    ggplot(
      plot_dat,
      aes(
        x = reorder(Genotype, Integrated_MET_Index),
        y = Integrated_MET_Index,
        fill = Rank_group
      )
    ) +
      geom_col(width = 0.62) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.7) +
      geom_text(
        aes(
          label = paste0("#", Integrated_rank, " | ", sprintf("%.2f", Integrated_MET_Index)),
          hjust = ifelse(Integrated_MET_Index >= 0, -0.08, 1.08)
        ),
        size = 3.0
      ) +
      scale_fill_manual(
        values = c(
          "Top 20%" = "#1D9E75",
          "Above average" = "#3C78D8",
          "Below average" = "#999999"
        )
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.12, 0.25))) +
      coord_flip() +
      labs(
        title = "Integrated MET ranking",
        subtitle = paste0(
          length(successful_traits),
          " parameter(s), Mean/FW/ASV weights = ",
          round(component_weights[["mean"]] * 100, 1), "/",
          round(component_weights[["fw"]] * 100, 1), "/",
          round(component_weights[["asv"]] * 100, 1)
        ),
        x = "Genotype",
        y = "Weighted standardized MET performance index",
        fill = NULL
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 15, margin = margin(b = 2)),
        plot.subtitle = element_text(color = "gray40", margin = margin(b = 6)),
        axis.title.x = element_text(size = 12, margin = margin(t = 6)),
        axis.title.y = element_text(size = 12, margin = margin(r = 6)),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        legend.position = "bottom"
      )
  }
  list(
    ranking = ranking,
    trait_weights = trait_weights,
    adjusted_performance = adjusted_output,
    standardized_scores = standardized_output,
    plot = p_integrated
  )
}
met_gge_chart_choices <- function() {
  c(
    "Mean vs Stability" = "mean_stability",
    "Which-Won-Where" = "which_won",
    "Discriminativeness vs Representativeness" = "environment",
    "Genotype Ranking" = "genotype_ranking",
    "Environment Ranking" = "environment_ranking"
  )
}
met_fw_chart_choices <- function() {
  c(
    "Mean vs Sensitivity" = "mean_sensitivity",
    "Response Correlation" = "response_correlation"
  )
}
met_fw_selected_view <- function(input_value) {
  choices <- met_fw_chart_choices()
  if (!is.null(input_value) && input_value %in% unname(choices)) input_value else unname(choices)[1]
}
met_plot_fw_selected <- function(result, fw_view) {
  switch(
    met_fw_selected_view(fw_view),
    mean_sensitivity = result$p_fw_mean_sens,
    response_correlation = result$p_fw_regression,
    result$p_fw_mean_sens
  )
}
met_fw_chart_file_slug <- function(fw_view) {
  switch(
    met_fw_selected_view(fw_view),
    mean_sensitivity = "MET_FW_mean_sensitivity",
    response_correlation = "MET_FW_response_correlation",
    "MET_FW"
  )
}
met_performance_heatmap_choices <- function(result = NULL) {
  estimate_label <- if (!is.null(result) && met_result_uses_anova(result)) {
    "Adjusted Means (BLUEs)"
  } else if (!is.null(result)) {
    "Model Predictions (BLUPs)"
  } else {
    "Estimated Values (BLUPs/BLUEs)"
  }
  stats::setNames(c("observed", "estimated"), c("Observed Means", estimate_label))
}
met_performance_heatmap_selected_view <- function(input_value) {
  choices <- met_performance_heatmap_choices()
  if (!is.null(input_value) && input_value %in% unname(choices)) input_value else unname(choices)[1]
}
met_plot_performance_heatmap <- function(result, heatmap_view) {
  switch(
    met_performance_heatmap_selected_view(heatmap_view),
    observed = result$p_perf_heatmap_observed %||% result$p_perf_heatmap,
    estimated = result$p_perf_heatmap,
    result$p_perf_heatmap_observed %||% result$p_perf_heatmap
  )
}
met_performance_heatmap_file_slug <- function(heatmap_view) {
  switch(
    met_performance_heatmap_selected_view(heatmap_view),
    observed = "MET_genotype_by_location_observed",
    estimated = "MET_genotype_by_location_estimated",
    "MET_genotype_by_location"
  )
}
met_gge_selected_view <- function(input_value) {
  choices <- met_gge_chart_choices()
  if (!is.null(input_value) && input_value %in% unname(choices)) input_value else unname(choices)[1]
}
met_plot_gge_selected <- function(result, gge_view) {
  gge_view <- met_gge_selected_view(gge_view)
  switch(
    gge_view,
    mean_stability = result$p_gge_mean_stability,
    which_won = result$p_gge_which_won,
    environment = result$p_gge_environment,
    genotype_ranking = result$p_gge_genotype_ranking %||% result$p_gge_mean_stability,
    environment_ranking = result$p_gge_environment_ranking %||% result$p_gge_environment,
    result$p_gge_mean_stability
  )
}
met_gge_chart_file_slug <- function(gge_view) {
  switch(
    met_gge_selected_view(gge_view),
    mean_stability = "MET_GGE_mean_stability",
    which_won = "MET_GGE_which_won_where",
    environment = "MET_GGE_discriminativeness_representativeness",
    genotype_ranking = "MET_GGE_genotype_ranking",
    environment_ranking = "MET_GGE_environment_ranking",
    "MET_GGE"
  )
}
plot_met_mgidi_ranking <- function(methods) {
  mgidi <- as.data.frame(methods$mgidi %||% data.frame())
  if (nrow(mgidi) == 0 || !all(c("Genotype", "MGIDI_distance", "MGIDI_rank") %in% names(mgidi))) {
    return(empty_plot("MGIDI needs at least two successful MET traits."))
  }
  membership <- as.data.frame(methods$method_membership %||% data.frame())
  if (all(c("Genotype", "MGIDI") %in% names(membership))) {
    mgidi <- mgidi %>%
      left_join(membership %>% dplyr::select(Genotype, MGIDI), by = "Genotype") %>%
      mutate(Selected = replace_na(as.logical(MGIDI), FALSE))
  } else {
    selected_n <- max(1, ceiling(nrow(mgidi) * 0.20))
    mgidi <- mgidi %>% mutate(Selected = MGIDI_rank <= selected_n)
  }
  mgidi <- mgidi %>%
    filter(is.finite(MGIDI_distance), !is.na(Genotype)) %>%
    arrange(MGIDI_rank, Genotype)
  if (nrow(mgidi) == 0) {
    return(empty_plot("MGIDI distances are not available."))
  }
  distance_ceiling <- max(mgidi$MGIDI_distance, na.rm = TRUE) * 1.08
  if (!is.finite(distance_ceiling) || distance_ceiling <= 0) distance_ceiling <- 1
  mgidi <- mgidi %>%
    mutate(
      Angle_index = row_number(),
      Radial_score = distance_ceiling - MGIDI_distance,
      Selection = ifelse(Selected, "Selected", "Nonselected")
    )
  cutoff_distance <- if (any(mgidi$Selected)) {
    max(mgidi$MGIDI_distance[mgidi$Selected], na.rm = TRUE)
  } else {
    stats::quantile(mgidi$MGIDI_distance, 0.20, na.rm = TRUE, names = FALSE)
  }
  cutoff_radius <- distance_ceiling - cutoff_distance
  radial_breaks <- pretty(c(0, max(mgidi$Radial_score, na.rm = TRUE)), n = 5)
  radial_breaks <- radial_breaks[radial_breaks >= 0 & radial_breaks <= distance_ceiling]
  ggplot(mgidi, aes(x = Angle_index, y = Radial_score)) +
    geom_hline(yintercept = cutoff_radius, color = "#E31A1C", linewidth = 0.9) +
    geom_path(aes(group = 1), color = "#242424", linewidth = 0.8) +
    geom_point(aes(fill = Selection), shape = 21, size = 3.4, color = "#555555", stroke = 0.7) +
    scale_fill_manual(values = c("Nonselected" = "#BDBDBD", "Selected" = "#E31A1C")) +
    scale_x_continuous(
      breaks = mgidi$Angle_index,
      labels = mgidi$Genotype,
      limits = c(0.5, nrow(mgidi) + 0.5),
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      breaks = radial_breaks,
      labels = function(value) round(distance_ceiling - value, 2),
      limits = c(0, distance_ceiling),
      expand = expansion(mult = c(0, 0.04))
    ) +
    coord_polar(start = -pi / 2) +
    labs(
      title = "MGIDI multi-trait selection",
      subtitle = "Lower ideotype distance is plotted farther from the center; the red circle is the selection cutoff.",
      x = NULL,
      y = "Multi-trait genotype-ideotype distance index",
      fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_line(color = "#E5E5E5", linewidth = 0.5),
      panel.grid.major.y = element_line(color = "#E5E5E5", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(color = "#303030", size = 9),
      axis.text.y = element_text(color = "#404040", size = 8),
      axis.title.y = element_text(margin = margin(r = 10)),
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )
}
run_met_pipeline <- function(
    df_raw,
    trait_used = NULL,
    check_varieties = NULL,
    replication_col = NULL,
    block_col = NULL,
    min_envs_for_biplot = NULL,
    model_type = "LMM") {
  input <- make_met_data(
    df_raw,
    trait_used,
    replication_col = replication_col,
    block_col = block_col
  )
  dat <- input$data
  trait_used <- input$trait_used
  trait_direction <- input$trait_direction %||% "Higher better"
  target_value <- suppressWarnings(as.numeric(input$target_value %||% NA_real_))
  trait_weight <- suppressWarnings(as.numeric(input$trait_weight %||% 1))
  has_replication <- !is.null(input$replication_col)
  has_block <- !is.null(input$block_col)
  model_type <- toupper(trimws(as.character(model_type %||% "LMM")))[1]
  if (model_type %in% c("ANOVA", "RCBD", "ANOVA (RCBD)", "ANOVA_RCBD")) {
    model_type <- "ANOVA_RCBD"
  }
  if (!model_type %in% c("LMM", "ANOVA_RCBD")) {
    stop("MET model must be LMM or ANOVA (RCBD).")
  }
  model_label <- if (identical(model_type, "LMM")) "LMM" else "ANOVA (RCBD)"
  estimate_label <- if (identical(model_type, "LMM")) "BLUP" else "adjusted mean (BLUE)"
  estimate_label_plural <- if (identical(model_type, "LMM")) "BLUPs" else "adjusted means (BLUEs)"
  estimate_matrix_source <- if (identical(model_type, "LMM")) {
    "LMM predicted BLUP matrix"
  } else {
    "ANOVA adjusted-mean (BLUE) matrix"
  }
  n_envs_total <- n_distinct(dat$Environment)
  min_envs_for_biplot <- suppressWarnings(as.integer(min_envs_for_biplot))
  if (length(min_envs_for_biplot) == 0 || is.na(min_envs_for_biplot) || min_envs_for_biplot < MET_MIN_ENVS_FOR_BIPLOT) {
    min_envs_for_biplot <- max(MET_MIN_ENVS_FOR_BIPLOT, ceiling(0.5 * n_envs_total))
  }
  min_envs_for_biplot <- max(MET_MIN_ENVS_FOR_BIPLOT, min_envs_for_biplot)
  if (n_envs_total >= MET_MIN_ENVS_FOR_BIPLOT) {
    min_envs_for_biplot <- min(min_envs_for_biplot, n_envs_total)
  }
  ammi_gge_available <- n_envs_total >= MET_MIN_ENVS_FOR_BIPLOT
  notes <- c()
  cell_bounds <- dat %>% group_by(Genotype, Environment) %>% summarise(Lower_c = quantile(Weight, 0.25, na.rm = TRUE) - 1.5 * IQR(Weight, na.rm = TRUE), Upper_c = quantile(Weight, 0.75, na.rm = TRUE) + 1.5 * IQR(Weight, na.rm = TRUE), .groups = "drop")
  dat_stepA <- dat %>% left_join(cell_bounds, by = c("Genotype", "Environment")) %>% filter(!(Weight < Lower_c | Weight > Upper_c)) %>% dplyr::select(-Lower_c, -Upper_c)
  if (nrow(dat_stepA) == 0) stop("All MET rows were removed by cell-level outlier filtering.")
  env_bounds <- dat_stepA %>% group_by(Environment) %>% summarise(Lower = quantile(Weight, 0.25, na.rm = TRUE) - 1.5 * IQR(Weight, na.rm = TRUE), Upper = quantile(Weight, 0.75, na.rm = TRUE) + 1.5 * IQR(Weight, na.rm = TRUE), .groups = "drop")
  p_before <- ggplot(dat_stepA %>% left_join(env_bounds, by = "Environment") %>% mutate(outlier_env = Weight < Lower | Weight > Upper), aes(x = Weight, fill = outlier_env)) + geom_histogram(bins = 30, color = "white", linewidth = 0.2) + geom_vline(data = env_bounds, aes(xintercept = Lower), linetype = "dashed", color = "#E74C3C", linewidth = 0.7) + geom_vline(data = env_bounds, aes(xintercept = Upper), linetype = "dashed", color = "#E74C3C", linewidth = 0.7) + facet_wrap(~Environment, scales = "free", ncol = 2) + scale_fill_manual(values = c("FALSE" = "#3498DB", "TRUE" = "#E74C3C"), labels = c("Kept", "Outlier")) + labs(title = paste0(trait_used, " distribution after cell-clean, before env-clean"), x = trait_used, y = "Count", fill = NULL) + theme_bw()
  dat_clean <- dat_stepA %>% left_join(env_bounds, by = "Environment") %>% filter(!(Weight < Lower | Weight > Upper)) %>% dplyr::select(-Lower, -Upper)
  if (nrow(dat_clean) == 0) stop("All MET rows were removed by environment-level outlier filtering.")
  p_after <- ggplot(dat_clean, aes(x = Weight)) + geom_histogram(bins = 30, fill = "#2ECC71", color = "white", linewidth = 0.2) + facet_wrap(~Environment, scales = "free", ncol = 2) + labs(title = paste0(trait_used, " distribution after both outlier steps"), x = trait_used, y = "Count") + theme_bw()
  outlier_summary <- data.frame(Trait_used = trait_used, Rows_original = nrow(dat), Rows_after_cell_clean = nrow(dat_stepA), Rows_after_environment_clean = nrow(dat_clean), Removed_total = nrow(dat) - nrow(dat_clean), stringsAsFactors = FALSE)
  presence <- dat_clean %>% distinct(Genotype, Environment) %>% mutate(Genotype = as.character(Genotype), Environment = as.character(Environment), present = 1) %>% pivot_wider(names_from = Environment, values_from = present, values_fill = 0) %>% mutate(n_envs = rowSums(across(-Genotype))) %>% arrange(desc(n_envs), Genotype)
  selected_controls <- unique(clean_text(check_varieties))
  selected_controls <- selected_controls[!is.na(selected_controls) & selected_controls != ""]
  controls_present <- intersect(selected_controls, as.character(unique(dat_clean$Genotype)))
  if (length(selected_controls) > 0 && length(controls_present) > 0) {
    controls_used <- controls_present
    missing_selected_controls <- setdiff(selected_controls, controls_present)
    if (length(missing_selected_controls) > 0) {
      notes <- c(
        notes,
        paste0(
          "Selected controls not found in cleaned MET data: ",
          paste(missing_selected_controls, collapse = ", "),
          "."
        )
      )
    }
  } else {
    controls_used <- presence %>% filter(n_envs == n_envs_total) %>% pull(Genotype)
    if (length(selected_controls) > 0) {
      notes <- c(
        notes,
        paste0(
          "No selected controls were found in cleaned MET data. ",
          "Genotypes present in all environments were used as the control set."
        )
      )
    } else {
      notes <- c(
        notes,
        "No controls were selected. Genotypes present in all environments were used as the control set."
      )
    }
  }
  if (length(controls_used) == 0) {
    controls_used <- as.character(unique(dat_clean$Genotype))
    notes <- c(notes, "No genotype was present in all environments. All genotypes were used as the control set for direction-adjusted check comparison.")
  }
  low_conf_genos <- presence %>% filter(n_envs == 1) %>% pull(Genotype) %>% as.character()
  if (length(low_conf_genos) > 0) notes <- c(notes, paste0("Single-environment genotypes have high uncertainty: ", paste(low_conf_genos, collapse = ", ")))
  n_r <- dat_clean %>% group_by(Genotype, Environment) %>% summarise(n = n(), .groups = "drop") %>% summarise(harmonic_r = 1 / mean(1 / n)) %>% pull(harmonic_r)
  raw_genotype_summary <- dat_clean %>%
    mutate(Genotype = as.character(Genotype)) %>%
    group_by(Genotype) %>%
    summarise(
      Raw_Mean = mean(Weight, na.rm = TRUE),
      N_Raw = sum(!is.na(Weight)),
      SE_Raw_Mean = ifelse(N_Raw > 1, stats::sd(Weight, na.rm = TRUE) / sqrt(N_Raw), NA_real_),
      .groups = "drop"
    ) %>%
    dplyr::select(Genotype, Raw_Mean, SE_Raw_Mean)
  if (identical(model_type, "LMM")) {
  met_model_formula <- function(include_environment = TRUE,
                                include_genotype = TRUE,
                                include_gxe = TRUE,
                                include_design = TRUE) {
    terms <- character(0)
    if (isTRUE(include_environment)) {
      terms <- c(terms, "(1|Environment)")
    }
    if (isTRUE(include_design) && has_replication) {
      terms <- c(terms, "(1|Environment:Rep)")
    }
    if (isTRUE(include_design) && has_block) {
      terms <- c(
        terms,
        if (has_replication) "(1|Environment:Rep:Block)" else "(1|Environment:Block)"
      )
    }
    if (isTRUE(include_genotype)) {
      terms <- c(terms, "(1|Genotype)")
    }
    if (isTRUE(include_gxe)) {
      terms <- c(terms, "(1|Genotype:Environment)")
    }
    if (length(terms) == 0) {
      return(Weight ~ 1)
    }
    as.formula(paste("Weight ~", paste(terms, collapse = " + ")))
  }
  full_formula <- met_model_formula()
  m_full <- lmer(full_formula, data = dat_clean, control = lmerControl(optimizer = "bobyqa"))
  vc <- as.data.frame(VarCorr(m_full))
  total_var <- sum(vc$vcov, na.rm = TRUE)
  vc$percent <- round((vc$vcov / total_var) * 100, 2)
  vc$label <- recode(
    vc$grp,
    "Environment" = "Environment (E)",
    "Environment:Rep" = "Replication within environment",
    "Environment:Block" = "Block within environment",
    "Environment:Rep:Block" = "Block within replication/environment",
    "Genotype" = "Genotype (G)",
    "Genotype:Environment" = "GxE Interaction",
    "Residual" = "Residual",
    .default = vc$grp
  )
  variance_components <- vc %>% transmute(Source = label, Variance = round(vcov, 5), Percent = percent)
  get_vc <- function(name) {
    value <- vc$vcov[vc$grp == name]
    if (length(value) == 0 || all(is.na(value))) 0 else value[1]
  }
  g_var <- get_vc("Genotype")
  gxe_var <- get_vc("Genotype:Environment")
  res_var <- get_vc("Residual")
  n_r <- dat_clean %>% group_by(Genotype, Environment) %>% summarise(n = n(), .groups = "drop") %>% summarise(harmonic_r = 1 / mean(1 / n)) %>% pull(harmonic_r)
  H2 <- g_var / (g_var + gxe_var / n_envs_total + res_var / (n_r * n_envs_total))
  stability_ratio <- g_var / (g_var + gxe_var)
  genotype_summary <- dat_clean %>%
    mutate(
      Genotype = as.character(Genotype),
      Environment = as.character(Environment),
      Rep = as.character(Rep),
      Block = as.character(Block)
    ) %>%
    group_by(Genotype) %>%
    summarise(
      Mean = round(mean(Weight, na.rm = TRUE), 4),
      Locations_tested = paste(sort(unique(Environment)), collapse = ", "),
      N_observations = n(),
      .groups = "drop"
    )
  replications_by_location <- dat_clean %>%
    mutate(
      Genotype = as.character(Genotype),
      Environment = as.character(Environment),
      Rep = as.character(Rep)
    ) %>%
    distinct(Genotype, Environment, Rep) %>%
    arrange(Genotype, Environment, Rep) %>%
    group_by(Genotype, Environment) %>%
    summarise(Reps = n_distinct(Rep), .groups = "drop") %>%
    group_by(Genotype) %>%
    summarise(
      Replications_by_location = paste(paste0(Environment, ": ", Reps), collapse = " | "),
      .groups = "drop"
    )
  genotype_summary <- genotype_summary %>%
    left_join(replications_by_location, by = "Genotype") %>%
    dplyr::select(Genotype, Mean, Locations_tested, Replications_by_location, N_observations)
  if (has_block) {
    blocks_by_location <- dat_clean %>%
      mutate(
        Genotype = as.character(Genotype),
        Environment = as.character(Environment),
        Rep = as.character(Rep),
        Block = as.character(Block)
      ) %>%
      distinct(Genotype, Environment, Rep, Block) %>%
      arrange(Genotype, Environment, Rep, Block) %>%
      group_by(Genotype, Environment, Rep) %>%
      summarise(Blocks = paste(sort(unique(Block)), collapse = ", "), .groups = "drop") %>%
      group_by(Genotype, Environment) %>%
      summarise(
        Rep_blocks = paste(paste0("Rep ", Rep, " = ", Blocks), collapse = "; "),
        .groups = "drop"
      ) %>%
      group_by(Genotype) %>%
      summarise(
        Blocks_by_location = paste(paste0(Environment, ": ", Rep_blocks), collapse = " | "),
        .groups = "drop"
      )
    genotype_summary <- genotype_summary %>%
      left_join(blocks_by_location, by = "Genotype") %>%
      dplyr::select(Genotype, Mean, Locations_tested, Replications_by_location, Blocks_by_location, N_observations)
  }
  model_summary <- data.frame(Trait_used = trait_used, Model = model_label, N_rows_clean = nrow(dat_clean), N_genotypes = n_distinct(dat_clean$Genotype), N_environments = n_distinct(dat_clean$Environment), Replication_column = input$replication_col %||% "", Block_column = input$block_col %||% "", AMMI_GGE_min_observed_locations = min_envs_for_biplot, Trait_direction = trait_direction, Target_value = ifelse(identical(trait_direction, "Target trait"), target_value, NA_real_), Input_trait_weight = trait_weight, N_replications = if (has_replication) n_distinct(dat_clean$Rep) else NA_integer_, N_blocks = if (has_block) n_distinct(dat_clean$Block) else NA_integer_, Model_formula = paste(deparse(full_formula), collapse = " "), Harmonic_replication = round(n_r, 3), Stability_ratio_G_over_G_plus_GxE = round(stability_ratio, 4), Broad_sense_H2 = round(H2, 4), Controls_used = paste(controls_used, collapse = ", "), Notes = paste(notes, collapse = " | "), stringsAsFactors = FALSE)
  p_variance <- ggplot(variance_components %>% mutate(Source = fct_reorder(Source, -Percent)), aes(x = Source, y = Percent, fill = Source)) + geom_bar(stat = "identity", width = 0.6) + geom_text(aes(label = paste0(Percent, "%")), vjust = -0.5, size = 4) + scale_fill_manual(values = c("#9B59B6", "#3498DB", "#E67E22", "#2ECC71")) + labs(title = paste0("Variance partitioning â€” ", trait_used), x = NULL, y = "% of Total Variance") + theme_bw() + theme(legend.position = "none")
  res_df <- data.frame(fitted = fitted(m_full), residual = residuals(m_full))
  p_qq <- ggplot(res_df, aes(sample = residual)) + stat_qq(color = "#3498DB", alpha = 0.6) + stat_qq_line(color = "#E74C3C", linewidth = 0.8) + labs(title = "Normal Q-Q", x = "Theoretical", y = "Sample") + theme_bw()
  p_rvf <- ggplot(res_df, aes(x = fitted, y = residual)) + geom_point(color = "#3498DB", alpha = 0.5, size = 1.5) + geom_hline(yintercept = 0, linetype = "dashed", color = "#E74C3C", linewidth = 0.8) + geom_smooth(method = "loess", se = FALSE, color = "#F39C12", linewidth = 0.8, span = 0.8) + labs(title = "Residuals vs Fitted", x = "Fitted", y = "Residuals") + theme_bw()
  p_residual <- p_qq + p_rvf
  met_fit_reduced <- function(formula) {
    tryCatch(
      lmer(formula, data = dat_clean, control = lmerControl(optimizer = "bobyqa")),
      error = function(e) structure(list(error = e$message), class = "met_model_error")
    )
  }
  m_no_e <- met_fit_reduced(met_model_formula(include_environment = FALSE))
  m_no_gxe <- met_fit_reduced(met_model_formula(include_gxe = FALSE))
  m_no_g <- met_fit_reduced(met_model_formula(include_genotype = FALSE))
  m_null <- met_fit_reduced(met_model_formula(include_genotype = FALSE, include_gxe = FALSE))
  lrt_table <- bind_rows(met_safe_lrt(m_full, m_no_e, "E"), met_safe_lrt(m_full, m_no_gxe, "GxE"), met_safe_lrt(m_full, m_no_g, "G"), met_safe_lrt(m_full, m_null, "G+GxE"))
  get_lrt_p <- function(test_name) {
    p_col <- intersect(c("Pr(>Chisq)", "Pr..Chisq."), names(lrt_table))
    if (length(p_col) == 0) return(NA_real_)
    p_value <- lrt_table %>%
      filter(Test == test_name) %>%
      `[[`(p_col[1])
    p_value <- as.numeric(p_value)
    p_value <- p_value[!is.na(p_value)]
    if (length(p_value) == 0) NA_real_ else p_value[1]
  }
  lmm_p_values <- c(
    "Environment (E)" = get_lrt_p("E"),
    "Genotype (G)" = get_lrt_p("G"),
    "GxE Interaction" = get_lrt_p("GxE")
  )
  variance_components <- variance_components %>%
    mutate(
      P_value = round(as.numeric(lmm_p_values[Source]), 5),
      Sig = sig_label(P_value)
    )
  grand_mean <- fixef(m_full)[["(Intercept)"]]
  ranef_full <- ranef(m_full, condVar = TRUE)
  g_re <- ranef_full$Genotype
  post_var <- attr(ranef_full$Genotype, "postVar")
  se_g <- if (!is.null(post_var)) sqrt(post_var[1, 1, ]) else rep(NA_real_, nrow(g_re))
  reliability_g <- if (is.finite(g_var) && g_var > 0) {
    pmin(pmax(1 - (se_g^2 / g_var), 0), 1)
  } else {
    rep(NA_real_, length(se_g))
  }
  BLUPs_main <- data.frame(Genotype = rownames(g_re), BLUP_G = grand_mean + g_re[, 1], SE_G = se_g, Reliability = round(reliability_g, 4), CI_lower = grand_mean + g_re[, 1] - 1.96 * se_g, CI_upper = grand_mean + g_re[, 1] + 1.96 * se_g) %>%
    left_join(raw_genotype_summary, by = "Genotype") %>%
    mutate(
      .Favorable_BLUP_score = case_when(
        identical(trait_direction, "Lower better") ~ -BLUP_G,
        identical(trait_direction, "Target trait") & is.finite(target_value) ~ -abs(BLUP_G - target_value),
        TRUE ~ BLUP_G
      ),
      Trait_direction = trait_direction
    ) %>%
    arrange(desc(.Favorable_BLUP_score)) %>%
    mutate(Rank_BLUP = row_number()) %>%
    dplyr::select(Genotype, Raw_Mean, SE_Raw_Mean, BLUP_G, SE_G, Reliability, CI_lower, CI_upper, Trait_direction, Rank_BLUP)
  p_blup <- ggplot(BLUPs_main, aes(x = reorder(Genotype, -Rank_BLUP), y = BLUP_G)) + geom_point(color = "#2C3E50", size = 3) + geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.4, color = "#3498DB", linewidth = 0.7) + coord_flip() + labs(title = "Genotype BLUPs with 95% CI", subtitle = paste0("Ranked by decision direction: ", trait_direction), x = "Genotype", y = paste0("BLUP for ", trait_used)) + theme_bw()
  fixed_intercept_se <- tryCatch(sqrt(as.numeric(stats::vcov(m_full)[1, 1])), error = function(e) NA_real_)
  env_re <- ranef_full$Environment
  env_post_var <- attr(env_re, "postVar")
  env_effects <- env_re %>%
    rownames_to_column("Environment") %>%
    rename(BLUP_E = `(Intercept)`) %>%
    mutate(SE_E = if (!is.null(env_post_var)) sqrt(env_post_var[1, 1, ]) else NA_real_)
  gxe_re <- ranef_full$`Genotype:Environment`
  gxe_post_var <- attr(gxe_re, "postVar")
  BLUPs_GxE <- gxe_re %>%
    rownames_to_column("Geno_Env") %>%
    rename(BLUP_GxE = `(Intercept)`) %>%
    mutate(SE_GxE = if (!is.null(gxe_post_var)) sqrt(gxe_post_var[1, 1, ]) else NA_real_) %>%
    separate(Geno_Env, into = c("Genotype", "Environment"), sep = ":", extra = "merge", fill = "right")
  observed_cells <- dat_clean %>%
    group_by(Genotype, Environment) %>%
    summarise(Observed_Mean = mean(Weight, na.rm = TRUE), N_Replications = n_distinct(Rep), .groups = "drop") %>%
    mutate(Genotype = as.character(Genotype), Environment = as.character(Environment))
  BLUPs_env_obs <- BLUPs_GxE %>%
    left_join(BLUPs_main[, c("Genotype", "BLUP_G", "SE_G")], by = "Genotype") %>%
    left_join(env_effects, by = "Environment") %>%
    left_join(observed_cells, by = c("Genotype", "Environment")) %>%
    mutate(
      BLUP_env = grand_mean + (BLUP_G - grand_mean) + BLUP_E + BLUP_GxE,
      Prediction_SE = sqrt(fixed_intercept_se^2 + SE_G^2 + SE_E^2 + SE_GxE^2),
      Cell_Status = "TESTED",
      Source = "Observed",
      Evidence_Flag = "Observed-supported model estimate",
      Uncertainty_Method = "Approximate conditional SE; cross-effect covariance omitted",
      Reportable_Estimate = BLUP_env
    ) %>%
    arrange(Environment, desc(BLUP_env))
  all_combos <- expand.grid(Genotype = as.character(unique(dat_clean$Genotype)), Environment = as.character(unique(dat_clean$Environment)), stringsAsFactors = FALSE)
  imputed_cells <- all_combos %>%
    anti_join(observed_cells %>% dplyr::select(Genotype, Environment), by = c("Genotype", "Environment")) %>%
    left_join(BLUPs_main[, c("Genotype", "BLUP_G", "SE_G")], by = "Genotype") %>%
    left_join(env_effects, by = "Environment") %>%
    mutate(
      BLUP_GxE = 0,
      SE_GxE = sqrt(pmax(gxe_var, 0)),
      BLUP_env = grand_mean + (BLUP_G - grand_mean) + BLUP_E,
      Prediction_SE = sqrt(fixed_intercept_se^2 + SE_G^2 + SE_E^2 + SE_GxE^2),
      Observed_Mean = NA_real_,
      N_Replications = 0L,
      Cell_Status = "UNTESTED",
      Source = dplyr::if_else(Genotype %in% low_conf_genos, "Imputed_low_confidence", "Imputed"),
      Evidence_Flag = dplyr::if_else(
        Genotype %in% low_conf_genos,
        "Predicted-untested; low genotype coverage",
        "Predicted-untested; GxE assumed zero"
      ),
      Uncertainty_Method = "Approximate prediction SE including unobserved GxE variance",
      Reportable_Estimate = BLUP_env
    )
  BLUPs_env_full <- bind_rows(
    BLUPs_env_obs %>% dplyr::select(Genotype, Environment, BLUP_G, BLUP_E, BLUP_GxE, BLUP_env, Reportable_Estimate, Prediction_SE, Observed_Mean, N_Replications, Cell_Status, Source, Evidence_Flag, Uncertainty_Method),
    imputed_cells %>% dplyr::select(Genotype, Environment, BLUP_G, BLUP_E, BLUP_GxE, BLUP_env, Reportable_Estimate, Prediction_SE, Observed_Mean, N_Replications, Cell_Status, Source, Evidence_Flag, Uncertainty_Method)
  ) %>%
    mutate(
      Prediction_CI_Lower = BLUP_env - 1.96 * Prediction_SE,
      Prediction_CI_Upper = BLUP_env + 1.96 * Prediction_SE
    ) %>%
    arrange(Environment, desc(BLUP_env))
  acc_dat <- dat_clean %>% group_by(Genotype, Environment) %>% summarise(obs = mean(Weight), .groups = "drop") %>% mutate(Genotype = as.character(Genotype), Environment = as.character(Environment)) %>% left_join(BLUPs_env_full %>% dplyr::select(Genotype, Environment, BLUP_env), by = c("Genotype", "Environment")) %>% filter(!is.na(BLUP_env))
  r_by_env <- acc_dat %>% group_by(Environment) %>% summarise(r_val = ifelse(n() >= 2, round(cor(obs, BLUP_env, use = "complete.obs"), 4), NA_real_), .groups = "drop") %>% mutate(label = paste0("r = ", r_val))
  p_accuracy <- ggplot(acc_dat, aes(x = obs, y = BLUP_env, color = Genotype)) + geom_point(size = 3, alpha = 0.85) + geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7, alpha = 0.15) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40", linewidth = 0.7) + geom_text(data = r_by_env, aes(label = label, x = -Inf, y = Inf), hjust = -0.15, vjust = 1.4, size = 3.5, fontface = "bold", color = "#2C3E50", inherit.aes = FALSE) + facet_wrap(~Environment, scales = "free", ncol = 2) + labs(title = "Prediction accuracy per hybrid by location", x = paste0("Observed mean ", trait_used), y = "Predicted BLUP", color = "Genotype") + theme_bw() + theme(legend.position = "bottom")
  } else {
    met_anova_formula <- function(include_environment = TRUE,
                                  include_genotype = TRUE,
                                  include_gxe = TRUE,
                                  include_design = TRUE) {
      terms <- character(0)
      if (isTRUE(include_environment)) terms <- c(terms, "Environment")
      if (isTRUE(include_design) && has_replication) terms <- c(terms, "Environment:Rep")
      if (isTRUE(include_design) && has_block) {
        terms <- c(terms, if (has_replication) "Environment:Rep:Block" else "Environment:Block")
      }
      if (isTRUE(include_genotype)) terms <- c(terms, "Genotype")
      if (isTRUE(include_gxe) && isTRUE(include_environment) && isTRUE(include_genotype)) {
        terms <- c(terms, "Genotype:Environment")
      }
      if (length(terms) == 0) return(Weight ~ 1)
      as.formula(paste("Weight ~", paste(terms, collapse = " + ")))
    }
    full_formula <- met_anova_formula()
    m_full <- stats::lm(full_formula, data = dat_clean)
    anova_table <- as.data.frame(stats::anova(m_full)) %>% tibble::rownames_to_column("Source")
    names(anova_table) <- make.names(names(anova_table))
    anova_table$Display_source <- dplyr::recode(
      anova_table$Source,
      "Environment" = "Environment (E)",
      "Environment:Rep" = "Replication within environment",
      "Environment:Block" = "Block within environment",
      "Environment:Rep:Block" = "Block within replication/environment",
      "Genotype" = "Genotype (G)",
      "Genotype:Environment" = "GxE Interaction",
      "Residuals" = "Residual",
      .default = anova_table$Source
    )
    ss_total <- sum(anova_table$Sum.Sq, na.rm = TRUE)
    variance_components <- anova_table %>%
      transmute(
        Source = Display_source,
        Df = Df,
        Sum_sq = round(Sum.Sq, 5),
        Mean_square = round(Mean.Sq, 5),
        F_value = round(F.value, 5),
        Percent = ifelse(ss_total > 0, round(100 * Sum.Sq / ss_total, 2), NA_real_),
        P_value = round(`Pr..F.`, 5),
        Sig = sig_label(P_value)
      )
    get_ms <- function(source_name) {
      value <- anova_table$Mean.Sq[anova_table$Source == source_name]
      if (length(value) == 0 || all(is.na(value))) NA_real_ else as.numeric(value[1])
    }
    ms_g <- get_ms("Genotype")
    ms_gxe <- get_ms("Genotype:Environment")
    ms_error <- get_ms("Residuals")
    g_var <- if (is.finite(ms_g) && is.finite(ms_gxe)) max((ms_g - ms_gxe) / (n_r * n_envs_total), 0) else NA_real_
    gxe_var <- if (is.finite(ms_gxe) && is.finite(ms_error)) max((ms_gxe - ms_error) / n_r, 0) else NA_real_
    res_var <- ms_error
    H2 <- if (is.finite(g_var) && is.finite(gxe_var) && is.finite(res_var)) {
      g_var / (g_var + gxe_var / n_envs_total + res_var / (n_r * n_envs_total))
    } else {
      NA_real_
    }
    stability_ratio <- if (is.finite(g_var) && is.finite(gxe_var) && (g_var + gxe_var) > 0) {
      g_var / (g_var + gxe_var)
    } else {
      NA_real_
    }
    model_summary <- data.frame(Trait_used = trait_used, Model = model_label, N_rows_clean = nrow(dat_clean), N_genotypes = n_distinct(dat_clean$Genotype), N_environments = n_distinct(dat_clean$Environment), Replication_column = input$replication_col %||% "", Block_column = input$block_col %||% "", AMMI_GGE_min_observed_locations = min_envs_for_biplot, Trait_direction = trait_direction, Target_value = ifelse(identical(trait_direction, "Target trait"), target_value, NA_real_), Input_trait_weight = trait_weight, N_replications = if (has_replication) n_distinct(dat_clean$Rep) else NA_integer_, N_blocks = if (has_block) n_distinct(dat_clean$Block) else NA_integer_, Model_formula = paste(deparse(full_formula), collapse = " "), Harmonic_replication = round(n_r, 3), Stability_ratio_G_over_G_plus_GxE = round(stability_ratio, 4), Broad_sense_H2 = round(H2, 4), Controls_used = paste(controls_used, collapse = ", "), Notes = paste(notes, collapse = " | "), stringsAsFactors = FALSE)
    p_variance <- ggplot(variance_components %>% mutate(Source = fct_reorder(Source, -Percent)), aes(x = Source, y = Percent, fill = Source)) + geom_bar(stat = "identity", width = 0.6) + geom_text(aes(label = paste0(Percent, "%")), vjust = -0.5, size = 4) + scale_fill_manual(values = c("#9B59B6", "#3498DB", "#E67E22", "#2ECC71", "#7F8C8D", "#1ABC9C"), na.value = "#95A5A6") + labs(title = paste0("ANOVA partitioning - ", trait_used), subtitle = model_label, x = NULL, y = "% of Total Sum of Squares") + theme_bw() + theme(legend.position = "none")
    res_df <- data.frame(fitted = fitted(m_full), residual = residuals(m_full))
    p_qq <- ggplot(res_df, aes(sample = residual)) + stat_qq(color = "#3498DB", alpha = 0.6) + stat_qq_line(color = "#E74C3C", linewidth = 0.8) + labs(title = "Normal Q-Q", x = "Theoretical", y = "Sample") + theme_bw()
    p_rvf <- ggplot(res_df, aes(x = fitted, y = residual)) + geom_point(color = "#3498DB", alpha = 0.5, size = 1.5) + geom_hline(yintercept = 0, linetype = "dashed", color = "#E74C3C", linewidth = 0.8) + geom_smooth(method = "loess", se = FALSE, color = "#F39C12", linewidth = 0.8, span = 0.8) + labs(title = "Residuals vs Fitted", x = "Fitted", y = "Residuals") + theme_bw()
    p_residual <- p_qq + p_rvf
    lrt_table <- variance_components %>%
      transmute(Test = Source, Model = "Fixed-effect ANOVA", Df, Sum_sq, Mean_sq = Mean_square, F_value, `Pr(>F)` = P_value, Sig, check.names = FALSE)
    genotype_emm <- tryCatch(as.data.frame(emmeans::emmeans(m_full, ~ Genotype)), error = function(e) data.frame())
    if (nrow(genotype_emm) > 0 && all(c("Genotype", "emmean") %in% names(genotype_emm))) {
      main_values <- genotype_emm$emmean
      names(main_values) <- as.character(genotype_emm$Genotype)
      se_g <- genotype_emm$SE
      names(se_g) <- as.character(genotype_emm$Genotype)
    } else {
      main_values <- dat_clean %>% group_by(Genotype) %>% summarise(Value = mean(Weight, na.rm = TRUE), .groups = "drop") %>% tibble::deframe()
      se_g <- setNames(rep(NA_real_, length(main_values)), names(main_values))
    }
    BLUPs_main <- data.frame(Genotype = names(main_values), BLUP_G = as.numeric(main_values), SE_G = as.numeric(se_g), Reliability = NA_real_, stringsAsFactors = FALSE) %>%
      left_join(raw_genotype_summary, by = "Genotype") %>%
      mutate(
        CI_lower = ifelse(is.na(SE_G), NA_real_, BLUP_G - 1.96 * SE_G),
        CI_upper = ifelse(is.na(SE_G), NA_real_, BLUP_G + 1.96 * SE_G),
        .Favorable_BLUP_score = case_when(
          identical(trait_direction, "Lower better") ~ -BLUP_G,
          identical(trait_direction, "Target trait") & is.finite(target_value) ~ -abs(BLUP_G - target_value),
          TRUE ~ BLUP_G
        ),
        Trait_direction = trait_direction
      ) %>%
      arrange(desc(.Favorable_BLUP_score)) %>%
      mutate(Rank_BLUP = row_number()) %>%
      dplyr::select(Genotype, Raw_Mean, SE_Raw_Mean, BLUP_G, SE_G, Reliability, CI_lower, CI_upper, Trait_direction, Rank_BLUP)
    p_blup <- ggplot(BLUPs_main, aes(x = reorder(Genotype, -Rank_BLUP), y = BLUP_G)) + geom_point(color = "#2C3E50", size = 3) + geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.4, color = "#3498DB", linewidth = 0.7, na.rm = TRUE) + coord_flip() + labs(title = "Genotype adjusted means with 95% CI", subtitle = paste0("Ranked by decision direction: ", trait_direction), x = "Genotype", y = paste0("Adjusted mean for ", trait_used)) + theme_bw()
    observed_cells <- dat_clean %>%
      group_by(Genotype, Environment) %>%
      summarise(Observed_Mean = mean(Weight, na.rm = TRUE), N_Replications = n_distinct(Rep), .groups = "drop") %>%
      mutate(Genotype = as.character(Genotype), Environment = as.character(Environment))
    all_combos <- expand.grid(Genotype = as.character(unique(dat_clean$Genotype)), Environment = as.character(unique(dat_clean$Environment)), stringsAsFactors = FALSE)
    cell_emm <- tryCatch(as.data.frame(emmeans::emmeans(m_full, ~ Genotype:Environment)), error = function(e) data.frame())
    cell_values <- if (nrow(cell_emm) > 0 && all(c("Genotype", "Environment", "emmean") %in% names(cell_emm))) {
      cell_emm %>%
        transmute(
          Genotype = as.character(Genotype),
          Environment = as.character(Environment),
          Estimated_cell = suppressWarnings(as.numeric(emmean)),
          Estimated_SE = if ("SE" %in% names(cell_emm)) suppressWarnings(as.numeric(SE)) else NA_real_
        )
    } else {
      data.frame(Genotype = character(), Environment = character(), Estimated_cell = numeric(), Estimated_SE = numeric())
    }
    env_means <- dat_clean %>% group_by(Environment) %>% summarise(Env_mean = mean(Weight, na.rm = TRUE), .groups = "drop") %>% mutate(Environment = as.character(Environment))
    geno_means <- dat_clean %>% group_by(Genotype) %>% summarise(Gen_mean = mean(Weight, na.rm = TRUE), .groups = "drop") %>% mutate(Genotype = as.character(Genotype))
    grand_mean <- mean(dat_clean$Weight, na.rm = TRUE)
    BLUPs_env_full <- all_combos %>%
      left_join(observed_cells, by = c("Genotype", "Environment")) %>%
      left_join(cell_values, by = c("Genotype", "Environment")) %>%
      left_join(BLUPs_main[, c("Genotype", "BLUP_G")], by = "Genotype") %>%
      left_join(env_means, by = "Environment") %>%
      left_join(geno_means, by = "Genotype") %>%
      mutate(
        Additive_fallback = Gen_mean + Env_mean - grand_mean,
        BLUP_env = dplyr::coalesce(Estimated_cell, Observed_Mean, Additive_fallback),
        Prediction_SE = Estimated_SE,
        Prediction_CI_Lower = ifelse(is.na(Prediction_SE), NA_real_, BLUP_env - 1.96 * Prediction_SE),
        Prediction_CI_Upper = ifelse(is.na(Prediction_SE), NA_real_, BLUP_env + 1.96 * Prediction_SE),
        BLUP_E = Env_mean - grand_mean,
        BLUP_GxE = BLUP_env - grand_mean - (BLUP_G - grand_mean) - BLUP_E,
        Cell_Status = ifelse(!is.na(Observed_Mean), "TESTED", "UNTESTED"),
        Source = case_when(
          !is.na(Observed_Mean) ~ "Observed",
          !is.na(Estimated_cell) ~ "ANOVA_estimated",
          Genotype %in% low_conf_genos ~ "Imputed_low_confidence",
          TRUE ~ "Imputed_additive"
        ),
        Evidence_Flag = case_when(
          !is.na(Observed_Mean) & !is.na(Estimated_cell) ~ "Observed-supported adjusted estimate",
          !is.na(Observed_Mean) ~ "Tested; BLUE not estimable",
          !is.na(Estimated_cell) ~ "Predicted-untested adjusted estimate",
          Genotype %in% low_conf_genos ~ "Predicted-untested additive fallback; low genotype coverage",
          TRUE ~ "Predicted-untested additive fallback"
        ),
        Uncertainty_Method = case_when(
          !is.na(Estimated_SE) ~ "emmeans SE",
          TRUE ~ "Not estimable for additive fallback"
        ),
        Reportable_Estimate = Estimated_cell
      ) %>%
      dplyr::select(
        Genotype, Environment, BLUP_G, BLUP_E, BLUP_GxE, BLUP_env, Reportable_Estimate,
        Prediction_SE, Prediction_CI_Lower, Prediction_CI_Upper,
        Observed_Mean, N_Replications, Cell_Status, Source, Evidence_Flag, Uncertainty_Method
      ) %>%
      arrange(Environment, desc(BLUP_env))
    acc_dat <- dat_clean %>% group_by(Genotype, Environment) %>% summarise(obs = mean(Weight), .groups = "drop") %>% mutate(Genotype = as.character(Genotype), Environment = as.character(Environment)) %>% left_join(BLUPs_env_full %>% dplyr::select(Genotype, Environment, BLUP_env), by = c("Genotype", "Environment")) %>% filter(!is.na(BLUP_env))
    r_by_env <- acc_dat %>% group_by(Environment) %>% summarise(r_val = ifelse(n() >= 2, round(cor(obs, BLUP_env, use = "complete.obs"), 4), NA_real_), .groups = "drop") %>% mutate(label = paste0("r = ", r_val))
    p_accuracy <- ggplot(acc_dat, aes(x = obs, y = BLUP_env, color = Genotype)) + geom_point(size = 3, alpha = 0.85) + geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7, alpha = 0.15) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40", linewidth = 0.7) + geom_text(data = r_by_env, aes(label = label, x = -Inf, y = Inf), hjust = -0.15, vjust = 1.4, size = 3.5, fontface = "bold", color = "#2C3E50", inherit.aes = FALSE) + facet_wrap(~Environment, scales = "free", ncol = 2) + labs(title = "Observed vs adjusted mean by location", x = paste0("Observed mean ", trait_used), y = "Adjusted mean", color = "Genotype") + theme_bw() + theme(legend.position = "bottom")
  }
  if (!exists("genotype_summary", inherits = FALSE)) {
    genotype_summary <- dat_clean %>%
      mutate(
        Genotype = as.character(Genotype),
        Environment = as.character(Environment),
        Rep = as.character(Rep),
        Block = as.character(Block)
      ) %>%
      group_by(Genotype) %>%
      summarise(
        Mean = round(mean(Weight, na.rm = TRUE), 4),
        Locations_tested = paste(sort(unique(Environment)), collapse = ", "),
        N_observations = n(),
        .groups = "drop"
      )
    replications_by_location <- dat_clean %>%
      mutate(
        Genotype = as.character(Genotype),
        Environment = as.character(Environment),
        Rep = as.character(Rep)
      ) %>%
      distinct(Genotype, Environment, Rep) %>%
      arrange(Genotype, Environment, Rep) %>%
      group_by(Genotype, Environment) %>%
      summarise(Reps = n_distinct(Rep), .groups = "drop") %>%
      group_by(Genotype) %>%
      summarise(
        Replications_by_location = paste(paste0(Environment, ": ", Reps), collapse = " | "),
        .groups = "drop"
      )
    genotype_summary <- genotype_summary %>%
      left_join(replications_by_location, by = "Genotype") %>%
      dplyr::select(Genotype, Mean, Locations_tested, Replications_by_location, N_observations)
  }
  BLUPs_env_full <- BLUPs_env_full %>%
    mutate(
      Decision_score_env = met_ds_direction_score(BLUP_env, trait_direction, target_value),
      Trait_Direction = trait_direction
    )
  best_check_by_env <- BLUPs_env_full %>%
    filter(Genotype %in% controls_used) %>%
    group_by(Environment) %>%
    arrange(desc(Decision_score_env), Genotype, .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    transmute(
      Environment,
      Benchmark_Check = Genotype,
      Benchmark_Estimate = BLUP_env,
      Benchmark_SE = Prediction_SE,
      Benchmark_Decision_Score = Decision_score_env
    )
  BLUPs_env_full <- BLUPs_env_full %>%
    left_join(best_check_by_env, by = "Environment") %>%
    mutate(
      Check_Advantage = Decision_score_env - Benchmark_Decision_Score,
      Check_Advantage_SE = ifelse(
        Genotype == Benchmark_Check,
        0,
        sqrt(Prediction_SE^2 + Benchmark_SE^2)
      ),
      Check_Advantage_CI_Lower = ifelse(
        is.na(Check_Advantage_SE), NA_real_, Check_Advantage - 1.96 * Check_Advantage_SE
      ),
      Check_Advantage_CI_Upper = ifelse(
        is.na(Check_Advantage_SE), NA_real_, Check_Advantage + 1.96 * Check_Advantage_SE
      ),
      Probability_Superior = case_when(
        Genotype == Benchmark_Check ~ 0.5,
        is.na(Check_Advantage_SE) ~ NA_real_,
        Check_Advantage_SE == 0 ~ as.numeric(Check_Advantage > 0),
        TRUE ~ stats::pnorm(Check_Advantage / Check_Advantage_SE)
      ),
      Direction_adjusted_control_advantage = Check_Advantage,
      Control_mean = Benchmark_Estimate,
      Pct_of_controls = case_when(
        identical(trait_direction, "Lower better") ~ round((Benchmark_Estimate / BLUP_env) * 100, 1),
        identical(trait_direction, "Target trait") & is.finite(target_value) ~ round(Check_Advantage, 3),
        TRUE ~ round((BLUP_env / Benchmark_Estimate) * 100, 1)
      )
    )
  best_observed_check_by_env <- BLUPs_env_full %>%
    filter(Genotype %in% controls_used, Cell_Status == "TESTED", !is.na(Observed_Mean)) %>%
    mutate(Observed_Decision_Score = met_ds_direction_score(Observed_Mean, trait_direction, target_value)) %>%
    group_by(Environment) %>%
    arrange(desc(Observed_Decision_Score), Genotype, .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    transmute(
      Environment,
      Observed_Benchmark_Check = Genotype,
      Observed_Benchmark_Mean = Observed_Mean,
      Observed_Benchmark_Score = Observed_Decision_Score
    )
  heatmap_dat <- BLUPs_env_full %>%
    left_join(best_observed_check_by_env, by = "Environment") %>%
    mutate(
      Observed_Decision_Score = met_ds_direction_score(Observed_Mean, trait_direction, target_value),
      Observed_Check_Advantage = Observed_Decision_Score - Observed_Benchmark_Score,
      Estimated_Check_Advantage = ifelse(is.na(Reportable_Estimate), NA_real_, Check_Advantage),
      label = ifelse(
        is.na(Reportable_Estimate),
        "",
        paste0(round(Reportable_Estimate, 1), "\nAdv ", round(Estimated_Check_Advantage, 2))
      ),
      observed_label = ifelse(is.na(Observed_Mean), "", paste0(round(Observed_Mean, 1), "\nAdv ", round(Observed_Check_Advantage, 2))),
      alpha_val = ifelse(Cell_Status == "UNTESTED", 0.45, 1.0)
    )
  heatmap_limit <- max(abs(c(heatmap_dat$Estimated_Check_Advantage, heatmap_dat$Observed_Check_Advantage)), na.rm = TRUE)
  if (!is.finite(heatmap_limit) || heatmap_limit <= 0) heatmap_limit <- 1
  p_perf_heatmap_observed <- ggplot(heatmap_dat, aes(x = Environment, y = reorder(Genotype, Decision_score_env), fill = Observed_Check_Advantage)) + geom_tile(color = "white", linewidth = 0.5) + geom_text(aes(label = observed_label), size = 2.5, lineheight = 0.9) + scale_fill_gradient2(low = "#E74C3C", mid = "white", high = "#2ECC71", midpoint = 0, limits = c(-heatmap_limit, heatmap_limit), na.value = "#D9D9D9", oob = scales::squish) + labs(title = "Genotype x location observed performance", subtitle = paste0("Raw cell means relative to the best observed check. Direction: ", trait_direction, ". Gray cells were not tested."), x = "Location", y = "Genotype", fill = "Check advantage") + theme_bw()
  estimated_heatmap_note <- if (identical(model_type, "LMM")) {
    "Faded cells were not tested and are model predictions."
  } else {
    "Gray cells were not tested and have no estimable BLUE under the fixed GxE model."
  }
  p_perf_heatmap <- ggplot(heatmap_dat, aes(x = Environment, y = reorder(Genotype, Decision_score_env), fill = Estimated_Check_Advantage)) + geom_tile(aes(alpha = alpha_val), color = "white", linewidth = 0.5) + geom_text(aes(label = label), size = 2.5, lineheight = 0.9) + scale_fill_gradient2(low = "#E74C3C", mid = "white", high = "#2ECC71", midpoint = 0, limits = c(-heatmap_limit, heatmap_limit), na.value = "#D9D9D9", oob = scales::squish) + scale_alpha_identity() + labs(title = paste0("Genotype x location ", ifelse(identical(model_type, "LMM"), "BLUPs", "BLUEs")), subtitle = paste0("Estimated values relative to the best check. Direction: ", trait_direction, ". ", estimated_heatmap_note), x = "Location", y = "Genotype", fill = "Check advantage") + theme_bw()
  GxE_matrix_wide <- BLUPs_env_full %>% dplyr::select(Genotype, Environment, BLUP_env) %>% pivot_wider(names_from = Environment, values_from = BLUP_env) %>% column_to_rownames("Genotype")
  GxE_long_complete <- GxE_matrix_wide %>% rownames_to_column("Genotype") %>% pivot_longer(-Genotype, names_to = "Environment", values_to = "BLUP_env")
  n_genos <- nrow(GxE_matrix_wide)
  FW_dat_loo <- GxE_long_complete %>% group_by(Environment) %>% mutate(EnvIndex = (sum(BLUP_env) - BLUP_env) / (n_genos - 1)) %>% ungroup()
  FW_results <- FW_dat_loo %>%
    group_by(Genotype) %>%
    summarise(GenMean = mean(BLUP_env), Sens = coef(lm(BLUP_env ~ EnvIndex))[2], .groups = "drop") %>%
    mutate(
      Favorable_GenMean = case_when(
        identical(trait_direction, "Lower better") ~ -GenMean,
        identical(trait_direction, "Target trait") & is.finite(target_value) ~ -abs(GenMean - target_value),
        TRUE ~ GenMean
      ),
      Trait_direction = trait_direction,
      Target_value = ifelse(identical(trait_direction, "Target trait"), target_value, NA_real_),
      b_interp = case_when(Sens > 1.1 ~ "Responsive", Sens < 0.9 ~ "Stable", TRUE ~ "Average")
    ) %>%
    arrange(desc(Favorable_GenMean))
  p_fw_mean_sens <- ggplot(FW_results, aes(x = GenMean, y = Sens, color = b_interp, label = Genotype)) + geom_point(size = 3) + geom_text(vjust = -0.8, size = 3) + geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") + scale_color_manual(values = c("Responsive" = "#E74C3C", "Average" = "#F39C12", "Stable" = "#2ECC71")) + labs(title = "Finlay-Wilkinson: mean vs sensitivity", subtitle = paste0("Decision direction: ", trait_direction), x = paste("Genotype", estimate_label, "mean"), y = "Sensitivity (b)", color = "Stability") + theme_bw()
  p_fw_regression <- ggplot(FW_dat_loo, aes(x = EnvIndex, y = BLUP_env, color = Genotype, group = Genotype)) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
    geom_point(shape = 18, size = 3) +
    labs(
      title = paste0("Finlay-Wilkinson response correlation - ", trait_used),
      subtitle = "Genotype response across the leave-one-out environmental index.",
      x = "Leave-one-out environment index",
      y = paste0(estimate_label, " for ", trait_used),
      color = "Genotype"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  metan_single <- if (ammi_gge_available) {
    run_met_single_trait_extensions(dat_clean)
  } else {
    structure(
      list(error = paste0("AMMI/GGE requires at least ", MET_MIN_ENVS_FOR_BIPLOT, " locations.")),
      class = "si_ext_error"
    )
  }
  biplot_confidence <- presence %>%
    transmute(
      Genotype = as.character(Genotype),
      Observed_locations = as.integer(n_envs),
      Total_locations = n_envs_total,
      Min_locations_required = min_envs_for_biplot,
      Coverage_pct = round(Observed_locations / Total_locations * 100, 1),
      Confidence_flag = ifelse(Observed_locations >= Min_locations_required, "OK", "Low confidence")
    )
  low_conf_biplot <- biplot_confidence %>%
    filter(Confidence_flag == "Low confidence") %>%
    pull(Genotype)
  GxE_obs_wide <- GxE_matrix_wide
  ammi_notes <- data.frame(
    Status = ifelse(ammi_gge_available, "AVAILABLE", "NOT AVAILABLE"),
    Availability_reason = ifelse(
      ammi_gge_available,
      paste0("At least ", MET_MIN_ENVS_FOR_BIPLOT, " locations are available."),
      paste0("AMMI and GGE require at least ", MET_MIN_ENVS_FOR_BIPLOT, " locations; the core MET analysis remains available.")
    ),
    Data_source = estimate_matrix_source,
    Imputed_cell_warning = ifelse(
      ammi_gge_available,
      paste0("AMMI/GGE uses a full ", estimate_label, " GxE matrix; unobserved genotype-location cells are model-estimated."),
      "AMMI/GGE was not calculated."
    ),
    N_genotypes_enter_AMMI_GGE = ifelse(ammi_gge_available, nrow(GxE_obs_wide), 0L),
    N_genotypes_total = nrow(GxE_matrix_wide),
    Min_locations_required = min_envs_for_biplot,
    N_low_confidence_genotypes = length(low_conf_biplot),
    Low_confidence_genotypes = ifelse(length(low_conf_biplot) > 0, paste(low_conf_biplot, collapse = ", "), "none"),
    Excluded_genotypes = "none",
    stringsAsFactors = FALSE
  )
  AMMI_geno <- data.frame()
  AMMI_env <- data.frame()
  GGE_geno <- data.frame()
  GGE_env <- data.frame()
  gge_decision_notes <- data.frame(
    Decision_item = "Availability",
    Value = ifelse(ammi_gge_available, "Pending calculation", "Not available"),
    Breeder_interpretation = ifelse(
      ammi_gge_available,
      "GGE will be calculated when the matrix checks pass.",
      paste0("GGE requires at least ", MET_MIN_ENVS_FOR_BIPLOT, " locations; use the core MET estimates and genotype x location views.")
    ),
    stringsAsFactors = FALSE
  )
  gge_mean_stability_decision <- data.frame()
  gge_winner_decision <- data.frame()
  gge_environment_decision <- data.frame()
  gge_pc_variance <- data.frame()
  ammi_gge_unavailable_message <- paste0(
    "AMMI and GGE require at least ", MET_MIN_ENVS_FOR_BIPLOT,
    " locations. Core MET estimates are still available."
  )
  p_ammi1 <- empty_plot(ammi_gge_unavailable_message)
  p_ammi2 <- empty_plot(ammi_gge_unavailable_message)
  p_gge <- empty_plot(ammi_gge_unavailable_message)
  p_gge_mean_stability <- empty_plot(ammi_gge_unavailable_message)
  p_gge_which_won <- empty_plot(ammi_gge_unavailable_message)
  p_gge_environment <- empty_plot(ammi_gge_unavailable_message)
  p_gge_genotype_ranking <- empty_plot(ammi_gge_unavailable_message)
  p_gge_environment_ranking <- empty_plot(ammi_gge_unavailable_message)
  p_env_cor <- empty_plot("Environment correlation needs at least 2 environments.")
  if (ammi_gge_available && nrow(GxE_obs_wide) >= 2 && ncol(GxE_obs_wide) >= MET_MIN_ENVS_FOR_BIPLOT) {
    row_means <- rowMeans(GxE_obs_wide)
    col_means <- colMeans(GxE_obs_wide)
    grand_mn <- mean(as.matrix(GxE_obs_wide))
    GxE_centered <- as.matrix(GxE_obs_wide) - outer(row_means, rep(1, ncol(GxE_obs_wide))) - outer(rep(1, nrow(GxE_obs_wide)), col_means) + grand_mn
    svd_result <- svd(GxE_centered)
    n_pc <- min(nrow(GxE_centered) - 1, ncol(GxE_centered) - 1)
    SS_pc <- svd_result$d[1:n_pc]^2
    PC_pct <- SS_pc / sum(SS_pc)
    AMMI_geno <- as.data.frame(svd_result$u[, 1:min(2, n_pc), drop = FALSE]) %>% setNames(paste0("PC", 1:min(2, n_pc))) %>% mutate(Genotype = rownames(GxE_obs_wide), GenMean = row_means, .before = 1)
    AMMI_env <- as.data.frame(svd_result$v[, 1:min(2, n_pc), drop = FALSE]) %>% setNames(paste0("PC", 1:min(2, n_pc))) %>% mutate(Environment = colnames(GxE_obs_wide), EnvMean = col_means, .before = 1)
    for (i in 1:min(2, n_pc)) {
      AMMI_geno[[paste0("PC", i)]] <- AMMI_geno[[paste0("PC", i)]] * svd_result$d[i]
      AMMI_env[[paste0("PC", i)]] <- AMMI_env[[paste0("PC", i)]] * svd_result$d[i]
    }
    AMMI_geno$ASV <- if (n_pc >= 2) sqrt((PC_pct[1] / PC_pct[2] * AMMI_geno$PC1)^2 + AMMI_geno$PC2^2) else abs(AMMI_geno$PC1)
    AMMI_geno <- AMMI_geno %>%
      left_join(biplot_confidence, by = "Genotype") %>%
      arrange(ASV)
    p_ammi1 <- ggplot() +
      geom_point(data = AMMI_geno, aes(x = PC1, y = GenMean, color = Confidence_flag, shape = Confidence_flag), size = 3) +
      geom_text(data = AMMI_geno, aes(x = PC1, y = GenMean, label = Genotype, color = Confidence_flag), vjust = -0.8, size = 3, show.legend = FALSE) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
      geom_point(data = AMMI_env, aes(x = PC1, y = EnvMean), size = 4, shape = 17, color = "#E74C3C") +
      geom_text(data = AMMI_env, aes(x = PC1, y = EnvMean, label = Environment), vjust = -0.9, size = 3.5, color = "#E74C3C", fontface = "bold") +
      scale_color_manual(values = c("OK" = "#2C3E50", "Low confidence" = "#D35400"), breaks = "Low confidence", na.value = "gray60") +
      scale_shape_manual(values = c("OK" = 16, "Low confidence" = 1), breaks = "Low confidence", na.value = 16) +
      labs(title = paste0("AMMI1 biplot - ", trait_used, " ", estimate_label_plural), x = paste0("IPCA1 (", round(PC_pct[1] * 100, 1), "%)"), y = paste0("Mean ", estimate_label, " for ", trait_used), color = "Coverage", shape = "Coverage") +
      theme_bw()
    if (n_pc >= 2) {
      scale_ammi <- ifelse(max(abs(AMMI_env$PC1), na.rm = TRUE) == 0, 1, max(abs(AMMI_geno$PC1), na.rm = TRUE) / max(abs(AMMI_env$PC1), na.rm = TRUE) * 0.7)
      AMMI_env_sc <- AMMI_env %>% mutate(PC1 = PC1 * scale_ammi, PC2 = PC2 * scale_ammi)
      p_ammi2 <- ggplot() +
        geom_segment(data = AMMI_env_sc, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.25, "cm"), type = "closed"), color = "#E74C3C", linewidth = 0.8) +
        geom_text(data = AMMI_env_sc, aes(x = PC1 * 1.12, y = PC2 * 1.12, label = Environment), color = "#E74C3C", size = 3.5, fontface = "bold") +
        geom_point(data = AMMI_geno, aes(x = PC1, y = PC2, color = Confidence_flag, shape = Confidence_flag), size = 2.5) +
        geom_text(data = AMMI_geno, aes(x = PC1, y = PC2, label = Genotype, color = Confidence_flag), vjust = -0.8, size = 2.8, show.legend = FALSE) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
        scale_color_manual(values = c("OK" = "#2C3E50", "Low confidence" = "#D35400"), breaks = "Low confidence", na.value = "gray60") +
        scale_shape_manual(values = c("OK" = 16, "Low confidence" = 1), breaks = "Low confidence", na.value = 16) +
        labs(title = paste0("AMMI2 biplot - ", trait_used, " ", estimate_label_plural), x = paste0("IPCA1 (", round(PC_pct[1] * 100, 1), "%)"), y = paste0("IPCA2 (", round(PC_pct[2] * 100, 1), "%)"), color = "Coverage", shape = "Coverage") +
        theme_bw()
    }
    gge_views <- build_met_gge_decision_views(
      as.matrix(GxE_obs_wide),
      direction = trait_direction,
      target_value = target_value,
      trait_name = trait_used,
      trait_weight = trait_weight,
      confidence = biplot_confidence,
      cell_source = BLUPs_env_full %>% dplyr::select(Genotype, Environment, Source),
      check_genotypes = controls_used,
      estimate_matrix_source = estimate_matrix_source
    )
    gge_decision_notes <- gge_views$notes
    gge_mean_stability_decision <- gge_views$mean_stability
    gge_winner_decision <- gge_views$winners
    gge_environment_decision <- gge_views$environments
    gge_pc_variance <- gge_views$pc_variance
    GGE_geno <- gge_views$overview_genotype
    GGE_env <- gge_views$overview_environment
    p_gge_mean_stability <- gge_views$p_mean_stability
    p_gge_which_won <- gge_views$p_which_won
    p_gge_environment <- gge_views$p_environment
    p_gge_genotype_ranking <- gge_views$p_genotype_ranking
    p_gge_environment_ranking <- gge_views$p_environment_ranking
    p_gge <- gge_views$p_overview
  }
  if (nrow(GxE_obs_wide) >= 2 && ncol(GxE_obs_wide) >= 2) {
    cor_mat <- cor(GxE_obs_wide, use = "pairwise.complete.obs")
    cor_long <- as.data.frame(cor_mat) %>% rownames_to_column("Env1") %>% pivot_longer(-Env1, names_to = "Env2", values_to = "r")
    p_env_cor <- ggplot(cor_long, aes(x = Env1, y = Env2, fill = r)) + geom_tile(color = "white") + geom_text(aes(label = round(r, 2)), size = 4.5) + scale_fill_gradient2(low = "#E74C3C", mid = "white", high = "#2ECC71", midpoint = 0, limits = c(-1, 1)) + labs(title = paste("Genotype", estimate_label, "correlation across environments"), x = NULL, y = NULL, fill = "r") + theme_bw()
  }
  ammi_support <- build_met_ammi_support(
    if (ammi_gge_available) GxE_obs_wide else NULL,
    ammi_genotype = AMMI_geno,
    metan_single = if (ammi_gge_available) metan_single else NULL,
    trait_name = trait_used
  )
  selection <- build_met_selection_ranking(
    list(
      blups_main = BLUPs_main,
      fw_results = FW_results,
      ammi_genotype = AMMI_geno,
      trait_direction = trait_direction,
      target_value = target_value
    ),
    component_weights = c(mean = MET_W_YIELD, fw = MET_W_FW, asv = MET_W_ASV)
  )
  p_met_selection <- plot_met_selection_ranking(selection, trait_used, estimate_label = estimate_label)
  return(list(raw_data = df_raw, met_data = dat, met_cleaned_data = dat_clean, trait_direction = trait_direction, target_value = target_value, trait_weight = trait_weight, outlier_summary = outlier_summary, presence = presence, genotype_summary = genotype_summary, model_summary = model_summary, variance_components = variance_components, lrt_table = lrt_table, blups_main = BLUPs_main, blups_environment = BLUPs_env_full, gxe_matrix = GxE_matrix_wide %>% rownames_to_column("Genotype"), fw_results = FW_results, ammi_notes = ammi_notes, ammi_support = ammi_support, ammi_genotype = AMMI_geno, ammi_environment = AMMI_env, gge_decision_notes = gge_decision_notes, gge_mean_stability_decision = gge_mean_stability_decision, gge_winner_decision = gge_winner_decision, gge_environment_decision = gge_environment_decision, gge_pc_variance = gge_pc_variance, gge_genotype = GGE_geno, gge_environment = GGE_env, met_selection = selection, p_before = p_before, p_after = p_after, p_variance = p_variance, p_residual = p_residual, p_blup = p_blup, p_accuracy = p_accuracy, p_perf_heatmap_observed = p_perf_heatmap_observed, p_perf_heatmap = p_perf_heatmap, p_fw_mean_sens = p_fw_mean_sens, p_fw_regression = p_fw_regression, p_ammi1 = p_ammi1, p_ammi2 = p_ammi2, p_gge_mean_stability = p_gge_mean_stability, p_gge_which_won = p_gge_which_won, p_gge_environment = p_gge_environment, p_gge_genotype_ranking = p_gge_genotype_ranking, p_gge_environment_ranking = p_gge_environment_ranking, p_gge = p_gge, p_env_cor = p_env_cor, p_met_selection = p_met_selection))
}
run_met_all_traits <- function(
    df_raw,
    check_varieties = NULL,
    trait_cols = NULL,
    replication_col = NULL,
    block_col = NULL,
    min_envs_for_biplot = NULL,
    model_type = "LMM") {
  available_traits <- get_met_trait_cols(df_raw)
  if (is.null(trait_cols) || length(trait_cols) == 0) {
    trait_cols <- available_traits
  } else {
    trait_cols <- intersect(as.character(trait_cols), available_traits)
  }
  if (length(trait_cols) == 0) {
    stop("Choose at least one valid MET trait to run.")
  }
  results <- list()
  failures <- data.frame(Trait = character(), Error = character(), stringsAsFactors = FALSE)
  for (trait in trait_cols) {
    result <- tryCatch({
      run_met_pipeline(
        df_raw,
        trait,
        check_varieties = check_varieties,
        replication_col = replication_col,
        block_col = block_col,
        min_envs_for_biplot = min_envs_for_biplot,
        model_type = model_type
      )
    }, error = function(e) {
      failures <<- bind_rows(
        failures,
        data.frame(Trait = trait, Error = e$message, stringsAsFactors = FALSE)
      )
      NULL
    })
    if (!is.null(result)) {
      results[[trait]] <- result
    }
  }
  if (length(results) == 0) {
    stop("MET pipeline failed for all numeric traits: ", paste(failures$Error, collapse = " | "))
  }
  integrated <- build_met_integrated_ranking(df_raw, results)
  list(
    settings = data.frame(
      Trait = trait_cols,
      Ran_MET = trait_cols %in% names(results),
      Replication_column = replication_col %||% "",
      Block_column = block_col %||% "",
      Model = ifelse(toupper(as.character(model_type %||% "LMM")) %in% c("ANOVA", "RCBD", "ANOVA (RCBD)", "ANOVA_RCBD"), "ANOVA (RCBD)", "LMM"),
      AMMI_GGE_min_observed_locations = min_envs_for_biplot %||% "",
      stringsAsFactors = FALSE
    ),
    met_by_trait = results,
    met_trait_names = names(results),
    met_failed_traits = failures,
    met_integrated_ranking = integrated$ranking,
    met_integrated_trait_weights = integrated$trait_weights,
    met_integrated_adjusted = integrated$adjusted_performance,
    met_integrated_standardized = integrated$standardized_scores,
    p_met_integrated_ranking = integrated$plot
  )
}
