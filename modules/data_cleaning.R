# Reviewable, trait-level cleaning. Row numbers refer to uploaded data rows.
si_clean_value_flag <- function(flags) {
  grepl("IQR flag|Standardized residual flag|Below trait minimum|Above trait maximum|Outside declared score range|Outside percent range", flags)
}

si_clean_empty <- function() {
  data.frame(Row = integer(), Trait = character(), Value = character(),
             Flag = character(), Action = character(), Reason = character(),
             stringsAsFactors = FALSE)
}

si_clean_rows <- function(data, identity_col) {
  if (!identity_col %in% names(data)) stop("Choose an identifier column.")
  labels <- toupper(trimws(as.character(data[[identity_col]])))
  which(!labels %in% metadata_row_labels)
}

# One fit service is used by cleaning flags and raw/proposed diagnostic plots.
si_trial_diagnostic_fit <- function(data, trait, context) {
  identity <- context$identity %||% id_col
  rows <- si_clean_rows(data, identity)
  d <- data[rows, , drop = FALSE]
  d$.response <- to_number(d[[trait]])
  rownames(d) <- as.character(rows)
  kind <- context$kind %||% "generic"
  add_factor <- function(name, column) {
    if (is.null(column) || !column %in% names(d)) stop("Diagnostic design column is missing: ", column %||% name)
    values <- trimws(as.character(d[[column]])); values[!nzchar(values)] <- NA_character_
    d[[name]] <<- factor(values)
  }
  if (kind == "multifactor") {
    prepared <- mf_prepare_data(data, trait, context$factors, context$rep,
      context$design, context$randomization, identity, metadata_row_labels)
    fit <- mf_fit_model(prepared$data, prepared$factors, context$design,
      context$randomization, context$model)$model
    # mf_prepare_data preserves the order of complete observations.
    required <- c(trait, context$factors,
      if (context$randomization == "RCBD" || context$design == "split_plot") context$rep)
    complete <- is.finite(to_number(d[[trait]]))
    for (col in setdiff(required, trait)) complete <- complete & !is.na(d[[col]]) & nzchar(trimws(as.character(d[[col]])))
    used <- rows[complete]
  } else {
    add_factor(".G", identity)
    terms <- character()
    if (kind == "single") {
      terms <- ".G"
      if (context$model != "CRD") {
        add_factor(".Rep", context$rep)
        terms <- c(terms, if (context$model == "LMM") "(1 | .Rep)" else ".Rep")
      }
    } else if (kind == "met") {
      add_factor(".Env", context$environment); add_factor(".Rep", context$rep)
      if (!is.null(context$block) && nzchar(context$block)) add_factor(".Block", context$block)
      terms <- if (context$model == "ANOVA_RCBD") {
        c(".Env", ".Env:.Rep", ".G", ".G:.Env", if (".Block" %in% names(d)) ".Env:.Rep:.Block")
      } else {
        c("(1 | .Env)", "(1 | .Env:.Rep)", "(1 | .G)", "(1 | .G:.Env)",
          if (".Block" %in% names(d)) "(1 | .Env:.Rep:.Block)")
      }
    } else {
      predictors <- intersect(context$predictors %||% c(identity, rep_col), names(d))
      for (j in seq_along(predictors)) {
        name <- paste0(".p", j); add_factor(name, predictors[j])
        if (nlevels(d[[name]]) > 1) terms <- c(terms, name)
      }
    }
    form <- as.formula(paste(".response ~", if (length(terms)) paste(terms, collapse = " + ") else "1"))
    required <- unique(all.vars(form))
    valid <- is.finite(d$.response) & complete.cases(d[, required, drop = FALSE])
    d <- droplevels(d[valid, , drop = FALSE])
    if (nrow(d) < 3) stop("Too few complete observations for this trial model.")
    fit <- if (any(grepl("|", terms, fixed = TRUE))) {
      lme4::lmer(form, d, REML = TRUE, control = lme4::lmerControl(optimizer = "bobyqa"))
    } else stats::lm(form, d)
    used <- as.integer(rownames(stats::model.frame(fit)))
  }
  if (inherits(fit, "merMod")) {
    messages <- fit@optinfo$conv$lme4$messages
    messages <- messages[!grepl("boundary.*singular", messages)]
    if (length(messages)) stop("Diagnostic model did not converge: ", paste(messages, collapse = "; "))
  }
  residual <- as.numeric(residuals(fit)); fitted <- as.numeric(fitted(fit))
  h <- tryCatch(as.numeric(hatvalues(fit)), error = function(e) rep(NA_real_, length(residual)))
  scale <- stats::sigma(fit)
  standardized <- residual / (scale * sqrt(pmax(1 - h, 0)))
  standardized[!is.finite(standardized)] <- NA_real_
  if (df.residual(fit) < 3) standardized[] <- NA_real_
  observations <- data.frame(Row = used, Observed = as.numeric(model.response(model.frame(fit))),
    Fitted = fitted, Residual = residual, Standardized = standardized)
  label <- paste(deparse(formula(fit)), collapse = " ")
  for (name in c(".response", ".G", ".Env", ".Rep", ".Block")) {
    replacement <- switch(name, .response = trait, .G = "Genotype", .Env = "Environment", .Rep = "Replication", .Block = "Block")
    label <- gsub(name, replacement, label, fixed = TRUE)
  }
  if (identical(context$model, "SOMMER")) label <- paste("Preliminary lme4 CS approximation for SOMMER:", label)
  attr(fit, "si_diagnostic_observations") <- observations
  attr(fit, "si_diagnostic_label") <- label
  list(model = fit, observations = observations, label = label)
}

si_clean_screen <- function(data, traits, identity_col, predictors = character(),
                            method = "residual", threshold = 3,
                            lower = NA_real_, upper = NA_real_, context = NULL) {
  rows <- si_clean_rows(data, identity_col)
  result <- si_clean_empty()
  notes <- character()
  predictors <- intersect(predictors, names(data))
  for (trait in traits) {
    value <- to_number(data[[trait]][rows])
    flags <- rep("", length(rows))
    add_flag <- function(mask, label) {
      mask[is.na(mask)] <- FALSE
      flags[mask] <<- ifelse(nzchar(flags[mask]), paste(flags[mask], label, sep = "; "), label)
    }
    add_flag(!is.finite(value), "Missing or non-numeric/non-finite value")
    add_flag(is.finite(value) & is.finite(lower) & value < lower, "Below trait minimum")
    add_flag(is.finite(value) & is.finite(upper) & value > upper, "Above trait maximum")
    metadata <- si_trait_metadata(data, trait, identity_col)
    if (metadata$Type[1] == "Score 1-5" && metadata$Metadata_source[1] != "Inferred (legacy)") {
      add_flag(is.finite(value) & (value < 1 | value > 5), "Outside declared score range 1-5")
    }
    if (metadata$Unit[1] == "percent") add_flag(is.finite(value) & (value < 0 | value > 100), "Outside percent range 0-100")
    good <- is.finite(value)
    if (method == "iqr" && sum(good) >= 4) {
      q <- quantile(value[good], c(.25, .75), names = FALSE)
      add_flag(good & (value < q[1] - threshold * diff(q) | value > q[2] + threshold * diff(q)), "IQR flag")
    } else if (method == "residual") {
      active_context <- context %||% list(kind = "generic", identity = identity_col, predictors = predictors, model = "LM")
      score <- tryCatch({
        diagnostic <- si_trial_diagnostic_fit(data, trait, active_context)
        notes <- c(notes, paste(trait, "screening model:", diagnostic$label))
        if (nrow(diagnostic$observations) < length(rows)) notes <- c(notes,
          paste(trait, length(rows) - nrow(diagnostic$observations), "observations unavailable to the diagnostic model because of missing response/design values."))
        diagnostic$observations$Standardized[match(rows, diagnostic$observations$Row)]
      }, error = function(e) {
        notes <<- c(notes, paste(trait, "residual screening unavailable:", conditionMessage(e)))
        rep(NA_real_, length(rows))
      })
      add_flag(is.finite(score) & abs(score) > threshold, "Standardized residual flag")
      if (!any(is.finite(score))) notes <- c(notes, paste(trait, "has no usable residual scores; inspect replication and model terms."))
    }
    key_cols <- unique(c(identity_col, predictors))
    key <- data[rows, key_cols, drop = FALSE]
    duplicate <- duplicated(key) | duplicated(key, fromLast = TRUE)
    if (length(predictors)) add_flag(duplicate, "Repeated identifier/design combination: verify plot identity")
    missing_id <- is.na(data[[identity_col]][rows]) | !nzchar(trimws(as.character(data[[identity_col]][rows])))
    add_flag(missing_id, "Missing identifier")
    flagged <- seq_along(rows)
    if (length(flagged)) {
      tab <- data.frame(Row = rows[flagged], Trait = trait,
        Value = as.character(data[[trait]][rows[flagged]]), Flag = flags[flagged],
        Action = "Keep", Reason = "", stringsAsFactors = FALSE)
      result <- rbind(result, tab)
    }
  }
  list(review = result, notes = unique(notes))
}

si_clean_apply <- function(data, review) {
  out <- data
  exclusions <- review[review$Action == "Exclude", , drop = FALSE]
  if (any(!review$Action %in% c("Keep", "Exclude"))) stop("Choose Keep or Exclude.")
  if (nrow(exclusions) && any(!nzchar(trimws(exclusions$Reason)))) stop("Enter a reason for every exclusion.")
  for (i in seq_len(nrow(exclusions))) {
    row <- exclusions$Row[i]
    trait <- exclusions$Trait[i]
    if (!trait %in% names(out) || row < 1 || row > nrow(out)) stop("Invalid exclusion reference.")
    out[[trait]][row] <- NA
  }
  out
}
