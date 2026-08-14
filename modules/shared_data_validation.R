# Shared data validation and pipeline-safety helpers for Modules 1-5.
# These functions are intentionally independent of the Shiny server so they can
# be reused by LPSI, MET, breeding, and mating workflows.

si_read_excel_upload <- function(path, name) {
  ext <- tolower(tools::file_ext(name))
  if (!ext %in% c("xlsx", "xls")) {
    stop("Please upload Excel file only: .xlsx or .xls", call. = FALSE)
  }
  as.data.frame(readxl::read_excel(path, .name_repair = "minimal"), check.names = FALSE)
}

si_to_number <- function(x) {
  suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
}

si_blank <- function(x) {
  is.na(x) | trimws(as.character(x)) == ""
}

si_detect_numeric_like_traits <- function(df, protected_cols = character(0), min_numeric = 1) {
  candidates <- setdiff(names(df), protected_cols)
  candidates[vapply(candidates, function(col) {
    sum(!is.na(si_to_number(df[[col]]))) >= min_numeric
  }, logical(1))]
}

si_validate_uploaded_table <- function(df, id_col = "Variety", rep_col = "Rep",
                                       remove_cols = character(0),
                                       metadata_labels = c("WEIGHT", "WEIGHTS", "IMPORTANCE",
                                                           "IMPORTANT", "DIRECTION",
                                                           "DIRECTIONS", "TRAIT_DIRECTION")) {
  errors <- character(0)
  warnings <- character(0)

  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0 || ncol(df) == 0) {
    return(list(ok = FALSE, errors = "The uploaded file appears to be empty.", warnings = warnings,
                trait_cols = character(0)))
  }

  nm <- names(df)
  if (any(is.na(nm) | trimws(nm) == "")) {
    errors <- c(errors, "One or more columns have no name. Please give every column a header.")
  }
  if (any(duplicated(nm))) {
    errors <- c(errors, paste0("Duplicate column name(s): ",
                               paste(unique(nm[duplicated(nm)]), collapse = ", "),
                               ". Column names must be unique."))
  }
  if (!id_col %in% nm) {
    errors <- c(errors, paste0("Required ID column is missing: ", id_col, "."))
  }
  if (!rep_col %in% nm) {
    warnings <- c(warnings, paste0("Replication column was not found: ", rep_col,
                                   ". LPSI diagnostics and RCBD-adjusted means may be limited."))
  }

  data_rows <- seq_len(nrow(df))
  if (id_col %in% nm) {
    id_upper <- toupper(trimws(as.character(df[[id_col]])))
    metadata_rows <- which(id_upper %in% metadata_labels)
    if (length(metadata_rows) > 0) {
      data_rows <- setdiff(data_rows, metadata_rows)
    }
    blank_ids <- sum(si_blank(df[[id_col]][data_rows]))
    if (blank_ids > 0) {
      warnings <- c(warnings, paste0(blank_ids, " data row(s) have blank ", id_col,
                                     " and will be ignored by most pipelines."))
    }
    id_trim <- trimws(as.character(df[[id_col]][data_rows]))
    canon <- tolower(gsub("[[:space:]]+", "", id_trim))
    nonblank_id <- !si_blank(id_trim)
    if (any(nonblank_id)) {
      by_canon <- tapply(id_trim[nonblank_id], canon[nonblank_id],
                         function(x) length(unique(x)))
      inconsistent <- names(by_canon)[by_canon > 1]
      if (length(inconsistent) > 0) {
        examples <- unique(unlist(lapply(head(inconsistent, 3), function(k) {
          unique(id_trim[canon == k])
        })))
        warnings <- c(warnings, paste0("Some genotype names differ only by spaces/capitals: ",
                                       paste(head(examples, 6), collapse = ", "),
                                       ". They may be treated as different genotypes."))
      }
    }
  }

  protected <- unique(c(id_col, rep_col, remove_cols))
  trait_cols <- si_detect_numeric_like_traits(df[data_rows, , drop = FALSE], protected)
  if (length(trait_cols) == 0) {
    errors <- c(errors, "No numeric-like trait columns were detected after excluding ID/replication columns.")
  } else {
    bad_cells <- vapply(trait_cols, function(tr) {
      raw <- df[[tr]][data_rows]
      num <- si_to_number(raw)
      sum(is.na(num) & !si_blank(raw))
    }, integer(1))
    bad_traits <- names(bad_cells)[bad_cells > 0]
    if (length(bad_traits) > 0) {
      warnings <- c(warnings, paste0(
        "Some trait cells contain text and will become missing: ",
        paste(paste0(bad_traits, " (", bad_cells[bad_traits], ")"), collapse = ", "),
        "."
      ))
    }
  }

  list(ok = length(errors) == 0, errors = errors, warnings = unique(warnings),
       trait_cols = trait_cols)
}

si_validate_uploaded_file <- function(df, remove_cols = character(0)) {
  errors <- character(0)
  warnings <- character(0)

  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0 || ncol(df) == 0) {
    return(list(ok = FALSE, errors = "The uploaded file appears to be empty.",
                warnings = warnings, trait_cols = character(0)))
  }

  nm <- names(df)
  if (any(is.na(nm) | trimws(nm) == "")) {
    errors <- c(errors, "One or more columns have no name. Please give every column a header.")
  }
  if (any(duplicated(nm))) {
    errors <- c(errors, paste0("Duplicate column name(s): ",
                               paste(unique(nm[duplicated(nm)]), collapse = ", "),
                               ". Column names must be unique."))
  }

  trait_cols <- si_detect_numeric_like_traits(df, protected_cols = remove_cols)
  if (length(trait_cols) == 0) {
    warnings <- c(warnings, "No numeric-like trait columns were detected yet. Some analyses may not run.")
  }

  list(ok = length(errors) == 0, errors = errors, warnings = unique(warnings),
       trait_cols = trait_cols)
}

si_format_validation_report <- function(report) {
  if (is.null(report)) return("Validation report is not available.")
  lines <- c(
    paste0("Validation: ", if (isTRUE(report$ok)) "passed" else "needs attention"),
    paste0("Detected numeric-like traits: ",
           if (length(report$trait_cols) > 0) paste(report$trait_cols, collapse = ", ") else "none")
  )
  if (length(report$errors) > 0) {
    lines <- c(lines, "Errors:", paste0("- ", report$errors))
  }
  if (length(report$warnings) > 0) {
    lines <- c(lines, "Warnings:", paste0("- ", report$warnings))
  }
  paste(lines, collapse = "\n")
}

si_reset_analysis_state <- function(analysis_results, analysis_used, analysis_message,
                                    saved_results, reason = "Input changed.") {
  analysis_results(NULL)
  analysis_used(NULL)
  analysis_message(paste(reason, "Run an analysis again."))
  for (nm in c("MATING", "BREEDING", "LPSI", "MET", "DIVERSITY")) {
    saved_results[[nm]] <- NULL
  }
  invisible(TRUE)
}

si_safe_table <- function(x) {
  if (is.null(x)) return(data.frame())
  tryCatch(as.data.frame(x, stringsAsFactors = FALSE), error = function(e) data.frame())
}

si_first_existing <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) NA_character_ else hit[1]
}

si_lpsi_pipeline_review <- function(results) {
  trait_info <- si_safe_table(results$trait_info)
  final_decision <- si_safe_table(results$final_decision)
  checks <- results$selected_checks
  if (is.null(checks)) checks <- character(0)

  weight_col <- si_first_existing(trait_info, c("Normalized_weight", "Raw_weight", "Weight"))
  max_weight <- if (!is.na(weight_col) && nrow(trait_info) > 0) {
    suppressWarnings(max(as.numeric(trait_info[[weight_col]]), na.rm = TRUE))
  } else {
    NA_real_
  }
  if (!is.finite(max_weight)) max_weight <- NA_real_

  decision_col <- si_first_existing(final_decision, c("Decision", "Final_decision", "Status", "Recommendation"))
  decision_text <- if (!is.na(decision_col) && nrow(final_decision) > 0) {
    paste(capture.output(print(table(final_decision[[decision_col]], useNA = "ifany"))), collapse = " | ")
  } else {
    "No explicit decision/status column found."
  }

  severe_col <- si_first_existing(final_decision, c("n_priority_severe_weak", "Severe_weakness_count"))
  severe_count <- if (!is.na(severe_col) && nrow(final_decision) > 0) {
    sum(suppressWarnings(as.numeric(final_decision[[severe_col]])) > 0, na.rm = TRUE)
  } else {
    NA_integer_
  }

  data.frame(
    Module = "LPSI",
    Check = c(
      "Trait coverage",
      "Benchmark checks",
      "Weight dominance",
      "Decision distribution",
      "Priority weakness screen"
    ),
    Result = c(
      paste0(nrow(trait_info), " trait(s) included"),
      if (length(checks) > 0) paste(checks, collapse = ", ") else "No selected check recorded",
      if (is.na(max_weight)) "Weight information not found" else paste0("Largest trait weight = ", round(max_weight, 4)),
      decision_text,
      if (is.na(severe_count)) "No severe-weakness column found" else paste0(severe_count, " genotype(s) with severe priority weakness")
    ),
    Why_this_matters = c(
      "LPSI is a multi-trait decision tool; missing traits change the final index.",
      "The benchmark defines above-check and below-check interpretation.",
      "A single dominant weight can make the index behave like direct selection on one trait.",
      "A compact count helps verify whether the pipeline is separating advance/retest/drop groups.",
      "Priority weaknesses protect against selecting genotypes with unacceptable key-trait defects."
    ),
    stringsAsFactors = FALSE
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

si_selection_bridge <- function(lpsi_results = NULL, met_results = NULL) {
  lpsi <- si_safe_table(if (!is.null(lpsi_results)) lpsi_results$final_decision else NULL)
  met <- si_safe_table(if (!is.null(met_results)) met_results$met_integrated_ranking else NULL)
  if (nrow(lpsi) == 0 && nrow(met) == 0) return(data.frame())

  lpsi_gen_col <- si_first_existing(lpsi, c("Original_ID", "Genotype", "GEN", "ID"))
  met_gen_col <- si_first_existing(met, c("Genotype", "Original_ID", "GEN", "ID"))

  lpsi_bridge <- if (!is.na(lpsi_gen_col) && nrow(lpsi) > 0) {
    rank_col <- si_first_existing(lpsi, c("Final_rank", "Rank", "Selection_rank"))
    decision_col <- si_first_existing(lpsi, c("Decision", "Final_decision", "Status", "Recommendation"))
    data.frame(
      Genotype = as.character(lpsi[[lpsi_gen_col]]),
      LPSI_rank = if (!is.na(rank_col)) suppressWarnings(as.numeric(lpsi[[rank_col]])) else seq_len(nrow(lpsi)),
      LPSI_decision = if (!is.na(decision_col)) as.character(lpsi[[decision_col]]) else "",
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(Genotype = character(0), LPSI_rank = numeric(0), LPSI_decision = character(0))
  }

  met_bridge <- if (!is.na(met_gen_col) && nrow(met) > 0) {
    rank_col <- si_first_existing(met, c("Integrated_rank", "Rank", "MET_rank"))
    score_col <- si_first_existing(met, c("Integrated_MET_Index", "MET_Index", "Selection_score"))
    data.frame(
      Genotype = as.character(met[[met_gen_col]]),
      MET_rank = if (!is.na(rank_col)) suppressWarnings(as.numeric(met[[rank_col]])) else seq_len(nrow(met)),
      MET_index = if (!is.na(score_col)) suppressWarnings(as.numeric(met[[score_col]])) else NA_real_,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(Genotype = character(0), MET_rank = numeric(0), MET_index = numeric(0))
  }

  all_genotypes <- sort(unique(c(lpsi_bridge$Genotype, met_bridge$Genotype)))
  out <- merge(data.frame(Genotype = all_genotypes, stringsAsFactors = FALSE), lpsi_bridge,
               by = "Genotype", all.x = TRUE)
  out <- merge(out, met_bridge, by = "Genotype", all.x = TRUE)
  out$Selected_by_LPSI <- !is.na(out$LPSI_rank)
  out$Selected_by_MET <- !is.na(out$MET_rank)
  out$Support_count <- as.integer(out$Selected_by_LPSI) + as.integer(out$Selected_by_MET)
  out$Bridge_recommendation <- ifelse(
    out$Support_count == 2, "Strong candidate: supported by LPSI and MET",
    ifelse(out$Selected_by_LPSI, "Trait-index candidate: confirm stability",
           ifelse(out$Selected_by_MET, "Stable MET candidate: check trait-index fit",
                  "Not selected by current summaries"))
  )
  out[order(-out$Support_count, out$LPSI_rank, out$MET_rank, na.last = TRUE), ]
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

si_lpsi_trait_direction <- function(results, trait) {
  info <- si_safe_table(results$trait_info)
  if (!"Trait" %in% names(info) || !"Direction" %in% names(info)) return("Higher better")
  hit <- info[as.character(info$Trait) == trait, , drop = FALSE]
  if (nrow(hit) == 0 || is.na(hit$Direction[1])) "Higher better" else as.character(hit$Direction[1])
}

si_lpsi_direct_selection <- function(results, primary_trait = NULL, selection_pct = 15) {
  means <- si_safe_table(results$actual_adjusted_means)
  final <- si_safe_table(results$final_decision)
  info <- si_safe_table(results$trait_info)
  if (nrow(means) == 0) return(data.frame())

  trait_candidates <- if ("Trait" %in% names(info)) as.character(info$Trait) else
    setdiff(names(means), c("ID", "Original_ID", "Check_Benchmark"))
  trait_candidates <- trait_candidates[trait_candidates %in% names(means)]
  if (length(trait_candidates) == 0) return(data.frame())
  if (is.null(primary_trait) || !primary_trait %in% trait_candidates) primary_trait <- trait_candidates[1]

  check_ids <- character(0)
  if (!is.null(results$selected_checks)) check_ids <- as.character(results$selected_checks)
  candidate_means <- means
  if ("Original_ID" %in% names(candidate_means) && length(check_ids) > 0) {
    candidate_means <- candidate_means[!as.character(candidate_means$Original_ID) %in% check_ids, , drop = FALSE]
  }
  if (nrow(candidate_means) == 0 && nrow(final) > 0) candidate_means <- means

  direction <- si_lpsi_trait_direction(results, primary_trait)
  values <- suppressWarnings(as.numeric(candidate_means[[primary_trait]]))
  target <- NA_real_
  if ("Trait" %in% names(info) && "Target_value" %in% names(info)) {
    tv <- suppressWarnings(as.numeric(info$Target_value[match(primary_trait, info$Trait)]))
    if (length(tv) > 0 && is.finite(tv[1])) target <- tv[1]
  }
  order_value <- if (direction == "Lower better") values else -values
  if (direction == "Target trait" && is.finite(target)) order_value <- abs(values - target)
  ord <- order(order_value, na.last = TRUE)
  candidate_means <- candidate_means[ord, , drop = FALSE]
  n_sel <- max(1, round(nrow(candidate_means) * suppressWarnings(as.numeric(selection_pct)) / 100))
  n_sel <- min(n_sel, nrow(candidate_means))

  out <- data.frame(
    Rank_Direct = seq_len(nrow(candidate_means)),
    ID = if ("ID" %in% names(candidate_means)) as.character(candidate_means$ID) else seq_len(nrow(candidate_means)),
    Original_ID = if ("Original_ID" %in% names(candidate_means)) as.character(candidate_means$Original_ID) else as.character(candidate_means[[1]]),
    Primary_trait = primary_trait,
    Direction = direction,
    Trait_value = values[ord],
    Selected_Direct = seq_len(nrow(candidate_means)) <= n_sel,
    stringsAsFactors = FALSE
  )
  out
}

si_lpsi_method_comparison <- function(results, primary_trait = NULL, selection_pct = 15) {
  final <- si_safe_table(results$final_decision)
  direct <- si_lpsi_direct_selection(results, primary_trait, selection_pct)
  if (nrow(final) == 0 && nrow(direct) == 0) return(data.frame())

  n_sel <- if (nrow(direct) > 0) sum(direct$Selected_Direct, na.rm = TRUE) else
    max(1, round(nrow(final) * suppressWarnings(as.numeric(selection_pct)) / 100))

  lpsi_gen <- si_first_existing(final, c("Original_ID", "Genotype", "GEN", "ID"))
  lpsi_rank_col <- si_first_existing(final, c("Final_rank", "Rank", "Selection_rank"))
  lpsi_decision_col <- si_first_existing(final, c("Decision", "Final_decision", "Status", "Recommendation"))
  lpsi_tbl <- if (!is.na(lpsi_gen) && nrow(final) > 0) {
    lpsi_rank <- if (!is.na(lpsi_rank_col)) suppressWarnings(as.numeric(final[[lpsi_rank_col]])) else seq_len(nrow(final))
    selected_by_decision <- if (!is.na(lpsi_decision_col)) {
      grepl("advance|select", as.character(final[[lpsi_decision_col]]), ignore.case = TRUE)
    } else rep(FALSE, nrow(final))
    if (!any(selected_by_decision, na.rm = TRUE)) selected_by_decision <- lpsi_rank <= n_sel
    data.frame(
      Genotype = as.character(final[[lpsi_gen]]),
      LPSI_rank = lpsi_rank,
      LPSI_decision = if (!is.na(lpsi_decision_col)) as.character(final[[lpsi_decision_col]]) else "",
      Selected_LPSI = selected_by_decision,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(Genotype = character(0), LPSI_rank = numeric(0), LPSI_decision = character(0), Selected_LPSI = logical(0))
  }

  direct_tbl <- if (nrow(direct) > 0) {
    data.frame(
      Genotype = as.character(direct$Original_ID),
      Direct_rank = direct$Rank_Direct,
      Direct_trait = direct$Primary_trait,
      Direct_value = direct$Trait_value,
      Selected_Direct = direct$Selected_Direct,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(Genotype = character(0), Direct_rank = numeric(0), Direct_trait = character(0), Direct_value = numeric(0), Selected_Direct = logical(0))
  }

  all_gen <- sort(unique(c(lpsi_tbl$Genotype, direct_tbl$Genotype)))
  out <- merge(data.frame(Genotype = all_gen, stringsAsFactors = FALSE), lpsi_tbl, by = "Genotype", all.x = TRUE)
  out <- merge(out, direct_tbl, by = "Genotype", all.x = TRUE)
  out$Selected_LPSI[is.na(out$Selected_LPSI)] <- FALSE
  out$Selected_Direct[is.na(out$Selected_Direct)] <- FALSE
  out$Support_count <- as.integer(out$Selected_LPSI) + as.integer(out$Selected_Direct)
  out$Recommendation <- ifelse(out$Support_count == 2, "Supported by LPSI and Direct Selection",
                               ifelse(out$Selected_LPSI, "LPSI candidate: check primary trait",
                                      ifelse(out$Selected_Direct, "Direct-selection candidate: check multi-trait balance",
                                             "Not selected by current methods")))
  out[order(-out$Support_count, out$LPSI_rank, out$Direct_rank, na.last = TRUE), , drop = FALSE]
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

si_lpsi_recommendation_table <- function(results) {
  final <- si_safe_table(if (!is.null(results)) results$final_decision else NULL)
  if (nrow(final) == 0) return(data.frame())

  genotype_col <- si_first_existing(final, c("Original_ID", "Genotype", "GEN", "ID"))
  decision_col <- si_first_existing(final, c("Decision", "Final_decision", "Status", "Recommendation"))
  index_col <- si_first_existing(final, c("Selection_Index", "Index", "Selection_score"))
  advantage_col <- si_first_existing(final, c("Index_Advantage", "SI_pct_index", "Advantage"))
  weakness_col <- si_first_existing(final, c("Weakness_trait", "Weakness", "Weak_traits"))
  reason_col <- si_first_existing(final, c("Decision_reason", "Reason", "Recommendation_reason"))
  rank_col <- si_first_existing(final, c("Final_rank", "Rank", "Selection_rank"))

  out <- data.frame(
    Genotype = if (!is.na(genotype_col)) as.character(final[[genotype_col]]) else as.character(seq_len(nrow(final))),
    Source = "LPSI",
    Final_action = if (!is.na(decision_col)) as.character(final[[decision_col]]) else "REVIEW",
    Rank = if (!is.na(rank_col)) suppressWarnings(as.numeric(final[[rank_col]])) else seq_len(nrow(final)),
    Index = if (!is.na(index_col)) suppressWarnings(as.numeric(final[[index_col]])) else NA_real_,
    Index_advantage = if (!is.na(advantage_col)) suppressWarnings(as.numeric(final[[advantage_col]])) else NA_real_,
    Weakness_trait = if (!is.na(weakness_col)) as.character(final[[weakness_col]]) else "",
    Reason = if (!is.na(reason_col)) as.character(final[[reason_col]]) else "Single-trial LPSI decision.",
    stringsAsFactors = FALSE
  )
  out[order(match(out$Final_action, c("ADVANCE", "RETEST", "DISCARD")), out$Rank, na.last = TRUE), , drop = FALSE]
}

si_met_action <- function(recommendation) {
  recommendation <- as.character(recommendation)
  ifelse(grepl("^Advance", recommendation, ignore.case = TRUE), "ADVANCE",
         ifelse(grepl("^Retest|review", recommendation, ignore.case = TRUE), "RETEST", "LOW_PRIORITY"))
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
  out[order(match(out$Final_action, c("ADVANCE", "RETEST", "LOW_PRIORITY")), out$Rank, na.last = TRUE), , drop = FALSE]
}

si_breeder_recommendation_table <- function(lpsi_results = NULL, met_results = NULL, top_n = 20) {
  lpsi <- si_lpsi_recommendation_table(lpsi_results)
  met <- si_met_breeder_recommendation(met_results, top_n = top_n)

  if (nrow(lpsi) == 0 && nrow(met) == 0) return(data.frame())
  if (nrow(lpsi) == 0) return(met)
  if (nrow(met) == 0) return(lpsi)

  lpsi_small <- data.frame(
    Genotype = lpsi$Genotype,
    LPSI_action = lpsi$Final_action,
    LPSI_rank = lpsi$Rank,
    LPSI_index = lpsi$Index,
    LPSI_advantage = lpsi$Index_advantage,
    Weakness_trait = lpsi$Weakness_trait,
    LPSI_reason = lpsi$Reason,
    stringsAsFactors = FALSE
  )
  met_small <- data.frame(
    Genotype = met$Genotype,
    MET_action = met$Final_action,
    MET_rank = met$Rank,
    MET_index = met$Index,
    MET_reason = met$Reason,
    stringsAsFactors = FALSE
  )

  all_genotypes <- sort(unique(c(lpsi_small$Genotype, met_small$Genotype)))
  out <- merge(data.frame(Genotype = all_genotypes, stringsAsFactors = FALSE), lpsi_small,
               by = "Genotype", all.x = TRUE)
  out <- merge(out, met_small, by = "Genotype", all.x = TRUE)
  out$LPSI_action[is.na(out$LPSI_action)] <- ""
  out$MET_action[is.na(out$MET_action)] <- ""
  out$Weakness_trait[is.na(out$Weakness_trait)] <- ""
  out$LPSI_reason[is.na(out$LPSI_reason)] <- ""
  out$MET_reason[is.na(out$MET_reason)] <- ""
  out$Support_count <- as.integer(out$LPSI_action == "ADVANCE") + as.integer(out$MET_action == "ADVANCE")
  out$Final_action <- ifelse(
    out$LPSI_action == "ADVANCE" & out$MET_action == "ADVANCE", "ADVANCE",
    ifelse(out$LPSI_action == "ADVANCE" | out$MET_action == "ADVANCE" |
             out$LPSI_action == "RETEST" | out$MET_action == "RETEST", "RETEST",
           "DISCARD")
  )
  out$Source <- "LPSI + MET"
  out$Reason <- ifelse(
    out$Support_count == 2,
    "Advance: supported by both single-trial trait index and MET stability/performance.",
    ifelse(out$LPSI_action == "ADVANCE",
           "Retest: strong LPSI candidate; confirm across-environment performance.",
           ifelse(out$MET_action == "ADVANCE",
                  "Retest: strong MET candidate; check single-trial trait balance and weaknesses.",
                  ifelse(out$Final_action == "RETEST",
                         "Retest/review: partial support from LPSI or MET.",
                         "Discard/lower priority: not supported by current LPSI or MET recommendation.")))
  )

  out <- out[, c("Genotype", "Source", "Final_action", "Support_count",
                 "LPSI_rank", "MET_rank", "LPSI_index", "MET_index",
                 "LPSI_advantage", "Weakness_trait", "LPSI_action",
                 "MET_action", "Reason", "LPSI_reason", "MET_reason")]
  out[order(match(out$Final_action, c("ADVANCE", "RETEST", "DISCARD")),
            -out$Support_count, out$LPSI_rank, out$MET_rank, na.last = TRUE), , drop = FALSE]
}
