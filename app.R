# Clear environment
rm(list = ls())
# Load packages
library(shiny)
library(bslib)
library(readxl)
library(tidyverse)
library(emmeans)
library(multcomp)
library(multcompView)
library(pheatmap)
library(ggplot2)
library(DT)
library(writexl)
library(grid)
library(lme4)
library(lmerTest)
library(patchwork)
options(shiny.maxRequestSize = 50 * 1024^2)

# Shared infrastructure used across the complete analysis pipeline.
source(file.path("modules", "shared_data_validation.R"), local = TRUE)
source(file.path("modules", "advanced_analysis_extensions.R"), local = TRUE)

# Numbered analysis modules. Keep this list as the pipeline backbone.
source(file.path("modules", "module_1_breeding.R"), local = TRUE)
source(file.path("modules", "module_2_genetic_diversity.R"), local = TRUE)
source(file.path("modules", "module_3_mating.R"), local = TRUE)
source(file.path("modules", "module_4_selection_index.R"), local = TRUE)
source(file.path("modules", "module_5_met.R"), local = TRUE)


# User settings
id_col  <- "Variety"
rep_col <- "Rep"
remove_cols <- c("Sex", "FsC", "Pr", "Gs", "DM")
weight_row_labels <- c("WEIGHT", "WEIGHTS", "IMPORTANCE", "IMPORTANT")
direction_row_labels <- c("DIRECTION", "DIRECTIONS", "TRAIT_DIRECTION")
type_col_candidates <- c("Type", "TYPE", "Entry_Type", "Entry type", "EntryType")
check_type_labels <- c(
  "CHECK", "CHECKS", "CONTROL", "CONTROLS", "BENCHMARK",
  "COMMERCIAL", "COMMERCIAL CHECK", "STANDARD CHECK"
)
priority_weight_cutoff <- 4
advance_index_cutoff <- 0
retest_index_cutoff  <- 0
priority_advance_cutoff_pct <- 0
priority_severe_weak_pct <- -10
lpsi_selection_intensity <- 0.10
run_simple_anova <- TRUE
run_lsd_test     <- TRUE
lsd_significance_alpha <- 0.05



# Run one mating-design analysis using columns selected in the app.
run_mating_pipeline <- function(
    df,
    design,
    trait_col,
    replication_col,
    parent1_col = NULL,
    parent2_col = NULL,
    line_col = NULL,
    tester_col = NULL,
    type_col = NULL) {
  design <- match.arg(
    design,
    c(
      "griffing_m1", "griffing_m2", "griffing_m3",
      "griffing_m4", "diallel_partial", "line_tester"
    )
  )
  
  required_cols <- if (design == "line_tester") {
    c(line_col, tester_col, replication_col, trait_col, type_col)
  } else {
    c(parent1_col, parent2_col, replication_col, trait_col)
  }
  missing_cols <- required_cols[
    is.na(required_cols) | required_cols == "" | !required_cols %in% names(df)
  ]
  if (length(missing_cols) > 0) {
    stop(
      "Select valid columns for every mating-analysis field. Missing: ",
      paste(unique(missing_cols), collapse = ", ")
    )
  }
  
  if (design == "line_tester") {
    result <- line_tester_manual(
      df,
      line_col = line_col,
      tester_col = tester_col,
      rep_col = replication_col,
      trait_col = trait_col,
      type_col = type_col
    )
    if (!is.null(result$error)) {
      stop(result$error)
    }
    return(result)
  }
  
  if (design == "diallel_partial") {
    return(diallel_partial_manual(
      df,
      trait = trait_col,
      p1 = parent1_col,
      p2 = parent2_col,
      rep = replication_col
    ))
  }
  
  if (design == "griffing_m2") {
    df$Parent1_std <- pmin(
      as.character(df[[parent1_col]]),
      as.character(df[[parent2_col]])
    )
    df$Parent2_std <- pmax(
      as.character(df[[parent1_col]]),
      as.character(df[[parent2_col]])
    )
    return(griffing_method2(
      df,
      rep_col = replication_col,
      male_col = "Parent1_std",
      female_col = "Parent2_std",
      trait_col = trait_col
    ))
  }
  
  mating_function <- switch(
    design,
    griffing_m1 = griffing_method1,
    griffing_m3 = griffing_method3,
    griffing_m4 = griffing_method4
  )
  mating_function(
    df,
    rep_col = replication_col,
    male_col = parent1_col,
    female_col = parent2_col,
    trait_col = trait_col
  )
}

get_mating_result_table <- function(results, table_type) {
  switch(
    table_type,
    anova = {
      table <- if (!is.null(results$anova_full)) results$anova_full else results$anova
      if (!is.null(table) && !"Source" %in% names(table)) {
        table <- tibble::rownames_to_column(as.data.frame(table), "Source")
      }
      table
    },
    gca = {
      if (!is.null(results$gca)) {
        results$gca
      } else {
        gca_tables <- list()
        if (!is.null(results$gca_lines)) {
          gca_tables$Line <- results$gca_lines
        }
        if (!is.null(results$gca_testers)) {
          gca_tables$Tester <- results$gca_testers
        }
        if (length(gca_tables) == 0) NULL else
          dplyr::bind_rows(gca_tables, .id = "Combiner")
      }
    },
    sca = results$sca,
    variance = {
      if (!is.null(results$var)) results$var else results$variance_components
    },
    classical_anova = results$griffing_anova,
    NULL
  )
}




build_export_tables <- function(analysis_type, results) {
  if (is.null(results)) {
    return(list())
  }
  
  output_list <- list()
  make_sheet_name <- function(prefix, trait, existing_names) {
    base <- gsub("[\\[\\]\\:\\*\\?/\\\\]", "_", paste(prefix, trait, sep = "_"))
    base <- substr(base, 1, 31)
    candidate <- base
    counter <- 2
    while (candidate %in% existing_names) {
      suffix <- paste0("_", counter)
      candidate <- paste0(substr(base, 1, 31 - nchar(suffix)), suffix)
      counter <- counter + 1
    }
    candidate
  }
  add_sheet <- function(prefix, trait, table) {
    if (is.null(table)) {
      return(invisible(NULL))
    }
    sheet_name <- make_sheet_name(prefix, trait, names(output_list))
    output_list[[sheet_name]] <<- as.data.frame(table)
  }
  
  if (analysis_type == "MATING") {
    add_sheet("01", "anova", get_mating_result_table(results, "anova"))
    add_sheet("02", "gca", get_mating_result_table(results, "gca"))
    add_sheet("03", "sca", get_mating_result_table(results, "sca"))
    add_sheet("04", "variance", get_mating_result_table(results, "variance"))
  } else if (analysis_type == "BREEDING") {
    add_sheet("00", "settings", results$settings)
    add_sheet("01", "genetic_stats", results$genetic_stats)
    add_sheet("02", "response_year", results$response_per_year)
    if (!is.null(results$realized_gain) && nrow(results$realized_gain) > 0) {
      add_sheet("03", "realized_gain", results$realized_gain)
    }
    if (!is.null(results$generation_stats) && nrow(results$generation_stats) > 0) {
      add_sheet("04", "generation", results$generation_stats)
    }
  } else if (analysis_type == "LPSI") {
    add_sheet("00", "decision_settings", results$decision_settings)
    add_sheet("00b", "breeder_recommendation", si_breeder_recommendation_table(lpsi_results = results))
    add_sheet("01", "summary", results$trait_info)
    add_sheet("02", "anova", results$anova_full)
    add_sheet("03", "mean_comparison", results$lsd_wide)
    add_sheet("04", "superiority_mean", results$superiority_index)
    add_sheet("05", "selection_index", results$index_ranking)
    add_sheet("06", "decision", results$final_decision)
    add_sheet("07", "heritability_gain", results$heritability_gain)
    add_sheet("08", "pipeline_review", si_lpsi_pipeline_review(results))
    default_trait <- tryCatch(as.character(results$trait_info$Trait[1]), error = function(e) NULL)
    add_sheet("09", "direct_selection", si_lpsi_direct_selection(results, default_trait, 15))
    add_sheet("10", "method_compare", si_lpsi_method_comparison(results, default_trait, 15))
  } else if (analysis_type == "MET") {
    add_sheet("00", "settings", results$settings)
    for (trait in results$met_trait_names) {
      result <- results$met_by_trait[[trait]]
      if (!is.null(result)) {
        add_sheet("01_summary", trait, result$genotype_summary)
        add_sheet("01b_qc", trait, build_met_qc_table(result))
        add_sheet("02_lmm", trait, result$model_summary)
        add_sheet("03_variance", trait, result$variance_components)
        add_sheet("04_blup", trait, result$blups_main)
        add_sheet("05_fw", trait, result$fw_results)
        add_sheet("06_ammi", trait, result$ammi_genotype)
        add_sheet("07_gge", trait, result$gge_genotype)
        add_sheet("08_selection", trait, result$met_selection)
        add_sheet("08b_ammi_gge_notes", trait, result$ammi_notes)
      }
    }
    add_sheet("09_overall", "selection", results$met_integrated_ranking)
    add_sheet("10", "pipeline_review", si_met_pipeline_review(results))
    add_sheet("11", "trait_quality", si_met_trait_quality(results))
  } else if (analysis_type == "DIVERSITY") {
    add_sheet("00", "benchmark_checks", data.frame(
      Benchmark_check = results$selected_checks %||% character(0),
      stringsAsFactors = FALSE
    ))
    add_sheet("01", "genotype_values", results$genotype_values)
    add_sheet("02", "clusters", results$clusters)
    add_sheet("03", "superiority", results$superiority)
    add_sheet("04", "correlation", tibble::rownames_to_column(as.data.frame(results$correlation), "Trait"))
  }
  
  output_list
}



# Helper functions
to_number <- function(x) {
  suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
}
clean_text <- function(x) {
  trimws(as.character(x))
}
natural_id_order <- function(x) {
  values <- unique(clean_text(x))
  values <- values[!is.na(values) & values != ""]
  if (length(values) == 0) {
    return(values)
  }
  numeric_key <- suppressWarnings(as.numeric(stringr::str_extract(values, "\\d+")))
  prefix_key <- stringr::str_replace(values, "\\d.*$", "")
  values[order(prefix_key, is.na(numeric_key), numeric_key, values)]
}
normalize_type_label <- function(x) {
  toupper(gsub("\\s+", " ", trimws(as.character(x))))
}
backtick_name <- function(x) {
  paste0("`", gsub("`", "``", x), "`")
}
sig_label <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    p < 0.10  ~ ".",
    TRUE ~ "ns"
  )
}
parse_trait_direction <- function(x) {
  value <- tolower(trimws(as.character(x)))
  if (is.na(value) || value == "" || value == "na") {
    return(list(direction = "Higher better", target_value = NA_real_))
  }
  value <- gsub("\\s+", " ", value)
  if (value %in% c("high", "higher", "higher better", "higher_better",
                   "more", "max", "maximize", "+", "positive")) {
    return(list(direction = "Higher better", target_value = NA_real_))
  }
  if (value %in% c("low", "lower", "lower better", "lower_better",
                   "less", "min", "minimize", "-", "negative")) {
    return(list(direction = "Lower better", target_value = NA_real_))
  }
  if (grepl("target", value)) {
    target_value <- stringr::str_extract(value, "-?\\d+\\.?\\d*")
    target_value <- to_number(target_value)
    if (!is.na(target_value)) {
      return(list(direction = "Target trait", target_value = target_value))
    }
  }
  return(list(direction = "Higher better", target_value = NA_real_))
}
standardize_trait <- function(x) {
  if (all(is.na(x))) {
    return(rep(0, length(x)))
  }
  if (sd(x, na.rm = TRUE) == 0) {
    return(rep(0, length(x)))
  }
  as.numeric(scale(x))
}
make_model_formula <- function(trait, data) {
  if (n_distinct(data$Rep[!is.na(data$Rep)]) > 1) {
    as.formula(paste0(backtick_name(trait), " ~ ID + Rep"))
  } else {
    as.formula(paste0(backtick_name(trait), " ~ ID"))
  }
}
# Prepare uploaded Excel data
prepare_excel_input <- function(df_raw) {
  names(df_raw) <- make.unique(trimws(names(df_raw)), sep = "_")
  if (!id_col %in% names(df_raw)) {
    stop("The ID column is missing. Expected column name: ", id_col)
  }
  if (!rep_col %in% names(df_raw)) {
    stop("The replication column is missing. Expected column name: ", rep_col)
  }
  notes <- c()
  id_values_upper <- toupper(clean_text(df_raw[[id_col]]))
  weight_rows <- which(id_values_upper %in% weight_row_labels)
  direction_rows <- which(id_values_upper %in% direction_row_labels)
  weight_row <- if (length(weight_rows) > 0) tail(weight_rows, 1) else integer(0)
  direction_row <- if (length(direction_rows) > 0) tail(direction_rows, 1) else integer(0)
  metadata_rows <- sort(unique(c(weight_rows, direction_rows)))
  df_data <- df_raw
  if (length(metadata_rows) > 0) {
    df_data <- df_data[-metadata_rows, , drop = FALSE]
  }
  empty_rows <- apply(
    df_data,
    1,
    function(z) all(is.na(z) | trimws(as.character(z)) == "")
  )
  df_data <- df_data[!empty_rows, , drop = FALSE]
  if (nrow(df_data) == 0) {
    stop("No real data rows found after removing DIRECTION/WEIGHT rows.")
  }
  df_data <- df_data %>%
    filter(!is.na(.data[[id_col]]) & trimws(as.character(.data[[id_col]])) != "")
  if (nrow(df_data) == 0) {
    stop("No rows with valid Variety/ID were found.")
  }
  type_col <- type_col_candidates[type_col_candidates %in% names(df_data)][1]
  if (length(type_col) == 0 || is.na(type_col)) {
    type_col <- NULL
  }
  detected_check_varieties <- character(0)
  if (!is.null(type_col)) {
    type_values <- normalize_type_label(df_data[[type_col]])
    detected_check_varieties <- unique(clean_text(df_data[[id_col]][type_values %in% check_type_labels]))
    detected_check_varieties <- detected_check_varieties[
      !is.na(detected_check_varieties) & detected_check_varieties != ""
    ]
  }
  check_original_name <- as.character(df_data[[id_col]][1])
  candidate_cols <- setdiff(
    names(df_data),
    c(id_col, rep_col, type_col, remove_cols)
  )
  trait_numeric_count <- sapply(candidate_cols, function(tr) {
    sum(!is.na(to_number(df_data[[tr]])))
  })
  trait_cols <- candidate_cols[trait_numeric_count > 0]
  if (length(trait_cols) == 0) {
    stop("No usable numeric trait columns found.")
  }
  weights_raw_used <- setNames(rep(1, length(trait_cols)), trait_cols)
  if (length(weight_row) > 0) {
    weight_values <- to_number(unlist(
      df_raw[weight_row, trait_cols, drop = FALSE],
      use.names = FALSE
    ))
    names(weight_values) <- trait_cols
    weights_raw_used <- weight_values
    missing_weight_traits <- names(weights_raw_used)[is.na(weights_raw_used)]
    if (length(missing_weight_traits) > 0) {
      notes <- c(
        notes,
        paste0(
          "Missing weight was detected for: ",
          paste(missing_weight_traits, collapse = ", "),
          ". These traits were given weight = 1."
        )
      )
      weights_raw_used[missing_weight_traits] <- 1
    }
    weights_raw_used[weights_raw_used < 0] <- 0
    if (sum(weights_raw_used, na.rm = TRUE) <= 0) {
      weights_raw_used <- setNames(rep(1, length(trait_cols)), trait_cols)
      notes <- c(
        notes,
        "All trait weights were zero or invalid, so all traits were given weight = 1."
      )
    }
  } else {
    notes <- c(
      notes,
      "No WEIGHT row was found. All detected traits were given weight = 1."
    )
  }
  trait_direction <- setNames(rep("Higher better", length(trait_cols)), trait_cols)
  target_traits <- numeric(0)
  if (length(direction_row) > 0) {
    direction_values <- unlist(
      df_raw[direction_row, trait_cols, drop = FALSE],
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
  } else {
    notes <- c(
      notes,
      "No DIRECTION row was found. All detected traits were treated as higher-better."
    )
  }
  direction_table <- data.frame(
    Trait = trait_cols,
    Direction = as.character(trait_direction[trait_cols]),
    Target_value = ifelse(
      trait_cols %in% names(target_traits),
      as.numeric(target_traits[trait_cols]),
      NA_real_
    ),
    Raw_weight = as.numeric(weights_raw_used[trait_cols])
  )
  list(
    data = df_data,
    trait_cols = trait_cols,
    weights_raw_used = weights_raw_used,
    trait_direction = trait_direction,
    target_traits = target_traits,
    direction_table = direction_table,
    check_original_name = check_original_name,
    type_col = type_col,
    detected_check_varieties = detected_check_varieties,
    notes = notes
  )
}
make_diversity_trait_info <- function(df_raw, genotype_col, trait_cols) {
  if (is.null(genotype_col) || !genotype_col %in% names(df_raw)) {
    return(data.frame(
      Trait = trait_cols,
      Direction = rep("Higher better", length(trait_cols)),
      Target_value = rep(NA_real_, length(trait_cols)),
      Raw_weight = rep(1, length(trait_cols)),
      stringsAsFactors = FALSE
    ))
  }
  id_values_upper <- toupper(clean_text(df_raw[[genotype_col]]))
  weight_rows <- which(id_values_upper %in% weight_row_labels)
  direction_rows <- which(id_values_upper %in% direction_row_labels)
  weight_row <- if (length(weight_rows) > 0) tail(weight_rows, 1) else integer(0)
  direction_row <- if (length(direction_rows) > 0) tail(direction_rows, 1) else integer(0)
  trait_cols <- trait_cols[trait_cols %in% names(df_raw)]
  trait_direction <- setNames(rep("Higher better", length(trait_cols)), trait_cols)
  target_values <- setNames(rep(NA_real_, length(trait_cols)), trait_cols)
  weights <- setNames(rep(1, length(trait_cols)), trait_cols)
  if (length(direction_row) > 0) {
    direction_values <- unlist(
      df_raw[direction_row, trait_cols, drop = FALSE],
      use.names = FALSE
    )
    names(direction_values) <- trait_cols
    for (tr in trait_cols) {
      parsed <- parse_trait_direction(direction_values[[tr]])
      trait_direction[[tr]] <- parsed$direction
      target_values[[tr]] <- parsed$target_value
    }
  }
  if (length(weight_row) > 0) {
    weight_values <- to_number(unlist(
      df_raw[weight_row, trait_cols, drop = FALSE],
      use.names = FALSE
    ))
    names(weight_values) <- trait_cols
    weights[!is.na(weight_values)] <- weight_values[!is.na(weight_values)]
  }
  data.frame(
    Trait = trait_cols,
    Direction = as.character(trait_direction[trait_cols]),
    Target_value = as.numeric(target_values[trait_cols]),
    Raw_weight = as.numeric(weights[trait_cols]),
    stringsAsFactors = FALSE
  )
}
# Diagnostic data preparation
make_diagnostic_data <- function(df_raw) {
  prepared <- prepare_excel_input(df_raw)
  df <- prepared$data %>%
    dplyr::select(all_of(c(id_col, rep_col, prepared$trait_cols))) %>%
    mutate(
      across(all_of(c(id_col, rep_col)), as.character),
      across(all_of(prepared$trait_cols), to_number)
    ) %>%
    rename(
      ID = all_of(id_col),
      Rep = all_of(rep_col)
    ) %>%
    filter(!is.na(ID) & trimws(ID) != "") %>%
    filter(!is.na(Rep) & trimws(Rep) != "")
  df$ID <- as.factor(df$ID)
  df$Rep <- as.factor(df$Rep)
  list(
    data = df,
    trait_cols = prepared$trait_cols,
    direction_table = prepared$direction_table,
    notes = prepared$notes
  )
}
fit_diagnostic_model <- function(df, trait) {
  model_data <- df %>%
    filter(!is.na(.data[[trait]]))
  if (nrow(model_data) < 3) {
    return(NULL)
  }
  rhs <- c()
  if (n_distinct(model_data$ID) > 1) {
    rhs <- c(rhs, "ID")
  }
  if (n_distinct(model_data$Rep) > 1) {
    rhs <- c(rhs, "Rep")
  }
  if (length(rhs) == 0) {
    form <- as.formula(paste0(backtick_name(trait), " ~ 1"))
  } else {
    form <- as.formula(paste0(backtick_name(trait), " ~ ", paste(rhs, collapse = " + ")))
  }
  lm(form, data = model_data)
}
make_shapiro_table <- function(diag) {
  map_dfr(diag$trait_cols, function(tr) {
    model <- tryCatch({
      fit_diagnostic_model(diag$data, tr)
    }, error = function(e) {
      NULL
    })
    if (is.null(model)) {
      return(data.frame(
        Trait = tr,
        N_residual = NA,
        Shapiro_W = NA,
        p_value = NA,
        Normality_note = "No",
        stringsAsFactors = FALSE
      ))
    }
    res <- residuals(model)
    res <- res[is.finite(res)]
    if (length(res) < 3) {
      return(data.frame(
        Trait = tr,
        N_residual = length(res),
        Shapiro_W = NA,
        p_value = NA,
        Normality_note = "No",
        stringsAsFactors = FALSE
      ))
    }
    if (length(res) > 5000) {
      return(data.frame(
        Trait = tr,
        N_residual = length(res),
        Shapiro_W = NA,
        p_value = NA,
        Normality_note = "No",
        stringsAsFactors = FALSE
      ))
    }
    if (length(unique(res)) < 2) {
      return(data.frame(
        Trait = tr,
        N_residual = length(res),
        Shapiro_W = NA,
        p_value = NA,
        Normality_note = "No",
        stringsAsFactors = FALSE
      ))
    }
    st <- tryCatch({
      shapiro.test(res)
    }, error = function(e) {
      NULL
    })
    if (is.null(st)) {
      return(data.frame(
        Trait = tr,
        N_residual = length(res),
        Shapiro_W = NA,
        p_value = NA,
        Normality_note = "No",
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      Trait = tr,
      N_residual = length(res),
      Shapiro_W = round(as.numeric(st$statistic), 4),
      p_value = round(st$p.value, 5),
      Normality_note = ifelse(st$p.value >= 0.05, "Yes", "No"),
      stringsAsFactors = FALSE
    )
  })
}
get_lpsi_normality_decision <- function(df, trait) {
  model <- tryCatch({
    fit_diagnostic_model(df, trait)
  }, error = function(e) {
    NULL
  })
  if (is.null(model)) {
    return(list(
      use_anova = TRUE,
      p_value = NA_real_,
      note = "Model could not be fitted; ANOVA used by default"
    ))
  }
  res <- residuals(model)
  res <- res[is.finite(res)]
  if (length(res) < 3) {
    return(list(
      use_anova = TRUE,
      p_value = NA_real_,
      note = "Too few residuals for Shapiro-Wilk test; ANOVA used by default"
    ))
  }
  if (length(res) > 5000) {
    return(list(
      use_anova = TRUE,
      p_value = NA_real_,
      note = "More than 5000 residuals; ANOVA used by default"
    ))
  }
  if (length(unique(res)) < 2) {
    return(list(
      use_anova = TRUE,
      p_value = NA_real_,
      note = "Residuals are constant; ANOVA used by default"
    ))
  }
  st <- tryCatch({
    shapiro.test(res)
  }, error = function(e) {
    NULL
  })
  if (is.null(st)) {
    return(list(
      use_anova = TRUE,
      p_value = NA_real_,
      note = "Shapiro-Wilk test failed; ANOVA used by default"
    ))
  }
  list(
    use_anova = st$p.value >= 0.05,
    p_value = as.numeric(st$p.value),
    note = ifelse(
      st$p.value >= 0.05,
      "Residuals passed Shapiro-Wilk; ANOVA used",
      "Residuals failed Shapiro-Wilk; Kruskal-Wallis used"
    )
  )
}
empty_plot <- function(message_text) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = message_text, size = 5) +
    theme_void()
}

format_superiority_pct_label <- function(x) {
  ifelse(is.na(x), "", paste0(sprintf("%.1f", x), "%"))
}

diversity_cluster_lookup <- function(result) {
  clusters <- as.data.frame(result$clusters)
  if (nrow(clusters) == 0 || !"GEN" %in% names(clusters) || !"Cluster" %in% names(clusters)) {
    return(data.frame(GEN = character(0), Cluster = factor()))
  }
  clusters %>%
    mutate(
      GEN = as.character(GEN),
      Cluster = factor(Cluster)
    )
}

diversity_auto_subcluster_k <- function(hc, main_k) {
  n <- length(hc$labels)
  main_k <- suppressWarnings(as.integer(main_k))
  if (length(main_k) == 0 || is.na(main_k) || main_k < 2) main_k <- 2
  if (n <= main_k) return(n)
  min_k <- min(n, main_k + 1)
  max_k <- min(n, max(main_k * 3, main_k + 3, 6))
  heights <- as.numeric(hc$height)
  candidates <- data.frame(k = min_k:max_k)
  candidates$i <- n - candidates$k
  candidates <- candidates[candidates$i >= 1 & candidates$i < length(heights), , drop = FALSE]
  if (nrow(candidates) == 0) return(min(n, max(main_k + 2, 2)))
  candidates$gap <- heights[candidates$i + 1] - heights[candidates$i]
  candidates$gap[!is.finite(candidates$gap)] <- -Inf
  candidates$k[which.max(candidates$gap)]
}

plot_diversity_dendrogram <- function(result, subcluster_k = NULL) {
  if (is.null(result$hclust) || length(result$hclust$labels) < 2) {
    return(empty_plot("Run Genetic Diversity to view dendrogram."))
  }
  hc <- result$hclust
  dend <- as.dendrogram(hc)
  leaf_order <- hc$labels[hc$order]
  x_pos <- setNames(seq_along(leaf_order), leaf_order)
  main_lookup <- diversity_cluster_lookup(result)
  n_leaves <- length(leaf_order)
  main_k <- if (nrow(main_lookup) > 0) length(unique(main_lookup$Cluster)) else 2
  subcluster_k <- suppressWarnings(as.integer(subcluster_k))
  if (length(subcluster_k) == 0 || is.na(subcluster_k) || subcluster_k < 2) {
    subcluster_k <- diversity_auto_subcluster_k(hc, main_k)
  }
  subcluster_k <- max(2, min(subcluster_k, n_leaves))
  subclusters <- stats::cutree(hc, k = subcluster_k)
  sub_lookup <- data.frame(
    GEN = names(subclusters),
    Subcluster = paste0("S", as.integer(subclusters)),
    stringsAsFactors = FALSE
  )
  leaf_subcluster <- setNames(sub_lookup$Subcluster, sub_lookup$GEN)
  branch_cluster <- function(labels) {
    groups <- unique(leaf_subcluster[labels])
    groups <- groups[!is.na(groups)]
    if (length(groups) == 1) groups else "Main split"
  }
  segments <- data.frame(x = numeric(0), y = numeric(0), xend = numeric(0), yend = numeric(0), Branch_group = character(0))
  
  walk_dend <- function(node) {
    if (is.leaf(node)) {
      label <- attr(node, "label")
      return(list(x = unname(x_pos[[label]]), y = 0, labels = label))
    }
    child_info <- lapply(node, walk_dend)
    y <- attr(node, "height")
    for (child in child_info) {
      segments <<- rbind(
        segments,
        data.frame(x = child$x, y = child$y, xend = child$x, yend = y, Branch_group = branch_cluster(child$labels))
      )
    }
    child_x <- vapply(child_info, `[[`, numeric(1), "x")
    child_labels <- unlist(lapply(child_info, `[[`, "labels"), use.names = FALSE)
    segments <<- rbind(
      segments,
      data.frame(x = min(child_x), y = y, xend = max(child_x), yend = y, Branch_group = branch_cluster(child_labels))
    )
    list(x = mean(child_x), y = y, labels = child_labels)
  }
  walk_dend(dend)
  
  max_height <- max(segments$y, segments$yend, na.rm = TRUE)
  if (!is.finite(max_height) || max_height <= 0) {
    max_height <- 1
  }
  y_step <- if (max_height <= 100) {
    10
  } else if (max_height <= 1000) {
    100
  } else if (max_height <= 5000) {
    500
  } else {
    1000
  }
  y_top <- ceiling(max_height / y_step) * y_step
  label_y <- -max_height * 0.08
  label_df <- data.frame(GEN = leaf_order, x = seq_along(leaf_order), y = 0, Label_y = label_y, stringsAsFactors = FALSE) %>%
    left_join(main_lookup, by = "GEN") %>%
    left_join(sub_lookup, by = "GEN") %>%
    mutate(
      Cluster = ifelse(is.na(Cluster), "Unclustered", as.character(Cluster)),
      Subcluster = ifelse(is.na(Subcluster), "Unclustered", as.character(Subcluster))
    )
  subcluster_levels <- unique(label_df$Subcluster)
  subcluster_palette <- setNames(
    grDevices::hcl.colors(length(subcluster_levels), palette = "Dark 3"),
    subcluster_levels
  )
  branch_palette <- c(subcluster_palette, "Main split" = "black")
  cluster_bands <- label_df %>%
    group_by(Cluster) %>%
    summarise(
      xmin = min(x) - 0.5,
      xmax = max(x) + 0.5,
      .groups = "drop"
    ) %>%
    arrange(xmin) %>%
    mutate(Cluster_band = paste0("Cluster ", Cluster))
  
  ggplot() +
    geom_rect(
      data = cluster_bands,
      aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = max_height * 0.52, fill = Cluster_band),
      alpha = 0.18,
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = segments,
      aes(x = x, y = y, xend = xend, yend = yend, color = Branch_group),
      linewidth = 0.75,
      show.legend = FALSE
    ) +
    geom_text(
      data = label_df,
      aes(x = x, y = Label_y, label = GEN, color = Subcluster),
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 3,
      show.legend = FALSE
    ) +
    scale_color_manual(values = branch_palette) +
    scale_y_continuous(
      breaks = seq(0, y_top, by = y_step),
      expand = expansion(mult = c(0, 0.03))
    ) +
    labs(
      title = "Ward's linkage cluster dendrogram",
      subtitle = "Genotypes are colored by cluster membership",
      x = NULL,
      y = "Mahalanobis D2 distance"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(color = "gray35", size = 13),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = 12, r = 20, b = 128, l = 16)
    ) +
    coord_cartesian(ylim = c(0, y_top), clip = "off")
}

plot_diversity_correlation_heatmap <- function(result) {
  corr <- as.data.frame(as.table(result$correlation), stringsAsFactors = FALSE)
  names(corr) <- c("Trait1", "Trait2", "Correlation")
  if (nrow(corr) == 0) {
    return(empty_plot("Trait correlation is not available."))
  }
  traits <- colnames(result$correlation)
  corr <- corr %>%
    mutate(
      Trait1 = factor(Trait1, levels = traits),
      Trait2 = factor(Trait2, levels = rev(traits)),
      Label = sprintf("%.2f", Correlation)
    )
  ggplot(corr, aes(x = Trait1, y = Trait2, fill = Correlation)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = Label), size = 3.1, color = "gray10") +
    scale_fill_gradient2(low = "#1F77B4", mid = "#F7F7F7", high = "#E68613", midpoint = 0, limits = c(-1, 1)) +
    labs(
      title = "Pearson correlation heatmap",
      subtitle = "Trait association among genotype values",
      x = NULL,
      y = NULL,
      fill = "r"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(color = "gray35"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
}

build_diversity_superiority_data <- function(result, check_genotypes = NULL) {
  values <- as.data.frame(result$genotype_values)
  if (nrow(values) == 0 || !"GEN" %in% names(values)) {
    return(data.frame())
  }
  traits <- setdiff(names(values), "GEN")
  traits <- traits[vapply(values[traits], is.numeric, logical(1))]
  if (length(traits) == 0) {
    return(data.frame())
  }
  trait_info <- result$trait_info
  trait_direction <- setNames(rep("Higher better", length(traits)), traits)
  target_values <- setNames(rep(NA_real_, length(traits)), traits)
  if (!is.null(trait_info) && nrow(as.data.frame(trait_info)) > 0) {
    trait_info <- as.data.frame(trait_info)
    if (all(c("Trait", "Direction") %in% names(trait_info))) {
      matched <- trait_info$Trait %in% traits
      trait_direction[trait_info$Trait[matched]] <- as.character(trait_info$Direction[matched])
    }
    if (all(c("Trait", "Target_value") %in% names(trait_info))) {
      matched <- trait_info$Trait %in% traits
      target_values[trait_info$Trait[matched]] <- suppressWarnings(as.numeric(trait_info$Target_value[matched]))
    }
  }
  check_genotypes <- unique(clean_text(check_genotypes %||% result$selected_checks))
  check_genotypes <- check_genotypes[!is.na(check_genotypes) & check_genotypes != ""]
  check_rows <- values[clean_text(values$GEN) %in% check_genotypes, , drop = FALSE]
  benchmark_label <- "Population mean"
  if (nrow(check_rows) > 0) {
    benchmark_label <- if (length(check_genotypes) > 1) {
      paste0("Mean of selected checks (", paste(check_genotypes, collapse = ", "), ")")
    } else {
      paste0("Check ", check_genotypes[1])
    }
  } else {
    check_rows <- values
  }
  benchmarks <- vapply(traits, function(tr) mean(check_rows[[tr]], na.rm = TRUE), numeric(1))
  calc_superiority <- function(candidate, benchmark, trait_name) {
    denom <- ifelse(abs(benchmark) < 0.0001, 1, abs(benchmark))
    target <- target_values[[trait_name]]
    if (is.finite(target)) {
      candidate_distance <- abs(candidate - target)
      benchmark_distance <- abs(benchmark - target)
      target_denom <- ifelse(abs(benchmark_distance) < 0.0001, 1, abs(benchmark_distance))
      return((benchmark_distance - candidate_distance) / target_denom * 100)
    }
    if (identical(trait_direction[[trait_name]], "Lower better")) {
      return((benchmark - candidate) / denom * 100)
    }
    (candidate - benchmark) / denom * 100
  }
  superiority <- values[, c("GEN", traits), drop = FALSE]
  for (tr in traits) {
    superiority[[tr]] <- calc_superiority(
      candidate = superiority[[tr]],
      benchmark = benchmarks[[tr]],
      trait_name = tr
    )
  }
  clusters <- diversity_cluster_lookup(result)
  superiority %>%
    left_join(clusters, by = "GEN") %>%
    mutate(Cluster = ifelse(is.na(Cluster), "Unclustered", as.character(Cluster))) %>%
    pivot_longer(
      cols = all_of(traits),
      names_to = "Trait",
      values_to = "Superiority_pct"
    ) %>%
    group_by(GEN) %>%
    mutate(Mean_superiority_pct = mean(Superiority_pct, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      Benchmark = benchmark_label,
      Benchmark_check = ifelse(clean_text(GEN) %in% check_genotypes, "YES", ""),
      Direction = unname(trait_direction[as.character(Trait)]),
      Superiority_pct = round(Superiority_pct, 1),
      Mean_superiority_pct = round(Mean_superiority_pct, 1)
    )
}

plot_diversity_superiority_heatmap <- function(result, check_genotypes = NULL) {
  superiority <- build_diversity_superiority_data(result, check_genotypes)
  if (nrow(superiority) == 0) {
    return(empty_plot("GDA superiority needs genotypic values for at least one trait."))
  }
  plot_data <- superiority %>%
    filter(Benchmark_check != "YES", !is.na(Superiority_pct))
  if (nrow(plot_data) == 0) {
    return(empty_plot("GDA superiority needs at least one non-check genotype."))
  }
  max_abs_sup <- max(abs(plot_data$Superiority_pct), na.rm = TRUE)
  if (!is.finite(max_abs_sup) || max_abs_sup < 0.0001) {
    max_abs_sup <- 1
  }
  trait_order <- rev(unique(plot_data$Trait))
  genotype_order <- superiority %>%
    distinct(GEN, Mean_superiority_pct) %>%
    filter(GEN %in% plot_data$GEN) %>%
    arrange(Mean_superiority_pct, GEN) %>%
    pull(GEN)
  dend_segments <- NULL
  cluster_wide <- plot_data %>%
    dplyr::select(GEN, Trait, Superiority_pct) %>%
    tidyr::pivot_wider(names_from = Trait, values_from = Superiority_pct)
  if (nrow(cluster_wide) >= 2) {
    cluster_mat <- as.matrix(cluster_wide[, setdiff(names(cluster_wide), "GEN"), drop = FALSE])
    storage.mode(cluster_mat) <- "double"
    cluster_mat[!is.finite(cluster_mat)] <- 0
    rownames(cluster_mat) <- cluster_wide$GEN
    dend_result <- tryCatch({
      hc <- stats::hclust(stats::dist(cluster_mat), method = "ward.D2")
      leaf_order <- hc$labels[hc$order]
      x_pos <- setNames(seq_along(leaf_order), leaf_order)
      segments <- data.frame(x = numeric(0), y = numeric(0), xend = numeric(0), yend = numeric(0))
      walk_dend <- function(node) {
        if (is.leaf(node)) {
          label <- attr(node, "label")
          return(list(x = unname(x_pos[[label]]), y = 0, labels = label))
        }
        child_info <- lapply(node, walk_dend)
        y <- attr(node, "height")
        for (child in child_info) {
          segments <<- rbind(
            segments,
            data.frame(x = child$x, y = child$y, xend = child$x, yend = y)
          )
        }
        child_x <- vapply(child_info, `[[`, numeric(1), "x")
        segments <<- rbind(
          segments,
          data.frame(x = min(child_x), y = y, xend = max(child_x), yend = y)
        )
        list(x = mean(child_x), y = y)
      }
      walk_dend(as.dendrogram(hc))
      list(order = leaf_order, segments = segments)
    }, error = function(e) {
      NULL
    })
    if (!is.null(dend_result)) {
      genotype_order <- dend_result$order
      dend_segments <- dend_result$segments
    }
  }
  genotype_positions <- data.frame(
    GEN = genotype_order,
    Genotype_x = seq_along(genotype_order),
    stringsAsFactors = FALSE
  )
  trait_positions <- data.frame(
    Trait = trait_order,
    Trait_y = seq_along(trait_order),
    stringsAsFactors = FALSE
  )
  plot_data <- plot_data %>%
    left_join(genotype_positions, by = "GEN") %>%
    left_join(trait_positions, by = "Trait") %>%
    filter(!is.na(Genotype_x), !is.na(Trait_y)) %>%
    mutate(
      Label = format_superiority_pct_label(Superiority_pct)
    )
  dend_plot_data <- NULL
  y_top <- length(trait_order) + 0.5
  if (!is.null(dend_segments) && nrow(dend_segments) > 0) {
    dend_height <- max(dend_segments$y, dend_segments$yend, na.rm = TRUE)
    if (!is.finite(dend_height) || dend_height <= 0) {
      dend_height <- 1
    }
    dend_base_y <- length(trait_order) + 0.52
    dend_span <- max(0.65, length(trait_order) * 0.16)
    dend_plot_data <- dend_segments %>%
      mutate(
        y = dend_base_y + (y / dend_height) * dend_span,
        yend = dend_base_y + (yend / dend_height) * dend_span
      )
    y_top <- dend_base_y + dend_span * 1.05
  }
  ggplot(plot_data, aes(x = Genotype_x, y = Trait_y, fill = Superiority_pct)) +
    geom_tile(width = 0.95, height = 0.95, color = "white", linewidth = 0.25) +
    geom_text(aes(label = Label), size = 2.4, color = "gray10") +
    {
      if (!is.null(dend_plot_data)) {
        geom_segment(
          data = dend_plot_data,
          aes(x = x, y = y, xend = xend, yend = yend),
          inherit.aes = FALSE,
          linewidth = 0.45,
          color = "gray25"
        )
      }
    } +
    scale_x_continuous(
      breaks = genotype_positions$Genotype_x,
      labels = genotype_positions$GEN,
      limits = c(0.5, nrow(genotype_positions) + 0.5),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = trait_positions$Trait_y,
      labels = trait_positions$Trait,
      limits = c(0.5, y_top),
      expand = c(0, 0)
    ) +
    scale_fill_gradient2(
      low = "#D85A30",
      mid = "white",
      high = "#1D9E75",
      midpoint = 0,
      limits = c(-max_abs_sup, max_abs_sup),
      na.value = "gray90"
    ) +
    labs(
      title = "GDA trait superiority",
      subtitle = paste0("Percent advantage versus ", unique(superiority$Benchmark)[1], "; positive values are better"),
      x = "Genotype",
      y = NULL,
      fill = "Superiority (%)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(color = "gray35"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = 9),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right"
    )
}

plot_genetic_gain_curve <- function(heritability_gain, trait_info, trait_name = NULL, selection_pct = NULL) {
  if (is.null(heritability_gain) || nrow(heritability_gain) == 0) {
    return(empty_plot("Run LPSI analysis to view genetic gain."))
  }
  valid_gain <- heritability_gain %>%
    filter(
      is.finite(Mean),
      is.finite(Phenotypic_variance),
      Phenotypic_variance > 0,
      is.finite(Genetic_advance)
    )
  if (nrow(valid_gain) == 0) {
    return(empty_plot("Genetic gain curve needs valid phenotypic variance and genetic advance."))
  }
  if (is.null(trait_name) || !trait_name %in% valid_gain$Trait) {
    trait_name <- valid_gain$Trait[1]
  }
  
  gain_row <- valid_gain %>% filter(Trait == trait_name) %>% slice(1)
  direction <- trait_info %>%
    filter(Trait == trait_name) %>%
    pull(Direction)
  target_value <- trait_info %>%
    filter(Trait == trait_name) %>%
    pull(Target_value)
  target_value <- if (length(target_value) == 0) NA_real_ else as.numeric(target_value[1])
  direction <- if (length(direction) == 0 || is.na(direction[1])) {
    "Higher better"
  } else {
    direction[1]
  }
  
  original_mean <- as.numeric(gain_row$Mean)
  sd_p <- sqrt(as.numeric(gain_row$Phenotypic_variance))
  h2 <- as.numeric(gain_row$Broad_sense_H2)
  select_pct <- suppressWarnings(as.numeric(selection_pct))
  if (length(select_pct) == 0 || !is.finite(select_pct[1])) {
    select_pct <- as.numeric(gain_row$Selection_intensity_pct)
  } else {
    select_pct <- select_pct[1]
  }
  select_prop <- select_pct / 100
  if (!is.finite(select_prop) || select_prop <= 0 || select_prop >= 1) {
    select_prop <- lpsi_selection_intensity
    select_pct <- select_prop * 100
  }
  selection_k <- dnorm(qnorm(1 - select_prop)) / select_prop
  ga <- abs(selection_k * sd_p * h2)
  sign_direction <- if (direction == "Target trait" && is.finite(target_value)) {
    target_delta <- target_value - original_mean
    ifelse(abs(target_delta) < 0.0001, 1, sign(target_delta))
  } else if (direction == "Lower better") {
    -1
  } else {
    1
  }
  selected_mean <- original_mean + sign_direction * ga
  threshold <- original_mean + sign_direction * qnorm(1 - select_prop) * sd_p
  
  x_min <- min(original_mean, selected_mean, threshold) - 3.6 * sd_p
  x_max <- max(original_mean, selected_mean, threshold) + 3.6 * sd_p
  curve_x <- seq(x_min, x_max, length.out = 500)
  curves <- bind_rows(
    data.frame(
      x = curve_x,
      density = dnorm(curve_x, original_mean, sd_p),
      Population = "Original population"
    ),
    data.frame(
      x = curve_x,
      density = dnorm(curve_x, selected_mean, sd_p),
      Population = "Expected selected population"
    )
  )
  selected_tail <- data.frame(
    x = curve_x,
    density = dnorm(curve_x, original_mean, sd_p)
  ) %>%
    filter(if (sign_direction > 0) x >= threshold else x <= threshold)
  
  y_max <- max(curves$density, na.rm = TRUE)
  original_peak_y <- dnorm(original_mean, original_mean, sd_p)
  selected_peak_y <- dnorm(selected_mean, selected_mean, sd_p)
  threshold_y <- dnorm(threshold, original_mean, sd_p)
  original_label_y <- original_peak_y * 1.05
  selected_label_y <- selected_peak_y * 1.05
  gain_y <- -y_max * 0.12
  gain_label_y <- gain_y * 1.80
  intensity_y <- gain_y * 0.65
  gain_label_x <- mean(c(original_mean, selected_mean))
  intensity_label_x <- threshold + sign_direction * sd_p * 0.16
  intensity_label_hjust <- ifelse(sign_direction > 0, 0, 1)
  
  ggplot(curves, aes(x = x, y = density, fill = Population, color = Population)) +
    geom_area(alpha = 0.50, position = "identity", linewidth = 0.4) +
    geom_line(linewidth = 0.9) +
    geom_area(
      data = selected_tail,
      aes(x = x, y = density),
      inherit.aes = FALSE,
      fill = "#F28E2B",
      alpha = 0.45
    ) +
    geom_hline(yintercept = 0, color = "gray10", linewidth = 0.45) +
    annotate(
      "segment",
      x = original_mean,
      xend = original_mean,
      y = 0,
      yend = original_peak_y,
      linetype = "longdash",
      color = "gray35",
      linewidth = 0.8
    ) +
    annotate(
      "segment",
      x = selected_mean,
      xend = selected_mean,
      y = 0,
      yend = selected_peak_y,
      linetype = "longdash",
      color = "gray35",
      linewidth = 0.8
    ) +
    annotate(
      "segment",
      x = threshold,
      xend = threshold,
      y = 0,
      yend = threshold_y,
      linetype = "dotted",
      color = "#C95F18",
      linewidth = 0.9
    ) +
    annotate(
      "segment",
      x = original_mean,
      xend = selected_mean,
      y = gain_y,
      yend = gain_y,
      arrow = arrow(ends = "both", length = unit(0.16, "cm")),
      color = "gray30",
      linewidth = 0.8
    ) +
    annotate(
      "text",
      x = gain_label_x,
      y = gain_label_y,
      label = paste0("Expected genetic gain = ", round(ga, 3)),
      color = "gray20",
      size = 3.7,
      hjust = 0.5
    ) +
    annotate(
      "text",
      x = original_mean,
      y = original_label_y,
      label = "Average\nperformance\n(original)",
      color = "gray15",
      size = 3.8,
      lineheight = 0.95,
      hjust = 0.5,
      vjust = 0
    ) +
    annotate(
      "text",
      x = selected_mean,
      y = selected_label_y,
      label = "Average\nperformance\n(selected)",
      color = "gray15",
      size = 3.8,
      lineheight = 0.95,
      hjust = 0.5,
      vjust = 0
    ) +
    annotate(
      "text",
      x = intensity_label_x,
      y = intensity_y,
      label = paste0("Selection intensity = ", round(select_pct, 1), "%"),
      color = "#8B3F0F",
      size = 3.4,
      lineheight = 0.95,
      hjust = intensity_label_hjust
    ) +
    scale_fill_manual(values = c(
      "Original population" = "#CFE8C5",
      "Expected selected population" = "#58B947"
    )) +
    scale_color_manual(values = c(
      "Original population" = "#557A50",
      "Expected selected population" = "#1F6F2A"
    )) +
    coord_cartesian(ylim = c(gain_y * 2.2, y_max * 1.40), clip = "off") +
    labs(
      title = paste("Genetic gain response curve -", trait_name),
      subtitle = paste0(
        "H2 = ", round(h2, 3),
        " | Phenotypic SD = ", round(sd_p, 3),
        " | Direction: ", direction
      ),
      x = "Trait value",
      y = NULL,
      fill = NULL,
      color = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(color = "gray35"),
      legend.position = "top",
      legend.justification = "right",
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(t = 12, r = 22, b = 28, l = 16)
    )
}
plot_lpsi_direct_selection <- function(direct_table, selection_pct = 15, lpsi_results = NULL) {
  direct_table <- as.data.frame(direct_table)
  if (nrow(direct_table) == 0) {
    return(empty_plot("Run LPSI and choose a trait for direct selection."))
  }
  selection_pct <- suppressWarnings(as.numeric(selection_pct))
  if (length(selection_pct) == 0 || !is.finite(selection_pct[1])) selection_pct <- 15
  selection_pct <- min(max(selection_pct[1], 1), 100)
  plot_data <- direct_table %>%
    mutate(
      Trait_value = suppressWarnings(as.numeric(Trait_value)),
      Selection_status = ifelse(Selected_Direct, "Selected", "Not selected"),
      Plot_ID = ifelse(is.na(Original_ID) | Original_ID == "", as.character(ID), as.character(Original_ID))
    ) %>%
    filter(is.finite(Trait_value)) %>%
    arrange(Rank_Direct)
  if (nrow(plot_data) == 0) {
    return(empty_plot("Direct selection needs numeric trait values."))
  }
  primary_trait <- as.character(plot_data$Primary_trait[1])
  check_refs <- data.frame()
  if (!is.null(lpsi_results) && !is.null(lpsi_results$actual_adjusted_means) && !is.null(lpsi_results$selected_checks)) {
    means <- as.data.frame(lpsi_results$actual_adjusted_means)
    selected_checks <- as.character(lpsi_results$selected_checks)
    if (length(selected_checks) > 0 && "Original_ID" %in% names(means) && primary_trait %in% names(means)) {
      check_refs <- means %>%
        filter(as.character(Original_ID) %in% selected_checks) %>%
        transmute(
          Check = as.character(Original_ID),
          Trait_value = suppressWarnings(as.numeric(.data[[primary_trait]]))
        ) %>%
        arrange(match(Check, selected_checks)) %>%
        filter(is.finite(Trait_value))
    }
  }
  if (nrow(check_refs) > 0) {
    check_line_palette <- data.frame(
      Color_name = c(
        "Gray", "Blue", "Purple", "Teal", "Brown",
        "Orange", "Pink", "Dark green"
      ),
      Line_color = c(
        "#6F6F6F", "#2F6FDB", "#8E44AD", "#00897B", "#795548",
        "#F39C12", "#C2185B", "#2E7D32"
      ),
      stringsAsFactors = FALSE
    )
    check_refs <- check_refs %>%
      mutate(
        Check_order = row_number(),
        Color_name = check_line_palette$Color_name[
          ((Check_order - 1) %% nrow(check_line_palette)) + 1
        ],
        Line_color = check_line_palette$Line_color[
          ((Check_order - 1) %% nrow(check_line_palette)) + 1
        ]
      )
  }
  subtitle <- paste0("Selection intensity = ", round(selection_pct, 1), "%")
  if (nrow(check_refs) > 0) {
    check_line_text <- paste0(
      check_refs$Color_name,
      " dot-dash = ",
      check_refs$Check
    )
    subtitle <- paste(subtitle, paste(check_line_text, collapse = " | "), sep = " | ")
  }
  ggplot(
    plot_data,
    aes(x = reorder(Plot_ID, -Rank_Direct), y = Trait_value, fill = Selection_status)
  ) +
    geom_col(width = 0.62) +
    {
      if (nrow(check_refs) > 0) {
        geom_hline(
          data = check_refs,
          aes(yintercept = Trait_value, color = Line_color),
          inherit.aes = FALSE,
          linetype = "dotdash",
          linewidth = 0.7
        )
      }
    } +
    geom_text(
      aes(label = paste0("#", Rank_Direct)),
      hjust = -0.12,
      size = 3.0,
      color = "gray20"
    ) +
    coord_flip() +
    scale_fill_manual(values = c("Selected" = "#1D9E75", "Not selected" = "#B8B8B8")) +
    scale_color_identity() +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.16))) +
    labs(
      title = paste("Single trait selection -", primary_trait),
      subtitle = subtitle,
      x = "Genotype",
      y = primary_trait,
      fill = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(color = "gray35"),
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
}
make_diagnostic_plot <- function(diag, trait, plot_type) {
  df <- diag$data
  if (!trait %in% names(df)) {
    return(empty_plot("Trait not found."))
  }
  model <- tryCatch({
    fit_diagnostic_model(df, trait)
  }, error = function(e) {
    NULL
  })
  if (plot_type == "Raw trait histogram") {
    return(
      ggplot(df, aes(x = .data[[trait]])) +
        geom_histogram(bins = 20, fill = "#2D89C8", color = "white") +
        theme_minimal(base_size = 14) +
        labs(
          title = paste("Raw data histogram:", trait),
          x = trait,
          y = "Frequency"
        )
    )
  }
  if (is.null(model)) {
    return(empty_plot("Model could not be fitted for this trait."))
  }
  plot_df <- data.frame(
    Observed = model$model[[1]],
    Fitted = fitted(model),
    Residual = residuals(model)
  )
  if (plot_type == "Residual histogram") {
    return(
      ggplot(plot_df, aes(x = Residual)) +
        geom_histogram(bins = 20, fill = "#2D89C8", color = "white") +
        theme_minimal(base_size = 14) +
        labs(
          title = paste("Residual histogram:", trait),
          x = "Residual",
          y = "Frequency"
        )
    )
  }
  if (plot_type == "QQ plot") {
    return(
      ggplot(plot_df, aes(sample = Residual)) +
        stat_qq(size = 2) +
        stat_qq_line(color = "#E85D24", linewidth = 1) +
        theme_minimal(base_size = 14) +
        labs(
          title = paste("QQ plot of residuals:", trait),
          x = "Theoretical quantiles",
          y = "Sample quantiles"
        )
    )
  }
  if (plot_type == "Residuals vs fitted") {
    return(
      ggplot(plot_df, aes(x = Fitted, y = Residual)) +
        geom_point(size = 2, alpha = 0.75) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "#E85D24") +
        theme_minimal(base_size = 14) +
        labs(
          title = paste("Residuals vs fitted:", trait),
          x = "Fitted value",
          y = "Residual"
        )
    )
  }
  if (plot_type == "Observed vs fitted") {
    return(
      ggplot(plot_df, aes(x = Fitted, y = Observed)) +
        geom_point(size = 2, alpha = 0.75) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#E85D24") +
        theme_minimal(base_size = 14) +
        labs(
          title = paste("Observed data against fitted value:", trait),
          x = "Fitted value",
          y = "Observed value"
        )
    )
  }
  empty_plot("Unknown plot type.")
}

plot_lpsi_mean_comparison <- function(results, trait) {
  letters_dat <- as.data.frame(results$lsd_long %||% data.frame())
  raw_dat <- as.data.frame(results$cleaned_data %||% data.frame())
  if (nrow(letters_dat) == 0 || !trait %in% letters_dat$Trait ||
      nrow(raw_dat) == 0 || !trait %in% names(raw_dat)) {
    return(empty_plot("Mean comparison is not available for the selected trait."))
  }
  letters_dat <- letters_dat[
    letters_dat$Trait == trait & is.finite(suppressWarnings(as.numeric(letters_dat$emmean))),
    , drop = FALSE
  ]
  if (nrow(letters_dat) == 0) return(empty_plot("No adjusted means are available for the selected trait."))
  letters_dat$emmean <- suppressWarnings(as.numeric(letters_dat$emmean))
  letters_dat$Label <- as.character(letters_dat$Original_ID %||% letters_dat$ID)
  letters_dat$Group <- trimws(as.character(letters_dat$LSD_group %||% ""))
  anova_trait <- as.data.frame(results$anova_full %||% data.frame())
  has_significant_difference <- nrow(anova_trait) > 0 &&
    all(c("Trait", "p_value") %in% names(anova_trait)) &&
    any(anova_trait$Trait == trait & suppressWarnings(as.numeric(anova_trait$p_value)) < 0.05, na.rm = TRUE)
  if (!has_significant_difference) letters_dat$Group <- ""

  plot_dat <- raw_dat[, c("ID", "Original_ID", trait), drop = FALSE]
  names(plot_dat)[3] <- "Value"
  plot_dat$Value <- suppressWarnings(as.numeric(plot_dat$Value))
  plot_dat$Label <- as.character(plot_dat$Original_ID %||% plot_dat$ID)
  plot_dat <- plot_dat[is.finite(plot_dat$Value) & !is.na(plot_dat$Label) & plot_dat$Label != "", , drop = FALSE]
  if (nrow(plot_dat) == 0) return(empty_plot("No phenotypic observations are available for the selected trait."))

  label_levels <- letters_dat$Label[order(letters_dat$emmean, decreasing = TRUE)]
  label_levels <- unique(c(label_levels, setdiff(unique(plot_dat$Label), label_levels)))
  plot_dat$Label <- factor(plot_dat$Label, levels = label_levels)
  letters_dat$Label <- factor(letters_dat$Label, levels = label_levels)
  value_range <- diff(range(plot_dat$Value, na.rm = TRUE))
  label_gap <- max(value_range * 0.06, max(abs(plot_dat$Value), na.rm = TRUE) * 0.015, 0.03)
  label_y <- aggregate(Value ~ Label, plot_dat, max, na.rm = TRUE)
  letters_dat <- merge(letters_dat, label_y, by = "Label", all.x = TRUE, sort = FALSE)
  letters_dat$Label <- factor(letters_dat$Label, levels = label_levels)

  ggplot(plot_dat, aes(x = Label, y = Value, fill = Label)) +
    geom_boxplot(width = 0.68, linewidth = 0.45, outlier.size = 1.5, show.legend = FALSE) +
    geom_text(
      data = letters_dat[letters_dat$Group != "" & !is.na(letters_dat$Group), , drop = FALSE],
      aes(x = Label, y = Value + label_gap, label = Group),
      fontface = "bold", vjust = 0, size = 4, inherit.aes = FALSE
    ) +
    scale_fill_hue(c = 85, l = 62) +
    labs(
      title = "Distribution of Phenotypic Values",
      subtitle = paste("Mean comparison for", trait),
      x = "Entry",
      y = trait
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", color = "#111111", size = 12, margin = margin(5, 8, 5, 8)),
      plot.title.position = "plot",
      plot.background = element_rect(fill = "#FFFFFF", color = "#9FC5E8", linewidth = 0.8),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      panel.grid.major.x = element_blank(),
      legend.position = "none"
    ) +
    expand_limits(y = max(plot_dat$Value, na.rm = TRUE) + 2.5 * label_gap)
}



# Main LPSI function
run_selection_pipeline <- function(df_raw, check_varieties = NULL,
                                   advance_cutoff = advance_index_cutoff,
                                   retest_cutoff = retest_index_cutoff,
                                   priority_cutoff_pct = priority_advance_cutoff_pct,
                                   severe_weak_pct = priority_severe_weak_pct) {
  advance_index_cutoff <- suppressWarnings(as.numeric(advance_cutoff))
  retest_index_cutoff <- suppressWarnings(as.numeric(retest_cutoff))
  priority_advance_cutoff_pct <- suppressWarnings(as.numeric(priority_cutoff_pct))
  priority_severe_weak_pct <- suppressWarnings(as.numeric(severe_weak_pct))
  if (!is.finite(advance_index_cutoff)) advance_index_cutoff <- 0
  if (!is.finite(retest_index_cutoff)) retest_index_cutoff <- 0
  if (!is.finite(priority_advance_cutoff_pct)) priority_advance_cutoff_pct <- 0
  if (!is.finite(priority_severe_weak_pct)) priority_severe_weak_pct <- -10
  
  prepared <- prepare_excel_input(df_raw)
  df_data <- prepared$data
  trait_cols <- prepared$trait_cols
  weights_raw_used <- prepared$weights_raw_used
  trait_direction <- prepared$trait_direction
  target_traits <- prepared$target_traits
  check_original_name <- prepared$check_original_name
  df <- df_data %>%
    dplyr::select(all_of(c(id_col, rep_col, trait_cols))) %>%
    mutate(
      across(all_of(c(id_col, rep_col)), as.character),
      across(all_of(trait_cols), to_number)
    ) %>%
    rename(
      ID  = all_of(id_col),
      Rep = all_of(rep_col)
    ) %>%
    filter(!is.na(ID) & trimws(ID) != "") %>%
    filter(!is.na(Rep) & trimws(Rep) != "")
  if (nrow(df) == 0) {
    stop("No valid rows remained after cleaning ID and Rep.")
  }
  trait_cols <- trait_cols[
    sapply(df[trait_cols], function(x) !all(is.na(x)))
  ]
  if (length(trait_cols) == 0) {
    stop("No usable trait columns found after numeric cleaning.")
  }
  weights_raw_used <- weights_raw_used[trait_cols]
  trait_direction <- trait_direction[trait_cols]
  target_traits <- target_traits[names(target_traits) %in% trait_cols]
  weights_raw_used[is.na(weights_raw_used)] <- 1
  weights_raw_used[weights_raw_used < 0] <- 0
  if (sum(weights_raw_used) <= 0) {
    weights_raw_used <- setNames(rep(1, length(trait_cols)), trait_cols)
  }
  weights <- weights_raw_used / sum(weights_raw_used)
  variety_levels <- unique(df$ID)
  selected_checks <- unique(clean_text(check_varieties))
  selected_checks <- selected_checks[
    !is.na(selected_checks) & selected_checks != "" & selected_checks %in% variety_levels
  ]
  check_selection_source <- "Manual selection"
  if (length(selected_checks) == 0) {
    selected_checks <- prepared$detected_check_varieties
    selected_checks <- selected_checks[
      !is.na(selected_checks) & selected_checks != "" & selected_checks %in% variety_levels
    ]
    check_selection_source <- if (length(selected_checks) > 0) {
      paste0("Detected from ", prepared$type_col, " column")
    } else {
      "First genotype fallback"
    }
  }
  if (length(selected_checks) == 0) {
    selected_checks <- check_original_name
  }
  variety_map <- setNames(as.character(seq_along(variety_levels)), variety_levels)
  if (!all(selected_checks %in% names(variety_map))) {
    stop(
      "Selected check genotype was not found after cleaning: ",
      paste(setdiff(selected_checks, names(variety_map)), collapse = ", ")
    )
  }
  check_labels <- as.character(variety_map[selected_checks])
  check_label <- check_labels[1]
  mapping_table <- data.frame(
    ID_number = as.character(variety_map),
    Original_variety = names(variety_map),
    Check = ifelse(names(variety_map) %in% selected_checks, "YES", "")
  )
  df$Original_ID <- df$ID
  df$ID <- as.character(variety_map[df$ID])
  df$ID  <- factor(df$ID, levels = as.character(seq_along(variety_levels)))
  df$Rep <- as.factor(df$Rep)
  df$ID <- relevel(df$ID, ref = check_label)
  id_lookup <- mapping_table %>%
    rename(
      ID = ID_number,
      Original_ID = Original_variety
    ) %>%
    dplyr::select(ID, Original_ID)
  check_original_label <- paste(selected_checks, collapse = ", ")
  score_cols <- trait_cols[
    sapply(trait_cols, function(tr) {
      vals <- na.omit(as.numeric(df[[tr]]))
      length(vals) > 0 &&
        all(vals == floor(vals)) &&
        min(vals, na.rm = TRUE) >= 1 &&
        max(vals, na.rm = TRUE) <= 5
    })
  ]
  df_index <- df
  for (tr in trait_cols) {
    if (tr %in% names(target_traits)) {
      target_value <- target_traits[[tr]]
      df_index[[tr]] <- -abs(df_index[[tr]] - target_value)
    } else if (trait_direction[[tr]] == "Lower better") {
      df_index[[tr]] <- -df_index[[tr]]
    } else {
      df_index[[tr]] <- df_index[[tr]]
    }
  }
  get_adjusted_means <- function(data, trait) {
    model_data <- data %>%
      filter(!is.na(.data[[trait]]))
    if (nrow(model_data) == 0) {
      return(
        data.frame(
          ID = levels(data$ID),
          Value = NA_real_,
          Trait = trait
        )
      )
    }
    out <- tryCatch({
      model <- aov(make_model_formula(trait, model_data), data = model_data)
      em <- emmeans(model, ~ ID)
      as.data.frame(em) %>%
        dplyr::select(ID, emmean) %>%
        rename(Value = emmean) %>%
        mutate(
          ID = as.character(ID),
          Trait = trait
        )
    }, error = function(e) {
      model_data %>%
        group_by(ID) %>%
        summarise(
          Value = mean(.data[[trait]], na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(
          ID = as.character(ID),
          Trait = trait
        )
    })
    out
  }
  adj_means_long <- map_dfr(
    trait_cols,
    ~ get_adjusted_means(df_index, .x)
  )
  adj_means <- adj_means_long %>%
    pivot_wider(
      names_from = Trait,
      values_from = Value
    ) %>%
    arrange(as.numeric(as.character(ID))) %>%
    left_join(id_lookup, by = "ID") %>%
    relocate(Original_ID, .after = ID)
  std_means <- adj_means %>%
    mutate(across(
      all_of(trait_cols),
      standardize_trait
    ))
  weight_table <- data.frame(
    Trait = names(weights),
    Direction = as.character(trait_direction[names(weights)]),
    Raw_weight = as.numeric(weights_raw_used[names(weights)]),
    Normalized_weight = round(as.numeric(weights), 4),
    Priority_trait = ifelse(
      as.numeric(weights_raw_used[names(weights)]) >= priority_weight_cutoff,
      "YES",
      ""
    )
  )
  priority_traits <- names(weights_raw_used[weights_raw_used >= priority_weight_cutoff])
  index_df <- std_means %>%
    rowwise() %>%
    mutate(
      Selection_Index = sum(
        c_across(all_of(trait_cols)) * weights[trait_cols],
        na.rm = TRUE
      )
    ) %>%
    ungroup() %>%
    dplyr::select(ID, Original_ID, Selection_Index) %>%
    arrange(desc(Selection_Index))
  check_index_details <- index_df %>%
    filter(ID %in% check_labels) %>%
    arrange(match(ID, check_labels)) %>%
    dplyr::select(ID, Original_ID, Selection_Index)
  check_index <- check_index_details %>%
    summarise(Check_index = mean(Selection_Index, na.rm = TRUE)) %>%
    pull(Check_index)
  if (length(check_index) == 0 || !is.finite(check_index)) {
    stop("Check index not found.")
  }
  index_df <- index_df %>%
    mutate(
      Index_Advantage = Selection_Index - check_index,
      SI_pct_index = round(
        Index_Advantage /
          ifelse(abs(check_index) < 0.0001, 1, abs(check_index)) * 100,
        2
      ),
      Status = case_when(
        ID %in% check_labels ~ "Check",
        Index_Advantage > 0 ~ "Above check",
        Index_Advantage == 0 ~ "Equal to check",
        TRUE ~ "Below check"
      )
    )
  actual_means_long <- map_dfr(
    trait_cols,
    ~ get_adjusted_means(df, .x)
  )
  actual_means <- actual_means_long %>%
    pivot_wider(
      names_from = Trait,
      values_from = Value
    ) %>%
    arrange(as.numeric(as.character(ID))) %>%
    left_join(id_lookup, by = "ID") %>%
    relocate(Original_ID, .after = ID)
  check_actual <- actual_means %>%
    filter(ID %in% check_labels)
  if (nrow(check_actual) == 0) {
    stop("Check not found in actual adjusted means.")
  }
  check_actual_values <- vapply(
    trait_cols,
    function(tr) mean(as.numeric(check_actual[[tr]]), na.rm = TRUE),
    numeric(1)
  )
  superiority_benchmark_label <- if (length(check_labels) > 1) {
    paste0("Mean of selected checks (", check_original_label, ")")
  } else {
    paste0("Check ", check_original_label)
  }
  calc_superiority_pct <- function(candidate, check, trait_name) {
    if (is.na(check)) {
      return(NA_real_)
    }
    denom <- ifelse(abs(check) < 0.0001, 1, abs(check))
    if (trait_name %in% names(target_traits)) {
      target_value <- target_traits[[trait_name]]
      cand_dist  <- abs(candidate - target_value)
      check_dist <- abs(check - target_value)
      denom2 <- ifelse(abs(check_dist) < 0.0001, 1, abs(check_dist))
      return((check_dist - cand_dist) / denom2 * 100)
    }
    if (trait_direction[[trait_name]] == "Lower better") {
      return((check - candidate) / denom * 100)
    }
    return((candidate - check) / denom * 100)
  }
  superiority_df <- actual_means %>%
    filter(!ID %in% check_labels)
  for (tr in trait_cols) {
    superiority_df[[tr]] <- calc_superiority_pct(
      candidate  = superiority_df[[tr]],
      check      = check_actual_values[[tr]],
      trait_name = tr
    )
  }
  superiority_df <- superiority_df %>%
    mutate(across(
      all_of(trait_cols),
      ~ round(.x, 1)
    )) %>%
    mutate(
      Check_Benchmark = superiority_benchmark_label,
      .after = Original_ID
    )
  candidate_actual_long <- actual_means %>%
    filter(!ID %in% check_labels) %>%
    dplyr::select(ID, Original_ID, all_of(trait_cols)) %>%
    pivot_longer(
      cols = all_of(trait_cols),
      names_to = "Trait",
      values_to = "Candidate_Value"
    )
  check_actual_long <- check_actual %>%
    dplyr::select(ID, Original_ID, all_of(trait_cols)) %>%
    rename(
      Check_ID = ID,
      Check_Original_ID = Original_ID
    ) %>%
    pivot_longer(
      cols = all_of(trait_cols),
      names_to = "Trait",
      values_to = "Check_Value"
    )
  superiority_by_check_long <- candidate_actual_long %>%
    left_join(
      check_actual_long,
      by = "Trait",
      relationship = "many-to-many"
    ) %>%
    rowwise() %>%
    mutate(
      Superiority_pct = calc_superiority_pct(
        Candidate_Value,
        Check_Value,
        Trait
      )
    ) %>%
    ungroup() %>%
    mutate(Superiority_pct = round(Superiority_pct, 1)) %>%
    dplyr::select(
      ID,
      Original_ID,
      Check_ID,
      Check_Original_ID,
      Trait,
      Candidate_Value,
      Check_Value,
      Superiority_pct
    )
  if (length(priority_traits) > 0 && nrow(superiority_df) > 0) {
    priority_summary <- superiority_df %>%
      dplyr::select(ID, Original_ID, all_of(priority_traits)) %>%
      pivot_longer(
        cols = all_of(priority_traits),
        names_to = "Trait",
        values_to = "Superiority_pct"
      ) %>%
      group_by(ID, Original_ID) %>%
      summarise(
        n_priority_traits = n(),
        n_priority_above_or_equal_check = sum(
          Superiority_pct >= priority_advance_cutoff_pct,
          na.rm = TRUE
        ),
        n_priority_below_check = sum(
          Superiority_pct < priority_advance_cutoff_pct,
          na.rm = TRUE
        ),
        n_priority_severe_weak = sum(
          Superiority_pct < priority_severe_weak_pct,
          na.rm = TRUE
        ),
        priority_pass_rate = n_priority_above_or_equal_check / n_priority_traits,
        .groups = "drop"
      )
  } else {
    priority_summary <- superiority_df %>%
      dplyr::select(ID, Original_ID) %>%
      mutate(
        n_priority_traits = 0,
        n_priority_above_or_equal_check = 0,
        n_priority_below_check = 0,
        n_priority_severe_weak = 0,
        priority_pass_rate = NA_real_
      )
  }
  weakness_traits <- if (length(priority_traits) > 0) priority_traits else trait_cols
  weakness_table <- superiority_df %>%
    dplyr::select(ID, Original_ID, all_of(weakness_traits)) %>%
    pivot_longer(
      cols = all_of(weakness_traits),
      names_to = "Trait",
      values_to = "Superiority_pct"
    ) %>%
    group_by(ID, Original_ID) %>%
    summarise(
      Weakness_trait = {
        weak <- Trait[!is.na(Superiority_pct) & Superiority_pct < priority_advance_cutoff_pct]
        if (length(weak) == 0) "None" else paste(unique(weak), collapse = ", ")
      },
      .groups = "drop"
    )
  final_decision <- index_df %>%
    filter(!ID %in% check_labels) %>%
    left_join(
      priority_summary %>%
        dplyr::select(
          ID,
          n_priority_traits,
          n_priority_above_or_equal_check,
          n_priority_below_check,
          n_priority_severe_weak,
          priority_pass_rate
        ),
      by = "ID"
    ) %>%
    left_join(
      weakness_table %>% dplyr::select(ID, Weakness_trait),
      by = "ID"
    ) %>%
    mutate(
      across(
        c(
          n_priority_traits,
          n_priority_above_or_equal_check,
          n_priority_below_check,
          n_priority_severe_weak
        ),
        ~ replace_na(.x, 0)
      ),
      Weakness_trait = replace_na(Weakness_trait, "None"),
      Plot_ID = Original_ID,
      Decision = case_when(
        Index_Advantage >= advance_index_cutoff &
          n_priority_below_check == 0 ~ "ADVANCE",
        Index_Advantage >= advance_index_cutoff &
          n_priority_below_check > 0 ~ "RETEST",
        Index_Advantage >= retest_index_cutoff &
          n_priority_severe_weak == 0 ~ "RETEST",
        TRUE ~ "DISCARD"
      ),
      Decision_reason = case_when(
        Decision == "ADVANCE" ~ paste0(
          "Selection index is equal to or higher than check, ",
          "and all priority traits are equal to or better than check."
        ),
        Decision == "RETEST" &
          Index_Advantage >= advance_index_cutoff ~ paste0(
            "Selection index is good, but one or more priority traits ",
            "are below check; retesting is recommended."
          ),
        Decision == "RETEST" ~ paste0(
          "Selection index is slightly below check, but priority traits ",
          "are not severely weak."
        ),
        TRUE ~ paste0(
          "Selection index is clearly below check or has severe weakness ",
          "in priority traits."
        )
      )
    ) %>%
    arrange(desc(Selection_Index))
  decision_colors <- c(
    "ADVANCE" = "#1D9E75",
    "RETEST"  = "#3c45e6",
    "DISCARD" = "#999999"
  )
  if (nrow(final_decision) > 0) {
    check_line_palette <- data.frame(
      Color_name = c(
        "Gray", "Blue", "Purple", "Teal", "Brown",
        "Orange", "Pink", "Dark green"
      ),
      Line_color = c(
        "#6F6F6F", "#2F6FDB", "#8E44AD", "#00897B", "#795548",
        "#F39C12", "#C2185B", "#2E7D32"
      ),
      stringsAsFactors = FALSE
    )
    check_index_details <- check_index_details %>%
      mutate(
        Check_order = row_number(),
        Color_name = check_line_palette$Color_name[
          ((Check_order - 1) %% nrow(check_line_palette)) + 1
        ],
        Line_color = check_line_palette$Line_color[
          ((Check_order - 1) %% nrow(check_line_palette)) + 1
        ]
      )
    check_line_text <- paste0(
      check_index_details$Color_name,
      " dotted = ",
      check_index_details$Original_ID
    )
    check_line_note <- if (nrow(check_index_details) > 1) {
      paste0(
        "Red dashed = mean checks | ",
        paste(check_line_text, collapse = " | ")
      )
    } else {
      paste0(
        "Red dashed = check ",
        check_original_label
      )
    }
    p_index <- ggplot(
      final_decision,
      aes(
        x = reorder(Plot_ID, Selection_Index),
        y = Selection_Index,
        fill = Decision
      )
    ) +
      geom_col(width = 0.58) +
      {
        if (nrow(check_index_details) > 1) {
          geom_hline(
            data = check_index_details,
            aes(yintercept = Selection_Index, color = Line_color),
            inherit.aes = FALSE,
            linetype = "dotted",
            linewidth = 0.55
          )
        }
      } +
      geom_hline(
        yintercept = check_index,
        linetype = "dashed",
        color = "#E85D24",
        linewidth = 0.8
      ) +
      geom_text(
        aes(
          label = paste0(
            sprintf("%.2f", Selection_Index),
            " | ",
            sprintf("%.2f", Index_Advantage)
          ),
          hjust = ifelse(Selection_Index >= 0, -0.10, 1.10)
        ),
        size = 3.0
      ) +
      scale_fill_manual(values = decision_colors) +
      scale_color_identity() +
      scale_y_continuous(
        expand = expansion(mult = c(0.12, 0.25))
      ) +
      coord_flip() +
      labs(
        title = "Hybrid selection index ranking",
        subtitle = check_line_note,
        x = "Hybrid ID",
        y = "Weighted standardized selection index",
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
        legend.position = "bottom",
        legend.margin = margin(t = -4),
        legend.box.margin = margin(t = -6),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 8, r = 12, b = 6, l = 8)
      )
  } else {
    p_index <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No candidate genotype found.") +
      theme_void()
  }
  heatmap_plot <- NULL
  heatmap_data <- superiority_by_check_long %>%
    filter(!is.na(Superiority_pct))
  if (nrow(heatmap_data) > 0) {
    max_abs_sup <- max(abs(heatmap_data$Superiority_pct), na.rm = TRUE)
    if (!is.finite(max_abs_sup) || max_abs_sup < 0.0001) {
      max_abs_sup <- 1
    }
    hybrid_levels <- natural_id_order(heatmap_data$Original_ID)
    if (length(hybrid_levels) == 0) {
      hybrid_levels <- unique(heatmap_data$Original_ID)
    }
    check_levels <- check_index_details$Original_ID
    heatmap_data <- heatmap_data %>%
      mutate(
        Original_ID = factor(Original_ID, levels = hybrid_levels),
        Trait = factor(Trait, levels = rev(trait_cols)),
        Check_Original_ID = factor(Check_Original_ID, levels = check_levels),
        Cell_label = format_superiority_pct_label(Superiority_pct)
      )
    heatmap_plot <- ggplot(
      heatmap_data,
      aes(x = Original_ID, y = Trait, fill = Superiority_pct)
    ) +
      geom_tile(color = "white", linewidth = 0.25) +
      geom_text(aes(label = Cell_label), size = 2.4, color = "gray10") +
      facet_grid(Check_Original_ID ~ .) +
      scale_fill_gradient2(
        low = "#D85A30",
        mid = "white",
        high = "#1D9E75",
        midpoint = 0,
        limits = c(-max_abs_sup, max_abs_sup)
      ) +
      labs(
        title = "Superiority index (%) by selected check",
        subtitle = paste0(
          "Each panel uses its own check as the reference. ",
          "Green = better than that panel's check; red = worse."
        ),
        x = "Hybrid ID",
        y = NULL,
        fill = "% vs check"
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
        plot.subtitle = element_text(face = "bold", size = 10, hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 9),
        strip.text.y = element_text(face = "bold"),
        panel.grid = element_blank(),
        legend.position = "right",
        plot.margin = margin(t = 8, r = 10, b = 8, l = 8)
      )
  }
  anova_full <- data.frame()
  if (run_simple_anova) {
    for (tr in trait_cols) {
      model_data <- df %>%
        filter(!is.na(.data[[tr]]))
      temp <- tryCatch({
        normality <- get_lpsi_normality_decision(df, tr)
        if (isTRUE(normality$use_anova)) {
          model <- aov(make_model_formula(tr, model_data), data = model_data)
          sm <- as.data.frame(summary(model)[[1]])
          data.frame(
            Trait = tr,
            Test = "ANOVA",
            Normality_p = round(normality$p_value, 5),
            Normality_note = "ANOVA used",
            Source = rownames(sm),
            Df = sm$Df,
            Sum_Sq = round(sm$`Sum Sq`, 4),
            Mean_Sq = round(sm$`Mean Sq`, 4),
            F_value = ifelse(is.na(sm$`F value`), NA, round(sm$`F value`, 4)),
            p_value = ifelse(is.na(sm$`Pr(>F)`), NA, round(sm$`Pr(>F)`, 5)),
            Sig = sig_label(sm$`Pr(>F)`),
            stringsAsFactors = FALSE
          )
        } else {
          kw_data <- model_data %>%
            filter(!is.na(ID), !is.na(.data[[tr]]))
          if (n_distinct(kw_data$ID) < 2) {
            stop("Kruskal-Wallis requires at least 2 ID groups.")
          }
          kw <- kruskal.test(
            as.formula(paste0(backtick_name(tr), " ~ ID")),
            data = kw_data
          )
          data.frame(
            Trait = tr,
            Test = "Kruskal-Wallis",
            Normality_p = round(normality$p_value, 5),
            Normality_note = "Kruskal-Wallis used",
            Source = "ID",
            Df = as.numeric(kw$parameter),
            Sum_Sq = NA_real_,
            Mean_Sq = NA_real_,
            F_value = round(as.numeric(kw$statistic), 4),
            p_value = round(kw$p.value, 5),
            Sig = sig_label(kw$p.value),
            stringsAsFactors = FALSE
          )
        }
      }, error = function(e) {
        data.frame(
          Trait = tr,
          Test = "Failed",
          Normality_p = NA,
          Normality_note = "",
          Source = "ANOVA failed",
          Df = NA,
          Sum_Sq = NA,
          Mean_Sq = NA,
          F_value = NA,
          p_value = NA,
          Sig = e$message,
          stringsAsFactors = FALSE
        )
      })
      anova_full <- bind_rows(anova_full, temp)
    }
  }
  selection_k <- dnorm(qnorm(1 - lpsi_selection_intensity)) / lpsi_selection_intensity
  heritability_gain <- map_dfr(trait_cols, function(tr) {
    model_data <- df %>%
      filter(!is.na(.data[[tr]]))
    tryCatch({
      if (n_distinct(model_data$ID) < 2) {
        stop("At least two varieties are required.")
      }
      model <- aov(make_model_formula(tr, model_data), data = model_data)
      sm <- as.data.frame(summary(model)[[1]])
      sm$Source <- trimws(rownames(sm))
      id_row <- sm[sm$Source == "ID", , drop = FALSE]
      residual_row <- sm[grepl("Residual", sm$Source), , drop = FALSE]
      if (nrow(id_row) == 0 || nrow(residual_row) == 0) {
        stop("ID or residual mean square was not available.")
      }
      ms_genotype <- as.numeric(id_row$`Mean Sq`[1])
      ms_error <- as.numeric(residual_row$`Mean Sq`[1])
      rep_harmonic <- model_data %>%
        group_by(ID) %>%
        summarise(n_rep = n_distinct(Rep), .groups = "drop") %>%
        summarise(value = 1 / mean(1 / n_rep)) %>%
        pull(value)
      genotypic_var <- max((ms_genotype - ms_error) / rep_harmonic, 0)
      error_var <- ms_error
      phenotypic_var <- genotypic_var + (error_var / rep_harmonic)
      h2 <- ifelse(phenotypic_var > 0, genotypic_var / phenotypic_var, NA_real_)
      trait_mean <- mean(model_data[[tr]], na.rm = TRUE)
      genetic_advance <- selection_k * sqrt(phenotypic_var) * h2
      genetic_advance_pct <- ifelse(
        is.finite(trait_mean) && abs(trait_mean) > 0.0001,
        (genetic_advance / trait_mean) * 100,
        NA_real_
      )
      cv_pct <- ifelse(
        is.finite(trait_mean) && abs(trait_mean) > 0.0001 && is.finite(error_var) && error_var >= 0,
        (sqrt(error_var) / abs(trait_mean)) * 100,
        NA_real_
      )
      data.frame(
        Trait = tr,
        Mean = round(trait_mean, 4),
        CV_pct = round(cv_pct, 2),
        MS_genotype = round(ms_genotype, 4),
        MS_error = round(ms_error, 4),
        Harmonic_replication = round(rep_harmonic, 3),
        Genotypic_variance = round(genotypic_var, 5),
        Error_variance = round(error_var, 5),
        Phenotypic_variance = round(phenotypic_var, 5),
        Broad_sense_H2 = round(h2, 4),
        Selection_intensity_pct = lpsi_selection_intensity * 100,
        Selection_intensity_k = round(selection_k, 4),
        Genetic_advance = round(genetic_advance, 4),
        Genetic_advance_pct_mean = round(genetic_advance_pct, 2),
        Note = "Single-environment replicated estimate",
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      data.frame(
        Trait = tr,
        Mean = NA_real_,
        CV_pct = NA_real_,
        MS_genotype = NA_real_,
        MS_error = NA_real_,
        Harmonic_replication = NA_real_,
        Genotypic_variance = NA_real_,
        Error_variance = NA_real_,
        Phenotypic_variance = NA_real_,
        Broad_sense_H2 = NA_real_,
        Selection_intensity_pct = lpsi_selection_intensity * 100,
        Selection_intensity_k = round(selection_k, 4),
        Genetic_advance = NA_real_,
        Genetic_advance_pct_mean = NA_real_,
        Note = paste("Not calculated:", e$message),
        stringsAsFactors = FALSE
      )
    })
  })
  lsd_all <- data.frame()
  lsd_wide <- data.frame()
  if (run_lsd_test) {
    for (tr in trait_cols) {
      model_data <- df %>%
        filter(!is.na(.data[[tr]]))
      temp <- tryCatch({
        normality <- get_lpsi_normality_decision(df, tr)
        model <- aov(make_model_formula(tr, model_data), data = model_data)
        sm <- as.data.frame(summary(model)[[1]])
        sm$Source <- trimws(rownames(sm))
        id_p_value <- sm$`Pr(>F)`[sm$Source == "ID"]
        anova_id_p_value <- anova_full %>%
          filter(
            Trait == tr,
            trimws(as.character(Test)) == "ANOVA",
            trimws(as.character(Source)) == "ID"
          ) %>%
          pull(p_value)
        use_lsd_letters <- isTRUE(normality$use_anova) &&
          (
            (length(id_p_value) > 0 && !is.na(id_p_value[1]) && id_p_value[1] < lsd_significance_alpha) ||
              (length(anova_id_p_value) > 0 && !is.na(anova_id_p_value[1]) && anova_id_p_value[1] < lsd_significance_alpha)
          )
        em <- emmeans(model, ~ ID)
        if (use_lsd_letters) {
          if (df.residual(model) <= 0) {
            stop("Residual degrees of freedom is zero. LSD cannot be calculated.")
          }
          cld_tbl <- multcomp::cld(
            em,
            Letters = c(letters, LETTERS),
            adjust = "none",
            alpha = lsd_significance_alpha,
            sort = TRUE,
            reversed = trait_direction[[tr]] == "Higher better"
          ) %>%
            as.data.frame()
          cld_tbl$.group <- trimws(as.character(cld_tbl$.group))
          cld_tbl %>%
            mutate(
              Trait = tr,
              ID = as.character(ID),
              emmean = round(emmean, 3),
              SE = round(SE, 3),
              LSD_group = .group
            ) %>%
            left_join(id_lookup, by = "ID") %>%
            dplyr::select(
              Trait,
              ID,
              Original_ID,
              emmean,
              SE,
              LSD_group
            )
        } else {
          as.data.frame(em) %>%
            mutate(
              Trait = tr,
              ID = as.character(ID),
              emmean = round(emmean, 3),
              SE = round(SE, 3),
              LSD_group = ""
            ) %>%
            left_join(id_lookup, by = "ID") %>%
            dplyr::select(
              Trait,
              ID,
              Original_ID,
              emmean,
              SE,
              LSD_group
            )
        }
      }, error = function(e) {
        data.frame(
          Trait = tr,
          ID = NA,
          Original_ID = NA,
          emmean = NA,
          SE = NA,
          LSD_group = paste("LSD failed:", e$message),
          stringsAsFactors = FALSE
        )
      })
      lsd_all <- bind_rows(lsd_all, temp)
    }
    if (nrow(lsd_all) > 0) {
      lsd_wide <- lsd_all %>%
        mutate(
          Mean_Group = ifelse(
            is.na(LSD_group) | trimws(LSD_group) == "",
            as.character(emmean),
            paste0(emmean, " ", LSD_group)
          )
        ) %>%
        dplyr::select(ID, Original_ID, Trait, Mean_Group) %>%
        pivot_wider(
          names_from = Trait,
          values_from = Mean_Group
        ) %>%
        arrange(as.numeric(as.character(ID)))
    }
  }
  trait_info <- data.frame(
    Trait = trait_cols,
    Type = ifelse(trait_cols %in% score_cols, "Score 1-5", "Quantitative"),
    Direction = as.character(trait_direction[trait_cols]),
    Target_value = ifelse(
      trait_cols %in% names(target_traits),
      as.numeric(target_traits[trait_cols]),
      NA_real_
    ),
    Weight = as.numeric(weights_raw_used[trait_cols]),
    Priority_trait = ifelse(
      weights_raw_used[trait_cols] >= priority_weight_cutoff,
      "YES",
      ""
    )
  )
  return(list(
    raw_data = df_raw,
    cleaned_data = df,
    trait_info = trait_info,
    adjusted_means_index = adj_means,
    standardized_scores = std_means,
    weight_table = weight_table,
    index_ranking = index_df,
    check_index_details = check_index_details,
    actual_adjusted_means = actual_means,
    superiority_index = superiority_df,
    priority_summary = priority_summary,
    final_decision = final_decision,
    heritability_gain = heritability_gain,
    anova_full = anova_full,
    lsd_long = lsd_all,
    lsd_wide = lsd_wide,
    ranking_plot = p_index,
    heatmap_plot = heatmap_plot,
    check_original_label = check_original_label,
    selected_checks = selected_checks,
    decision_settings = data.frame(
      Advance_index_cutoff = advance_index_cutoff,
      Retest_index_cutoff = retest_index_cutoff,
      Priority_trait_cutoff_pct = priority_advance_cutoff_pct,
      Severe_weakness_cutoff_pct = priority_severe_weak_pct,
      Check_selection_source = check_selection_source,
      stringsAsFactors = FALSE
    )
  ))
}



# MET function
MET_W_YIELD <- 0.50
MET_W_FW <- 0.25
MET_W_ASV <- 0.25
MET_MIN_ENVS_FOR_BIPLOT <- 2
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
    c("Genotype", "Environment", met_rep_col_candidates, met_block_col_candidates, remove_cols)
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
      "Missing cells are predicted in the BLUP matrix used by AMMI/GGE."
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
      "High removal changes may affect BLUP and ranking decisions."
    )
  }
  if (nrow(genotype_summary) > 0 && "N_observations" %in% names(genotype_summary)) {
    n_obs <- suppressWarnings(as.numeric(genotype_summary$N_observations))
    add_qc(
      "Replication balance",
      "Observations per genotype",
      paste0(min(n_obs, na.rm = TRUE), " to ", max(n_obs, na.rm = TRUE)),
      if (length(unique(stats::na.omit(n_obs))) > 1) "Review" else "OK",
      "Unequal observations can increase shrinkage differences among genotype BLUPs."
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
  rank_blup <- rank(-selection_base$BLUP_G, ties.method = "average")
  rank_stability <- met_fill_rank(rank(abs(selection_base$Sens - 1), ties.method = "average", na.last = "keep"))
  rank_asv <- met_fill_rank(rank(selection_base$ASV, ties.method = "average", na.last = "keep"))
  selection_base %>%
    mutate(
      Rank_BLUP = rank_blup,
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
plot_met_selection_ranking <- function(selection, trait_used) {
  ggplot(selection, aes(x = reorder(Genotype, -BLUP_G), y = BLUP_G, fill = b_interp)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.3, color = "#2C3E50", linewidth = 0.6) +
    geom_text(aes(y = BLUP_G / 2, label = paste0("#", Final_rank)), vjust = 0.5, size = 3.5, fontface = "bold", color = "white") +
    scale_fill_manual(values = c("Responsive" = "#E74C3C", "Average" = "#F39C12", "Stable" = "#2ECC71"), na.value = "gray70") +
    labs(
      title = paste0("Hybrid selection - ", trait_used, " performance & stability"),
      subtitle = paste0(
        "Weights: Mean=", round(unique(selection$Mean_weight)[1], 1),
        "%, FW=", round(unique(selection$FW_weight)[1], 1),
        "%, ASV=", round(unique(selection$ASV_weight)[1], 1), "%"
      ),
      x = "Genotype",
      y = paste0("BLUP for ", trait_used),
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
run_met_pipeline <- function(
    df_raw,
    trait_used = NULL,
    check_varieties = NULL,
    replication_col = NULL,
    block_col = NULL,
    min_envs_for_biplot = NULL) {
  input <- make_met_data(
    df_raw,
    trait_used,
    replication_col = replication_col,
    block_col = block_col
  )
  dat <- input$data
  trait_used <- input$trait_used
  has_replication <- !is.null(input$replication_col)
  has_block <- !is.null(input$block_col)
  n_envs_total <- n_distinct(dat$Environment)
  min_envs_for_biplot <- suppressWarnings(as.integer(min_envs_for_biplot))
  if (length(min_envs_for_biplot) == 0 || is.na(min_envs_for_biplot) || min_envs_for_biplot < 1) {
    min_envs_for_biplot <- max(MET_MIN_ENVS_FOR_BIPLOT, ceiling(0.5 * n_envs_total))
  }
  min_envs_for_biplot <- min(max(1, min_envs_for_biplot), n_envs_total)
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
    notes <- c(notes, "No genotype was present in all environments. All genotypes were used as the control set for percent-of-control calculation.")
  }
  low_conf_genos <- presence %>% filter(n_envs == 1) %>% pull(Genotype) %>% as.character()
  if (length(low_conf_genos) > 0) notes <- c(notes, paste0("Single-environment genotypes have high uncertainty: ", paste(low_conf_genos, collapse = ", ")))
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
  model_summary <- data.frame(Trait_used = trait_used, N_rows_clean = nrow(dat_clean), N_genotypes = n_distinct(dat_clean$Genotype), N_environments = n_distinct(dat_clean$Environment), Replication_column = input$replication_col %||% "", Block_column = input$block_col %||% "", AMMI_GGE_min_observed_locations = min_envs_for_biplot, N_replications = if (has_replication) n_distinct(dat_clean$Rep) else NA_integer_, N_blocks = if (has_block) n_distinct(dat_clean$Block) else NA_integer_, Model_formula = paste(deparse(full_formula), collapse = " "), Harmonic_replication = round(n_r, 3), Stability_ratio_G_over_G_plus_GxE = round(stability_ratio, 4), Broad_sense_H2 = round(H2, 4), Controls_used = paste(controls_used, collapse = ", "), Notes = paste(notes, collapse = " | "), stringsAsFactors = FALSE)
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
  BLUPs_main <- data.frame(Genotype = rownames(g_re), BLUP_G = grand_mean + g_re[, 1], SE_G = se_g, Reliability = round(reliability_g, 4), CI_lower = grand_mean + g_re[, 1] - 1.96 * se_g, CI_upper = grand_mean + g_re[, 1] + 1.96 * se_g) %>% arrange(desc(BLUP_G)) %>% mutate(Rank_BLUP = row_number())
  p_blup <- ggplot(BLUPs_main, aes(x = reorder(Genotype, BLUP_G), y = BLUP_G)) + geom_point(color = "#2C3E50", size = 3) + geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.4, color = "#3498DB", linewidth = 0.7) + coord_flip() + labs(title = "Genotype BLUPs with 95% CI", x = "Genotype", y = paste0("BLUP for ", trait_used)) + theme_bw()
  env_effects <- ranef_full$Environment %>% rownames_to_column("Environment") %>% rename(BLUP_E = `(Intercept)`)
  BLUPs_GxE <- ranef_full$`Genotype:Environment` %>% rownames_to_column("Geno_Env") %>% rename(BLUP_GxE = `(Intercept)`) %>% separate(Geno_Env, into = c("Genotype", "Environment"), sep = ":", extra = "merge", fill = "right")
  BLUPs_env_obs <- BLUPs_GxE %>% left_join(BLUPs_main[, c("Genotype", "BLUP_G")], by = "Genotype") %>% left_join(env_effects, by = "Environment") %>% mutate(BLUP_env = grand_mean + (BLUP_G - grand_mean) + BLUP_E + BLUP_GxE, Source = "Observed") %>% arrange(Environment, desc(BLUP_env))
  all_combos <- expand.grid(Genotype = as.character(unique(dat_clean$Genotype)), Environment = as.character(unique(dat_clean$Environment)), stringsAsFactors = FALSE)
  imputed_cells <- all_combos %>% anti_join(dat_clean %>% distinct(Genotype = as.character(Genotype), Environment = as.character(Environment)), by = c("Genotype", "Environment")) %>% left_join(BLUPs_main[, c("Genotype", "BLUP_G")], by = "Genotype") %>% left_join(env_effects, by = "Environment") %>% mutate(BLUP_GxE = 0, BLUP_env = grand_mean + (BLUP_G - grand_mean) + BLUP_E, Source = dplyr::if_else(Genotype %in% low_conf_genos, "Imputed_low_confidence", "Imputed"))
  BLUPs_env_full <- bind_rows(BLUPs_env_obs %>% dplyr::select(Genotype, Environment, BLUP_G, BLUP_E, BLUP_GxE, BLUP_env, Source), imputed_cells %>% dplyr::select(Genotype, Environment, BLUP_G, BLUP_E, BLUP_GxE, BLUP_env, Source)) %>% arrange(Environment, desc(BLUP_env))
  acc_dat <- dat_clean %>% group_by(Genotype, Environment) %>% summarise(obs = mean(Weight), .groups = "drop") %>% mutate(Genotype = as.character(Genotype), Environment = as.character(Environment)) %>% left_join(BLUPs_env_full %>% dplyr::select(Genotype, Environment, BLUP_env), by = c("Genotype", "Environment")) %>% filter(!is.na(BLUP_env))
  r_by_env <- acc_dat %>% group_by(Environment) %>% summarise(r_val = ifelse(n() >= 2, round(cor(obs, BLUP_env, use = "complete.obs"), 4), NA_real_), .groups = "drop") %>% mutate(label = paste0("r = ", r_val))
  p_accuracy <- ggplot(acc_dat, aes(x = obs, y = BLUP_env, color = Genotype)) + geom_point(size = 3, alpha = 0.85) + geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7, alpha = 0.15) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40", linewidth = 0.7) + geom_text(data = r_by_env, aes(label = label, x = -Inf, y = Inf), hjust = -0.15, vjust = 1.4, size = 3.5, fontface = "bold", color = "#2C3E50", inherit.aes = FALSE) + facet_wrap(~Environment, scales = "free", ncol = 2) + labs(title = "Prediction accuracy per hybrid by location", x = paste0("Observed mean ", trait_used), y = "Predicted BLUP", color = "Genotype") + theme_bw() + theme(legend.position = "bottom")
  ctrl_means_env <- BLUPs_env_full %>% filter(Genotype %in% controls_used) %>% group_by(Environment) %>% summarise(Control_mean = mean(BLUP_env, na.rm = TRUE), .groups = "drop")
  BLUPs_env_full <- BLUPs_env_full %>% left_join(ctrl_means_env, by = "Environment") %>% mutate(Pct_of_controls = round((BLUP_env / Control_mean) * 100, 1))
  heatmap_dat <- BLUPs_env_full %>% mutate(label = paste0(round(BLUP_env, 1), "\n(", Pct_of_controls, "%)"), alpha_val = ifelse(grepl("Imputed", Source), 0.45, 1.0))
  p_perf_heatmap <- ggplot(heatmap_dat, aes(x = Environment, y = reorder(Genotype, BLUP_env), fill = Pct_of_controls)) + geom_tile(aes(alpha = alpha_val), color = "white", linewidth = 0.5) + geom_text(aes(label = label), size = 2.5, lineheight = 0.9) + scale_fill_gradient2(low = "#E74C3C", mid = "white", high = "#2ECC71", midpoint = 100) + scale_alpha_identity() + labs(title = "Hybrid performance relative to controls", subtitle = "Faded = LMM-imputed", x = "Environment", y = "Genotype", fill = "% of controls") + theme_bw()
  GxE_matrix_wide <- BLUPs_env_full %>% dplyr::select(Genotype, Environment, BLUP_env) %>% pivot_wider(names_from = Environment, values_from = BLUP_env) %>% column_to_rownames("Genotype")
  GxE_long_complete <- GxE_matrix_wide %>% rownames_to_column("Genotype") %>% pivot_longer(-Genotype, names_to = "Environment", values_to = "BLUP_env")
  n_genos <- nrow(GxE_matrix_wide)
  FW_dat_loo <- GxE_long_complete %>% group_by(Environment) %>% mutate(EnvIndex = (sum(BLUP_env) - BLUP_env) / (n_genos - 1)) %>% ungroup()
  FW_results <- FW_dat_loo %>% group_by(Genotype) %>% summarise(GenMean = mean(BLUP_env), Sens = coef(lm(BLUP_env ~ EnvIndex))[2], .groups = "drop") %>% mutate(b_interp = case_when(Sens > 1.1 ~ "Responsive", Sens < 0.9 ~ "Stable", TRUE ~ "Average")) %>% arrange(desc(GenMean))
  p_fw_mean_sens <- ggplot(FW_results, aes(x = GenMean, y = Sens, color = b_interp, label = Genotype)) + geom_point(size = 3) + geom_text(vjust = -0.8, size = 3) + geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") + scale_color_manual(values = c("Responsive" = "#E74C3C", "Average" = "#F39C12", "Stable" = "#2ECC71")) + labs(title = "Finlay-Wilkinson: mean vs sensitivity", x = "Genotype BLUP mean", y = "Sensitivity (b)", color = "Stability") + theme_bw()
  env_index_plot <- FW_dat_loo %>% group_by(Environment) %>% summarise(env_mean = mean(EnvIndex), .groups = "drop")
  p_fw_regression <- ggplot(FW_dat_loo, aes(x = EnvIndex, y = BLUP_env, color = Genotype, group = Genotype)) + geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) + geom_point(shape = 18, size = 3) + geom_vline(xintercept = mean(env_index_plot$env_mean), linetype = "dashed", color = "red", linewidth = 0.7) + geom_vline(data = env_index_plot, aes(xintercept = env_mean), linetype = "solid", color = "gray40", linewidth = 0.3, alpha = 0.5, inherit.aes = FALSE) + geom_text(data = env_index_plot, aes(x = env_mean, y = -Inf, label = Environment), angle = 0, vjust = -1.0, hjust = 0.5, size = 3, color = "gray20", inherit.aes = FALSE) + labs(title = paste0("Finlay & Wilkinson analysis for ", trait_used, " (LOO index)"), x = "LOO environment index", y = "BLUP", color = "Genotype") + scale_x_continuous(expand = expansion(mult = c(0.05, 0.1))) + coord_cartesian(clip = "off") + theme_bw() + theme(plot.margin = margin(t = 8, r = 8, b = 28, l = 8))
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
    Data_source = "LMM predicted BLUP matrix",
    Imputed_cell_warning = "AMMI/GGE uses a full BLUP GxE matrix; unobserved genotype-location cells are model-predicted.",
    N_genotypes_enter_AMMI_GGE = nrow(GxE_obs_wide),
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
  p_ammi1 <- empty_plot("AMMI needs at least 2 genotypes and 2 environments.")
  p_ammi2 <- empty_plot("AMMI2 needs at least 2 PCs.")
  p_gge <- empty_plot("GGE needs at least 2 genotypes and 2 environments.")
  p_env_cor <- empty_plot("Environment correlation needs at least 2 environments.")
  if (nrow(GxE_obs_wide) >= 2 && ncol(GxE_obs_wide) >= 2) {
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
      labs(title = paste0("AMMI1 biplot - ", trait_used, " BLUPs"), x = paste0("IPCA1 (", round(PC_pct[1] * 100, 1), "%)"), y = paste0("Mean BLUP for ", trait_used), color = "Coverage", shape = "Coverage") +
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
        labs(title = paste0("AMMI2 biplot - ", trait_used, " BLUPs"), x = paste0("IPCA1 (", round(PC_pct[1] * 100, 1), "%)"), y = paste0("IPCA2 (", round(PC_pct[2] * 100, 1), "%)"), color = "Coverage", shape = "Coverage") +
        theme_bw()
    }
    GGE_centered <- sweep(as.matrix(GxE_obs_wide), 2, col_means)
    svd_gge <- svd(GGE_centered)
    n_pc_gge <- min(nrow(GGE_centered) - 1, ncol(GGE_centered) - 1)
    GGE_pct <- svd_gge$d[1:n_pc_gge]^2 / sum(svd_gge$d[1:n_pc_gge]^2)
    GGE_geno <- as.data.frame(svd_gge$u[, 1:min(2, n_pc_gge), drop = FALSE]) %>% setNames(paste0("PC", 1:min(2, n_pc_gge))) %>% mutate(Genotype = rownames(GxE_obs_wide), .before = 1)
    GGE_env <- as.data.frame(svd_gge$v[, 1:min(2, n_pc_gge), drop = FALSE]) %>% setNames(paste0("PC", 1:min(2, n_pc_gge))) %>% mutate(Environment = colnames(GxE_obs_wide), .before = 1)
    for (i in 1:min(2, n_pc_gge)) {
      GGE_geno[[paste0("PC", i)]] <- GGE_geno[[paste0("PC", i)]] * svd_gge$d[i]
      GGE_env[[paste0("PC", i)]] <- GGE_env[[paste0("PC", i)]] * svd_gge$d[i]
    }
    gge_pc1_cor <- suppressWarnings(cor(GGE_geno$PC1, row_means[match(GGE_geno$Genotype, names(row_means))], use = "complete.obs"))
    if (n_pc_gge >= 1 && is.finite(gge_pc1_cor) && gge_pc1_cor < 0) {
      GGE_geno$PC1 <- -GGE_geno$PC1
      GGE_env$PC1 <- -GGE_env$PC1
    }
    GGE_geno <- GGE_geno %>% left_join(biplot_confidence, by = "Genotype")
    if (n_pc_gge >= 2) {
      scale_gge <- ifelse(max(abs(GGE_env$PC1), na.rm = TRUE) == 0, 1, max(abs(GGE_geno$PC1), na.rm = TRUE) / max(abs(GGE_env$PC1), na.rm = TRUE) * 0.7)
      GGE_env_sc <- GGE_env %>% mutate(PC1 = PC1 * scale_gge, PC2 = PC2 * scale_gge, label_x = PC1 * 1.15, label_y = PC2 * 1.15)
      p_gge <- ggplot() +
        geom_segment(data = GGE_env_sc, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.25, "cm"), type = "closed"), color = "#E74C3C", linewidth = 0.8) +
        geom_text(data = GGE_env_sc, aes(x = label_x, y = label_y, label = Environment), color = "#E74C3C", size = 3.5, fontface = "bold") +
        geom_point(data = GGE_geno, aes(x = PC1, y = PC2, color = Confidence_flag, shape = Confidence_flag), size = 2.5) +
        geom_text(data = GGE_geno, aes(x = PC1, y = PC2, label = Genotype, color = Confidence_flag), vjust = -0.8, size = 2.8, show.legend = FALSE) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
        scale_color_manual(values = c("OK" = "#2C3E50", "Low confidence" = "#D35400"), breaks = "Low confidence", na.value = "gray60") +
        scale_shape_manual(values = c("OK" = 16, "Low confidence" = 1), breaks = "Low confidence", na.value = 16) +
        labs(title = paste0("GGE biplot - ", trait_used, " BLUPs"), x = paste0("PC1 (", round(GGE_pct[1] * 100, 1), "%)"), y = paste0("PC2 (", round(GGE_pct[2] * 100, 1), "%)"), color = "Coverage", shape = "Coverage") +
        theme_bw()
    }
    cor_mat <- cor(GxE_obs_wide, use = "pairwise.complete.obs")
    cor_long <- as.data.frame(cor_mat) %>% rownames_to_column("Env1") %>% pivot_longer(-Env1, names_to = "Env2", values_to = "r")
    p_env_cor <- ggplot(cor_long, aes(x = Env1, y = Env2, fill = r)) + geom_tile(color = "white") + geom_text(aes(label = round(r, 2)), size = 4.5) + scale_fill_gradient2(low = "#E74C3C", mid = "white", high = "#2ECC71", midpoint = 0, limits = c(-1, 1)) + labs(title = "Genotype BLUP correlation across environments", x = NULL, y = NULL, fill = "r") + theme_bw()
  }
  selection <- build_met_selection_ranking(
    list(
      blups_main = BLUPs_main,
      fw_results = FW_results,
      ammi_genotype = AMMI_geno
    ),
    component_weights = c(mean = MET_W_YIELD, fw = MET_W_FW, asv = MET_W_ASV)
  )
  p_met_selection <- plot_met_selection_ranking(selection, trait_used)
  return(list(raw_data = df_raw, met_data = dat, met_cleaned_data = dat_clean, outlier_summary = outlier_summary, presence = presence, genotype_summary = genotype_summary, model_summary = model_summary, variance_components = variance_components, lrt_table = lrt_table, blups_main = BLUPs_main, blups_environment = BLUPs_env_full, gxe_matrix = GxE_matrix_wide %>% rownames_to_column("Genotype"), fw_results = FW_results, ammi_notes = ammi_notes, ammi_genotype = AMMI_geno, ammi_environment = AMMI_env, gge_genotype = GGE_geno, gge_environment = GGE_env, met_selection = selection, p_before = p_before, p_after = p_after, p_variance = p_variance, p_residual = p_residual, p_blup = p_blup, p_accuracy = p_accuracy, p_perf_heatmap = p_perf_heatmap, p_fw_mean_sens = p_fw_mean_sens, p_fw_regression = p_fw_regression, p_ammi1 = p_ammi1, p_ammi2 = p_ammi2, p_gge = p_gge, p_env_cor = p_env_cor, p_met_selection = p_met_selection))
}
run_met_all_traits <- function(
    df_raw,
    check_varieties = NULL,
    trait_cols = NULL,
    replication_col = NULL,
    block_col = NULL,
    min_envs_for_biplot = NULL) {
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
        min_envs_for_biplot = min_envs_for_biplot
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



# UI helpers
sidebar_radio_menu <- function(input_id, groups, selected, group_controls = list()) {
  tags$div(
    id = input_id,
    class = "form-group shiny-input-radiogroup shiny-input-container split-radio-menu",
    role = "radiogroup",
    lapply(names(groups), function(group_name) {
      choices <- groups[[group_name]]
      tags$div(
        class = "side-subpanel",
        tags$div(class = "side-subpanel-title", group_name),
        group_controls[[group_name]],
        lapply(names(choices), function(label) {
          value <- unname(choices[[label]])
          radio_id <- paste(input_id, value, sep = "_")
          tags$div(
            class = paste(
              "form-check",
              if (value %in% c("met_integrated", "met_integrated_plot")) {
                "integrated-choice"
              } else {
                ""
              }
            ),
            tags$input(
              id = radio_id,
              type = "radio",
              name = input_id,
              value = value,
              class = "form-check-input",
              checked = if (identical(value, selected)) "checked" else NULL
            ),
            tags$label(
              class = "form-check-label",
              `for` = radio_id,
              label
            )
          )
        })
      )
    })
  )
}

sidebar_button_menu <- function(input_id, choices, selected = NULL, style = "module", children = list()) {
  if (is.null(selected) || !selected %in% unname(choices)) {
    selected <- unname(choices)[1]
  }
  tags$div(
    id = input_id,
    class = paste(
      "form-group shiny-input-radiogroup shiny-input-container button-radio-menu",
      paste0(style, "-button-menu")
    ),
    role = "radiogroup",
    lapply(names(choices), function(label) {
      value <- unname(choices[[label]])
      radio_id <- paste(input_id, value, sep = "_")
      tags$div(
        class = "button-menu-item",
        tags$div(
          class = "form-check",
          tags$input(
            id = radio_id,
            type = "radio",
            name = input_id,
            value = value,
            class = "form-check-input",
            checked = if (identical(value, selected)) "checked" else NULL
          ),
          tags$label(
            class = "form-check-label",
            `for` = radio_id,
            label
          )
        ),
        children[[value]]
      )
    })
  )
}

sidebar_detail_panel <- function(title, choices = NULL, input_id = NULL, selected = NULL, controls = NULL, style = "detail") {
  tags$div(
    class = "side-subpanel detail-subpanel",
    if (!is.null(choices) && !is.null(input_id)) {
      selectInput(
        inputId = input_id,
        label = title,
        choices = choices,
        selected = selected,
        width = "100%"
      )
    } else {
      tags$div(class = "side-subpanel-title", title)
    },
    controls
  )
}

nav_step_title <- function(number, label) {
  tags$span(
    class = "step-label",
    tags$span(class = "step-number", number),
    tags$span(class = "step-text", label)
  )
}

panel_header <- function(title, subtitle = NULL) {
  card_header(
    tags$div(
      class = "panel-heading",
      tags$div(class = "panel-title", title),
      if (!is.null(subtitle)) {
        tags$div(class = "panel-subtitle", subtitle)
      }
    )
  )
}

export_row <- function(status_class, filename, output_id) {
  tags$div(
    class = "export-row",
    tags$div(
      class = "export-file",
      tags$span(class = paste("export-dot", status_class)),
      tags$span(filename)
    ),
    downloadButton(
      outputId = output_id,
      label = "EXPORT",
      class = "export-link"
    )
  )
}

# UI
ui <- page_navbar(
  title = tags$span(
    class = "app-brand",
    tags$span(
      class = "brand-mark",
      tags$span(class = "brand-bar bar-one"),
      tags$span(class = "brand-bar bar-two"),
      tags$span(class = "brand-bar bar-three")
    ),
    tags$span("Selection Analysis Pipeline")
  ),
  navbar_options = navbar_options(
    bg = "#FFFFFF",
    theme = "light",
    underline = FALSE
  ),
  theme = bs_theme(
    version = 5,
    bg = "#F8F6F0",
    fg = "#263123",
    primary = "#315F28"
  ),
  header = tags$style(HTML("
    body, .bslib-page-navbar {
      background: #F8F6F0 !important;
      color: #263123;
    }
    .navbar {
      background: #FFFFFF !important;
      border: 0;
      border-bottom: 1px solid #D9DEE7;
      border-radius: 0;
      box-shadow: none;
      margin: 0;
      padding: 7px 10px;
      position: sticky;
      top: 0;
      z-index: 1030;
    }
    .navbar > .container-fluid {
      max-width: none;
      padding: 0 4px;
    }
    .navbar-nav {
      width: auto;
      display: flex;
      flex-direction: row;
    }
    .navbar-nav .nav-item {
      flex: 0 0 auto;
    }
    .navbar-nav .nav-link {
      color: #52634E !important;
      font-size: 16px;
      padding: 10px 14px !important;
    }
    .navbar-nav .nav-link.active,
    .navbar-nav .nav-link:hover {
      color: #315F28 !important;
      font-weight: 700;
    }
    .app-brand {
      align-items: center;
      color: #315F28;
      display: inline-flex;
      font-size: 22px;
      font-weight: 700;
      gap: 10px;
      white-space: nowrap;
    }
    .brand-mark {
      align-items: flex-end;
      background: #315F28;
      border-radius: 6px;
      display: inline-flex;
      gap: 2px;
      height: 25px;
      justify-content: center;
      padding: 6px;
      width: 25px;
    }
    .brand-bar {
      border: 1px solid #FFFFFF;
      border-radius: 1px;
      display: inline-block;
      width: 3px;
    }
    .bar-one { height: 6px; }
    .bar-two { height: 10px; }
    .bar-three { height: 14px; }
    .navbar-brand {
      margin-right: 28px;
      padding: 0;
    }
    .container-fluid {
      max-width: none;
      width: 100%;
    }
    .bslib-page-navbar > .container-fluid {
      padding-left: 16px;
      padding-right: 16px;
    }
    .card {
      margin-bottom: 18px;
      border: 1px solid #E1DCCD;
      border-radius: 12px;
      background: #FFFFFF;
      box-shadow: 0 6px 18px rgba(78, 68, 42, 0.06);
      overflow: hidden;
    }
    .card-header {
      background: #FFFFFF;
      border-bottom: 1px solid #E5DFD0;
      padding: 14px 18px 12px;
    }
    .card-body {
      padding: 18px;
    }
    .panel-title {
      color: #273122;
      font-family: Georgia, serif;
      font-size: 17px;
      line-height: 1.15;
    }
    .panel-subtitle {
      color: #929789;
      font-size: 12px;
      margin-top: 2px;
    }
    .btn-primary,
    .btn-primary:hover,
    .btn-primary:focus {
      background-color: #315F28;
      border-color: #315F28;
    }
    #run_analysis {
      width: 100%;
      margin-top: 8px;
      font-weight: 700;
      padding: 10px 14px;
    }
    .small-note {
      color: #8C9187;
      font-size: 12px;
    }
    .control-section {
      border-bottom: 1px solid rgba(49, 95, 40, 0.20);
      margin: 0;
      padding: 10px 0;
    }
    .control-section:last-child {
      border-bottom: 0;
      margin-bottom: 0;
    }
    .control-label {
      color: #34402F;
      font-size: 13px;
      font-weight: 700;
      margin-bottom: 7px;
    }
    .analysis-choice .form-check {
      border: 0;
      border-radius: 0;
      margin: 0 0 1px;
      padding: 4px 0 4px 28px;
      background: transparent;
    }
    .analysis-choice .form-check:has(input:checked) {
      background: transparent;
    }
    .analysis-choice .form-check-label {
      display: flex;
      flex-direction: column;
      line-height: 1.2;
    }
    .analysis-name {
      font-weight: 700;
      color: #283525;
      font-size: 13px;
    }
    .analysis-description {
      color: #8E9487;
      font-size: 11px;
      margin-top: 3px;
    }
    .diagnostic-choice .form-check {
      align-items: center;
      display: flex;
      margin: 0 0 2px;
      padding: 3px 0;
    }
    .diagnostic-choice .form-check-input {
      margin: 0 9px 0 0;
    }
    .workspace-controls .control-section .form-group,
    .workspace-controls .control-section .shiny-input-container,
    .workspace-controls .control-section .shiny-options-group {
      margin-bottom: 0;
    }
    .workspace-controls .control-label {
      margin-bottom: 4px;
    }
    .analysis-status {
      background: #F4F7EE;
      border-radius: 8px;
      color: #4B6543;
      font-size: 12px;
      margin-top: 12px;
      padding: 9px 10px;
    }
    .full-window-workspace {
      height: calc(100vh - 112px);
      height: calc(100dvh - 112px);
      overflow: hidden;
    }
    .full-window-workspace > .bslib-grid {
      column-gap: 12px !important;
      grid-template-columns: minmax(0, 1fr) minmax(0, 3fr) !important;
      row-gap: 0 !important;
      height: 100%;
    }
    .full-window-workspace > .bslib-grid > * {
      grid-column: auto !important;
      margin-bottom: 0;
      min-height: 0;
    }
    .analyze-pane {
      height: 100%;
      max-height: 100%;
      min-height: 0;
      overflow: hidden;
      padding: 0;
    }
    .analyze-pane > .card {
      height: 100%;
      margin-bottom: 0;
      min-height: 0;
    }
    .analyze-left-pane > .card,
    .analyze-right-pane > .card {
      overflow-y: auto;
    }
    .analyze-left-pane .workspace-controls > .card-body {
      gap: 0;
      justify-content: flex-start;
    }
    .analyze-left-pane .workspace-controls .shiny-panel-conditional {
      flex: 0 0 auto !important;
    }
    pre.shiny-text-output {
      white-space: pre-wrap;
      margin: 0;
    }
    .split-radio-menu {
      width: 100%;
    }
    .side-subpanel {
      border: 0;
      border-bottom: 1px solid #E5DFD0;
      border-radius: 0;
      padding: 12px 0 14px;
      margin-bottom: 0;
      background-color: #FFFFFF;
    }
    .side-subpanel:last-child {
      border-bottom: 0;
    }
    .side-subpanel-title {
      color: #92916F;
      font-size: 11px;
      font-weight: 700;
      letter-spacing: .04em;
      margin-bottom: 8px;
      text-transform: uppercase;
    }
    .side-subpanel .form-check {
      display: flex;
      align-items: center;
      gap: 10px;
      margin-bottom: 4px;
      padding: 7px 9px;
      border-radius: 7px;
    }
    .side-subpanel .form-check:has(input:checked) {
      background: #EAF2E3;
      color: #315F28;
      font-weight: 700;
    }
    .side-subpanel .form-check-input {
      flex: 0 0 auto;
      margin-top: 0;
      border-color: #B8C3AD;
    }
    .side-subpanel .form-check-input:checked {
      background-color: #47743D;
      border-color: #47743D;
    }
    .side-subpanel .form-check-label {
      line-height: 1.25;
      font-size: 13px;
    }
    .side-subpanel .form-check:last-child {
      margin-bottom: 0;
    }
    .side-subpanel .integrated-choice {
      border-top: 1px solid #DEE2E6;
      margin-top: 18px;
      padding-top: 14px;
    }
    .side-subpanel .form-group {
      margin-bottom: 14px;
    }
    .button-radio-menu {
      width: 100%;
    }
    .button-radio-menu > .button-menu-item > .form-check {
      display: block;
      margin-bottom: 8px;
      padding: 0;
    }
    .button-radio-menu > .button-menu-item > .form-check > .form-check-input {
      opacity: 0;
      position: absolute;
      pointer-events: none;
    }
    .button-radio-menu > .button-menu-item > .form-check > .form-check-label {
      border: 1px solid #D9D1BF;
      border-radius: 8px;
      cursor: pointer;
      display: block;
      font-size: 13px;
      font-weight: 700;
      line-height: 1.2;
      padding: 9px 12px;
      text-align: center;
      width: 100%;
    }
    .module-button-menu > .button-menu-item > .form-check > .form-check-label {
      border-color: #315F28;
      color: #315F28;
      background: #FFFFFF;
    }
    .module-button-menu > .button-menu-item > .form-check:has(> input:checked) > .form-check-label {
      background: #315F28;
      border-color: #315F28;
      color: #FFFFFF;
    }
    .sub-button-menu > .button-menu-item > .form-check > .form-check-label {
      border-color: #8FB284;
      color: #315F28;
      background: #FFFFFF;
    }
    .sub-button-menu > .button-menu-item > .form-check:has(> input:checked) > .form-check-label {
      background: #EAF2E3;
      border-color: #4E7A43;
      color: #174313;
    }
    .detail-button-menu > .button-menu-item > .form-check > .form-check-label {
      border-color: #C9D5BE;
      color: #3D5038;
      background: #FFFFFF;
      text-align: left;
    }
    .detail-button-menu > .button-menu-item > .form-check:has(> input:checked) > .form-check-label {
      background: #EAF2E3;
      border-color: #4E7A43;
      color: #315F28;
    }
    .analysis-output-workspace .analysis-output-controls,
    .analysis-output-workspace .analysis-output-main {
      height: 100%;
    }
    .analysis-output-workspace .analysis-output-controls,
    .analysis-output-workspace .analysis-output-main,
    .data-workspace .workspace-controls,
    .data-workspace .data-preview-card {
      min-height: 0;
      overflow-y: auto;
    }
    .analysis-output-controls .module-button-menu > .button-menu-item > .form-check {
      align-items: center;
      display: flex;
      margin-bottom: 1px;
      padding-left: 8px;
    }
    .analysis-output-controls .module-button-menu > .button-menu-item > .form-check > .form-check-input {
      border-color: #B8C3AD;
      flex: 0 0 auto;
      margin: 0 10px 0 0;
      opacity: 1;
      pointer-events: auto;
      position: static;
    }
    .analysis-output-controls .module-button-menu > .button-menu-item > .form-check > .form-check-label {
      background: transparent;
      border: 0;
      color: #4F594B;
      flex: 1 1 auto;
      font-weight: 600;
      padding: 4px;
      text-align: left;
    }
    .analysis-output-controls .module-button-menu > .button-menu-item > .form-check:has(> input:checked),
    .analysis-output-controls .module-button-menu > .button-menu-item > .form-check:has(> input:checked) > .form-check-label {
      background: transparent;
      color: #315F28;
      font-weight: 800;
    }
    .analyze-preview-header {
      align-items: flex-start;
      display: flex;
      gap: 16px;
      justify-content: space-between;
      width: 100%;
    }
    .analyze-preview-header .form-group {
      margin: 0;
      min-width: 190px;
    }
    .analysis-output-controls .side-subpanel {
      border-bottom: 1px solid #E1DACB;
      margin: 0 14px;
      padding: 9px 0;
    }
    .analysis-output-controls .side-subpanel-title {
      margin-bottom: 4px;
    }
    .analysis-output-controls .side-subpanel .form-group {
      margin-bottom: 7px;
    }
    .analysis-output-controls .side-subpanel .form-check {
      padding-bottom: 4px;
      padding-top: 4px;
    }
    .analysis-output-controls .side-subpanel:last-child {
      border-bottom: 0;
    }
    .analysis-output-controls,
    .analysis-output-controls > .card-body,
    .workspace-controls,
    .workspace-controls > .card-body {
      background: #D4E1CF;
    }
    .analysis-output-controls .card-header,
    .workspace-controls .card-header {
      background: #C4D6BD;
    }
    .analysis-output-controls .side-subpanel {
      background: transparent;
      border-bottom-color: #BAC9B3;
    }
    .analysis-output-controls .side-subpanel-title {
      color: #526A4B;
    }
    .output-header-controls > .shiny-html-output:last-child {
      flex: 0 0 min(320px, 38vw);
      min-width: 210px;
    }
    .output-header-controls > .shiny-html-output:last-child:empty {
      display: none;
    }
    .output-header-controls .form-group {
      margin: 0;
      width: 100%;
    }
    .result-analysis-section > .form-group:first-child > label,
    .result-analysis-section > .shiny-input-container:first-child > label {
      color: #979172;
      font-size: 11px;
      font-weight: 800;
      letter-spacing: 0.03em;
      text-transform: uppercase;
    }
    .result-analysis-section .selectize-control,
    .result-analysis-section .form-group {
      width: 100%;
    }
    .slicer-button-menu > .button-menu-item > .form-check {
      align-items: center;
      display: flex;
      margin-bottom: 8px;
      padding-left: 8px;
    }
    .slicer-button-menu > .button-menu-item > .form-check > .form-check-input {
      border-color: #B8C3AD;
      flex: 0 0 auto;
      margin: 0 10px 0 0;
      opacity: 1;
      pointer-events: auto;
      position: static;
    }
    .slicer-button-menu > .button-menu-item > .form-check > .form-check-label {
      border: 0;
      border-radius: 8px;
      color: #4F594B;
      flex: 1 1 auto;
      font-weight: 500;
      padding: 8px 10px;
      text-align: left;
    }
    .slicer-button-menu > .button-menu-item > .form-check:has(> input:checked) {
      background: #EAF2E3;
      border-radius: 8px;
    }
    .slicer-button-menu > .button-menu-item > .form-check:has(> input:checked) > .form-check-label {
      color: #315F28;
      font-weight: 700;
    }
    .button-menu-child {
      margin: 0 0 8px 42px;
      background: #F3F8EF;
      border-left: 3px solid #315F28;
      padding: 8px 8px 2px 12px;
    }
    .module-button-menu > .button-menu-item > .button-menu-child {
      margin-left: 42px;
    }
    .sub-button-menu > .button-menu-item > .button-menu-child {
      margin-left: 34px;
    }
    .button-menu-child .side-subpanel,
    .button-menu-child .detail-subpanel {
      border-bottom: 0;
      padding: 0 0 8px;
    }
    .sidebar-empty-note {
      background: #F7F9F3;
      border: 1px solid #DDE5D6;
      border-radius: 8px;
      color: #67715E;
      font-size: 12px;
      padding: 10px 12px;
    }
    .chart-download-panel {
      background: #F7F9F3;
      border: 1px solid #DDE5D6;
      border-radius: 9px;
      margin: 14px auto 18px;
      max-width: 420px;
      padding: 12px;
      width: calc(100% - 28px);
    }
    .chart-options-subpanel {
      border-top: 1px solid #E5DFD0;
      margin-top: 0;
      padding-top: 9px;
    }
    .chart-options-subpanel .form-group {
      margin-bottom: 10px;
    }
    .met-weight-panel {
      border-top: 1px solid #E6E0D2;
      margin-top: 18px;
      padding-top: 16px;
    }
    .met-weight-title {
      color: #34402F;
      font-size: 13px;
      font-weight: 700;
      margin-bottom: 10px;
    }
    .met-weight-panel .form-group {
      margin-bottom: 0;
    }
    .met-weight-panel .irs--shiny .irs-bar,
    .met-weight-panel .irs--shiny .irs-single {
      background: #315F28;
      border-color: #315F28;
    }
    .met-weight-panel .irs--shiny .irs-handle {
      border-color: #315F28;
    }
    .chart-download-title {
      color: #34402F;
      font-size: 13px;
      font-weight: 700;
      margin-bottom: 9px;
    }
    .chart-download-panel .form-group {
      margin-bottom: 9px;
    }
    .chart-download-details {
      position: relative;
    }
    .chart-download-details > summary {
      align-items: center;
      background: #315F28;
      border-radius: 6px;
      color: #FFFFFF;
      cursor: pointer;
      display: flex;
      font-size: 18px;
      height: 40px;
      justify-content: center;
      list-style: none;
      padding: 0;
      width: 44px;
    }
    .chart-download-details > summary::-webkit-details-marker { display: none; }
    .chart-download-details[open] > summary {
      background: #244A1E;
    }
    .chart-download-popover {
      background: #D4E1CF;
      border: 1px solid rgba(49, 95, 40, 0.28);
      border-radius: 8px;
      box-shadow: 0 8px 24px rgba(30, 45, 25, 0.18);
      padding: 14px;
      position: absolute;
      right: 0;
      top: 48px;
      width: 300px;
      z-index: 1000;
    }
    .chart-header-actions {
      align-items: center;
      display: flex;
      gap: 10px;
    }
    .chart-header-actions .shiny-html-output {
      flex: 0 0 auto;
      min-width: 0;
    }
    .chart-header-actions .form-group {
      align-items: center;
      display: flex;
      gap: 7px;
      margin: 0;
      width: auto;
    }
    .chart-header-actions .control-label {
      margin: 0;
      white-space: nowrap;
    }
    .chart-header-actions .selectize-control,
    .chart-header-actions select.form-select {
      margin: 0;
      width: 170px;
    }
    #download_chart {
      background: #315F28;
      border: 0;
      color: #FFFFFF;
      font-size: 12px;
      font-weight: 700;
      margin-top: 4px;
      padding: 9px 12px;
      width: 100%;
    }
    .chart-preview-container {
      align-items: center;
      display: flex;
      justify-content: center;
      min-height: 300px;
      overflow: auto;
      padding: 18px;
      width: 100%;
    }
    .chart-preview-frame {
      flex: 0 0 auto;
      max-width: 100%;
    }
    .chart-preview-frame .shiny-plot-output {
      height: 100% !important;
      width: 100% !important;
    }
    table.dataTable thead th {
      color: #8D8D70;
      font-size: 11px;
      text-transform: uppercase;
    }
    table.dataTable tbody td {
      border-color: #ECE6D8 !important;
      font-size: 13px;
    }
    .dataTables_wrapper .dataTables_paginate .paginate_button.current {
      background: #47743D !important;
      border-color: #47743D !important;
      color: #FFFFFF !important;
    }
    .upload-note {
      background: #EEF4E7;
      border-radius: 9px;
      color: #45613F;
      margin-top: 14px;
      padding: 12px;
    }
    .export-panel {
      max-width: none;
      width: 100%;
    }
    .export-row {
      align-items: center;
      background: #FCFBF7;
      border: 1px solid #E2DCCC;
      border-radius: 8px;
      display: flex;
      justify-content: space-between;
      margin-bottom: 10px;
      padding: 9px 12px;
    }
    .export-file {
      align-items: center;
      display: flex;
      font-family: monospace;
      font-size: 12px;
      gap: 9px;
    }
    .export-dot {
      background: #C7992D;
      border-radius: 50%;
      height: 7px;
      width: 7px;
    }
    .export-dot.ready {
      background: #4C7C40;
    }
    .export-link {
      background: transparent !important;
      border: 0 !important;
      box-shadow: none !important;
      color: #315F28 !important;
      font-size: 11px;
      font-weight: 700;
      padding: 3px !important;
      text-decoration: none;
    }
    #download_all {
      background: #315F28;
      border: 0;
      color: #FFFFFF;
      display: block;
      font-weight: 700;
      margin-top: 14px;
      padding: 11px;
      text-align: center;
      width: 100%;
    }
    @media (max-width: 768px) {
      .navbar {
        margin: 0;
      }
      .full-window-workspace {
        height: auto;
        overflow: visible;
      }
      .full-window-workspace > .bslib-grid {
        grid-template-columns: minmax(0, 1fr) !important;
        row-gap: 12px !important;
      }
      .analyze-pane {
        height: auto;
        max-height: none;
        overflow-y: visible;
        padding-right: 0;
      }
    }
  ")),
  # Upload data navbar
  nav_panel(
    title = "Data",
    tags$div(
      class = "data-workspace full-window-workspace",
      layout_columns(
        col_widths = c(3, 9),
        card(
          class = "workspace-controls",
          panel_header("Add your data", "One Excel file, one trial"),
          fileInput(
            inputId = "excel_file",
            label = "Drop or choose an Excel file",
            accept = c(".xlsx", ".xls")
          ),
          tags$div(
            class = "upload-note",
            verbatimTextOutput("upload_message")
          )
        ),
        card(
          class = "data-preview-card",
          panel_header("Preview", "First rows of the uploaded workbook"),
          DTOutput("raw_table")
        )
      )
    )
  ),
  # Analysis navbar
  nav_panel(
    title = "Analyze",
    tags$div(
      class = "analyze-workspace full-window-workspace",
      layout_columns(
        col_widths = c(3, 9),
        tags$div(
          class = "analyze-pane analyze-left-pane",
          card(
            class = "workspace-controls",
            panel_header("Analyze", "Choose a trait and an analysis"),
            tags$div(
              class = "control-section",
              tags$div(class = "control-label", "Which trait?"),
              selectInput(
                inputId = "eval_trait",
                label = NULL,
                choices = NULL
              )
            ),
            tags$div(
              class = "control-section diagnostic-choice",
              tags$div(class = "control-label", "Check the model fit"),
              radioButtons(
                inputId = "diagnostic_plot_type",
                label = NULL,
                choices = c(
                  "Histogram" = "Residual histogram",
                  "Q-Q plot" = "QQ plot",
                  "Fitted vs actual" = "Observed vs fitted"
                ),
                selected = "QQ plot",
                inline = FALSE
              )
            ),
            tags$div(
              class = "control-section",
              tags$div(class = "control-label", "Choose the analysis"),
              selectInput(
                inputId = "analysis_method",
                label = NULL,
                choices = c(
                  "Breeding" = "BREEDING",
                  "Genetic Diversity Analysis" = "DIVERSITY",
                  "Mating" = "MATING",
                  "Selection Index" = "LPSI",
                  "Multi-Environment Trial" = "MET"
                ),
                selected = "BREEDING"
              )
            ),
            conditionalPanel(
              condition = "input.analysis_method == 'MATING'",
              tags$div(
                class = "control-section",
                selectInput(
                  inputId = "mating_design",
                  label = "Design",
                  choices = c(
                    "Griffing Method I" = "griffing_m1",
                    "Griffing Method II" = "griffing_m2",
                    "Griffing Method III" = "griffing_m3",
                    "Griffing Method IV" = "griffing_m4",
                    "Partial Diallel" = "diallel_partial",
                    "Line x Tester" = "line_tester"
                  ),
                  selected = "griffing_m1"
                ),
                uiOutput("mating_column_inputs")
              )
            ),
            conditionalPanel(
              condition = "input.analysis_method == 'BREEDING'",
              tags$div(
                class = "control-section",
                uiOutput("breeding_column_inputs")
              )
            ),
            conditionalPanel(
              condition = "input.analysis_method == 'LPSI'",
              tags$div(
                class = "control-section",
                uiOutput("lpsi_benchmark_check_inputs")
              )
            ),
            conditionalPanel(
              condition = "input.analysis_method == 'MET'",
              tags$div(
                class = "control-section",
                uiOutput("check_variety_inputs")
              )
            ),
            conditionalPanel(
              condition = "input.analysis_method == 'DIVERSITY'",
              tags$div(
                class = "control-section",
                uiOutput("diversity_column_inputs")
              )
            ),
            actionButton(
              inputId = "run_analysis",
              label = tagList("Run analysis", tags$span(" ->")),
              class = "btn-primary"
            ),
            tags$div(
              class = "analysis-status",
              verbatimTextOutput("analysis_status")
            )
          )
        ),
        tags$div(
          class = "analyze-pane analyze-right-pane",
          card(
            card_header(
              tags$div(
                class = "analyze-preview-header",
                uiOutput("analyze_preview_header"),
                selectInput(
                  "analyze_preview_view",
                  label = NULL,
                  choices = c("Preview chart" = "chart", "Data summary" = "summary"),
                  selected = "chart"
                )
              )
            ),
            conditionalPanel("input.analyze_preview_view == 'chart'", plotOutput("diagnostic_plot", height = "540px")),
            conditionalPanel("input.analyze_preview_view == 'summary'", DTOutput("shapiro_table"))
          )
        )
      )
    )
  ),
  # Result navbar
  nav_panel(
    title = "Results",
    tags$div(
      class = "results-workspace analysis-output-workspace full-window-workspace",
      layout_columns(
        col_widths = c(3, 9),
        card(
          class = "analysis-output-controls",
        panel_header("Results", "Grouped by analysis question"),
        uiOutput("result_sidebar_menu"),
        conditionalPanel(
          "input.result_module == 'met' && (input.result_view == 'met_selection' || input.result_view == 'met_integrated')",
          tags$div(
            class = "side-subpanel met-weight-panel",
            tags$div(class = "side-subpanel-title", "Weight"),
            sliderInput(
              inputId = "met_weight_mean",
              label = "Mean",
              min = 0,
              max = 100,
              value = 50,
              step = 5,
              post = "%"
            ),
            sliderInput(
              inputId = "met_weight_fw",
              label = "FW",
              min = 0,
              max = 100,
              value = 25,
              step = 5,
              post = "%"
            ),
            sliderInput(
              inputId = "met_weight_asv",
              label = "ASV",
              min = 0,
              max = 100,
              value = 25,
              step = 5,
              post = "%"
            )
          )
        )
      ),
        card(
          class = "analysis-output-main",
        card_header(
          tags$div(
            class = "analyze-preview-header output-header-controls",
            uiOutput("result_header"),
            uiOutput("result_header_trait_control")
          )
        ),
        conditionalPanel("input.result_view == 'mating_anova'", DTOutput("mating_anova_table")),
        conditionalPanel("input.result_view == 'mating_gca'", DTOutput("mating_gca_table")),
        conditionalPanel("input.result_view == 'mating_sca'", DTOutput("mating_sca_table")),
        conditionalPanel("input.result_view == 'mating_variance'", DTOutput("mating_variance_table")),
        conditionalPanel("input.result_view == 'breeding_stats'", DTOutput("breeding_stats_table")),
        conditionalPanel("input.result_view == 'breeding_response'", DTOutput("breeding_response_table")),
        conditionalPanel("input.result_view == 'breeding_realized'", DTOutput("breeding_realized_table")),
        conditionalPanel("input.result_view == 'breeding_generation'", DTOutput("breeding_generation_table")),
        conditionalPanel("input.result_view == 'lpsi_trait'", DTOutput("trait_table")),
        conditionalPanel("input.result_view == 'lpsi_ranking'", DTOutput("index_table")),
        conditionalPanel("input.result_view == 'lpsi_superiority'", DTOutput("superiority_table")),
        conditionalPanel("input.result_view == 'lpsi_anova'", DTOutput("anova_full_table")),
        conditionalPanel("input.result_view == 'lpsi_lsd'", DTOutput("lsd_wide_table")),
        conditionalPanel("input.result_view == 'lpsi_heritability'", DTOutput("heritability_gain_table")),
        conditionalPanel("input.result_view == 'lpsi_direct'", DTOutput("lpsi_direct_table")),
        conditionalPanel("input.result_view == 'lpsi_compare'", DTOutput("lpsi_compare_table")),
        conditionalPanel("input.result_view == 'met_summary'", DTOutput("met_summary_table")),
        conditionalPanel("input.result_view == 'met_qc'", DTOutput("met_qc_table")),
        conditionalPanel("input.result_view == 'met_variance'", DTOutput("met_variance_table")),
        conditionalPanel("input.result_view == 'met_blup'", DTOutput("met_blup_table")),
        conditionalPanel("input.result_view == 'met_fw'", DTOutput("met_fw_table")),
        conditionalPanel(
          "input.result_view == 'met_ammi'",
          tagList(DTOutput("met_ammi_notes_table"), DTOutput("met_ammi_table"))
        ),
        conditionalPanel(
          "input.result_view == 'met_gge'",
          tagList(DTOutput("met_gge_notes_table"), DTOutput("met_gge_table"))
        ),
        conditionalPanel("input.result_view == 'met_selection'", DTOutput("met_selection_table")),
        conditionalPanel("input.result_view == 'met_integrated'", DTOutput("met_integrated_table")),
        conditionalPanel("input.result_view == 'diversity_values'", DTOutput("diversity_values_table")),
        conditionalPanel("input.result_view == 'diversity_clusters'", DTOutput("diversity_clusters_table")),
        conditionalPanel("input.result_view == 'diversity_superiority'", DTOutput("diversity_superiority_table")),
        conditionalPanel("input.result_view == 'diversity_corr'", DTOutput("diversity_corr_table"))
        )
      )
    )
  ),
  # Charts navbar
  nav_panel(
    title = "Charts",
    tags$div(
      class = "charts-workspace analysis-output-workspace full-window-workspace",
      layout_columns(
        col_widths = c(3, 9),
        card(
          class = "analysis-output-controls",
        panel_header("Charts", "Visual patterns from completed analyses"),
        uiOutput("chart_sidebar_menu")
      ),
        card(
          class = "analysis-output-main",
        card_header(
          tags$div(
            class = "analyze-preview-header output-header-controls",
            tags$div(
              class = "panel-heading",
              tags$div(class = "panel-title", "Analysis chart"),
              tags$div(class = "panel-subtitle", "Live preview using the selected width and height")
            ),
            tags$div(
              class = "chart-header-actions",
              uiOutput("chart_header_trait_control"),
              tags$details(
                class = "chart-download-details",
                tags$summary(
                  title = "Download chart",
                  `aria-label` = "Open chart download options",
                  icon("download")
                ),
                tags$div(
                  class = "chart-download-panel chart-download-popover",
                  numericInput("chart_width", "Width (in)", value = 12, min = 4, max = 30, step = 0.5),
                  numericInput("chart_height", "Height (in)", value = 7, min = 4, max = 30, step = 0.5),
                  selectInput(
                    "chart_format",
                    "File format",
                    choices = c(
                      "PNG (high resolution)" = "png",
                      "PDF (vector)" = "pdf"
                    ),
                    selected = "png"
                  ),
                  downloadButton("download_chart", "DOWNLOAD CHART"),
                  tags$div(
                    class = "small-note",
                    "PNG is exported at 300 DPI. PDF stays sharp at any size."
                  )
                )
              )
            )
          )
        ),
        uiOutput("chart_preview_ui")
        )
      )
    )
  ),
  # Export navbar
  nav_panel(
    title = "Export",
    layout_columns(
      col_widths = c(12),
      card(
        class = "export-panel",
        panel_header("Export", "Take your completed result tables with you"),
        uiOutput("export_panel_rows"),
        downloadButton(
          outputId = "download_all",
          label = "Download everything (.zip)"
        )
      )
    )
  )
)



# Server
server <- function(input, output, session) {
  uploaded_data <- reactive({
    req(input$excel_file)
    df <- tryCatch({
      si_read_excel_upload(
        input$excel_file$datapath,
        input$excel_file$name
      )
    }, error = function(e) {
      validate(need(FALSE, e$message))
    })
    report <- si_validate_uploaded_file(
      df,
      remove_cols = remove_cols
    )
    validate(
      need(report$ok, paste(report$errors, collapse = "\n"))
    )
    attr(df, "si_validation_report") <- report
    df
  })
  diagnostic_data <- reactive({
    req(uploaded_data())
    tryCatch({
      make_diagnostic_data(uploaded_data())
    }, error = function(e) {
      NULL
    })
  })
  analysis_results <- reactiveVal(NULL)
  analysis_used <- reactiveVal(NULL)
  analysis_message <- reactiveVal("No analysis has been run yet.")
  saved_results <- reactiveValues(
    MATING = NULL,
    BREEDING = NULL,
    DIVERSITY = NULL,
    LPSI = NULL,
    MET = NULL
  )
  observeEvent(input$excel_file, {
    si_reset_analysis_state(
      analysis_results,
      analysis_used,
      analysis_message,
      saved_results,
      reason = "New file uploaded."
    )
  }, ignoreInit = TRUE)
  mating_result_for_table <- reactive({
    req(analysis_results())
    validate(need(
      analysis_used() == "MATING",
      "Run Mating Analysis to view this table."
    ))
    table_type <- sub("^mating_", "", input$result_view)
    table <- get_mating_result_table(analysis_results(), table_type)
    validate(need(
      !is.null(table) && nrow(as.data.frame(table)) > 0,
      paste("No", gsub("_", " ", table_type), "table is available for this mating design.")
    ))
    as.data.frame(table)
  })
  breeding_result <- reactive({
    req(analysis_results())
    validate(need(
      analysis_used() == "BREEDING",
      "Run Breeding Analysis to view this result."
    ))
    analysis_results()
  })
  lpsi_settings <- reactive({
    checks <- if (identical(analysis_used(), "LPSI")) {
      input$lpsi_benchmark_checks
    } else {
      input$lpsi_run_benchmark_checks
    }
    if (is.null(checks) || length(checks) == 0) {
      checks <- input$lpsi_run_benchmark_checks
    }
    list(
      checks = checks,
      advance = input$lpsi_advance_cutoff %||% advance_index_cutoff,
      retest = input$lpsi_retest_cutoff %||% retest_index_cutoff,
      priority = input$lpsi_priority_cutoff_pct %||% priority_advance_cutoff_pct,
      severe = input$lpsi_severe_weak_pct %||% priority_severe_weak_pct
    )
  })
  run_lpsi_with_settings <- function(settings = lpsi_settings()) {
    run_selection_pipeline(
      uploaded_data(),
      check_varieties = settings$checks,
      advance_cutoff = settings$advance,
      retest_cutoff = settings$retest,
      priority_cutoff_pct = settings$priority,
      severe_weak_pct = settings$severe
    )
  }
  lpsi_result <- reactive({
    result <- if (identical(analysis_used(), "LPSI")) {
      analysis_results()
    } else {
      saved_results$LPSI
    }
    validate(need(!is.null(result), "Run LPSI analysis to view this result."))
    result
  })
  observeEvent({
    list(
      input$lpsi_benchmark_checks,
      input$lpsi_run_benchmark_checks,
      input$lpsi_advance_cutoff,
      input$lpsi_retest_cutoff,
      input$lpsi_priority_cutoff_pct,
      input$lpsi_severe_weak_pct
    )
  }, {
    if (!identical(analysis_used(), "LPSI")) return(invisible(NULL))
    req(uploaded_data())
    analysis_message("Updating LPSI decision settings...")
    res <- tryCatch({
      run_lpsi_with_settings()
    }, error = function(e) {
      showNotification(
        paste("LPSI settings update failed:", e$message),
        type = "error",
        duration = NULL
      )
      NULL
    })
    if (!is.null(res)) {
      analysis_results(res)
      saved_results$LPSI <- res
      analysis_message("LPSI decision settings updated.")
    }
  }, ignoreInit = TRUE)
  met_result_for_table <- reactive({
    req(analysis_results())
    validate(need(analysis_used() == "MET", "Run MET analysis to view this table."))
    req(input$met_result_trait)
    result <- analysis_results()$met_by_trait[[input$met_result_trait]]
    validate(need(!is.null(result), "This MET trait result is not available."))
    result
  })
  met_result_for_plot <- reactive({
    validate(need(
      !is.null(analysis_results()) && identical(analysis_used(), "MET"),
      "Run MET analysis to view this plot."
    ))
    validate(need(!is.null(input$met_plot_trait) && input$met_plot_trait != "", "Choose a trait."))
    result <- analysis_results()$met_by_trait[[input$met_plot_trait]]
    validate(need(!is.null(result), "This MET trait plot is not available."))
    result
  })
  met_component_weights <- reactive({
    normalize_met_component_weights(
      input$met_weight_mean,
      input$met_weight_fw,
      input$met_weight_asv
    )
  })
  weighted_met_selection_for_table <- reactive({
    build_met_selection_ranking(
      met_result_for_table(),
      component_weights = met_component_weights()
    )
  })
  weighted_met_selection_for_plot <- reactive({
    result <- met_result_for_plot()
    build_met_selection_ranking(
      result,
      component_weights = met_component_weights()
    )
  })
  weighted_met_integrated <- reactive({
    req(analysis_results(), uploaded_data())
    validate(need(analysis_used() == "MET", "Run MET analysis to view this result."))
    build_met_integrated_ranking(
      uploaded_data(),
      analysis_results()$met_by_trait,
      component_weights = met_component_weights()
    )
  })
  breeder_recommendation <- reactive({
    lpsi_res <- if (!is.null(saved_results$LPSI)) saved_results$LPSI else NULL
    met_res <- if (!is.null(saved_results$MET)) saved_results$MET else NULL
    
    if (identical(analysis_used(), "LPSI")) {
      lpsi_res <- analysis_results()
    }
    if (identical(analysis_used(), "MET")) {
      met_res <- analysis_results()
      integrated <- weighted_met_integrated()
      met_res$met_integrated_ranking <- integrated$ranking
      met_res$met_integrated_trait_weights <- integrated$trait_weights
      met_res$met_integrated_adjusted <- integrated$adjusted_performance
      met_res$met_integrated_standardized <- integrated$standardized_scores
    }
    
    si_breeder_recommendation_table(
      lpsi_results = lpsi_res,
      met_results = met_res
    )
  })
  lpsi_direct_selection_r <- reactive({
    result <- lpsi_result()
    si_lpsi_direct_selection(
      result,
      primary_trait = input$lpsi_direct_trait %||% input$eval_trait,
      selection_pct = input$lpsi_selection_pct %||% 15
    )
  })
  lpsi_direct_selection_chart_r <- reactive({
    result <- lpsi_result()
    si_lpsi_direct_selection(
      result,
      primary_trait = input$lpsi_chart_trait %||% input$lpsi_direct_trait %||% input$eval_trait,
      selection_pct = input$lpsi_selection_pct %||% 15
    )
  })
  lpsi_method_comparison_r <- reactive({
    result <- lpsi_result()
    si_lpsi_method_comparison(
      result,
      primary_trait = input$lpsi_direct_trait %||% input$eval_trait,
      selection_pct = input$lpsi_selection_pct %||% 15
    )
  })
  diversity_result <- reactive({
    validate(need(
      !is.null(analysis_results()) && identical(analysis_used(), "DIVERSITY"),
      "Run Genetic Diversity to view this result."
    ))
    analysis_results()
  })
  get_diversity_selected_checks <- function(input_id = NULL) {
    result <- if (identical(analysis_used(), "DIVERSITY")) {
      analysis_results()
    } else {
      saved_results$DIVERSITY
    }
    if (is.null(result)) {
      checks <- unique(clean_text(input$diversity_benchmark_checks))
      return(checks[!is.na(checks) & checks != ""])
    }
    choices <- unique(clean_text(as.data.frame(result$genotype_values)$GEN))
    choices <- choices[!is.na(choices) & choices != ""]
    checks <- if (!is.null(input_id)) unique(clean_text(input[[input_id]])) else character(0)
    checks <- checks[!is.na(checks) & checks != "" & checks %in% choices]
    if (length(checks) == 0) {
      checks <- unique(clean_text(input$diversity_benchmark_checks))
    }
    checks <- checks[!is.na(checks) & checks != "" & checks %in% choices]
    if (length(checks) == 0) {
      checks <- unique(clean_text(result$selected_checks))
      checks <- checks[!is.na(checks) & checks != "" & checks %in% choices]
    }
    if (length(checks) == 0 && length(choices) > 0) {
      checks <- choices[1]
    }
    checks
  }
  diversity_result_selected_checks <- reactive({
    get_diversity_selected_checks("diversity_result_benchmark_checks")
  })
  diversity_chart_selected_checks <- reactive({
    get_diversity_selected_checks("diversity_chart_benchmark_checks")
  })
  diversity_superiority_result <- reactive({
    build_diversity_superiority_data(diversity_result(), diversity_result_selected_checks())
  })
  selected_chart <- reactive({
    req(input$plot_view)
    chart_module <- input$chart_module %||% "selection_index"
    if (chart_module == "mating") {
      validate(need(FALSE, "No chart view is available for this module yet."))
    }
    view <- input$plot_view
    if (chart_module == "selection_index") {
      lpsi_mode <- input$chart_lpsi_mode %||% "single"
      lpsi_views <- if (identical(lpsi_mode, "multi")) {
        c("lpsi_ranking_plot", "lpsi_heatmap")
      } else {
        c("lpsi_mean_comparison_plot", "lpsi_gain_curve", "lpsi_direct_plot")
      }
      if (!(view %in% lpsi_views)) {
        view <- unname(lpsi_views)[1]
      }
    }
    if (chart_module == "breeding" && !(view %in% c("breeding_trend", "breeding_gam", "breeding_h2_heatmap", "breeding_distribution"))) {
      view <- "breeding_trend"
    }
    if (chart_module == "met" && !(view %in% c("met_env_cor", "met_fw_plot", "met_fw_regression", "met_ammi1", "met_ammi2", "met_gge", "met_selection_plot", "met_integrated_plot"))) {
      view <- "met_env_cor"
    }
    if (chart_module == "diversity" && !(view %in% c("diversity_superiority_plot", "diversity_dendrogram_plot", "diversity_corr_heatmap_plot"))) {
      view <- "diversity_superiority_plot"
    }
    
    if (view %in% c("lpsi_direct_plot", "lpsi_ranking_plot", "lpsi_heatmap", "lpsi_gain_curve", "lpsi_mean_comparison_plot")) {
      lpsi <- lpsi_result()
      if (view == "lpsi_mean_comparison_plot") {
        comparison_trait <- input$lpsi_chart_trait
        validate(need(!is.null(comparison_trait) && comparison_trait != "", "Choose a trait for mean comparison."))
        return(list(
          plot = plot_lpsi_mean_comparison(lpsi, comparison_trait),
          name = paste0("LPSI_mean_comparison_", gsub("[^A-Za-z0-9_-]+", "_", comparison_trait))
        ))
      }
      if (view == "lpsi_direct_plot") {
        return(list(
          plot = plot_lpsi_direct_selection(
            lpsi_direct_selection_chart_r(),
            selection_pct = input$lpsi_selection_pct %||% 15,
            lpsi_results = lpsi
          ),
          name = "LPSI_single_trait_selection"
        ))
      }
      if (view == "lpsi_ranking_plot") {
        return(list(
          plot = lpsi$ranking_plot,
          name = "LPSI_multi_trait_selection"
        ))
      }
      if (view == "lpsi_gain_curve") {
        gain_trait <- input$lpsi_chart_trait
        validate(need(!is.null(gain_trait) && gain_trait != "", "Choose a trait for the genetic gain chart."))
        return(list(
          plot = plot_genetic_gain_curve(
            lpsi$heritability_gain,
            lpsi$trait_info,
            gain_trait,
            selection_pct = input$lpsi_selection_pct %||% 15
          ),
          name = paste0("LPSI_genetic_gain_", gsub("[^A-Za-z0-9_-]+", "_", gain_trait))
        ))
      }
      heatmap <- lpsi$heatmap_plot
      validate(need(!is.null(heatmap), "The superiority heatmap is not available."))
      return(list(
        plot = heatmap,
        name = "LPSI_superiority_heatmap"
      ))
    }
    
    if (view %in% c("breeding_trend", "breeding_gam", "breeding_h2_heatmap", "breeding_distribution")) {
      req(analysis_results())
      validate(need(analysis_used() == "BREEDING", "Run Breeding Analysis before downloading this chart."))
      trait_name <- input$breeding_plot_trait
      generation_col <- input$breeding_generation_col
      generation_stats <- analysis_results()$generation_stats
      
      if (view == "breeding_trend") {
        validate(need(
          !is.null(generation_stats) && nrow(generation_stats) > 0,
          "Genetic trend needs a generation/stage column."
        ))
        validate(need(!is.null(trait_name) && trait_name != "", "Choose a trait."))
        trend_data <- generation_stats %>% filter(Trait == trait_name)
        return(list(
          plot = breeding_plot_genetic_trend(trend_data),
          name = paste0("Breeding_genetic_trend_", gsub("[^A-Za-z0-9_-]+", "_", trait_name))
        ))
      }
      
      if (view == "breeding_gam") {
        return(list(
          plot = breeding_plot_gam(
            analysis_results()$genetic_stats,
            x_col = "Trait",
            title = "Genetic advance as percent of mean"
          ),
          name = "Breeding_GAM"
        ))
      }
      
      if (view == "breeding_h2_heatmap") {
        h2_data <- if (!is.null(generation_stats) && nrow(generation_stats) > 0) {
          generation_stats
        } else {
          analysis_results()$genetic_stats %>% mutate(Generation = "Overall")
        }
        return(list(
          plot = breeding_plot_heritability_heatmap(h2_data),
          name = "Breeding_heritability_heatmap"
        ))
      }
      
      validate(need(
        !is.null(generation_col) && generation_col != "",
        "Distribution shift needs a generation/stage column."
      ))
      validate(need(!is.null(trait_name) && trait_name != "", "Choose a trait."))
      return(list(
        plot = breeding_plot_distribution_shift(
          uploaded_data(),
          trait = trait_name,
          generation_col = generation_col
        ),
        name = paste0("Breeding_distribution_", gsub("[^A-Za-z0-9_-]+", "_", trait_name))
      ))
    }
    
    if (view %in% c("diversity_dendrogram_plot", "diversity_superiority_plot", "diversity_corr_heatmap_plot")) {
      result <- diversity_result()
      chart_details <- switch(
        view,
        diversity_dendrogram_plot = list(plot_diversity_dendrogram(result), "GDA_dendrogram"),
        diversity_superiority_plot = list(plot_diversity_superiority_heatmap(result, diversity_chart_selected_checks()), "GDA_superiority"),
        diversity_corr_heatmap_plot = list(plot_diversity_correlation_heatmap(result), "GDA_trait_correlation"),
        NULL
      )
      validate(need(!is.null(chart_details), "Select a GDA chart to download."))
      return(list(
        plot = chart_details[[1]],
        name = chart_details[[2]]
      ))
    }
    
    req(analysis_results())
    validate(need(analysis_used() == "MET", "Run MET analysis before downloading this chart."))
    if (view == "met_integrated_plot") {
      return(list(
        plot = weighted_met_integrated()$plot,
        name = "MET_overall_ranking"
      ))
    }
    
    result <- met_result_for_plot()
    chart_details <- switch(
      view,
      met_selection_plot = list(
        plot_met_selection_ranking(
          weighted_met_selection_for_plot(),
          input$met_plot_trait
        ),
        "MET_ranking"
      ),
      met_fw_plot = list(result$p_fw_mean_sens, "MET_FW_sensitivity"),
      met_fw_regression = list(result$p_fw_regression, "MET_FW_regression"),
      met_ammi1 = list(result$p_ammi1, "MET_AMMI1"),
      met_ammi2 = list(result$p_ammi2, "MET_AMMI2"),
      met_gge = list(result$p_gge, "MET_GGE"),
      met_env_cor = list(result$p_env_cor, "MET_environment_correlation"),
      NULL
    )
    validate(need(!is.null(chart_details), "Select a chart to download."))
    validate(need(!is.null(chart_details[[1]]), "The selected chart is not available."))
    trait_name <- gsub("[^A-Za-z0-9_-]+", "_", input$met_plot_trait)
    list(
      plot = chart_details[[1]],
      name = paste(chart_details[[2]], trait_name, sep = "_")
    )
  })
  draw_chart <- function(chart) {
    grid::grid.newpage()
    if (inherits(chart, c("grob", "gTree", "gtable"))) {
      grid::grid.draw(chart)
    } else {
      print(chart)
    }
  }
  output$chart_preview_ui <- renderUI({
    chart_module <- input$chart_module %||% "selection_index"
    if (chart_module == "mating") {
      return(tags$div(
        class = "chart-download-panel",
        "No chart view is available for this module yet."
      ))
    }
    view <- if (is.null(input$plot_view)) "lpsi_gain_curve" else input$plot_view
    if (chart_module == "selection_index") {
      lpsi_mode <- input$chart_lpsi_mode %||% "single"
      lpsi_views <- if (identical(lpsi_mode, "multi")) {
        c("lpsi_ranking_plot", "lpsi_heatmap")
      } else {
        c("lpsi_mean_comparison_plot", "lpsi_gain_curve", "lpsi_direct_plot")
      }
      if (!(view %in% lpsi_views)) {
        view <- unname(lpsi_views)[1]
      }
    }
    if (chart_module == "breeding" && !(view %in% c("breeding_trend", "breeding_gam", "breeding_h2_heatmap", "breeding_distribution"))) {
      view <- "breeding_trend"
    }
    if (chart_module == "met" && !(view %in% c("met_env_cor", "met_fw_plot", "met_fw_regression", "met_ammi1", "met_ammi2", "met_gge", "met_selection_plot", "met_integrated_plot"))) {
      view <- "met_env_cor"
    }
    if (chart_module == "diversity" && !(view %in% c("diversity_superiority_plot", "diversity_dendrogram_plot", "diversity_corr_heatmap_plot"))) {
      view <- "diversity_superiority_plot"
    }
    output_id <- switch(
      view,
      lpsi_direct_plot = "lpsi_direct_plot",
      lpsi_mean_comparison_plot = "lpsi_mean_comparison_plot",
      lpsi_ranking_plot = "ranking_plot",
      lpsi_heatmap = "heatmap_plot",
      lpsi_gain_curve = "gain_curve_plot",
      breeding_trend = "breeding_trend_plot",
      breeding_gam = "breeding_gam_plot",
      breeding_h2_heatmap = "breeding_h2_heatmap_plot",
      breeding_distribution = "breeding_distribution_plot",
      met_selection_plot = "met_selection_plot",
      met_fw_plot = "met_fw_plot",
      met_fw_regression = "met_fw_regression_plot",
      met_ammi1 = "met_ammi1_plot",
      met_ammi2 = "met_ammi2_plot",
      met_gge = "met_gge_plot",
      met_env_cor = "met_env_cor_plot",
      met_integrated_plot = "met_integrated_plot",
      diversity_dendrogram_plot = "diversity_dendrogram_plot",
      diversity_superiority_plot = "diversity_superiority_plot",
      diversity_corr_heatmap_plot = "diversity_corr_heatmap_plot",
      "ranking_plot"
    )
    
    width_in <- suppressWarnings(as.numeric(input$chart_width))
    height_in <- suppressWarnings(as.numeric(input$chart_height))
    if (!is.finite(width_in)) width_in <- 12
    if (!is.finite(height_in)) height_in <- 7
    width_in <- min(max(width_in, 4), 30)
    height_in <- min(max(height_in, 4), 30)
    
    preview_scale <- min(70, 1100 / width_in, 580 / height_in)
    preview_width <- round(width_in * preview_scale)
    
    tags$div(
      class = "chart-preview-container",
      tags$div(
        class = "chart-preview-frame",
        style = sprintf(
          "width: %dpx; aspect-ratio: %s / %s;",
          preview_width,
          format(width_in, scientific = FALSE, trim = TRUE),
          format(height_in, scientific = FALSE, trim = TRUE)
        ),
        plotOutput(output_id, width = "100%", height = "100%")
      )
    )
  })
  output$upload_message <- renderPrint({
    req(uploaded_data())
    data <- uploaded_data()
    diag <- tryCatch({
      make_diagnostic_data(data)
    }, error = function(e) {
      NULL
    })
    cat("Upload successful\n")
    cat("Number of rows:", nrow(data), "\n")
    cat("Number of columns:", ncol(data), "\n")
    if (!is.null(diag)) {
      cat("Number of traits observed:", length(diag$trait_cols), "\n")
      cat("Replication number:", n_distinct(diag$data$Rep), "\n")
    } else {
      cat("Trait and replication information could not be detected.\n")
      cat("Choose an analysis to see the required columns for that workflow.\n")
    }
  })
  output$diagnostic_header <- renderUI({
    trait <- if (is.null(input$eval_trait) || input$eval_trait == "") {
      "selected trait"
    } else {
      input$eval_trait
    }
    tags$div(
      class = "panel-heading",
      tags$div(class = "panel-title", paste0("Is ", trait, " well behaved?")),
      tags$div(
        class = "panel-subtitle",
        "Residual diagnostic using variety and replication when available"
      )
    )
  })
  output$analyze_preview_header <- renderUI({
    if (identical(input$analyze_preview_view %||% "chart", "summary")) {
      return(tags$div(
        class = "panel-heading",
        tags$div(class = "panel-title", "Data summary by trait"),
        tags$div(class = "panel-subtitle", "Residual normality test for every measured trait")
      ))
    }
    trait <- input$eval_trait %||% "selected trait"
    if (trait == "") trait <- "selected trait"
    tags$div(
      class = "panel-heading",
      tags$div(class = "panel-title", paste0("Is ", trait, " well behaved?")),
      tags$div(class = "panel-subtitle", "Residual diagnostic using variety and replication when available")
    )
  })
  output$result_header <- renderUI({
    view <- input$result_view
    if (is.null(view) || view == "") {
      view <- "mating_anova"
    }
    trait <- if (startsWith(view, "mating_")) {
      input$mating_trait_col
    } else if (startsWith(view, "breeding_")) {
      input$breeding_plot_trait
    } else if (startsWith(view, "met_")) {
      input$met_result_trait
    } else {
      input$eval_trait
    }
    if (is.null(trait) || trait == "") {
      trait <- "selected trait"
    }
    
    title <- switch(
      view,
      mating_anova = "ANOVA",
      mating_gca = "GCA parent effects",
      mating_sca = "SCA cross effects",
      mating_variance = "Variance breakdown",
      breeding_stats = "Breeding genetic parameters",
      breeding_response = "Response per year",
      breeding_realized = "Realized gain",
      breeding_generation = "Generation summary",
      lpsi_trait = "Trait summary",
      lpsi_anova = "LPSI analysis of variance",
      lpsi_lsd = "Mean comparison",
      lpsi_superiority = "Trait superiority vs mean check benchmark",
      lpsi_heritability = "Heritability and genetic gain",
      lpsi_ranking = "Recommendation",
      lpsi_direct = "Direct selection on primary trait",
      lpsi_compare = "Selection method comparison",
      met_summary = paste("MET summary -", trait),
      met_qc = paste("MET quality control -", trait),
      met_variance = paste("Variance -", trait),
      met_blup = paste("BLUP -", trait),
      met_fw = paste("Finlay-Wilkinson stability -", trait),
      met_ammi = paste("AMMI -", trait),
      met_gge = paste("GGE -", trait),
      met_selection = paste("Selection -", trait),
      met_integrated = "Integrated MET selection",
      diversity_values = "Genotypic values for diversity",
      diversity_clusters = "Diversity clusters",
      diversity_superiority = "Trait superiority",
      diversity_corr = "Trait correlation",
      "Analysis results"
    )
    tags$div(
      class = "panel-heading",
      tags$div(class = "panel-title", title),
      tags$div(class = "panel-subtitle", "Results from the selected analysis view")
    )
  })
  output$result_header_trait_control <- renderUI({
    module <- input$result_module %||% "mating"
    if (identical(analysis_used(), "MET")) {
      result <- analysis_results()
      traits <- if (!is.null(result)) as.character(result$met_trait_names) else character(0)
      if (length(traits) == 0) return(NULL)
      return(tags$div(
        class = "chart-header-actions",
        selectInput("met_result_trait", "Trait", choices = traits, selected = traits[1])
      ))
    }
    if (identical(module, "selection_index") && identical(input$result_lpsi_mode %||% "single", "single")) {
      data <- safe_uploaded_data()
      prepared <- tryCatch(prepare_excel_input(data), error = function(e) NULL)
      if (is.null(prepared) || length(prepared$trait_cols) == 0) return(NULL)
      selected <- input$lpsi_direct_trait %||% input$eval_trait %||% prepared$trait_cols[1]
      if (!selected %in% prepared$trait_cols) selected <- prepared$trait_cols[1]
      return(selectInput("lpsi_direct_trait", NULL, choices = prepared$trait_cols, selected = selected))
    }
    NULL
  })
  output$chart_header_trait_control <- renderUI({
    module <- input$chart_module %||% "selection_index"
    if (identical(module, "breeding")) {
      traits <- safe_diagnostic_traits()
      if (length(traits) == 0) return(NULL)
      selected <- input$breeding_plot_trait
      if (is.null(selected) || !selected %in% traits) selected <- traits[1]
      return(selectInput("breeding_plot_trait", "Trait", choices = traits, selected = selected))
    }
    if (identical(module, "met")) {
      traits <- safe_met_traits()
      if (length(traits) == 0) return(NULL)
      selected <- input$met_plot_trait
      if (is.null(selected) || !selected %in% traits) selected <- traits[1]
      return(selectInput("met_plot_trait", "Trait", choices = traits, selected = selected))
    }
    if (identical(module, "selection_index") && identical(input$chart_lpsi_mode %||% "single", "single")) {
      return(uiOutput("lpsi_chart_trait_control"))
    }
    NULL
  })
  output$export_panel_rows <- renderUI({
    tagList(
      export_row(
        if (is.null(saved_results$MATING)) "pending" else "ready",
        "Mating_analysis_results.xlsx",
        "download_mating"
      ),
      export_row(
        if (is.null(saved_results$BREEDING)) "pending" else "ready",
        "Breeding_analysis_results.xlsx",
        "download_breeding"
      ),
      export_row(
        if (is.null(saved_results$DIVERSITY)) "pending" else "ready",
        "Genetic_diversity_results.xlsx",
        "download_diversity"
      ),
      export_row(
        if (is.null(saved_results$LPSI)) "pending" else "ready",
        "LPSI_selection_results.xlsx",
        "download_lpsi"
      ),
      export_row(
        if (is.null(saved_results$MET)) "pending" else "ready",
        "MET_across_locations_results.xlsx",
        "download_met"
      )
    )
  })
  output$raw_table <- renderDT({
    req(uploaded_data())
    datatable(
      uploaded_data(),
      options = list(pageLength = 10, scrollX = TRUE)
    )
  })
  safe_uploaded_data <- function() {
    tryCatch(uploaded_data(), error = function(e) NULL)
  }
  safe_diagnostic_traits <- function() {
    data <- safe_uploaded_data()
    if (is.null(data)) return(character(0))
    tryCatch(make_diagnostic_data(data)$trait_cols, error = function(e) character(0))
  }
  safe_met_traits <- function() {
    if (!is.null(analysis_results()) && identical(analysis_used(), "MET")) {
      completed_traits <- analysis_results()$met_trait_names
      if (!is.null(completed_traits) && length(completed_traits) > 0) {
        return(as.character(completed_traits))
      }
    }
    if (!is.null(saved_results$MET)) {
      completed_traits <- saved_results$MET$met_trait_names
      if (!is.null(completed_traits) && length(completed_traits) > 0) {
        return(as.character(completed_traits))
      }
    }
    selected_traits <- input$met_trait_cols
    if (!is.null(selected_traits) && length(selected_traits) > 0) {
      excluded_columns <- c(input$met_rep_col, input$met_block_col)
      selected_traits <- setdiff(as.character(selected_traits), excluded_columns)
      if (length(selected_traits) > 0) return(selected_traits)
    }
    data <- safe_uploaded_data()
    if (is.null(data)) return(character(0))
    tryCatch(get_met_trait_cols(data), error = function(e) character(0))
  }
  output$result_sidebar_menu <- renderUI({
    tagList(
      tags$div(
        style = "display:none;",
        selectInput(
          "result_module", NULL,
          choices = c("mating", "breeding", "diversity", "selection_index", "met"),
          selected = "mating"
        ),
        selectInput(
          "result_view", NULL,
          choices = c(
            "mating_anova", "mating_gca", "mating_sca", "mating_variance",
            "breeding_stats", "breeding_response", "breeding_realized", "breeding_generation",
            "diversity_values", "diversity_clusters", "diversity_superiority", "diversity_corr",
            "lpsi_direct", "lpsi_compare", "lpsi_trait", "lpsi_anova", "lpsi_lsd",
            "lpsi_superiority", "lpsi_heritability", "lpsi_ranking",
            "met_summary", "met_qc", "met_variance", "met_blup", "met_fw",
            "met_ammi", "met_gge", "met_selection", "met_integrated"
          ),
          selected = "mating_anova"
        )
      ),
      tags$div(
        class = "side-subpanel result-analysis-section",
        selectInput(
          "result_mating_view", "WHICH CROSS BEST",
          choices = c(
            "ANOVA" = "mating_anova",
            "GCA - parent effects" = "mating_gca",
            "SCA - cross effects" = "mating_sca",
            "Variance breakdown" = "mating_variance"
          ),
          selected = input$result_mating_view %||% "mating_anova"
        )
      ),
      tags$div(
        class = "side-subpanel result-analysis-section",
        selectInput(
          "result_breeding_view", "BREEDING PROCESS",
          choices = c(
            "Genetic parameters" = "breeding_stats",
            "Response per year" = "breeding_response",
            "Realized gain" = "breeding_realized",
            "Generation summary" = "breeding_generation"
          ),
          selected = input$result_breeding_view %||% "breeding_stats"
        )
      ),
      tags$div(
        class = "side-subpanel result-analysis-section",
        selectInput(
          "result_diversity_view", "GENETIC DIVERSITY",
          choices = c(
            "Genotypic values" = "diversity_values",
            "Clusters" = "diversity_clusters",
            "Superiority" = "diversity_superiority",
            "Trait correlation" = "diversity_corr"
          ),
          selected = input$result_diversity_view %||% "diversity_values"
        ),
        conditionalPanel(
          "input.result_diversity_view == 'diversity_superiority'",
          uiOutput("diversity_result_benchmark_controls")
        )
      ),
      tags$div(
        class = "side-subpanel result-analysis-section",
        selectInput(
          "result_lpsi_mode", "SELECTION INDEX (LPSI)",
          choices = c(
            "Single trait selection" = "single",
            "Multi trait selection" = "multi"
          ),
          selected = input$result_lpsi_mode %||% "single"
        ),
        uiOutput("result_lpsi_view_slicer"),
        uiOutput("result_lpsi_threshold_detail")
      ),
      tags$div(
        class = "side-subpanel result-analysis-section",
        selectInput(
          "result_met_view", "Choose result",
          choices = c(
            "Summary" = "met_summary", "QC" = "met_qc",
            "Variance" = "met_variance", "BLUP" = "met_blup",
            "Stability (FW)" = "met_fw", "AMMI" = "met_ammi",
            "GGE" = "met_gge", "Selection" = "met_selection",
            "Integrated Selection" = "met_integrated"
          ),
          selected = input$result_met_view %||% "met_summary"
        )
      )
    )
  })
  output$result_lpsi_view_slicer <- renderUI({
    mode <- input$result_lpsi_mode %||% "single"
    choices <- if (identical(mode, "multi")) {
      c(
        "Superiority" = "lpsi_superiority",
        "Heritability & gain" = "lpsi_heritability",
        "Recommendation" = "lpsi_ranking"
      )
    } else {
      c(
        "Summary" = "lpsi_trait",
        "ANOVA" = "lpsi_anova",
        "Mean comparison" = "lpsi_lsd",
        "Selection" = "lpsi_direct",
        "Compare with LPSI" = "lpsi_compare"
      )
    }
    current <- input$result_lpsi_view %||% ""
    selectInput(
      "result_lpsi_view", "Choose result",
      choices = choices,
      selected = if (current %in% choices) current else unname(choices)[1]
    )
  })
  observeEvent(input$result_mating_view, {
    updateSelectInput(session, "result_module", selected = "mating")
    updateSelectInput(session, "result_view", selected = input$result_mating_view)
  }, ignoreInit = TRUE)
  observeEvent(input$result_breeding_view, {
    updateSelectInput(session, "result_module", selected = "breeding")
    updateSelectInput(session, "result_view", selected = input$result_breeding_view)
  }, ignoreInit = TRUE)
  observeEvent(input$result_diversity_view, {
    updateSelectInput(session, "result_module", selected = "diversity")
    updateSelectInput(session, "result_view", selected = input$result_diversity_view)
  }, ignoreInit = TRUE)
  observeEvent(input$result_lpsi_view, {
    updateSelectInput(session, "result_module", selected = "selection_index")
    updateSelectInput(session, "result_view", selected = input$result_lpsi_view)
  }, ignoreInit = TRUE)
  observeEvent(input$result_met_view, {
    updateSelectInput(session, "result_module", selected = "met")
    updateSelectInput(session, "result_view", selected = input$result_met_view)
  }, ignoreInit = TRUE)
  # Hybrid result navigation: persistent analysis choices with contextual dropdowns.
  output$result_sidebar_menu <- renderUI({
    choices <- c(
      "Breeding" = "breeding",
      "Genetic Diversity Analysis" = "diversity",
      "Mating" = "mating",
      "Selection Index" = "selection_index",
      "Multi-Environment Trial" = "met"
    )
    module <- input$result_module %||% "mating"
    if (!(module %in% choices)) module <- "mating"

    tagList(
      tags$div(
        class = "side-subpanel main-analysis-slicer",
        tags$div(class = "side-subpanel-title", "Choose analysis"),
        sidebar_button_menu(
          input_id = "result_module",
          choices = choices,
          selected = module,
          style = "module"
        )
      ),
      switch(
        module,
        mating = uiOutput("result_mating_detail"),
        breeding = uiOutput("result_breeding_detail"),
        diversity = uiOutput("result_diversity_detail"),
        selection_index = tagList(
          tags$div(
            class = "side-subpanel result-analysis-section",
            tags$div(class = "side-subpanel-title", "Selection type"),
            sidebar_button_menu(
              input_id = "result_lpsi_mode",
              choices = c(
                "Single trait selection" = "single",
                "Multi-trait selection" = "multi"
              ),
              selected = input$result_lpsi_mode %||% "single",
              style = "module"
            )
          ),
          if (identical(input$result_lpsi_mode %||% "single", "multi")) {
            uiOutput("result_lpsi_multi_detail")
          } else {
            uiOutput("result_lpsi_single_detail")
          },
          uiOutput("result_lpsi_threshold_detail")
        ),
        met = uiOutput("result_met_detail")
      )
    )
  })
  output$chart_sidebar_menu <- renderUI({
    choices <- c(
      "Breeding" = "breeding",
      "Genetic Diversity Analysis" = "diversity",
      "Mating" = "mating",
      "Selection Index" = "selection_index",
      "Multi-Environment Trial" = "met"
    )
    module <- input$chart_module %||% "selection_index"
    if (!(module %in% choices)) {
      module <- "selection_index"
    }
    
    tagList(
      tags$div(
        class = "side-subpanel main-analysis-slicer",
        tags$div(class = "side-subpanel-title", "Choose analysis"),
        sidebar_button_menu(
          input_id = "chart_module",
          choices = choices,
          selected = module,
          style = "module"
        )
      ),
      switch(
        module,
        mating = uiOutput("chart_mating_detail"),
        breeding = uiOutput("chart_breeding_detail"),
        selection_index = uiOutput("chart_lpsi_detail"),
        met = uiOutput("chart_met_detail"),
        diversity = uiOutput("chart_diversity_detail")
      )
    )
  })
  output$result_mating_detail <- renderUI({
    choices <- c(
      "ANOVA" = "mating_anova",
      "GCA - parent effects" = "mating_gca",
      "SCA - cross effects" = "mating_sca",
      "Variance breakdown" = "mating_variance"
    )
    current <- input$result_view %||% ""
    sidebar_detail_panel(
      "Which cross best",
      choices = choices,
      input_id = "result_view",
      selected = if (current %in% choices) current else "mating_anova"
    )
  })
  output$result_breeding_detail <- renderUI({
    choices <- c(
      "Genetic parameters" = "breeding_stats",
      "Response per year" = "breeding_response",
      "Realized gain" = "breeding_realized",
      "Generation summary" = "breeding_generation"
    )
    current <- input$result_view %||% ""
    sidebar_detail_panel(
      "Breeding process",
      choices = choices,
      input_id = "result_view",
      selected = if (current %in% choices) current else "breeding_stats"
    )
  })
  output$result_lpsi_single_detail <- renderUI({
    choices <- c(
      "Summary" = "lpsi_trait",
      "ANOVA" = "lpsi_anova",
      "Mean comparison" = "lpsi_lsd",
      "Selection" = "lpsi_direct",
      "Compare with LPSI" = "lpsi_compare"
    )
    current <- input$result_view %||% ""
    sidebar_detail_panel(
      "Single trait selection",
      choices = choices,
      input_id = "result_view",
      selected = if (current %in% choices) current else "lpsi_direct",
      controls = uiOutput("lpsi_direct_controls")
    )
  })
  output$result_lpsi_multi_detail <- renderUI({
    choices <- c(
      "Superiority" = "lpsi_superiority",
      "Heritability & gain" = "lpsi_heritability",
      "Recommendation" = "lpsi_ranking"
    )
    current <- input$result_view %||% ""
    sidebar_detail_panel(
      "Multi trait selection",
      choices = choices,
      input_id = "result_view",
      selected = if (current %in% choices) current else "lpsi_superiority"
    )
  })
  output$result_lpsi_threshold_detail <- renderUI({
    sidebar_detail_panel(
      "Benchmark & thresholds",
      controls = uiOutput("lpsi_result_controls")
    )
  })
  output$result_met_detail <- renderUI({
    choices <- c(
      "Summary" = "met_summary",
      "QC" = "met_qc",
      "Variance" = "met_variance",
      "BLUP" = "met_blup",
      "Stability (FW)" = "met_fw",
      "AMMI" = "met_ammi",
      "GGE" = "met_gge",
      "Selection" = "met_selection",
      "Integrated Selection" = "met_integrated"
    )
    current <- input$result_view %||% ""
    met_traits <- safe_met_traits()
    met_selected <- input$met_result_trait
    if (is.null(met_selected) || !met_selected %in% met_traits) {
      met_selected <- if (length(met_traits) > 0) met_traits[1] else character(0)
    }
    sidebar_detail_panel(
      "Across Location",
      choices = choices,
      input_id = "result_view",
      selected = if (current %in% choices) current else "met_summary"
    )
  })
  output$result_diversity_detail <- renderUI({
    choices <- c(
      "Genotypic values" = "diversity_values",
      "Clusters" = "diversity_clusters",
      "Superiority" = "diversity_superiority",
      "Trait correlation" = "diversity_corr"
    )
    current <- input$result_view %||% ""
    sidebar_detail_panel(
      "Genetic Diversity",
      choices = choices,
      input_id = "result_view",
      selected = if (current %in% choices) current else "diversity_values",
      controls = conditionalPanel(
        "input.result_view == 'diversity_superiority'",
        uiOutput("diversity_result_benchmark_controls")
      )
    )
  })
  output$chart_mating_detail <- renderUI({
    tags$div(
      class = "side-subpanel detail-subpanel",
      tags$div(class = "side-subpanel-title", "Mating"),
      tags$div(class = "sidebar-empty-note", "No chart view is available for this module yet.")
    )
  })
  output$chart_breeding_detail <- renderUI({
    choices <- c(
      "Genetic trend" = "breeding_trend",
      "GAM" = "breeding_gam",
      "Heritability heatmap" = "breeding_h2_heatmap",
      "Distribution shift" = "breeding_distribution"
    )
    current <- input$plot_view %||% ""
    traits <- safe_diagnostic_traits()
    trait_selected <- input$breeding_plot_trait
    if (is.null(trait_selected) || !trait_selected %in% traits) {
      trait_selected <- if (length(traits) > 0) traits[1] else character(0)
    }
    sidebar_detail_panel(
      "Breeding",
      choices = choices,
      input_id = "plot_view",
      selected = if (current %in% choices) current else "breeding_trend"
    )
  })
  output$chart_lpsi_detail <- renderUI({
    mode <- input$chart_lpsi_mode %||% "single"
    mode_choices <- c(
      "Single trait selection" = "single",
      "Multi-trait selection" = "multi"
    )
    if (!(mode %in% mode_choices)) mode <- "single"
    current <- input$plot_view %||% ""
    choices <- if (identical(mode, "multi")) {
      c(
        "Ranking" = "lpsi_ranking_plot",
        "Superiority" = "lpsi_heatmap"
      )
    } else {
      c(
        "Mean comparison" = "lpsi_mean_comparison_plot",
        "Genetic gain" = "lpsi_gain_curve",
        "Ranking" = "lpsi_direct_plot"
      )
    }
    tagList(
      tags$div(
        class = "side-subpanel result-analysis-section",
        tags$div(class = "side-subpanel-title", "Selection type"),
        sidebar_button_menu(
          input_id = "chart_lpsi_mode",
          choices = mode_choices,
          selected = mode,
          style = "module"
        )
      ),
      sidebar_detail_panel(
        if (identical(mode, "multi")) "Multi-trait selection" else "Single trait selection",
        choices = choices,
        input_id = "plot_view",
        selected = if (current %in% choices) current else unname(choices)[1],
        style = "sub"
      )
    )
  })
  output$chart_met_detail <- renderUI({
    choices <- c(
      "Environment Correlation" = "met_env_cor",
      "FW Sensitivity" = "met_fw_plot",
      "FW Regression" = "met_fw_regression",
      "AMMI1" = "met_ammi1",
      "AMMI2" = "met_ammi2",
      "GGE" = "met_gge",
      "Ranking" = "met_selection_plot",
      "Integrated Ranking" = "met_integrated_plot"
    )
    current <- input$plot_view %||% ""
    traits <- safe_met_traits()
    trait_selected <- input$met_plot_trait
    if (is.null(trait_selected) || !trait_selected %in% traits) {
      trait_selected <- if (length(traits) > 0) traits[1] else character(0)
    }
    sidebar_detail_panel(
      "MET",
      choices = choices,
      input_id = "plot_view",
      selected = if (current %in% choices) current else "met_env_cor"
    )
  })
  output$chart_diversity_detail <- renderUI({
    choices <- c(
      "Superiority" = "diversity_superiority_plot",
      "Dendrogram" = "diversity_dendrogram_plot",
      "Trait correlation" = "diversity_corr_heatmap_plot"
    )
    current <- input$plot_view %||% ""
    sidebar_detail_panel(
      "GDA",
      choices = choices,
      input_id = "plot_view",
      selected = if (current %in% choices) current else "diversity_superiority_plot",
      controls = conditionalPanel(
        "input.plot_view == 'diversity_superiority_plot'",
        uiOutput("diversity_chart_benchmark_controls")
      ),
      style = "sub"
    )
  })
  render_diversity_benchmark_control <- function(input_id) {
    result <- if (identical(analysis_used(), "DIVERSITY")) {
      analysis_results()
    } else {
      saved_results$DIVERSITY
    }
    if (is.null(result) || is.null(result$genotype_values)) {
      return(NULL)
    }
    choices <- unique(clean_text(as.data.frame(result$genotype_values)$GEN))
    choices <- choices[!is.na(choices) & choices != ""]
    if (length(choices) == 0) {
      return(NULL)
    }
    selected <- unique(clean_text(input[[input_id]]))
    selected <- selected[!is.na(selected) & selected != "" & selected %in% choices]
    if (length(selected) == 0) {
      selected <- unique(clean_text(result$selected_checks))
      selected <- selected[!is.na(selected) & selected != "" & selected %in% choices]
    }
    if (length(selected) == 0) {
      selected <- choices[1]
    }
    selectizeInput(
      inputId = input_id,
      label = "Benchmark check",
      choices = choices,
      selected = selected,
      multiple = TRUE,
      options = list(placeholder = "Choose one or more check genotypes")
    )
  }
  output$diversity_result_benchmark_controls <- renderUI({
    render_diversity_benchmark_control("diversity_result_benchmark_checks")
  })
  output$diversity_chart_benchmark_controls <- renderUI({
    render_diversity_benchmark_control("diversity_chart_benchmark_checks")
  })
  observe({
    req(diagnostic_data())
    traits <- diagnostic_data()$trait_cols
    updateSelectInput(
      session = session,
      inputId = "eval_trait",
      choices = traits,
      selected = traits[1]
    )
    updateSelectInput(
      session = session,
      inputId = "lpsi_chart_trait",
      choices = traits,
      selected = traits[1]
    )
    updateSelectInput(
      session = session,
      inputId = "breeding_plot_trait",
      choices = traits,
      selected = traits[1]
    )
  })
  observe({
    req(uploaded_data())
    traits <- tryCatch({
      get_met_trait_cols(uploaded_data())
    }, error = function(e) {
      character(0)
    })
    updateSelectInput(
      session = session,
      inputId = "met_result_trait",
      choices = traits,
      selected = if (length(traits) > 0) traits[1] else character(0)
    )
    updateSelectInput(
      session = session,
      inputId = "met_plot_trait",
      choices = traits,
      selected = if (length(traits) > 0) traits[1] else character(0)
    )
  })
  output$check_variety_inputs <- renderUI({
    req(uploaded_data(), input$analysis_method)
    if (input$analysis_method != "MET") {
      return(NULL)
    }
    
    prepared <- tryCatch({
      prepare_met_trait_settings(uploaded_data())
    }, error = function(e) {
      NULL
    })
    validate(need(!is.null(prepared), "Upload valid MET data before choosing checks."))
    all_cols <- names(prepared$data)
    choose_optional_col <- function(current_value, candidates) {
      if (!is.null(current_value) && current_value %in% all_cols) {
        return(current_value)
      }
      matched <- candidates[candidates %in% all_cols]
      if (length(matched) > 0) matched[1] else ""
    }
    choose_required_col <- function(current_value, candidates) {
      selected <- choose_optional_col(current_value, candidates)
      if (selected != "") selected else all_cols[1]
    }
    check_choices <- unique(clean_text(prepared$data$Genotype))
    check_choices <- check_choices[!is.na(check_choices) & check_choices != ""]
    priority_traits <- prepared$trait_info$Trait[
      suppressWarnings(as.numeric(prepared$trait_info$Raw_weight)) >= priority_weight_cutoff
    ]
    selected_traits <- input$met_trait_cols
    if (is.null(selected_traits) || length(selected_traits) == 0) {
      selected_traits <- if (length(priority_traits) > 0) priority_traits else head(prepared$trait_cols, min(3, length(prepared$trait_cols)))
    }
    selected_traits <- selected_traits[selected_traits %in% prepared$trait_cols]
    if (length(selected_traits) == 0) selected_traits <- prepared$trait_cols[1]
    n_met_envs <- n_distinct(prepared$data$Environment)
    default_min_envs <- max(MET_MIN_ENVS_FOR_BIPLOT, ceiling(0.5 * n_met_envs))
    default_min_envs <- min(max(1, default_min_envs), max(1, n_met_envs))
    selected_min_envs <- suppressWarnings(as.integer(input$met_min_envs_for_biplot %||% default_min_envs))
    if (is.na(selected_min_envs) || selected_min_envs < 1 || selected_min_envs > n_met_envs) {
      selected_min_envs <- default_min_envs
    }
    
    tagList(
      selectizeInput(
        inputId = "met_trait_cols",
        label = "MET traits to run",
        choices = prepared$trait_cols,
        selected = selected_traits,
        multiple = TRUE,
        options = list(placeholder = "Choose priority traits for a faster decision run")
      ),
      selectInput(
        inputId = "met_rep_col",
        label = "Replication column",
        choices = all_cols,
        selected = choose_required_col(input$met_rep_col, met_rep_col_candidates)
      ),
      selectInput(
        inputId = "met_block_col",
        label = "Block column",
        choices = c("No block column" = "", all_cols),
        selected = choose_optional_col(input$met_block_col, met_block_col_candidates)
      ),
      numericInput(
        inputId = "met_min_envs_for_biplot",
        label = "Min observed locations for AMMI/GGE flag",
        value = selected_min_envs,
        min = 1,
        max = max(1, n_met_envs),
        step = 1
      ),
      selectizeInput(
        inputId = "met_reference_checks",
        label = "Reference check",
        choices = check_choices,
        selected = character(0),
        multiple = TRUE,
        options = list(placeholder = "Optional: choose one or more check varieties")
      )
    )
  })
  output$lpsi_benchmark_check_inputs <- renderUI({
    if (
      !identical(input$analysis_method, "LPSI") &&
      is.null(saved_results$LPSI) &&
      !identical(analysis_used(), "LPSI")
    ) {
      return(NULL)
    }
    req(uploaded_data())
    prepared <- tryCatch({
      prepare_excel_input(uploaded_data())
    }, error = function(e) {
      NULL
    })
    if (is.null(prepared)) {
      return(NULL)
    }
    
    choices <- unique(clean_text(prepared$data[[id_col]]))
    choices <- choices[!is.na(choices) & choices != ""]
    selected <- if (!is.null(input$lpsi_run_benchmark_checks) && length(input$lpsi_run_benchmark_checks) > 0) {
      input$lpsi_run_benchmark_checks
    } else if (!is.null(saved_results$LPSI) && length(saved_results$LPSI$selected_checks) > 0) {
      saved_results$LPSI$selected_checks
    } else if (length(prepared$detected_check_varieties) > 0) {
      prepared$detected_check_varieties
    } else {
      prepared$check_original_name
    }
    selected <- selected[selected %in% choices]
    if (length(selected) == 0) {
      selected <- prepared$check_original_name
    }
    
    tagList(
      selectizeInput(
        inputId = "lpsi_run_benchmark_checks",
        label = "Benchmark check",
        choices = choices,
        selected = selected,
        multiple = TRUE,
        options = list(placeholder = "Choose one or more check varieties")
      ),
      if (length(prepared$detected_check_varieties) > 0) {
        tags$small(
          class = "text-muted",
          paste0(
            "Default detected from ",
            prepared$type_col,
            ": ",
            paste(prepared$detected_check_varieties, collapse = ", ")
          )
        )
      }
    )
  })
  output$lpsi_direct_controls <- renderUI({
    if (is.null(saved_results$LPSI) && !identical(analysis_used(), "LPSI")) {
      return(NULL)
    }
    data <- safe_uploaded_data()
    if (is.null(data)) return(NULL)
    prepared <- tryCatch(prepare_excel_input(data), error = function(e) NULL)
    if (is.null(prepared)) return(NULL)
    
    tagList(
      numericInput(
        inputId = "lpsi_selection_pct",
        label = "Selection intensity (%)",
        value = input$lpsi_selection_pct %||% 15,
        min = 1,
        max = 50,
        step = 1
      )
    )
  })
  
  output$lpsi_chart_trait_control <- renderUI({
    if (is.null(saved_results$LPSI) && !identical(analysis_used(), "LPSI")) {
      return(NULL)
    }
    data <- safe_uploaded_data()
    if (is.null(data)) return(NULL)
    prepared <- tryCatch(prepare_excel_input(data), error = function(e) NULL)
    if (is.null(prepared)) return(NULL)
    
    selected <- input$lpsi_chart_trait %||% input$lpsi_direct_trait %||% input$eval_trait %||% prepared$trait_cols[1]
    if (!selected %in% prepared$trait_cols) {
      selected <- prepared$trait_cols[1]
    }
    selectInput(
      inputId = "lpsi_chart_trait",
      label = "Trait",
      choices = prepared$trait_cols,
      selected = selected
    )
  })
  
  output$lpsi_result_controls <- renderUI({
    if (is.null(saved_results$LPSI) && !identical(analysis_used(), "LPSI")) {
      return(NULL)
    }
    data <- safe_uploaded_data()
    if (is.null(data)) return(NULL)
    prepared <- tryCatch(prepare_excel_input(data), error = function(e) NULL)
    if (is.null(prepared)) return(NULL)
    
    choices <- unique(clean_text(prepared$data[[id_col]]))
    choices <- choices[!is.na(choices) & choices != ""]
    selected <- if (!is.null(saved_results$LPSI) && length(saved_results$LPSI$selected_checks) > 0) {
      saved_results$LPSI$selected_checks
    } else if (!is.null(input$lpsi_run_benchmark_checks) && length(input$lpsi_run_benchmark_checks) > 0) {
      input$lpsi_run_benchmark_checks
    } else if (length(prepared$detected_check_varieties) > 0) {
      prepared$detected_check_varieties
    } else {
      prepared$check_original_name
    }
    selected <- selected[selected %in% choices]
    if (length(selected) == 0) selected <- prepared$check_original_name
    decision_settings <- if (!is.null(saved_results$LPSI$decision_settings)) {
      saved_results$LPSI$decision_settings
    } else {
      data.frame(
        Advance_index_cutoff = advance_index_cutoff,
        Retest_index_cutoff = retest_index_cutoff,
        Priority_trait_cutoff_pct = priority_advance_cutoff_pct,
        Severe_weakness_cutoff_pct = priority_severe_weak_pct
      )
    }
    setting_value <- function(input_id, column, fallback) {
      current <- input[[input_id]]
      current <- suppressWarnings(as.numeric(current))
      if (length(current) > 0 && is.finite(current[1])) {
        return(current[1])
      }
      value <- suppressWarnings(as.numeric(decision_settings[[column]][1]))
      if (length(value) > 0 && is.finite(value[1])) value[1] else fallback
    }
    
    tagList(
      actionButton(
        inputId = "toggle_lpsi_thresholds",
        label = "Benchmark & threshold settings",
        class = "btn-outline-secondary"
      ),
      conditionalPanel(
        "input.toggle_lpsi_thresholds % 2 == 1",
        selectizeInput(
          inputId = "lpsi_benchmark_checks",
          label = "Benchmark check",
          choices = choices,
          selected = selected,
          multiple = TRUE,
          options = list(placeholder = "Choose one or more check varieties")
        ),
        numericInput(
          inputId = "lpsi_advance_cutoff",
          label = "Advance index cutoff",
          value = setting_value("lpsi_advance_cutoff", "Advance_index_cutoff", advance_index_cutoff),
          step = 0.05
        ),
        numericInput(
          inputId = "lpsi_retest_cutoff",
          label = "Retest index cutoff",
          value = setting_value("lpsi_retest_cutoff", "Retest_index_cutoff", retest_index_cutoff),
          step = 0.05
        ),
        numericInput(
          inputId = "lpsi_priority_cutoff_pct",
          label = "Priority trait pass cutoff (%)",
          value = setting_value("lpsi_priority_cutoff_pct", "Priority_trait_cutoff_pct", priority_advance_cutoff_pct),
          step = 1
        ),
        numericInput(
          inputId = "lpsi_severe_weak_pct",
          label = "Severe weakness cutoff (%)",
          value = setting_value("lpsi_severe_weak_pct", "Severe_weakness_cutoff_pct", priority_severe_weak_pct),
          step = 1
        )
      )
    )
  })
  output$diversity_column_inputs <- renderUI({
    req(uploaded_data())
    data <- uploaded_data()
    all_cols <- names(data)
    choose_default <- function(candidates, choices, fallback = 1, allow_empty = FALSE) {
      matched <- candidates[candidates %in% choices]
      if (length(matched) > 0) matched[1]
      else if (allow_empty) ""
      else choices[min(fallback, length(choices))]
    }
    genotype_default <- choose_default(c(id_col, "Genotype", "GEN", "genotype", "ID", "Hybrid", "Variety"), all_cols)
    genotype_col <- input$diversity_genotype_col
    if (is.null(genotype_col) || !genotype_col %in% all_cols) {
      genotype_col <- genotype_default
    }
    rep_default <- choose_default(c(rep_col, "Replication", "REP", "Block", "Replicate"), all_cols, allow_empty = TRUE)
    genotype_values_upper <- toupper(clean_text(data[[genotype_col]]))
    metadata_rows <- genotype_values_upper %in% c(weight_row_labels, direction_row_labels)
    trait_detection_data <- data[!metadata_rows, , drop = FALSE]
    numeric_cols <- all_cols[vapply(trait_detection_data, function(column) {
      sum(!is.na(to_number(column))) > 0
    }, logical(1))]
    trait_choices <- setdiff(numeric_cols, c(genotype_col, rep_default, remove_cols))
    if (length(trait_choices) == 0) trait_choices <- numeric_cols
    check_choices <- unique(clean_text(trait_detection_data[[genotype_col]]))
    check_choices <- check_choices[!is.na(check_choices) & check_choices != ""]
    type_col <- type_col_candidates[type_col_candidates %in% names(trait_detection_data)][1]
    detected_checks <- character(0)
    if (length(type_col) > 0 && !is.na(type_col)) {
      type_values <- normalize_type_label(trait_detection_data[[type_col]])
      detected_checks <- unique(clean_text(trait_detection_data[[genotype_col]][type_values %in% check_type_labels]))
      detected_checks <- detected_checks[!is.na(detected_checks) & detected_checks != ""]
    }
    selected_checks <- input$diversity_benchmark_checks
    if (is.null(selected_checks) || length(selected_checks) == 0) {
      selected_checks <- if (length(detected_checks) > 0) detected_checks else check_choices[1]
    }
    selected_checks <- selected_checks[selected_checks %in% check_choices]
    if (length(selected_checks) == 0 && length(check_choices) > 0) {
      selected_checks <- check_choices[1]
    }
    tagList(
      selectInput("diversity_genotype_col", "Genotype / variety column", all_cols, selected = genotype_col),
      selectInput("diversity_rep_col", "Replication column", c("No replication column" = "", all_cols), selected = rep_default),
      selectizeInput(
        "diversity_benchmark_checks",
        "Benchmark check",
        choices = check_choices,
        selected = selected_checks,
        multiple = TRUE,
        options = list(placeholder = "Choose one or more check genotypes")
      ),
      selectizeInput(
        "diversity_trait_cols",
        "Traits for GDA",
        choices = trait_choices,
        selected = head(trait_choices, min(4, length(trait_choices))),
        multiple = TRUE
      ),
      selectInput(
        "diversity_value_method",
        "Genotypic value method",
        choices = c("Arithmetic mean" = "mean", "BLUE" = "blue", "BLUP" = "blup"),
        selected = "mean"
      )
    )
  })
  output$mating_column_inputs <- renderUI({
    req(uploaded_data(), input$mating_design)
    data <- uploaded_data()
    all_cols <- names(data)
    numeric_cols <- all_cols[vapply(data, function(column) {
      sum(!is.na(to_number(column))) > 0
    }, logical(1))]
    trait_candidates <- setdiff(
      numeric_cols,
      c(
        "Rep", "Replication", "Block", "Male", "Female",
        "Parent1", "Parent2", "Parent 1", "Parent 2",
        "Line", "Tester", "Type", "type"
      )
    )
    if (length(trait_candidates) > 0) {
      numeric_cols <- trait_candidates
    }
    
    choose_default <- function(candidates, choices, fallback = 1) {
      matched <- candidates[candidates %in% choices]
      if (length(matched) > 0) matched[1] else choices[min(fallback, length(choices))]
    }
    
    rep_default <- choose_default(c("Rep", "Replication", "Block"), all_cols)
    trait_default <- if (length(numeric_cols) > 0) numeric_cols[1] else character(0)
    
    if (input$mating_design == "line_tester") {
      tagList(
        selectInput(
          "mating_line_col", "Line column", all_cols,
          selected = choose_default(c("Line", "line"), all_cols)
        ),
        selectInput(
          "mating_tester_col", "Tester column", all_cols,
          selected = choose_default(c("Tester", "tester"), all_cols, 2)
        ),
        selectInput(
          "mating_rep_col", "Replication column", all_cols,
          selected = rep_default
        ),
        selectInput(
          "mating_type_col", "Type column", all_cols,
          selected = choose_default(c("Type", "type"), all_cols)
        ),
        selectInput(
          "mating_trait_col", "Trait", numeric_cols,
          selected = trait_default
        ),
        tags$p(
          class = "small-note",
          "The Type column must label hybrid rows exactly as 'cross'; other rows are treated as parents."
        )
      )
    } else {
      tagList(
        selectInput(
          "mating_parent1_col", "Parent 1 / male column", all_cols,
          selected = choose_default(c("Male", "Parent1", "Parent 1"), all_cols)
        ),
        selectInput(
          "mating_parent2_col", "Parent 2 / female column", all_cols,
          selected = choose_default(c("Female", "Parent2", "Parent 2"), all_cols, 2)
        ),
        selectInput(
          "mating_rep_col", "Replication column", all_cols,
          selected = rep_default
        ),
        selectInput(
          "mating_trait_col", "Trait", numeric_cols,
          selected = trait_default
        )
      )
    }
  })
  output$breeding_column_inputs <- renderUI({
    req(uploaded_data())
    data <- uploaded_data()
    all_cols <- names(data)
    numeric_cols <- all_cols[vapply(data, function(column) {
      sum(!is.na(to_number(column))) > 0
    }, logical(1))]
    numeric_cols <- setdiff(
      numeric_cols,
      c("Rep", "Replication", "Block", "Cycle", "Generation", "Stage")
    )
    if (length(numeric_cols) == 0) {
      numeric_cols <- all_cols
    }
    
    choose_default <- function(candidates, choices, fallback = 1, allow_empty = FALSE) {
      matched <- candidates[candidates %in% choices]
      if (length(matched) > 0) {
        matched[1]
      } else if (allow_empty) {
        ""
      } else {
        choices[min(fallback, length(choices))]
      }
    }
    
    genotype_default <- choose_default(
      c(id_col, "Genotype", "genotype", "ID", "Hybrid", "Variety"),
      all_cols
    )
    rep_default <- choose_default(
      c(rep_col, "Replication", "Block", "Replicate"),
      all_cols,
      allow_empty = TRUE
    )
    generation_default <- choose_default(
      c("Generation", "Stage", "Cycle", "Selection_Cycle", "selection_cycle"),
      all_cols,
      allow_empty = TRUE
    )
    cycle_group_default <- choose_default(
      c("cycle_group", "Cycle_group", "CycleGroup", "Group", "Selection_group"),
      all_cols,
      allow_empty = TRUE
    )
    group_choices <- c("No realized-gain group" = "", all_cols)
    generation_choices <- c("No generation/stage column" = "", all_cols)
    rep_choices <- c("No replication column" = "", all_cols)
    
    cycle_label_ui <- NULL
    if (!is.null(input$breeding_cycle_group_col) && input$breeding_cycle_group_col != "") {
      labels <- unique(clean_text(data[[input$breeding_cycle_group_col]]))
      labels <- labels[!is.na(labels) & labels != ""]
      if (length(labels) > 0) {
        current_default <- if ("current_cycle" %in% labels) "current_cycle" else labels[1]
        check_default <- if ("check_prior_cycle" %in% labels) {
          "check_prior_cycle"
        } else if (length(labels) >= 2) {
          labels[2]
        } else {
          labels[1]
        }
        cycle_label_ui <- tagList(
          selectInput(
            "breeding_current_label", "Current-cycle label",
            choices = labels,
            selected = current_default
          ),
          selectInput(
            "breeding_check_label", "Prior/check label",
            choices = labels,
            selected = check_default
          ),
          checkboxInput(
            "breeding_higher_is_better",
            "Higher values are better for realized gain",
            value = TRUE
          )
        )
      }
    }
    
    tagList(
      selectInput(
        "breeding_genotype_col", "Genotype / variety column",
        all_cols,
        selected = genotype_default
      ),
      selectInput(
        "breeding_rep_col", "Replication column",
        rep_choices,
        selected = rep_default
      ),
      selectizeInput(
        "breeding_trait_cols", "Traits",
        choices = numeric_cols,
        selected = numeric_cols[1],
        multiple = TRUE,
        options = list(placeholder = "Choose one or more traits")
      ),
      selectInput(
        "breeding_generation_col", "Generation / stage column",
        generation_choices,
        selected = generation_default
      ),
      selectInput(
        "breeding_cycle_group_col", "Current vs check group",
        group_choices,
        selected = cycle_group_default
      ),
      cycle_label_ui,
      numericInput(
        "breeding_selection_pct",
        "Selection proportion (%)",
        value = 5,
        min = 0.1,
        max = 99,
        step = 0.5
      ),
      numericInput(
        "breeding_years_per_cycle",
        "Years per cycle",
        value = 3,
        min = 0.1,
        max = 20,
        step = 0.5
      )
    )
  })
  output$diagnostic_plot <- renderPlot({
    req(diagnostic_data())
    req(input$eval_trait)
    req(input$diagnostic_plot_type)
    make_diagnostic_plot(
      diag = diagnostic_data(),
      trait = input$eval_trait,
      plot_type = input$diagnostic_plot_type
    )
  })
  output$shapiro_table <- renderDT({
    req(diagnostic_data())
    summary_table <- make_shapiro_table(diagnostic_data()) %>%
      transmute(
        Trait,
        `Residual observations` = N_residual,
        `Shapiro W` = Shapiro_W,
        `P value` = p_value,
        `Shapiro result` = ifelse(
          Normality_note == "Yes",
          "Residuals are approximately normal",
          "Review residual distribution"
        )
      )
    datatable(
      summary_table,
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE)
    )
  })
  observeEvent(input$run_analysis, {
    req(uploaded_data())
    analysis_results(NULL)
    analysis_used(input$analysis_method)
    if (input$analysis_method == "MATING") {
      req(
        input$mating_design,
        input$mating_trait_col,
        input$mating_rep_col
      )
      if (input$mating_design == "line_tester") {
        req(input$mating_line_col, input$mating_tester_col, input$mating_type_col)
      } else {
        req(input$mating_parent1_col, input$mating_parent2_col)
      }
      
      analysis_message("Running mating analysis...")
      res <- tryCatch({
        run_mating_pipeline(
          df = uploaded_data(),
          design = input$mating_design,
          trait_col = input$mating_trait_col,
          replication_col = input$mating_rep_col,
          parent1_col = input$mating_parent1_col,
          parent2_col = input$mating_parent2_col,
          line_col = input$mating_line_col,
          tester_col = input$mating_tester_col,
          type_col = input$mating_type_col
        )
      }, error = function(e) {
        showNotification(
          paste("Mating analysis failed:", e$message),
          type = "error",
          duration = NULL
        )
        NULL
      })
      analysis_results(res)
      if (!is.null(res)) {
        saved_results$MATING <- res
        analysis_message("Mating analysis complete. Check the Result navbar.")
        showNotification("Mating analysis complete.", type = "message")
      } else {
        analysis_message("Mating analysis failed. Please check the error message.")
      }
    } else if (input$analysis_method == "BREEDING") {
      req(input$breeding_genotype_col, input$breeding_trait_cols)
      selection_proportion <- suppressWarnings(as.numeric(input$breeding_selection_pct) / 100)
      years_per_cycle <- suppressWarnings(as.numeric(input$breeding_years_per_cycle))
      if (!is.finite(selection_proportion) || selection_proportion <= 0 || selection_proportion >= 1) {
        selection_proportion <- 0.05
      }
      if (!is.finite(years_per_cycle) || years_per_cycle <= 0) {
        years_per_cycle <- 3
      }
      rep_col_used <- if (!is.null(input$breeding_rep_col) && input$breeding_rep_col != "") {
        input$breeding_rep_col
      } else {
        NULL
      }
      generation_col_used <- if (!is.null(input$breeding_generation_col) && input$breeding_generation_col != "") {
        input$breeding_generation_col
      } else {
        NULL
      }
      cycle_group_used <- if (!is.null(input$breeding_cycle_group_col) && input$breeding_cycle_group_col != "") {
        input$breeding_cycle_group_col
      } else {
        NULL
      }
      current_label <- if (!is.null(cycle_group_used)) input$breeding_current_label else NULL
      check_label <- if (!is.null(cycle_group_used)) input$breeding_check_label else NULL
      if (
        !is.null(cycle_group_used) &&
        (is.null(current_label) || current_label == "" || is.null(check_label) || check_label == "")
      ) {
        cycle_group_used <- NULL
        current_label <- NULL
        check_label <- NULL
      }
      
      analysis_message("Running Breeding Analysis...")
      res <- tryCatch({
        out <- breeding_run_gain_pipeline(
          data = uploaded_data(),
          traits = input$breeding_trait_cols,
          genotype = input$breeding_genotype_col,
          replication = rep_col_used,
          cycle_group_col = cycle_group_used,
          current_label = current_label,
          check_label = check_label,
          selection_proportion = selection_proportion,
          years_per_cycle = years_per_cycle,
          higher_is_better = isTRUE(input$breeding_higher_is_better)
        )
        if (!is.null(generation_col_used)) {
          out$generation_stats <- purrr::map_dfr(input$breeding_trait_cols, function(trait_name) {
            breeding_compute_generation_stats(
              data = uploaded_data(),
              trait = trait_name,
              genotype = input$breeding_genotype_col,
              replication = rep_col_used,
              generation_col = generation_col_used,
              selection_proportion = selection_proportion
            )
          })
        } else {
          out$generation_stats <- data.frame()
        }
        out$settings <- data.frame(
          Genotype_column = input$breeding_genotype_col,
          Replication_column = ifelse(is.null(rep_col_used), "", rep_col_used),
          Traits = paste(input$breeding_trait_cols, collapse = ", "),
          Generation_column = ifelse(is.null(generation_col_used), "", generation_col_used),
          Cycle_group_column = ifelse(is.null(cycle_group_used), "", cycle_group_used),
          Selection_proportion_pct = selection_proportion * 100,
          Years_per_cycle = years_per_cycle,
          stringsAsFactors = FALSE
        )
        out
      }, error = function(e) {
        showNotification(
          paste("Breeding Analysis failed:", e$message),
          type = "error",
          duration = NULL
        )
        NULL
      })
      analysis_results(res)
      if (!is.null(res)) {
        saved_results$BREEDING <- res
        updateSelectInput(
          session = session,
          inputId = "breeding_plot_trait",
          choices = input$breeding_trait_cols,
          selected = input$breeding_trait_cols[1]
        )
        analysis_message("Breeding Analysis complete. Check the Result and Chart navbars.")
        showNotification("Breeding Analysis complete.", type = "message")
      } else {
        analysis_message("Breeding Analysis failed. Please check the error message.")
      }
    } else if (input$analysis_method == "LPSI") {
      analysis_message("Running LPSI analysis...")
      res <- tryCatch({
        run_lpsi_with_settings()
      }, error = function(e) {
        showNotification(
          paste("LPSI pipeline failed:", e$message),
          type = "error",
          duration = NULL
        )
        NULL
      })
      analysis_results(res)
      if (!is.null(res)) {
        saved_results$LPSI <- res
        chart_traits <- res$trait_info$Trait
        updateSelectInput(
          session = session,
          inputId = "lpsi_chart_trait",
          choices = chart_traits,
          selected = chart_traits[1]
        )
        analysis_message("LPSI analysis complete. Check the Result and Chart navbars.")
        showNotification("LPSI analysis complete.", type = "message")
      } else {
        analysis_message("LPSI analysis failed. Please check the error message.")
      }
    } else if (input$analysis_method == "MET") {
      met_traits <- input$met_trait_cols
      if (is.null(met_traits) || length(met_traits) == 0) {
        met_traits <- tryCatch(head(get_met_trait_cols(uploaded_data()), 3), error = function(e) character(0))
      }
      validate(need(length(met_traits) > 0, "Choose at least one MET trait to run."))
      met_rep_col <- if (!is.null(input$met_rep_col) && input$met_rep_col != "") input$met_rep_col else NULL
      met_block_col <- if (!is.null(input$met_block_col) && input$met_block_col != "") input$met_block_col else NULL
      met_min_envs_for_biplot <- suppressWarnings(as.integer(input$met_min_envs_for_biplot %||% NA_integer_))
      validate(need(!is.null(met_rep_col), "MET requires a Rep/Replication column."))
      met_traits <- setdiff(met_traits, c(met_rep_col, met_block_col))
      validate(need(length(met_traits) > 0, "Choose at least one MET trait that is not the Rep or Block column."))
      analysis_message(paste0("Running MET analysis for ", length(met_traits), " selected trait(s)..."))
      res <- tryCatch({
        run_met_all_traits(
          uploaded_data(),
          check_varieties = input$met_reference_checks,
          trait_cols = met_traits,
          replication_col = met_rep_col,
          block_col = met_block_col,
          min_envs_for_biplot = met_min_envs_for_biplot
        )
      }, error = function(e) {
        showNotification(
          paste("MET pipeline failed:", e$message),
          type = "error",
          duration = NULL
        )
        NULL
      })
      analysis_results(res)
      if (!is.null(res)) {
        saved_results$MET <- res
        updateSelectInput(
          session = session,
          inputId = "met_result_trait",
          choices = res$met_trait_names,
          selected = res$met_trait_names[1]
        )
        updateSelectInput(
          session = session,
          inputId = "met_plot_trait",
          choices = res$met_trait_names,
          selected = res$met_trait_names[1]
        )
        failed_count <- nrow(res$met_failed_traits)
        failed_note <- if (failed_count > 0) {
          paste0(" ", failed_count, " trait(s) failed and were skipped.")
        } else {
          ""
        }
        analysis_message(paste0(
          "MET analysis complete for ",
          length(res$met_trait_names),
          " trait(s). Check the Result and Chart navbars.",
          failed_note
        ))
        showNotification("MET analysis complete.", type = "message")
      } else {
        analysis_message("MET analysis failed. Please check the error message.")
      }
    } else if (input$analysis_method == "DIVERSITY") {
      req(input$diversity_genotype_col, input$diversity_trait_cols)
      diversity_checks <- unique(clean_text(input$diversity_benchmark_checks))
      diversity_checks <- diversity_checks[!is.na(diversity_checks) & diversity_checks != ""]
      validate(need(length(diversity_checks) > 0, "Choose at least one benchmark check for GDA superiority."))
      rep_col_used <- if (!is.null(input$diversity_rep_col) && input$diversity_rep_col != "") {
        input$diversity_rep_col
      } else {
        NULL
      }
      diversity_data <- uploaded_data()
      diversity_trait_info <- make_diversity_trait_info(
        diversity_data,
        input$diversity_genotype_col,
        input$diversity_trait_cols
      )
      genotype_values_upper <- toupper(clean_text(diversity_data[[input$diversity_genotype_col]]))
      metadata_rows <- genotype_values_upper %in% c(weight_row_labels, direction_row_labels)
      if (any(metadata_rows, na.rm = TRUE)) {
        diversity_data <- diversity_data[!metadata_rows, , drop = FALSE]
      }
      available_genotypes <- unique(clean_text(diversity_data[[input$diversity_genotype_col]]))
      available_genotypes <- available_genotypes[!is.na(available_genotypes) & available_genotypes != ""]
      validate(need(
        all(diversity_checks %in% available_genotypes),
        paste("Selected GDA benchmark check not found:", paste(setdiff(diversity_checks, available_genotypes), collapse = ", "))
      ))
      analysis_message("Running Genetic Diversity analysis...")
      res <- tryCatch({
        out <- si_ext_run_genetic_diversity(
          data = diversity_data,
          genotype = input$diversity_genotype_col,
          traits = input$diversity_trait_cols,
          replication = rep_col_used,
          value_method = input$diversity_value_method %||% "mean",
          n_clusters = NULL
        )
        out$trait_info <- diversity_trait_info
        out$selected_checks <- diversity_checks
        out$superiority <- build_diversity_superiority_data(out, diversity_checks)
        out
      }, error = function(e) {
        showNotification(
          paste("Genetic Diversity failed:", e$message),
          type = "error",
          duration = NULL
        )
        NULL
      })
      analysis_results(res)
      if (!is.null(res)) {
        saved_results$DIVERSITY <- res
        analysis_message("Genetic Diversity complete. Check the Result navbar.")
        showNotification("Genetic Diversity complete.", type = "message")
      } else {
        analysis_message("Genetic Diversity failed. Please check the error message.")
      }
    }
  })
  output$analysis_status <- renderPrint({
    cat(analysis_message(), "\n")
    if (!is.null(analysis_used())) {
      cat("Selected analysis:", analysis_used(), "\n")
    }
  })
  output$mating_anova_table <- renderDT({
    datatable(
      mating_result_for_table(),
      options = list(pageLength = 100, scrollX = TRUE)
    )
  })
  output$mating_gca_table <- renderDT({
    datatable(
      mating_result_for_table(),
      options = list(pageLength = 50, scrollX = TRUE)
    )
  })
  output$mating_sca_table <- renderDT({
    datatable(
      mating_result_for_table(),
      options = list(pageLength = 100, scrollX = TRUE)
    )
  })
  output$mating_variance_table <- renderDT({
    datatable(
      mating_result_for_table(),
      options = list(pageLength = 50, scrollX = TRUE)
    )
  })
  output$breeding_stats_table <- renderDT({
    result <- breeding_result()
    datatable(
      result$genetic_stats,
      options = list(pageLength = 50, scrollX = TRUE)
    )
  })
  output$breeding_response_table <- renderDT({
    result <- breeding_result()
    datatable(
      result$response_per_year,
      options = list(pageLength = 50, scrollX = TRUE)
    )
  })
  output$breeding_realized_table <- renderDT({
    result <- breeding_result()
    validate(need(
      !is.null(result$realized_gain) && nrow(result$realized_gain) > 0,
      "Realized gain needs a current-vs-check group column and labels."
    ))
    datatable(
      result$realized_gain,
      options = list(pageLength = 50, scrollX = TRUE)
    )
  })
  output$breeding_generation_table <- renderDT({
    result <- breeding_result()
    validate(need(
      !is.null(result$generation_stats) && nrow(result$generation_stats) > 0,
      "Generation summary needs a generation/stage column."
    ))
    datatable(
      result$generation_stats,
      options = list(pageLength = 50, scrollX = TRUE)
    )
  })
  output$breeder_recommendation_table <- renderDT({
    rec <- breeder_recommendation()
    validate(need(nrow(rec) > 0, "Run LPSI or MET to create a breeder recommendation table."))
    datatable(
      rec,
      options = list(pageLength = 50, scrollX = TRUE)
    )
  })
  output$trait_table <- renderDT({
    result <- lpsi_result()
    datatable(
      result$trait_info,
      options = list(pageLength = 50, scrollX = TRUE)
    )
  })
  output$final_table <- renderDT({
    result <- lpsi_result()
    datatable(
      result$final_decision,
      options = list(pageLength = 50, scrollX = TRUE)
    )
  })
  output$index_table <- renderDT({
    result <- lpsi_result()
    datatable(
      result$final_decision,
      options = list(pageLength = 50, scrollX = TRUE)
    )
  })
  output$superiority_table <- renderDT({
    result <- lpsi_result()
    datatable(
      result$superiority_index,
      options = list(pageLength = 50, scrollX = TRUE)
    )
  })
  output$anova_full_table <- renderDT({
    result <- lpsi_result()
    validate(
      need(nrow(result$anova_full) > 0, "ANOVA table was not generated.")
    )
    datatable(
      result$anova_full,
      options = list(pageLength = 100, scrollX = TRUE)
    )
  })
  output$lsd_wide_table <- renderDT({
    result <- lpsi_result()
    validate(
      need(nrow(result$lsd_wide) > 0, "LSD table was not generated.")
    )
    datatable(
      result$lsd_wide,
      options = list(pageLength = 50, scrollX = TRUE)
    )
  })
  output$heritability_gain_table <- renderDT({
    result <- lpsi_result()
    validate(
      need(nrow(result$heritability_gain) > 0, "Heritability table was not generated.")
    )
    datatable(
      result$heritability_gain,
      options = list(pageLength = 50, scrollX = TRUE)
    )
  })
  output$lpsi_direct_table <- renderDT({
    direct <- lpsi_direct_selection_r()
    validate(need(nrow(direct) > 0, "Run LPSI and choose a primary trait to view direct selection."))
    datatable(direct, options = list(pageLength = 50, scrollX = TRUE))
  })
  output$lpsi_compare_table <- renderDT({
    comparison <- lpsi_method_comparison_r()
    validate(need(nrow(comparison) > 0, "Run LPSI to compare selection methods."))
    datatable(comparison, options = list(pageLength = 50, scrollX = TRUE))
  })
  output$adjusted_means_table <- renderDT({
    result <- lpsi_result()
    datatable(
      result$actual_adjusted_means,
      options = list(pageLength = 50, scrollX = TRUE)
    )
  })
  output$lpsi_direct_plot <- renderPlot({
    result <- lpsi_result()
    direct_selection <- lpsi_direct_selection_chart_r()
    print(plot_lpsi_direct_selection(
      direct_selection,
      selection_pct = input$lpsi_selection_pct %||% 15,
      lpsi_results = result
    ))
  }, res = 96, execOnResize = TRUE)
  output$lpsi_mean_comparison_plot <- renderPlot({
    result <- lpsi_result()
    validate(need(!is.null(input$lpsi_chart_trait) && input$lpsi_chart_trait != "", "Choose a trait."))
    print(plot_lpsi_mean_comparison(result, input$lpsi_chart_trait))
  }, res = 96, execOnResize = TRUE)
  output$ranking_plot <- renderPlot({
    result <- lpsi_result()
    print(result$ranking_plot)
  }, res = 96, execOnResize = TRUE)
  output$heatmap_plot <- renderPlot({
    result <- lpsi_result()
    if (is.null(result$heatmap_plot)) {
      plot.new()
      text(0.5, 0.5, "Heatmap could not be created.")
    } else {
      draw_chart(result$heatmap_plot)
    }
  }, res = 96, execOnResize = TRUE)
  output$gain_curve_plot <- renderPlot({
    result <- lpsi_result()
    validate(need(!is.null(input$lpsi_chart_trait) && input$lpsi_chart_trait != "", "Choose a trait for the genetic gain chart."))
    print(plot_genetic_gain_curve(
      result$heritability_gain,
      result$trait_info,
      input$lpsi_chart_trait,
      selection_pct = input$lpsi_selection_pct %||% 15
    ))
  }, res = 96, execOnResize = TRUE)
  output$breeding_trend_plot <- renderPlot({
    result <- breeding_result()
    validate(need(
      !is.null(result$generation_stats) && nrow(result$generation_stats) > 0,
      "Genetic trend needs a generation/stage column."
    ))
    validate(need(!is.null(input$breeding_plot_trait) && input$breeding_plot_trait != "", "Choose a trait."))
    trend_data <- result$generation_stats %>% filter(Trait == input$breeding_plot_trait)
    print(breeding_plot_genetic_trend(trend_data))
  }, res = 96, execOnResize = TRUE)
  output$breeding_gam_plot <- renderPlot({
    result <- breeding_result()
    print(breeding_plot_gam(
      result$genetic_stats,
      x_col = "Trait",
      title = "Genetic advance as percent of mean"
    ))
  }, res = 96, execOnResize = TRUE)
  output$breeding_h2_heatmap_plot <- renderPlot({
    result <- breeding_result()
    h2_data <- if (!is.null(result$generation_stats) && nrow(result$generation_stats) > 0) {
      result$generation_stats
    } else {
      result$genetic_stats %>% mutate(Generation = "Overall")
    }
    print(breeding_plot_heritability_heatmap(h2_data))
  }, res = 96, execOnResize = TRUE)
  output$breeding_distribution_plot <- renderPlot({
    result <- breeding_result()
    validate(need(
      !is.null(input$breeding_generation_col) && input$breeding_generation_col != "",
      "Distribution shift needs a generation/stage column."
    ))
    validate(need(!is.null(input$breeding_plot_trait) && input$breeding_plot_trait != "", "Choose a trait."))
    print(breeding_plot_distribution_shift(
      uploaded_data(),
      trait = input$breeding_plot_trait,
      generation_col = input$breeding_generation_col
    ))
  }, res = 96, execOnResize = TRUE)
  output$met_summary_table <- renderDT({
    datatable(met_result_for_table()$genotype_summary, options = list(pageLength = 50, scrollX = TRUE))
  })
  output$met_qc_table <- renderDT({
    qc <- build_met_qc_table(met_result_for_table())
    datatable(qc, options = list(pageLength = 25, scrollX = TRUE))
  })
  output$met_variance_table <- renderDT({
    datatable(met_result_for_table()$variance_components, options = list(pageLength = 20, scrollX = TRUE))
  })
  output$met_blup_table <- renderDT({
    datatable(met_result_for_table()$blups_main, options = list(pageLength = 50, scrollX = TRUE))
  })
  output$met_env_blup_table <- renderDT({
    datatable(met_result_for_table()$blups_environment, options = list(pageLength = 100, scrollX = TRUE))
  })
  output$met_fw_table <- renderDT({
    datatable(met_result_for_table()$fw_results, options = list(pageLength = 50, scrollX = TRUE))
  })
  output$met_ammi_notes_table <- renderDT({
    datatable(met_result_for_table()$ammi_notes, options = list(dom = "t", scrollX = TRUE))
  })
  output$met_gge_notes_table <- renderDT({
    datatable(met_result_for_table()$ammi_notes, options = list(dom = "t", scrollX = TRUE))
  })
  output$met_ammi_table <- renderDT({
    datatable(met_result_for_table()$ammi_genotype, options = list(pageLength = 50, scrollX = TRUE))
  })
  output$met_gge_table <- renderDT({
    datatable(met_result_for_table()$gge_genotype, options = list(pageLength = 50, scrollX = TRUE))
  })
  output$met_selection_table <- renderDT({
    datatable(weighted_met_selection_for_table(), options = list(pageLength = 50, scrollX = TRUE))
  })
  output$met_integrated_table <- renderDT({
    req(analysis_results())
    validate(need(analysis_used() == "MET", "Run MET analysis to view this table."))
    datatable(
      weighted_met_integrated()$ranking,
      options = list(pageLength = 50, scrollX = TRUE)
    )
  })
  output$diversity_values_table <- renderDT({
    datatable(diversity_result()$genotype_values, options = list(pageLength = 50, scrollX = TRUE))
  })
  output$diversity_clusters_table <- renderDT({
    datatable(diversity_result()$clusters, options = list(pageLength = 50, scrollX = TRUE))
  })
  output$diversity_superiority_table <- renderDT({
    superiority <- diversity_superiority_result()
    datatable(superiority, options = list(pageLength = 50, scrollX = TRUE))
  })
  output$diversity_corr_table <- renderDT({
    datatable(tibble::rownames_to_column(as.data.frame(diversity_result()$correlation), "Trait"),
              options = list(pageLength = 50, scrollX = TRUE))
  })
  output$diversity_dendrogram_plot <- renderPlot({
    result <- diversity_result()
    print(plot_diversity_dendrogram(result))
  }, res = 96, execOnResize = TRUE)
  output$diversity_superiority_plot <- renderPlot({
    result <- diversity_result()
    print(plot_diversity_superiority_heatmap(result, diversity_chart_selected_checks()))
  }, res = 96, execOnResize = TRUE)
  output$diversity_corr_heatmap_plot <- renderPlot({
    result <- diversity_result()
    print(plot_diversity_correlation_heatmap(result))
  }, res = 96, execOnResize = TRUE)
  output$met_selection_plot <- renderPlot({
    selection <- weighted_met_selection_for_plot()
    print(plot_met_selection_ranking(selection, input$met_plot_trait))
  }, res = 96, execOnResize = TRUE)
  output$met_integrated_plot <- renderPlot({
    validate(need(
      !is.null(analysis_results()) && identical(analysis_used(), "MET"),
      "Run MET analysis to view this plot."
    ))
    integrated <- weighted_met_integrated()
    print(integrated$plot)
  }, res = 96, execOnResize = TRUE)
  output$met_performance_heatmap <- renderPlot({
    result <- met_result_for_plot()
    print(result$p_perf_heatmap)
  })
  output$met_fw_plot <- renderPlot({
    result <- met_result_for_plot()
    print(result$p_fw_mean_sens)
  }, res = 96, execOnResize = TRUE)
  output$met_fw_regression_plot <- renderPlot({
    result <- met_result_for_plot()
    print(result$p_fw_regression)
  }, res = 96, execOnResize = TRUE)
  output$met_ammi1_plot <- renderPlot({
    result <- met_result_for_plot()
    print(result$p_ammi1)
  }, res = 96, execOnResize = TRUE)
  output$met_ammi2_plot <- renderPlot({
    result <- met_result_for_plot()
    print(result$p_ammi2)
  }, res = 96, execOnResize = TRUE)
  output$met_gge_plot <- renderPlot({
    result <- met_result_for_plot()
    print(result$p_gge)
  }, res = 96, execOnResize = TRUE)
  output$met_env_cor_plot <- renderPlot({
    result <- met_result_for_plot()
    print(result$p_env_cor)
  }, res = 96, execOnResize = TRUE)
  output$met_variance_plot <- renderPlot({
    result <- met_result_for_plot()
    print(result$p_variance)
  })
  
  output$download_chart <- downloadHandler(
    filename = function() {
      chart <- selected_chart()
      extension <- if (is.null(input$chart_format)) "png" else input$chart_format
      paste0(chart$name, ".", extension)
    },
    content = function(file) {
      chart <- selected_chart()
      format <- if (is.null(input$chart_format)) "png" else input$chart_format
      width <- as.numeric(input$chart_width)
      height <- as.numeric(input$chart_height)
      dpi <- 300
      
      validate(need(is.finite(width) && width >= 4 && width <= 30, "Width must be between 4 and 30 inches."))
      validate(need(is.finite(height) && height >= 4 && height <= 30, "Height must be between 4 and 30 inches."))
      if (format == "png") {
        validate(need(
          width * height * dpi^2 <= 100000000,
          "This PNG would be too large. Reduce its width or height."
        ))
        grDevices::png(
          filename = file,
          width = width,
          height = height,
          units = "in",
          res = dpi,
          bg = "white"
        )
      } else if (format == "pdf") {
        grDevices::pdf(
          file = file,
          width = width,
          height = height,
          onefile = TRUE,
          bg = "white"
        )
      } else {
        stop("Unsupported chart format.")
      }
      on.exit(grDevices::dev.off(), add = TRUE)
      draw_chart(chart$plot)
    }
  )
  
  write_analysis_workbook <- function(analysis_type, results, file) {
    tables <- build_export_tables(analysis_type, results)
    if (length(tables) == 0) {
      stop("No result tables are available for this analysis.")
    }
    writexl::write_xlsx(tables, path = file)
  }
  
  output$download_mating <- downloadHandler(
    filename = function() paste0("Mating_analysis_results_", Sys.Date(), ".xlsx"),
    content = function(file) {
      req(saved_results$MATING)
      write_analysis_workbook("MATING", saved_results$MATING, file)
    }
  )
  output$download_breeding <- downloadHandler(
    filename = function() paste0("Breeding_analysis_results_", Sys.Date(), ".xlsx"),
    content = function(file) {
      req(saved_results$BREEDING)
      write_analysis_workbook("BREEDING", saved_results$BREEDING, file)
    }
  )
  output$download_lpsi <- downloadHandler(
    filename = function() paste0("LPSI_selection_results_", Sys.Date(), ".xlsx"),
    content = function(file) {
      req(saved_results$LPSI)
      result <- if (identical(analysis_used(), "LPSI")) {
        lpsi_result()
      } else {
        saved_results$LPSI
      }
      write_analysis_workbook("LPSI", result, file)
    }
  )
  output$download_met <- downloadHandler(
    filename = function() paste0("MET_across_locations_results_", Sys.Date(), ".xlsx"),
    content = function(file) {
      req(saved_results$MET)
      result <- if (identical(analysis_used(), "MET")) {
        result <- analysis_results()
        integrated <- weighted_met_integrated()
        result$met_integrated_ranking <- integrated$ranking
        result$met_integrated_trait_weights <- integrated$trait_weights
        result$met_integrated_adjusted <- integrated$adjusted_performance
        result$met_integrated_standardized <- integrated$standardized_scores
        result
      } else {
        saved_results$MET
      }
      write_analysis_workbook("MET", result, file)
    }
  )
  output$download_diversity <- downloadHandler(
    filename = function() paste0("Genetic_diversity_results_", Sys.Date(), ".xlsx"),
    content = function(file) {
      req(saved_results$DIVERSITY)
      result <- if (identical(analysis_used(), "DIVERSITY")) {
        analysis_results()
      } else {
        saved_results$DIVERSITY
      }
      checks <- get_diversity_selected_checks("diversity_result_benchmark_checks")
      result$selected_checks <- checks
      result$superiority <- build_diversity_superiority_data(result, checks)
      write_analysis_workbook("DIVERSITY", result, file)
    }
  )
  output$download_all <- downloadHandler(
    filename = function() paste0("Selection_analysis_results_", Sys.Date(), ".zip"),
    content = function(file) {
      available <- c(
        MATING = !is.null(saved_results$MATING),
        BREEDING = !is.null(saved_results$BREEDING),
        DIVERSITY = !is.null(saved_results$DIVERSITY),
        LPSI = !is.null(saved_results$LPSI),
        MET = !is.null(saved_results$MET)
      )
      validate(need(any(available), "Run at least one analysis before exporting."))
      
      temp_dir <- tempfile("selection_exports_")
      dir.create(temp_dir)
      on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
      
      export_files <- character()
      export_names <- c(
        MATING = "Mating_analysis_results.xlsx",
        BREEDING = "Breeding_analysis_results.xlsx",
        DIVERSITY = "Genetic_diversity_results.xlsx",
        LPSI = "LPSI_selection_results.xlsx",
        MET = "MET_across_locations_results.xlsx"
      )
      for (analysis_type in names(available)[available]) {
        path <- file.path(temp_dir, export_names[[analysis_type]])
        result <- saved_results[[analysis_type]]
        if (analysis_type == "LPSI" && identical(analysis_used(), "LPSI")) {
          result <- lpsi_result()
        }
        if (analysis_type == "MET" && identical(analysis_used(), "MET")) {
          result <- analysis_results()
          integrated <- weighted_met_integrated()
          result$met_integrated_ranking <- integrated$ranking
          result$met_integrated_trait_weights <- integrated$trait_weights
          result$met_integrated_adjusted <- integrated$adjusted_performance
          result$met_integrated_standardized <- integrated$standardized_scores
        }
        if (analysis_type == "DIVERSITY") {
          checks <- get_diversity_selected_checks("diversity_result_benchmark_checks")
          result$selected_checks <- checks
          result$superiority <- build_diversity_superiority_data(result, checks)
        }
        write_analysis_workbook(
          analysis_type,
          result,
          path
        )
        export_files <- c(export_files, path)
      }
      if (!is.null(saved_results$LPSI)) {
        lpsi_bridge_result <- if (!is.null(saved_results$LPSI) && identical(analysis_used(), "LPSI")) {
          lpsi_result()
        } else {
          saved_results$LPSI
        }
        recommendation <- si_breeder_recommendation_table(lpsi_bridge_result, NULL)
        if (nrow(recommendation) > 0) {
          recommendation_path <- file.path(temp_dir, "Breeder_recommendation.xlsx")
          writexl::write_xlsx(
            list(Breeder_recommendation = recommendation),
            path = recommendation_path
          )
          export_files <- c(export_files, recommendation_path)
        }
      }
      zip::zipr(
        zipfile = file,
        files = export_files,
        compression_level = 9
      )
    }
  )
}


# Run app
shinyApp(ui = ui, server = server)
