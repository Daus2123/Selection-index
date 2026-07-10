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

# Standalone mating-design analysis functions
source("mating_design_module.R", local = TRUE)

# Source-safe breeding gain analysis functions
source("Breeding Modul.R", local = TRUE)


# User settings
id_col  <- "Variety"
rep_col <- "Rep"
remove_cols <- c("Sex", "FsC", "Pr", "Gs", "DM")
weight_row_labels <- c("WEIGHT", "WEIGHTS", "IMPORTANCE", "IMPORTANT")
direction_row_labels <- c("DIRECTION", "DIRECTIONS", "TRAIT_DIRECTION")
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
    add_sheet("01", "summary", results$trait_info)
    add_sheet("02", "anova", results$anova_full)
    add_sheet("03", "mean_comparison", results$lsd_wide)
    add_sheet("04", "superiority_mean", results$superiority_index)
    add_sheet("04b", "superiority_by_check", results$superiority_by_check)
    add_sheet("05", "selection_index", results$index_ranking)
    add_sheet("06", "decision", results$final_decision)
    add_sheet("07", "heritability_gain", results$heritability_gain)
  } else if (analysis_type == "MET") {
    for (trait in results$met_trait_names) {
      result <- results$met_by_trait[[trait]]
      if (!is.null(result)) {
        add_sheet("01_lmm", trait, result$variance_components)
        add_sheet("02_blup", trait, result$blups_main)
        add_sheet("03_fw", trait, result$fw_results)
        add_sheet("04_ammi", trait, result$ammi_genotype)
        add_sheet("05_gge", trait, result$gge_genotype)
        add_sheet("06_selection", trait, result$met_selection)
      }
    }
    add_sheet("07_overall", "selection", results$met_integrated_ranking)
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
  check_original_name <- as.character(df_data[[id_col]][1])
  candidate_cols <- setdiff(
    names(df_data),
    c(id_col, rep_col, remove_cols)
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
    notes = notes
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
plot_genetic_gain_curve <- function(heritability_gain, trait_info, trait_name = NULL) {
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
  ga <- abs(as.numeric(gain_row$Genetic_advance))
  h2 <- as.numeric(gain_row$Broad_sense_H2)
  select_pct <- as.numeric(gain_row$Selection_intensity_pct)
  select_prop <- select_pct / 100
  if (!is.finite(select_prop) || select_prop <= 0 || select_prop >= 1) {
    select_prop <- lpsi_selection_intensity
  }
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



# Main LPSI function
run_selection_pipeline <- function(df_raw, check_varieties = NULL) {
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
    hybrid_levels <- final_decision %>%
      arrange(desc(Selection_Index)) %>%
      pull(Original_ID)
    if (length(hybrid_levels) == 0) {
      hybrid_levels <- unique(heatmap_data$Original_ID)
    }
    check_levels <- check_index_details$Original_ID
    heatmap_data <- heatmap_data %>%
      mutate(
        Original_ID = factor(Original_ID, levels = hybrid_levels),
        Trait = factor(Trait, levels = rev(trait_cols)),
        Check_Original_ID = factor(Check_Original_ID, levels = check_levels),
        Cell_label = paste0(sprintf("%.1f", Superiority_pct), "%")
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
      data.frame(
        Trait = tr,
        Mean = round(trait_mean, 4),
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
    superiority_by_check = superiority_by_check_long,
    priority_summary = priority_summary,
    final_decision = final_decision,
    heritability_gain = heritability_gain,
    anova_full = anova_full,
    lsd_wide = lsd_wide,
    ranking_plot = p_index,
    heatmap_plot = heatmap_plot,
    check_original_label = check_original_label,
    selected_checks = selected_checks
  ))
}



# MET function
MET_W_YIELD <- 0.50
MET_W_FW <- 0.25
MET_W_ASV <- 0.25
MET_MIN_ENVS_FOR_BIPLOT <- 2
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
  trait_cols <- pick_met_numeric_traits(df_data, c("Genotype", "Environment", rep_col, remove_cols))
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
make_met_data <- function(df_raw, trait_used = NULL) {
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
  dat <- df %>% mutate(Genotype = trimws(as.character(Genotype)), Environment = trimws(as.character(Environment))) %>% filter(!is.na(Genotype), Genotype != "", !is.na(Environment), Environment != "", !is.na(Weight)) %>% mutate(Genotype = factor(Genotype), Environment = factor(Environment))
  if (nrow(dat) == 0) {
    stop("No valid MET rows remained after cleaning Genotype, Environment, and Weight.")
  }
  if (n_distinct(dat$Genotype) < 2) {
    stop("MET requires at least 2 genotypes.")
  }
  if (n_distinct(dat$Environment) < 2) {
    stop("MET requires at least 2 environments.")
  }
  list(data = dat, trait_used = trait_used)
}
met_safe_lrt <- function(model_a, model_b, test_name) {
  tryCatch({
    as.data.frame(anova(model_a, model_b)) %>% rownames_to_column("Model") %>% mutate(Test = test_name, .before = 1)
  }, error = function(e) {
    data.frame(Test = test_name, Model = "LRT failed", npar = NA, AIC = NA, BIC = NA, logLik = NA, deviance = NA, Chisq = NA, Df = NA, `Pr(>Chisq)` = NA, Note = e$message, check.names = FALSE)
  })
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
run_met_pipeline <- function(df_raw, trait_used = NULL, check_varieties = NULL) {
  input <- make_met_data(df_raw, trait_used)
  dat <- input$data
  trait_used <- input$trait_used
  n_envs_total <- n_distinct(dat$Environment)
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
  m_full <- lmer(Weight ~ (1|Environment) + (1|Genotype) + (1|Genotype:Environment), data = dat_clean, control = lmerControl(optimizer = "bobyqa"))
  vc <- as.data.frame(VarCorr(m_full))
  total_var <- sum(vc$vcov, na.rm = TRUE)
  vc$percent <- round((vc$vcov / total_var) * 100, 2)
  vc$label <- recode(vc$grp, "Environment" = "Environment (E)", "Genotype" = "Genotype (G)", "Genotype:Environment" = "GxE Interaction", "Residual" = "Residual")
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
  model_summary <- data.frame(Trait_used = trait_used, N_rows_clean = nrow(dat_clean), N_genotypes = n_distinct(dat_clean$Genotype), N_environments = n_distinct(dat_clean$Environment), Harmonic_replication = round(n_r, 3), Stability_ratio_G_over_G_plus_GxE = round(stability_ratio, 4), Broad_sense_H2 = round(H2, 4), Controls_used = paste(controls_used, collapse = ", "), Notes = paste(notes, collapse = " | "), stringsAsFactors = FALSE)
  p_variance <- ggplot(variance_components %>% mutate(Source = fct_reorder(Source, -Percent)), aes(x = Source, y = Percent, fill = Source)) + geom_bar(stat = "identity", width = 0.6) + geom_text(aes(label = paste0(Percent, "%")), vjust = -0.5, size = 4) + scale_fill_manual(values = c("#9B59B6", "#3498DB", "#E67E22", "#2ECC71")) + labs(title = paste0("Variance partitioning â€” ", trait_used), x = NULL, y = "% of Total Variance") + theme_bw() + theme(legend.position = "none")
  res_df <- data.frame(fitted = fitted(m_full), residual = residuals(m_full))
  p_qq <- ggplot(res_df, aes(sample = residual)) + stat_qq(color = "#3498DB", alpha = 0.6) + stat_qq_line(color = "#E74C3C", linewidth = 0.8) + labs(title = "Normal Q-Q", x = "Theoretical", y = "Sample") + theme_bw()
  p_rvf <- ggplot(res_df, aes(x = fitted, y = residual)) + geom_point(color = "#3498DB", alpha = 0.5, size = 1.5) + geom_hline(yintercept = 0, linetype = "dashed", color = "#E74C3C", linewidth = 0.8) + geom_smooth(method = "loess", se = FALSE, color = "#F39C12", linewidth = 0.8, span = 0.8) + labs(title = "Residuals vs Fitted", x = "Fitted", y = "Residuals") + theme_bw()
  p_residual <- p_qq + p_rvf
  m_no_e <- lmer(Weight ~ (1|Genotype) + (1|Genotype:Environment), data = dat_clean, control = lmerControl(optimizer = "bobyqa"))
  m_no_gxe <- lmer(Weight ~ (1|Environment) + (1|Genotype), data = dat_clean, control = lmerControl(optimizer = "bobyqa"))
  m_no_g <- lmer(Weight ~ (1|Environment) + (1|Genotype:Environment), data = dat_clean, control = lmerControl(optimizer = "bobyqa"))
  m_null <- lmer(Weight ~ (1|Environment), data = dat_clean, control = lmerControl(optimizer = "bobyqa"))
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
  genos_for_biplot <- presence %>% filter(n_envs >= MET_MIN_ENVS_FOR_BIPLOT) %>% pull(Genotype) %>% as.character()
  GxE_obs_wide <- GxE_matrix_wide[rownames(GxE_matrix_wide) %in% genos_for_biplot, , drop = FALSE]
  excluded_genos <- setdiff(rownames(GxE_matrix_wide), genos_for_biplot)
  ammi_notes <- data.frame(N_genotypes_enter_AMMI_GGE = nrow(GxE_obs_wide), N_genotypes_total = nrow(GxE_matrix_wide), Excluded_genotypes = ifelse(length(excluded_genos) > 0, paste(excluded_genos, collapse = ", "), "none"), stringsAsFactors = FALSE)
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
    AMMI_geno <- AMMI_geno[order(AMMI_geno$ASV), ]
    p_ammi1 <- ggplot() + geom_point(data = AMMI_geno, aes(x = PC1, y = GenMean), size = 3, color = "#2C3E50") + geom_text(data = AMMI_geno, aes(x = PC1, y = GenMean, label = Genotype), vjust = -0.8, size = 3, color = "#2C3E50") + geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") + geom_point(data = AMMI_env, aes(x = PC1, y = EnvMean), size = 4, shape = 17, color = "#E74C3C") + geom_text(data = AMMI_env, aes(x = PC1, y = EnvMean, label = Environment), vjust = -0.9, size = 3.5, color = "#E74C3C", fontface = "bold") + labs(title = paste0("AMMI1 biplot â€” ", trait_used, " BLUPs"), x = paste0("IPCA1 (", round(PC_pct[1] * 100, 1), "%)"), y = paste0("Mean BLUP for ", trait_used)) + theme_bw()
    if (n_pc >= 2) {
      scale_ammi <- ifelse(max(abs(AMMI_env$PC1), na.rm = TRUE) == 0, 1, max(abs(AMMI_geno$PC1), na.rm = TRUE) / max(abs(AMMI_env$PC1), na.rm = TRUE) * 0.7)
      AMMI_env_sc <- AMMI_env %>% mutate(PC1 = PC1 * scale_ammi, PC2 = PC2 * scale_ammi)
      p_ammi2 <- ggplot() + geom_segment(data = AMMI_env_sc, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.25, "cm"), type = "closed"), color = "#E74C3C", linewidth = 0.8) + geom_text(data = AMMI_env_sc, aes(x = PC1 * 1.12, y = PC2 * 1.12, label = Environment), color = "#E74C3C", size = 3.5, fontface = "bold") + geom_point(data = AMMI_geno, aes(x = PC1, y = PC2), color = "#2C3E50", size = 2.5) + geom_text(data = AMMI_geno, aes(x = PC1, y = PC2, label = Genotype), color = "#2C3E50", vjust = -0.8, size = 2.8) + geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") + geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") + labs(title = paste0("AMMI2 biplot â€” ", trait_used, " BLUPs"), x = paste0("IPCA1 (", round(PC_pct[1] * 100, 1), "%)"), y = paste0("IPCA2 (", round(PC_pct[2] * 100, 1), "%)")) + theme_bw()
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
    if (n_pc_gge >= 2) {
      scale_gge <- ifelse(max(abs(GGE_env$PC1), na.rm = TRUE) == 0, 1, max(abs(GGE_geno$PC1), na.rm = TRUE) / max(abs(GGE_env$PC1), na.rm = TRUE) * 0.7)
      GGE_env_sc <- GGE_env %>% mutate(PC1 = PC1 * scale_gge, PC2 = PC2 * scale_gge, label_x = PC1 * 1.15, label_y = PC2 * 1.15)
      p_gge <- ggplot() + geom_segment(data = GGE_env_sc, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.25, "cm"), type = "closed"), color = "#E74C3C", linewidth = 0.8) + geom_text(data = GGE_env_sc, aes(x = label_x, y = label_y, label = Environment), color = "#E74C3C", size = 3.5, fontface = "bold") + geom_point(data = GGE_geno, aes(x = PC1, y = PC2), color = "#2C3E50", size = 2.5) + geom_text(data = GGE_geno, aes(x = PC1, y = PC2, label = Genotype), color = "#2C3E50", vjust = -0.8, size = 2.8) + geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") + geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") + labs(title = paste0("GGE biplot â€” ", trait_used, " BLUPs"), x = paste0("PC1 (", round(GGE_pct[1] * 100, 1), "%)"), y = paste0("PC2 (", round(GGE_pct[2] * 100, 1), "%)")) + theme_bw()
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
  return(list(raw_data = df_raw, met_data = dat, met_cleaned_data = dat_clean, outlier_summary = outlier_summary, presence = presence, model_summary = model_summary, variance_components = variance_components, lrt_table = lrt_table, blups_main = BLUPs_main, blups_environment = BLUPs_env_full, gxe_matrix = GxE_matrix_wide %>% rownames_to_column("Genotype"), fw_results = FW_results, ammi_notes = ammi_notes, ammi_genotype = AMMI_geno, ammi_environment = AMMI_env, gge_genotype = GGE_geno, gge_environment = GGE_env, met_selection = selection, p_before = p_before, p_after = p_after, p_variance = p_variance, p_residual = p_residual, p_blup = p_blup, p_accuracy = p_accuracy, p_perf_heatmap = p_perf_heatmap, p_fw_mean_sens = p_fw_mean_sens, p_fw_regression = p_fw_regression, p_ammi1 = p_ammi1, p_ammi2 = p_ammi2, p_gge = p_gge, p_env_cor = p_env_cor, p_met_selection = p_met_selection))
}
run_met_all_traits <- function(df_raw, check_varieties = NULL) {
  trait_cols <- get_met_trait_cols(df_raw)
  results <- list()
  failures <- data.frame(Trait = character(), Error = character(), stringsAsFactors = FALSE)
  for (trait in trait_cols) {
    result <- tryCatch({
      run_met_pipeline(df_raw, trait, check_varieties = check_varieties)
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
      border-bottom: 1px solid #E6E0D2;
      margin-bottom: 16px;
      padding-bottom: 14px;
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
      border: 1px solid #E1DACA;
      border-radius: 9px;
      margin-bottom: 7px;
      padding: 10px 12px 10px 34px;
      background: #FFFFFF;
    }
    .analysis-choice .form-check:has(input:checked) {
      background: #EEF5E9;
      border-color: #4E7A43;
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
      display: inline-block;
      margin-right: 16px;
    }
    .analysis-status {
      background: #F4F7EE;
      border-radius: 8px;
      color: #4B6543;
      font-size: 12px;
      margin-top: 12px;
      padding: 9px 10px;
    }
    .analyze-workspace {
      height: calc(100vh - 112px);
      height: calc(100dvh - 112px);
      overflow: hidden;
    }
    .analyze-workspace > .bslib-grid {
      height: 100%;
    }
    .analyze-pane {
      max-height: 100%;
      min-height: 0;
      overflow-y: auto;
      padding-right: 6px;
      scrollbar-gutter: stable;
    }
    .analyze-pane > .card:last-child {
      margin-bottom: 0;
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
      padding-top: 16px;
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
      .analyze-workspace {
        height: auto;
        overflow: visible;
      }
      .analyze-pane {
        max-height: none;
        overflow-y: visible;
        padding-right: 0;
      }
    }
  ")),
  # Upload data navbar
  nav_panel(
    title = "Data",
    layout_columns(
      col_widths = c(4, 8),
      card(
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
        panel_header("Preview", "First rows of the uploaded workbook"),
        DTOutput("raw_table")
      )
    )
  ),
  # Analysis navbar
  nav_panel(
    title = "Analyze",
    tags$div(
      class = "analyze-workspace",
      layout_columns(
        col_widths = c(4, 8),
        tags$div(
          class = "analyze-pane analyze-left-pane",
          card(
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
            inline = TRUE
          )
        ),
        tags$div(
          class = "control-section analysis-choice",
          tags$div(class = "control-label", "Choose the analysis"),
          radioButtons(
            inputId = "analysis_method",
            label = NULL,
            choiceNames = list(
              tags$span(
                tags$span(class = "analysis-name", "Mating analysis"),
                tags$span(class = "analysis-description", "Compare parents and cross combinations")
              ),
              tags$span(
                tags$span(class = "analysis-name", "Breeding Analysis"),
                tags$span(class = "analysis-description", "Track heritability, genetic gain, and response per year")
              ),
              tags$span(
                tags$span(class = "analysis-name", "LPSI"),
                tags$span(class = "analysis-description", "Blend every trait into one selection score")
              ),
              tags$span(
                tags$span(class = "analysis-name", "MET"),
                tags$span(class = "analysis-description", "Compare performance and stability across locations")
              )
            ),
            choiceValues = c("MATING", "BREEDING", "LPSI", "MET"),
            selected = "MATING"
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
          condition = "input.analysis_method == 'MET'",
          tags$div(
            class = "control-section",
            uiOutput("check_variety_inputs")
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
            card_header(uiOutput("diagnostic_header")),
            plotOutput("diagnostic_plot", height = "470px")
          ),
          card(
            panel_header("Data summary by trait", "Residual normality test for every measured trait"),
            DTOutput("shapiro_table")
          )
        )
      )
    )
  ),
  # Result navbar
  nav_panel(
    title = "Results",
    layout_columns(
      col_widths = c(3, 9),
      card(
        panel_header("Results", "Grouped by analysis question"),
        sidebar_radio_menu(
          input_id = "result_view",
          selected = "mating_anova",
          group_controls = list(
            "One overall ranking (LPSI)" = uiOutput("lpsi_benchmark_check_inputs"),
            "Across locations (MET)" = selectInput(
              inputId = "met_result_trait",
              label = "Choose trait",
              choices = NULL
            )
          ),
          groups = list(
            "Which crosses work best" = c(
              "ANOVA" = "mating_anova",
              "GCA - parent effects" = "mating_gca",
              "SCA - cross effects" = "mating_sca",
              "Variance breakdown" = "mating_variance"
            ),
            "Breeding progress" = c(
              "Genetic parameters" = "breeding_stats",
              "Response per year" = "breeding_response",
              "Realized gain" = "breeding_realized",
              "Generation summary" = "breeding_generation"
            ),
            "One overall ranking (LPSI)" = c(
              "Summary" = "lpsi_trait",
              "ANOVA" = "lpsi_anova",
              "Mean comparison" = "lpsi_lsd",
              "Superiority" = "lpsi_superiority",
              "Superiority by check" = "lpsi_superiority_by_check",
              "Heritability & gain" = "lpsi_heritability",
              "Selected varieties" = "lpsi_ranking"
            ),
            "Across locations (MET)" = c(
              "Mixed model" = "met_variance",
              "BLUP" = "met_blup",
              "Stability (FW)" = "met_fw",
              "AMMI" = "met_ammi",
              "GGE" = "met_gge",
              "Selection" = "met_selection",
              "Overall Selection" = "met_integrated"
            )
          )
        ),
        conditionalPanel(
          "input.result_view == 'met_selection' || input.result_view == 'met_integrated'",
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
        card_header(uiOutput("result_header")),
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
        conditionalPanel("input.result_view == 'lpsi_superiority_by_check'", DTOutput("superiority_by_check_table")),
        conditionalPanel("input.result_view == 'lpsi_anova'", DTOutput("anova_full_table")),
        conditionalPanel("input.result_view == 'lpsi_lsd'", DTOutput("lsd_wide_table")),
        conditionalPanel("input.result_view == 'lpsi_heritability'", DTOutput("heritability_gain_table")),
        conditionalPanel("input.result_view == 'met_variance'", DTOutput("met_variance_table")),
        conditionalPanel("input.result_view == 'met_blup'", DTOutput("met_blup_table")),
        conditionalPanel("input.result_view == 'met_fw'", DTOutput("met_fw_table")),
        conditionalPanel("input.result_view == 'met_ammi'", DTOutput("met_ammi_table")),
        conditionalPanel("input.result_view == 'met_gge'", DTOutput("met_gge_table")),
        conditionalPanel("input.result_view == 'met_selection'", DTOutput("met_selection_table")),
        conditionalPanel("input.result_view == 'met_integrated'", DTOutput("met_integrated_table"))
      )
    )
  ),
  # Plot navbar
  nav_panel(
    title = "Charts",
    layout_columns(
      col_widths = c(3, 9),
      card(
        panel_header("Charts", "Visual patterns from completed analyses"),
        sidebar_radio_menu(
          input_id = "plot_view",
          selected = "lpsi_ranking_plot",
          group_controls = list(
            "LPSI" = conditionalPanel(
              "input.plot_view == 'lpsi_gain_curve'",
              selectInput(
                inputId = "lpsi_gain_trait",
                label = "Choose gain trait",
                choices = NULL
              )
            ),
            "Breeding Analysis" = conditionalPanel(
              "input.plot_view == 'breeding_trend' || input.plot_view == 'breeding_distribution'",
              selectInput(
                inputId = "breeding_plot_trait",
                label = "Choose trait",
                choices = NULL
              )
            ),
            "MET" = selectInput(
              inputId = "met_plot_trait",
              label = "Choose trait",
              choices = NULL
            )
          ),
          groups = list(
            "Breeding Analysis" = c(
              "Genetic trend" = "breeding_trend",
              "GAM" = "breeding_gam",
              "Heritability heatmap" = "breeding_h2_heatmap",
              "Distribution shift" = "breeding_distribution"
            ),
            "LPSI" = c(
              "Ranking" = "lpsi_ranking_plot",
              "Superiority" = "lpsi_heatmap",
              "Genetic gain" = "lpsi_gain_curve"
            ),
            "MET" = c(
              "Environment Correlation" = "met_env_cor",
              "FW Sensitivity" = "met_fw_plot",
              "FW Regression" = "met_fw_regression",
              "AMMI1" = "met_ammi1",
              "AMMI2" = "met_ammi2",
              "GGE" = "met_gge",
              "Ranking" = "met_selection_plot",
              "Overall Ranking" = "met_integrated_plot"
            )
          )
        ),
        tags$div(
          class = "side-subpanel chart-options-subpanel",
          tags$div(class = "side-subpanel-title", "Preview and download"),
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
      ),
      card(
        panel_header("Analysis chart", "Live preview using the selected width and height"),
        uiOutput("chart_preview_ui")
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
    ext <- tools::file_ext(input$excel_file$name)
    validate(
      need(ext %in% c("xlsx", "xls"), "Please upload Excel file only: .xlsx or .xls")
    )
    read_excel(
      input$excel_file$datapath,
      .name_repair = "unique"
    )
  })
  diagnostic_data <- reactive({
    req(uploaded_data())
    tryCatch({
      make_diagnostic_data(uploaded_data())
    }, error = function(e) {
      showNotification(
        paste("Data preparation failed:", e$message),
        type = "error",
        duration = NULL
      )
      NULL
    })
  })
  analysis_results <- reactiveVal(NULL)
  analysis_used <- reactiveVal(NULL)
  analysis_message <- reactiveVal("No analysis has been run yet.")
  saved_results <- reactiveValues(
    MATING = NULL,
    BREEDING = NULL,
    LPSI = NULL,
    MET = NULL
  )
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
  lpsi_result <- reactive({
    req(analysis_results(), uploaded_data())
    validate(need(analysis_used() == "LPSI", "Run LPSI analysis to view this result."))

    benchmark_checks <- input$lpsi_benchmark_checks
    if (is.null(benchmark_checks) || length(benchmark_checks) == 0) {
      benchmark_checks <- analysis_results()$selected_checks
    }
    run_selection_pipeline(
      uploaded_data(),
      check_varieties = benchmark_checks
    )
  })
  met_result_for_table <- reactive({
    req(analysis_results())
    validate(need(analysis_used() == "MET", "Run MET analysis to view this table."))
    req(input$met_result_trait)
    result <- analysis_results()$met_by_trait[[input$met_result_trait]]
    validate(need(!is.null(result), "This MET trait result is not available."))
    result
  })
  met_result_for_plot <- reactive({
    req(analysis_results())
    validate(need(analysis_used() == "MET", "Run MET analysis to view this plot."))
    req(input$met_plot_trait)
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
    build_met_selection_ranking(
      met_result_for_plot(),
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
  selected_chart <- reactive({
    req(input$plot_view)
    view <- input$plot_view

    if (view %in% c("lpsi_ranking_plot", "lpsi_heatmap", "lpsi_gain_curve")) {
      lpsi <- lpsi_result()
      validate(need(analysis_used() == "LPSI", "Run LPSI analysis before downloading this chart."))
      if (view == "lpsi_ranking_plot") {
        return(list(
          plot = lpsi$ranking_plot,
          name = "LPSI_ranking"
        ))
      }
      if (view == "lpsi_gain_curve") {
        gain_trait <- input$lpsi_gain_trait
        validate(need(!is.null(gain_trait) && gain_trait != "", "Choose a trait for the genetic gain chart."))
        return(list(
          plot = plot_genetic_gain_curve(
            lpsi$heritability_gain,
            lpsi$trait_info,
            gain_trait
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
    view <- if (is.null(input$plot_view)) "lpsi_ranking_plot" else input$plot_view
    output_id <- switch(
      view,
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
      cat("Detected traits:", paste(diag$trait_cols, collapse = ", "), "\n")
    } else {
      cat("Trait and replication information could not be detected.\n")
      cat("Please check that the file has columns named Variety and Rep.\n")
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
      lpsi_superiority_by_check = "Trait superiority by selected check",
      lpsi_heritability = "Heritability and genetic gain",
      lpsi_ranking = "Selected varieties",
      met_variance = paste("Mixed model -", trait),
      met_blup = paste("BLUP -", trait),
      met_fw = paste("Finlay-Wilkinson stability -", trait),
      met_ammi = paste("AMMI -", trait),
      met_gge = paste("GGE -", trait),
      met_selection = paste("Selection -", trait),
      met_integrated = "Overall MET selection",
      "Analysis results"
    )
    tags$div(
      class = "panel-heading",
      tags$div(class = "panel-title", title),
      tags$div(class = "panel-subtitle", "Results from the selected analysis view")
    )
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
      inputId = "lpsi_gain_trait",
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
    choices <- unique(clean_text(prepared$data$Genotype))
    choices <- choices[!is.na(choices) & choices != ""]
    selectizeInput(
      inputId = "met_reference_checks",
      label = "Reference check",
      choices = choices,
      selected = character(0),
      multiple = TRUE,
      options = list(placeholder = "Optional: choose one or more check varieties")
    )
  })
  output$lpsi_benchmark_check_inputs <- renderUI({
    if (is.null(saved_results$LPSI) && !identical(analysis_used(), "LPSI")) {
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
    selected <- if (!is.null(saved_results$LPSI) && length(saved_results$LPSI$selected_checks) > 0) {
      saved_results$LPSI$selected_checks
    } else {
      prepared$check_original_name
    }
    selected <- selected[selected %in% choices]
    if (length(selected) == 0) {
      selected <- prepared$check_original_name
    }

    selectizeInput(
      inputId = "lpsi_benchmark_checks",
      label = "Benchmark check",
      choices = choices,
      selected = selected,
      multiple = TRUE,
      options = list(placeholder = "Choose one or more check varieties")
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
        run_selection_pipeline(
          uploaded_data()
        )
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
        gain_traits <- res$heritability_gain %>%
          filter(
            is.finite(Mean),
            is.finite(Phenotypic_variance),
            Phenotypic_variance > 0,
            is.finite(Genetic_advance)
          ) %>%
          pull(Trait)
        if (length(gain_traits) == 0) {
          gain_traits <- res$trait_info$Trait
        }
        updateSelectInput(
          session = session,
          inputId = "lpsi_gain_trait",
          choices = gain_traits,
          selected = gain_traits[1]
        )
        analysis_message("LPSI analysis complete. Check the Result and Plot navbars.")
        showNotification("LPSI analysis complete.", type = "message")
      } else {
        analysis_message("LPSI analysis failed. Please check the error message.")
      }
    } else if (input$analysis_method == "MET") {
      analysis_message("Running MET analysis for all detected numeric traits...")
      res <- tryCatch({
        run_met_all_traits(
          uploaded_data(),
          check_varieties = input$met_reference_checks
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
          " trait(s). Check the Result and Plot navbars.",
          failed_note
        ))
        showNotification("MET analysis complete.", type = "message")
      } else {
        analysis_message("MET analysis failed. Please check the error message.")
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
  output$superiority_by_check_table <- renderDT({
    result <- lpsi_result()
    datatable(
      result$superiority_by_check,
      options = list(pageLength = 100, scrollX = TRUE)
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
  output$adjusted_means_table <- renderDT({
    result <- lpsi_result()
    datatable(
      result$actual_adjusted_means,
      options = list(pageLength = 50, scrollX = TRUE)
    )
  })
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
    validate(need(!is.null(input$lpsi_gain_trait) && input$lpsi_gain_trait != "", "Choose a trait for the genetic gain chart."))
    print(plot_genetic_gain_curve(
      result$heritability_gain,
      result$trait_info,
      input$lpsi_gain_trait
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
  output$met_model_summary_table <- renderDT({
    datatable(met_result_for_table()$model_summary, options = list(pageLength = 20, scrollX = TRUE))
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
  output$met_selection_plot <- renderPlot({
    print(plot_met_selection_ranking(weighted_met_selection_for_plot(), input$met_plot_trait))
  }, res = 96, execOnResize = TRUE)
  output$met_integrated_plot <- renderPlot({
    req(analysis_results())
    validate(need(analysis_used() == "MET", "Run MET analysis to view this plot."))
    print(weighted_met_integrated()$plot)
  }, res = 96, execOnResize = TRUE)
  output$met_performance_heatmap <- renderPlot({
    print(met_result_for_plot()$p_perf_heatmap)
  })
  output$met_fw_plot <- renderPlot({
    print(met_result_for_plot()$p_fw_mean_sens)
  }, res = 96, execOnResize = TRUE)
  output$met_fw_regression_plot <- renderPlot({
    print(met_result_for_plot()$p_fw_regression)
  }, res = 96, execOnResize = TRUE)
  output$met_ammi1_plot <- renderPlot({
    print(met_result_for_plot()$p_ammi1)
  }, res = 96, execOnResize = TRUE)
  output$met_ammi2_plot <- renderPlot({
    print(met_result_for_plot()$p_ammi2)
  }, res = 96, execOnResize = TRUE)
  output$met_gge_plot <- renderPlot({
    print(met_result_for_plot()$p_gge)
  }, res = 96, execOnResize = TRUE)
  output$met_env_cor_plot <- renderPlot({
    print(met_result_for_plot()$p_env_cor)
  }, res = 96, execOnResize = TRUE)
  output$met_variance_plot <- renderPlot({
    print(met_result_for_plot()$p_variance)
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
      write_analysis_workbook("MET", saved_results$MET, file)
    }
  )
  output$download_all <- downloadHandler(
    filename = function() paste0("Selection_analysis_results_", Sys.Date(), ".zip"),
    content = function(file) {
      available <- c(
        MATING = !is.null(saved_results$MATING),
        BREEDING = !is.null(saved_results$BREEDING),
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
        LPSI = "LPSI_selection_results.xlsx",
        MET = "MET_across_locations_results.xlsx"
      )
      for (analysis_type in names(available)[available]) {
        path <- file.path(temp_dir, export_names[[analysis_type]])
        result <- saved_results[[analysis_type]]
        if (analysis_type == "LPSI" && identical(analysis_used(), "LPSI")) {
          result <- lpsi_result()
        }
        write_analysis_workbook(
          analysis_type,
          result,
          path
        )
        export_files <- c(export_files, path)
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
