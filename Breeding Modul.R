if (FALSE) {
############################################################
# GENETIC GAIN TRACKING PIPELINE FOR RECURRENT SELECTION
# Tracks: variance components, heritability, GA, GAM,
#         realized gain per cycle, and response per year (R)
############################################################

# ---- 0. SETUP ----
required_pkgs <- c("dplyr", "ggplot2", "tidyr", "agricolae")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs)) install.packages(new_pkgs)

library(dplyr)
library(ggplot2)
library(tidyr)
library(agricolae)

set.seed(123)  # reproducibility for the simulated example data


############################################################
# PART 1: VARIANCE COMPONENTS, HERITABILITY, GA, GAM
# Used at Stage 1 replicated trials (RCBD: genotype x replication)
############################################################

compute_genetic_stats <- function(data, trait, genotype, replication, k = 2.06) {
  # data: data frame with columns for trait, genotype, replication
  # k: selection intensity constant (2.06 = top 5%, 1.4 = top 20%, 2.67 = top 1%)
  
  formula_str <- paste(trait, "~", genotype, "+", replication)
  model <- aov(as.formula(formula_str), data = data)
  anova_tab <- summary(model)[[1]]
  
  ms_genotype <- anova_tab[genotype, "Mean Sq"]
  ms_error    <- anova_tab["Residuals", "Mean Sq"]
  r           <- length(unique(data[[replication]]))   # number of reps
  
  # Genotypic and phenotypic variance
  GV <- (ms_genotype - ms_error) / r
  GV <- max(GV, 0)                # variance cannot be negative
  PV <- GV + ms_error
  
  # Coefficients of variation
  trait_mean <- mean(data[[trait]], na.rm = TRUE)
  GCV <- (sqrt(GV) / trait_mean) * 100
  PCV <- (sqrt(PV) / trait_mean) * 100
  
  # Heritability (broad sense)
  H2 <- (GV / PV) * 100
  
  # Genetic advance and GA as % of mean
  GA  <- k * sqrt(PV) * (H2 / 100)
  GAM <- (GA / trait_mean) * 100
  
  # Descriptive stats
  se    <- sd(data[[trait]], na.rm = TRUE) / sqrt(nrow(data))
  range_vals <- range(data[[trait]], na.rm = TRUE)
  
  data.frame(
    Mean       = round(trait_mean, 2),
    SE         = round(se, 2),
    Min        = round(range_vals[1], 2),
    Max        = round(range_vals[2], 2),
    GV         = round(GV, 3),
    PV         = round(PV, 3),
    `GCV_pct`  = round(GCV, 2),
    `PCV_pct`  = round(PCV, 2),
    `H2_pct`   = round(H2, 2),
    GA         = round(GA, 3),
    `GAM_pct`  = round(GAM, 2)
  )
}


############################################################
# PART 2: REALIZED GAIN PER CYCLE
# Compares this cycle's selected lines against prior-cycle
# check lines grown in the SAME trial and year
############################################################

compute_realized_gain <- function(data, trait, group_col,
                                  current_label, check_label) {
  # group_col: column identifying "current_cycle" vs "check_prior_cycle"
  
  means <- data %>%
    group_by(.data[[group_col]]) %>%
    summarise(mean_trait = mean(.data[[trait]], na.rm = TRUE), .groups = "drop")
  
  current_mean <- means$mean_trait[means[[group_col]] == current_label]
  check_mean   <- means$mean_trait[means[[group_col]] == check_label]
  
  realized_gain     <- current_mean - check_mean
  realized_gain_pct <- (realized_gain / check_mean) * 100
  
  data.frame(
    Current_Mean       = round(current_mean, 2),
    Check_Mean         = round(check_mean, 2),
    Realized_Gain      = round(realized_gain, 2),
    Realized_Gain_Pct  = round(realized_gain_pct, 2)
  )
}


############################################################
# PART 3: RESPONSE PER YEAR (R)
# The overall pipeline efficiency score
############################################################

compute_response_per_year <- function(GV, H2_pct, k = 2.06, Y = 3) {
  # GV: genotypic variance (from compute_genetic_stats)
  # H2_pct: heritability in percent
  # Y: years per cycle
  sigma_G <- sqrt(GV)
  r <- sqrt(H2_pct / 100)      # selection accuracy
  R <- (sigma_G * k * r) / Y
  round(R, 4)
}


############################################################
# PART 4: EXAMPLE — SIMULATED DATA FOR ONE TRAIT (YIELD)
# Replace this section with your real trial data
############################################################

simulate_stage1_trial <- function(n_genotypes = 20, n_reps = 3, base_mean = 40, gen_sd = 4, error_sd = 3) {
  genotype_effects <- rnorm(n_genotypes, mean = 0, sd = gen_sd)
  df <- expand.grid(genotype = paste0("G", 1:n_genotypes), rep = paste0("R", 1:n_reps))
  df$genotype_effect <- genotype_effects[as.numeric(factor(df$genotype))]
  df$yield <- base_mean + df$genotype_effect + rnorm(nrow(df), 0, error_sd)
  df
}

# Simulate Stage 1 trials for three successive generations within one cycle
F2_trial <- simulate_stage1_trial(gen_sd = 3, error_sd = 4)
F3_trial <- simulate_stage1_trial(gen_sd = 4, error_sd = 3)
F4_trial <- simulate_stage1_trial(gen_sd = 5, error_sd = 2)

# Run Part 1 for each generation
stats_F2 <- compute_genetic_stats(F2_trial, "yield", "genotype", "rep")
stats_F3 <- compute_genetic_stats(F3_trial, "yield", "genotype", "rep")
stats_F4 <- compute_genetic_stats(F4_trial, "yield", "genotype", "rep")

summary_table <- bind_rows(
  cbind(Generation = "F2", stats_F2),
  cbind(Generation = "F3", stats_F3),
  cbind(Generation = "F4", stats_F4)
)

print(summary_table)


# Simulate realized gain data: current cycle vs. prior-cycle checks
realized_gain_data <- data.frame(
  yield = c(rnorm(30, 44, 3), rnorm(10, 40, 3)),
  cycle_group = c(rep("current_cycle", 30), rep("check_prior_cycle", 10))
)

gain_result <- compute_realized_gain(realized_gain_data, "yield", "cycle_group",
                                     "current_cycle", "check_prior_cycle")
print(gain_result)


# Response per year using F4 stats (final selection stage of the cycle)
R_value <- compute_response_per_year(GV = stats_F4$GV, H2_pct = stats_F4$H2_pct, Y = 3)
cat("Response per year (R):", R_value, "\n")


############################################################
# PART 5: VISUALIZATIONS
############################################################

# 5a. Genetic trend plot (population mean by generation)
trend_plot <- ggplot(summary_table, aes(x = Generation, y = Mean, group = 1)) +
  geom_line(color = "#1D9E75", linewidth = 1) +
  geom_point(size = 3, color = "#0F6E56") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.1) +
  labs(title = "Genetic trend across generations", y = "Trait mean", x = "Generation") +
  theme_minimal()
print(trend_plot)

# 5b. Overlapping histograms (distribution shift)
hist_data <- bind_rows(
  cbind(F2_trial, Generation = "F2") %>% select(yield, Generation),
  cbind(F3_trial, Generation = "F3") %>% select(yield, Generation),
  cbind(F4_trial, Generation = "F4") %>% select(yield, Generation)
)

hist_plot <- ggplot(hist_data, aes(x = yield, fill = Generation)) +
  geom_density(alpha = 0.4) +
  labs(title = "Trait distribution shift across generations", x = "Yield") +
  theme_minimal()
print(hist_plot)

# 5c. GAM bar chart
gam_plot <- ggplot(summary_table, aes(x = Generation, y = GAM_pct, fill = Generation)) +
  geom_col() +
  labs(title = "Genetic advance as % of mean (GAM)", y = "GAM (%)") +
  theme_minimal() +
  theme(legend.position = "none")
print(gam_plot)

# 5d. Heritability heatmap across generations (single trait shown as example;
#     extend rows with more traits or environments as needed)
heatmap_data <- summary_table %>% select(Generation, H2_pct)
heatmap_plot <- ggplot(heatmap_data, aes(x = Generation, y = "Yield", fill = H2_pct)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#FAEEDA", high = "#854F0B") +
  labs(title = "Heritability (H2 %) across generations", y = "", fill = "H2 (%)") +
  theme_minimal()
print(heatmap_plot)


############################################################
# PART 6: EXPORT SUMMARY TABLE
############################################################
write.csv(summary_table, "genetic_gain_summary_table.csv", row.names = FALSE)
}


############################################################
# BREEDING GAIN MODULE
#
# Source-safe functions for recurrent-selection and genetic
# gain summaries. This section is the active implementation.
#
# It is safe to source from app.R because it:
# - does not install packages,
# - does not clear the environment,
# - does not run simulations automatically,
# - does not write files automatically,
# - uses the breeding_ prefix to avoid helper-name conflicts.
############################################################

BREEDING_MODULE_VERSION <- "1.0.0"


# ---- 1. Utility helpers ----

breeding_check_packages <- function(required) {
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing required package(s): ",
      paste(missing, collapse = ", "),
      ". Install them before running this function.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

breeding_as_number <- function(x) {
  suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
}

breeding_selection_intensity_k <- function(selection_proportion = 0.05) {
  selection_proportion <- suppressWarnings(as.numeric(selection_proportion))
  if (
    length(selection_proportion) != 1 ||
      !is.finite(selection_proportion) ||
      selection_proportion <= 0 ||
      selection_proportion >= 1
  ) {
    stop("selection_proportion must be one number between 0 and 1.", call. = FALSE)
  }

  stats::dnorm(stats::qnorm(1 - selection_proportion)) / selection_proportion
}

breeding_require_columns <- function(data, columns) {
  columns <- columns[!is.null(columns)]
  missing <- columns[is.na(columns) | columns == "" | !columns %in% names(data)]
  if (length(missing) > 0) {
    stop(
      "Missing required column(s): ",
      paste(unique(missing), collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

breeding_harmonic_replication <- function(genotype_values) {
  reps_per_genotype <- table(genotype_values)
  if (length(reps_per_genotype) == 0 || any(reps_per_genotype <= 0)) {
    return(NA_real_)
  }
  1 / mean(1 / as.numeric(reps_per_genotype))
}

breeding_empty_trait_result <- function(trait, note) {
  data.frame(
    Trait = trait,
    N_obs = NA_integer_,
    N_genotypes = NA_integer_,
    N_replications = NA_integer_,
    Mean = NA_real_,
    SE = NA_real_,
    Min = NA_real_,
    Max = NA_real_,
    MS_genotype = NA_real_,
    MS_error = NA_real_,
    Harmonic_replication = NA_real_,
    Genotypic_variance = NA_real_,
    Error_variance = NA_real_,
    Phenotypic_variance = NA_real_,
    GCV_pct = NA_real_,
    PCV_pct = NA_real_,
    Broad_sense_H2 = NA_real_,
    Broad_sense_H2_pct = NA_real_,
    Selection_proportion_pct = NA_real_,
    Selection_intensity_k = NA_real_,
    Genetic_advance = NA_real_,
    Genetic_advance_pct_mean = NA_real_,
    Note = note,
    stringsAsFactors = FALSE
  )
}


# ---- 2. Variance, heritability, genetic advance ----

breeding_compute_genetic_stats <- function(
    data,
    trait,
    genotype,
    replication = NULL,
    selection_proportion = 0.05,
    selection_k = NULL) {
  breeding_require_columns(data, c(trait, genotype, replication))

  if (is.null(selection_k)) {
    selection_k <- breeding_selection_intensity_k(selection_proportion)
  }

  model_data <- data.frame(
    TraitValue = breeding_as_number(data[[trait]]),
    Genotype = as.character(data[[genotype]]),
    stringsAsFactors = FALSE
  )
  model_data$Replication <- if (is.null(replication)) {
    "R1"
  } else {
    as.character(data[[replication]])
  }

  keep <- !is.na(model_data$TraitValue) &
    !is.na(model_data$Genotype) &
    trimws(model_data$Genotype) != "" &
    !is.na(model_data$Replication) &
    trimws(model_data$Replication) != ""
  model_data <- model_data[keep, , drop = FALSE]
  model_data$Genotype <- factor(model_data$Genotype)
  model_data$Replication <- factor(model_data$Replication)

  if (nrow(model_data) < 3) {
    stop("At least three valid observations are required.", call. = FALSE)
  }
  if (length(unique(model_data$Genotype)) < 2) {
    stop("At least two genotypes are required.", call. = FALSE)
  }

  use_replication <- length(unique(model_data$Replication)) > 1
  model_formula <- if (use_replication) {
    TraitValue ~ Genotype + Replication
  } else {
    TraitValue ~ Genotype
  }

  model <- stats::aov(model_formula, data = model_data)
  anova_tab <- as.data.frame(summary(model)[[1]])
  anova_tab$Source <- trimws(rownames(anova_tab))
  genotype_row <- anova_tab[anova_tab$Source == "Genotype", , drop = FALSE]
  residual_row <- anova_tab[grepl("Residual", anova_tab$Source), , drop = FALSE]

  if (nrow(genotype_row) == 0 || nrow(residual_row) == 0) {
    stop("Genotype or residual mean square could not be estimated.", call. = FALSE)
  }
  if (is.na(residual_row$Df[1]) || residual_row$Df[1] <= 0) {
    stop("Residual degrees of freedom must be greater than zero.", call. = FALSE)
  }

  ms_genotype <- as.numeric(genotype_row$`Mean Sq`[1])
  ms_error <- as.numeric(residual_row$`Mean Sq`[1])
  harmonic_replication <- breeding_harmonic_replication(model_data$Genotype)

  genotypic_variance <- max((ms_genotype - ms_error) / harmonic_replication, 0)
  error_variance <- ms_error

  # Genotype-mean phenotypic variance. This matches selection decisions
  # made from replicated genotype means and avoids inflating error variance.
  phenotypic_variance <- genotypic_variance + (error_variance / harmonic_replication)
  broad_sense_h2 <- ifelse(
    is.finite(phenotypic_variance) && phenotypic_variance > 0,
    genotypic_variance / phenotypic_variance,
    NA_real_
  )

  trait_mean <- mean(model_data$TraitValue, na.rm = TRUE)
  trait_se <- stats::sd(model_data$TraitValue, na.rm = TRUE) / sqrt(nrow(model_data))
  trait_range <- range(model_data$TraitValue, na.rm = TRUE)
  mean_denominator <- ifelse(abs(trait_mean) < 0.0001, NA_real_, abs(trait_mean))

  gcv <- sqrt(genotypic_variance) / mean_denominator * 100
  pcv <- sqrt(phenotypic_variance) / mean_denominator * 100
  genetic_advance <- selection_k * sqrt(phenotypic_variance) * broad_sense_h2
  genetic_advance_pct <- genetic_advance / mean_denominator * 100

  data.frame(
    Trait = trait,
    N_obs = nrow(model_data),
    N_genotypes = length(unique(model_data$Genotype)),
    N_replications = length(unique(model_data$Replication)),
    Mean = round(trait_mean, 4),
    SE = round(trait_se, 4),
    Min = round(trait_range[1], 4),
    Max = round(trait_range[2], 4),
    MS_genotype = round(ms_genotype, 4),
    MS_error = round(ms_error, 4),
    Harmonic_replication = round(harmonic_replication, 4),
    Genotypic_variance = round(genotypic_variance, 6),
    Error_variance = round(error_variance, 6),
    Phenotypic_variance = round(phenotypic_variance, 6),
    GCV_pct = round(gcv, 2),
    PCV_pct = round(pcv, 2),
    Broad_sense_H2 = round(broad_sense_h2, 4),
    Broad_sense_H2_pct = round(broad_sense_h2 * 100, 2),
    Selection_proportion_pct = round(selection_proportion * 100, 2),
    Selection_intensity_k = round(selection_k, 4),
    Genetic_advance = round(genetic_advance, 4),
    Genetic_advance_pct_mean = round(genetic_advance_pct, 2),
    Note = "Replicated genotype-mean estimate",
    stringsAsFactors = FALSE
  )
}

breeding_compute_multi_trait_stats <- function(
    data,
    traits,
    genotype,
    replication = NULL,
    selection_proportion = 0.05,
    selection_k = NULL,
    fail_silently = TRUE) {
  breeding_require_columns(data, c(traits, genotype, replication))

  rows <- lapply(traits, function(trait_name) {
    tryCatch(
      breeding_compute_genetic_stats(
        data = data,
        trait = trait_name,
        genotype = genotype,
        replication = replication,
        selection_proportion = selection_proportion,
        selection_k = selection_k
      ),
      error = function(e) {
        if (!isTRUE(fail_silently)) {
          stop(e, call. = FALSE)
        }
        breeding_empty_trait_result(trait_name, paste("Not calculated:", e$message))
      }
    )
  })

  do.call(rbind, rows)
}


# ---- 3. Realized gain against prior-cycle checks ----

breeding_compute_realized_gain <- function(
    data,
    trait,
    group_col,
    current_label,
    check_label,
    higher_is_better = TRUE) {
  breeding_require_columns(data, c(trait, group_col))

  gain_data <- data.frame(
    TraitValue = breeding_as_number(data[[trait]]),
    Group = as.character(data[[group_col]]),
    stringsAsFactors = FALSE
  )
  keep <- !is.na(gain_data$TraitValue) &
    !is.na(gain_data$Group) &
    trimws(gain_data$Group) != ""
  gain_data <- gain_data[keep, , drop = FALSE]

  if (!current_label %in% gain_data$Group) {
    stop("current_label was not found in group_col.", call. = FALSE)
  }
  if (!check_label %in% gain_data$Group) {
    stop("check_label was not found in group_col.", call. = FALSE)
  }

  group_mean <- stats::aggregate(TraitValue ~ Group, data = gain_data, FUN = mean)
  group_n <- stats::aggregate(TraitValue ~ Group, data = gain_data, FUN = length)
  names(group_n)[2] <- "N"
  means <- merge(group_mean, group_n, by = "Group", all.x = TRUE)

  current_row <- means[means$Group == current_label, , drop = FALSE]
  check_row <- means[means$Group == check_label, , drop = FALSE]
  raw_difference <- current_row$TraitValue[1] - check_row$TraitValue[1]
  improvement <- if (isTRUE(higher_is_better)) raw_difference else -raw_difference
  denominator <- ifelse(abs(check_row$TraitValue[1]) < 0.0001, 1, abs(check_row$TraitValue[1]))

  data.frame(
    Trait = trait,
    Current_label = current_label,
    Check_label = check_label,
    Current_N = current_row$N[1],
    Check_N = check_row$N[1],
    Current_Mean = round(current_row$TraitValue[1], 4),
    Check_Mean = round(check_row$TraitValue[1], 4),
    Raw_difference_current_minus_check = round(raw_difference, 4),
    Improvement = round(improvement, 4),
    Improvement_pct_of_check = round(improvement / denominator * 100, 2),
    Direction = ifelse(isTRUE(higher_is_better), "Higher better", "Lower better"),
    stringsAsFactors = FALSE
  )
}


# ---- 4. Expected response per year ----

breeding_compute_response_per_year <- function(
    genotypic_variance,
    h2,
    selection_proportion = 0.05,
    selection_k = NULL,
    years_per_cycle = 3) {
  genotypic_variance <- suppressWarnings(as.numeric(genotypic_variance))
  h2 <- suppressWarnings(as.numeric(h2))
  years_per_cycle <- suppressWarnings(as.numeric(years_per_cycle))

  if (is.null(selection_k)) {
    selection_k <- breeding_selection_intensity_k(selection_proportion)
  }
  if (!is.finite(genotypic_variance) || genotypic_variance < 0) {
    return(NA_real_)
  }
  if (!is.finite(h2)) {
    return(NA_real_)
  }
  if (h2 > 1) {
    h2 <- h2 / 100
  }
  if (h2 < 0 || h2 > 1) {
    return(NA_real_)
  }
  if (!is.finite(years_per_cycle) || years_per_cycle <= 0) {
    stop("years_per_cycle must be greater than zero.", call. = FALSE)
  }

  sigma_g <- sqrt(genotypic_variance)
  selection_accuracy <- sqrt(h2)
  round((sigma_g * selection_k * selection_accuracy) / years_per_cycle, 6)
}

breeding_compute_response_table <- function(
    genetic_stats,
    years_per_cycle = 3,
    h2_col = "Broad_sense_H2",
    variance_col = "Genotypic_variance") {
  breeding_require_columns(genetic_stats, c("Trait", h2_col, variance_col))

  selection_k <- if ("Selection_intensity_k" %in% names(genetic_stats)) {
    genetic_stats$Selection_intensity_k
  } else {
    rep(NA_real_, nrow(genetic_stats))
  }

  response <- mapply(
    function(genotypic_variance, h2, row_selection_k) {
      k_used <- if (is.finite(row_selection_k)) row_selection_k else NULL
      breeding_compute_response_per_year(
        genotypic_variance = genotypic_variance,
        h2 = h2,
        selection_k = k_used,
        years_per_cycle = years_per_cycle
      )
    },
    genotypic_variance = genetic_stats[[variance_col]],
    h2 = genetic_stats[[h2_col]],
    row_selection_k = selection_k
  )

  data.frame(
    Trait = genetic_stats$Trait,
    Genotypic_variance = genetic_stats[[variance_col]],
    Broad_sense_H2 = genetic_stats[[h2_col]],
    Selection_intensity_k = selection_k,
    Years_per_cycle = years_per_cycle,
    Response_per_year = as.numeric(response),
    stringsAsFactors = FALSE
  )
}


# ---- 5. Combined breeding gain pipeline ----

breeding_run_gain_pipeline <- function(
    data,
    traits,
    genotype,
    replication = NULL,
    cycle_group_col = NULL,
    current_label = NULL,
    check_label = NULL,
    selection_proportion = 0.05,
    selection_k = NULL,
    years_per_cycle = 3,
    higher_is_better = TRUE) {
  stats_table <- breeding_compute_multi_trait_stats(
    data = data,
    traits = traits,
    genotype = genotype,
    replication = replication,
    selection_proportion = selection_proportion,
    selection_k = selection_k
  )

  response_table <- breeding_compute_response_table(
    stats_table,
    years_per_cycle = years_per_cycle
  )

  realized_gain <- data.frame()
  if (
    !is.null(cycle_group_col) &&
      !is.null(current_label) &&
      !is.null(check_label)
  ) {
    realized_gain <- do.call(rbind, lapply(traits, function(trait_name) {
      tryCatch(
        breeding_compute_realized_gain(
          data = data,
          trait = trait_name,
          group_col = cycle_group_col,
          current_label = current_label,
          check_label = check_label,
          higher_is_better = higher_is_better
        ),
        error = function(e) {
          data.frame(
            Trait = trait_name,
            Current_label = current_label,
            Check_label = check_label,
            Current_N = NA_integer_,
            Check_N = NA_integer_,
            Current_Mean = NA_real_,
            Check_Mean = NA_real_,
            Raw_difference_current_minus_check = NA_real_,
            Improvement = NA_real_,
            Improvement_pct_of_check = NA_real_,
            Direction = ifelse(isTRUE(higher_is_better), "Higher better", "Lower better"),
            Note = paste("Not calculated:", e$message),
            stringsAsFactors = FALSE
          )
        }
      )
    }))
  }

  list(
    module_version = BREEDING_MODULE_VERSION,
    genetic_stats = stats_table,
    response_per_year = response_table,
    realized_gain = realized_gain
  )
}


# ---- 6. Optional generation-based summary ----

breeding_compute_generation_stats <- function(
    data,
    trait,
    genotype,
    replication,
    generation_col,
    selection_proportion = 0.05,
    selection_k = NULL) {
  breeding_require_columns(data, c(trait, genotype, replication, generation_col))

  generations <- unique(as.character(data[[generation_col]]))
  generations <- generations[!is.na(generations) & trimws(generations) != ""]

  rows <- lapply(generations, function(generation_name) {
    subset_data <- data[as.character(data[[generation_col]]) == generation_name, , drop = FALSE]
    out <- tryCatch(
      breeding_compute_genetic_stats(
        data = subset_data,
        trait = trait,
        genotype = genotype,
        replication = replication,
        selection_proportion = selection_proportion,
        selection_k = selection_k
      ),
      error = function(e) {
        breeding_empty_trait_result(trait, paste("Not calculated:", e$message))
      }
    )
    data.frame(Generation = generation_name, out, stringsAsFactors = FALSE)
  })

  do.call(rbind, rows)
}


# ---- 7. Plot helpers ----

breeding_plot_genetic_trend <- function(
    stats_table,
    x_col = "Generation",
    mean_col = "Mean",
    se_col = "SE",
    title = "Genetic trend across generations") {
  breeding_check_packages("ggplot2")
  breeding_require_columns(stats_table, c(x_col, mean_col, se_col))

  ggplot2::ggplot(
    stats_table,
    ggplot2::aes(x = .data[[x_col]], y = .data[[mean_col]], group = 1)
  ) +
    ggplot2::geom_line(color = "#1D9E75", linewidth = 1) +
    ggplot2::geom_point(size = 3, color = "#0F6E56") +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = .data[[mean_col]] - .data[[se_col]],
        ymax = .data[[mean_col]] + .data[[se_col]]
      ),
      width = 0.1
    ) +
    ggplot2::labs(title = title, x = x_col, y = "Trait mean") +
    ggplot2::theme_minimal()
}

breeding_plot_distribution_shift <- function(
    data,
    trait,
    generation_col,
    title = "Trait distribution shift and realized mean gain") {
  breeding_check_packages("ggplot2")
  breeding_require_columns(data, c(trait, generation_col))

  plot_data <- data.frame(
    TraitValue = breeding_as_number(data[[trait]]),
    Generation = as.character(data[[generation_col]]),
    stringsAsFactors = FALSE
  )
  plot_data <- plot_data[!is.na(plot_data$TraitValue) & !is.na(plot_data$Generation), ]
  if (nrow(plot_data) == 0) {
    stop("No valid values are available for the selected trait and generation column.", call. = FALSE)
  }

  generation_levels <- unique(plot_data$Generation)
  plot_data$Generation <- factor(plot_data$Generation, levels = generation_levels)

  mean_table <- stats::aggregate(
    TraitValue ~ Generation,
    data = plot_data,
    FUN = mean
  )
  names(mean_table)[2] <- "Mean"
  mean_table$Generation <- factor(mean_table$Generation, levels = generation_levels)
  mean_table$Mean_label <- paste0(
    "Average\n",
    mean_table$Generation,
    "\n",
    round(mean_table$Mean, 2)
  )

  density_values <- lapply(split(plot_data$TraitValue, plot_data$Generation), function(x) {
    x <- x[is.finite(x)]
    if (length(unique(x)) < 2) {
      return(0)
    }
    stats::density(x)$y
  })
  max_density <- max(unlist(density_values), na.rm = TRUE)
  if (!is.finite(max_density) || max_density <= 0) {
    max_density <- 1
  }

  gain_segments <- data.frame(
    From = character(),
    To = character(),
    Mean_from = numeric(),
    Mean_to = numeric(),
    Gain = numeric(),
    Label_x = numeric(),
    Arrow_y = numeric(),
    Label_y = numeric(),
    Gain_label = character(),
    stringsAsFactors = FALSE
  )
  if (nrow(mean_table) >= 2) {
    gain_segments <- data.frame(
      From = as.character(mean_table$Generation[-nrow(mean_table)]),
      To = as.character(mean_table$Generation[-1]),
      Mean_from = mean_table$Mean[-nrow(mean_table)],
      Mean_to = mean_table$Mean[-1],
      stringsAsFactors = FALSE
    )
    gain_segments$Gain <- gain_segments$Mean_to - gain_segments$Mean_from
    gain_segments$Label_x <- (gain_segments$Mean_from + gain_segments$Mean_to) / 2
    gain_segments$Arrow_y <- max_density * (0.10 + 0.10 * ((seq_len(nrow(gain_segments)) - 1) %% 2))
    gain_segments$Label_y <- gain_segments$Arrow_y + max_density * 0.06
    gain_segments$Gain_label <- paste0(
      gain_segments$From,
      " to ",
      gain_segments$To,
      "\nGain = ",
      round(gain_segments$Gain, 2)
    )
  }

  palette <- grDevices::colorRampPalette(c("#DDEFD4", "#8FCF76", "#43A947"))(
    max(length(generation_levels), 3)
  )[seq_along(generation_levels)]

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$TraitValue, fill = .data$Generation, color = .data$Generation)
  ) +
    ggplot2::geom_density(alpha = 0.30, linewidth = 0.9) +
    ggplot2::geom_vline(
      data = mean_table,
      ggplot2::aes(xintercept = .data$Mean, color = .data$Generation),
      linetype = "dashed",
      linewidth = 0.8,
      show.legend = FALSE
    ) +
    ggplot2::geom_label(
      data = mean_table,
      ggplot2::aes(x = .data$Mean, y = max_density * 1.12, label = .data$Mean_label),
      inherit.aes = FALSE,
      size = 3.4,
      label.size = 0,
      fill = "white",
      color = "#263123"
    ) +
    ggplot2::geom_segment(
      data = gain_segments,
      ggplot2::aes(
        x = .data$Mean_from,
        xend = .data$Mean_to,
        y = .data$Arrow_y,
        yend = .data$Arrow_y
      ),
      inherit.aes = FALSE,
      color = "#E85D24",
      linewidth = 1,
      arrow = grid::arrow(length = grid::unit(0.22, "cm"), type = "closed")
    ) +
    ggplot2::geom_label(
      data = gain_segments,
      ggplot2::aes(x = .data$Label_x, y = .data$Label_y, label = .data$Gain_label),
      inherit.aes = FALSE,
      size = 3.3,
      label.size = 0,
      fill = "#FFF3EA",
      color = "#7A2D12"
    ) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.28))) +
    ggplot2::labs(
      title = title,
      subtitle = "Dashed lines show cycle averages. Orange arrows show average gain from one cycle to the next.",
      x = trait,
      y = "Population density",
      fill = generation_col,
      color = generation_col
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle = ggplot2::element_text(color = "gray35")
    )
}

breeding_plot_gam <- function(
    stats_table,
    x_col = "Generation",
    gam_col = "Genetic_advance_pct_mean",
    title = "Genetic advance as percent of mean") {
  breeding_check_packages("ggplot2")
  breeding_require_columns(stats_table, c(x_col, gam_col))

  ggplot2::ggplot(
    stats_table,
    ggplot2::aes(x = .data[[x_col]], y = .data[[gam_col]], fill = .data[[x_col]])
  ) +
    ggplot2::geom_col() +
    ggplot2::labs(title = title, x = x_col, y = "Genetic advance (% of mean)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}

breeding_plot_heritability_heatmap <- function(
    stats_table,
    stage_col = "Generation",
    trait_col = "Trait",
    h2_col = "Broad_sense_H2_pct",
    title = "Broad-sense heritability") {
  breeding_check_packages("ggplot2")
  breeding_require_columns(stats_table, c(stage_col, trait_col, h2_col))

  ggplot2::ggplot(
    stats_table,
    ggplot2::aes(x = .data[[stage_col]], y = .data[[trait_col]], fill = .data[[h2_col]])
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = round(.data[[h2_col]], 1)), size = 3.5) +
    ggplot2::scale_fill_gradient(low = "#FAEEDA", high = "#854F0B", na.value = "grey90") +
    ggplot2::labs(title = title, x = stage_col, y = trait_col, fill = "H2 (%)") +
    ggplot2::theme_minimal()
}


# ---- 8. Optional simulation helper ----

breeding_simulate_stage_trial <- function(
    n_genotypes = 20,
    n_reps = 3,
    base_mean = 40,
    gen_sd = 4,
    error_sd = 3,
    seed = NULL,
    trait_name = "yield") {
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    on.exit({
      if (!is.null(old_seed)) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  genotype_effects <- stats::rnorm(n_genotypes, mean = 0, sd = gen_sd)
  out <- expand.grid(
    genotype = paste0("G", seq_len(n_genotypes)),
    rep = paste0("R", seq_len(n_reps)),
    stringsAsFactors = FALSE
  )
  out[[trait_name]] <- base_mean +
    genotype_effects[as.numeric(factor(out$genotype))] +
    stats::rnorm(nrow(out), 0, error_sd)
  out
}
