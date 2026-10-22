suppressMessages(invisible(source("app.R")))
dir.create("tmp/mgidi-fa", recursive = TRUE, showWarnings = FALSE)

set.seed(42)
genotypes <- sprintf("G%02d", 1:16)
traits <- c("Yield", "Fruit", "Brix", "Disease", "Maturity")
values <- matrix(rnorm(80), 16, 5, dimnames = list(genotypes, traits))
values[, 2] <- values[, 1] * 0.7 + values[, 2] * 0.3
trial <- data.frame(Variety = genotypes, Rep = 1, Environment = "E1", values)
metadata <- trial[1:2, ]
metadata$Variety <- c("WEIGHT", "DIRECTION")
metadata[1, traits] <- c(5, 2, 1, 4, 2)
metadata[2, traits] <- c("Higher better", "Higher better", "Higher better", "Lower better", "Target 0")
trial <- rbind(metadata, trial)
results <- setNames(lapply(seq_along(traits), function(i) list(
  blups_main = data.frame(Genotype = genotypes, BLUP_G = values[, i], Reliability = 0.9),
  trait_direction = if (i == 4) "Lower better" else if (i == 5) "Target trait" else "Higher better",
  target_value = if (i == 5) 0 else NA_real_)), traits)
info <- met_ds_trait_matrix(trial, results)
fa <- met_mgidi_fa(info, 20)
stopifnot(nrow(fa$ranking) == 16, all(is.finite(fa$ranking$MGIDI_distance)),
          sum(fa$ranking$Selected) == 4, nrow(fa$loadings) == 5,
          any(fa$variance$Retained), nrow(fa$scores) == 16,
          isTRUE(all.equal(unname(info$adjusted_wide$Disease), unname(-values[, 4]))),
          isTRUE(all.equal(unname(info$adjusted_wide$Maturity), unname(-abs(values[, 5])))))
reference <- metan::mgidi(info$adjusted_wide, ideotype = rep("h", 5), SI = 100, verbose = FALSE)
stopifnot(isTRUE(all.equal(fa$ranking$MGIDI_distance,
  reference$MGIDI$MGIDI[match(fa$ranking$Genotype, reference$MGIDI$Genotype)])))
plot <- plot_met_mgidi_ranking(list(mgidi_fa = fa), use_fa = TRUE)
stopifnot(nrow(ggplot2::ggplot_build(plot)$data[[3]]) == 16)
ggplot2::ggsave("tmp/mgidi-fa/chart.png", plot, width = 11, height = 8, dpi = 110)
info$adjusted_wide$Disease[1] <- NA_real_
sparse <- met_mgidi_fa(info)
stopifnot(nrow(sparse$ranking) == 16, is.na(sparse$ranking$MGIDI_distance[sparse$ranking$Genotype == "G01"]),
  any(grepl("G01", sparse$notes$Note)))
info$adjusted_wide[, traits] <- 1
unavailable <- met_mgidi_fa(info)
stopifnot(all(is.na(unavailable$ranking$MGIDI_rank)), any(unavailable$notes$Status == "Not available"))

# All supported model choices feed the same FA calculation and export contract.
set.seed(19)
trial <- expand.grid(Variety = genotypes, Rep = as.character(1:2), Environment = c("E1", "E2", "E3"))
for (i in seq_along(traits)) trial[[traits[i]]] <- 20 + values[match(trial$Variety, genotypes), i] * 3 + rnorm(nrow(trial), sd = 0.5)
# The optional battery of external stability methods is unrelated to FA ranking.
si_ext_run_metan_met <- function(...) list()
for (model_type in c("LMM", "ANOVA_RCBD", "SOMMER")) {
  result <- run_met_all_traits(trial, check_varieties = "G01", trait_cols = traits,
    replication_col = "Rep", min_envs_for_biplot = 4, model_type = model_type)
  fa <- result$met_multitrait_methods$mgidi_fa
  stopifnot(nrow(fa$ranking) == 16, sum(is.finite(fa$ranking$MGIDI_distance)) >= 3,
    all(fa$notes$Status == "OK"))
  for (genotype in fa$ranking$Genotype[!is.finite(fa$ranking$MGIDI_distance)]) {
    stopifnot(any(grepl(genotype, fa$notes$Note, fixed = TRUE)))
  }
  exports <- build_export_tables("MET", result)
  stopifnot(any(grepl("FA.*ranking", names(exports))), any(grepl("FA.*loadings", names(exports))),
    any(grepl("FA.*notes", names(exports))))
}
cat("MGIDI FA reference, direction, missing-data, chart and all-model export checks passed.\n")
