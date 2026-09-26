suppressMessages(invisible(source("app.R")))
trial <- expand.grid(Variety = c("A", "B", "C"), Rep = 1:2,
                     Environment = c("E1", "E2"))
trial$Yield <- rep(c(10, 20, 30), 4)
fake_result <- list(blups_main = data.frame(Genotype = c("A", "B", "C"), BLUP_G = c(10, 20, 30)),
  fw_results = data.frame(Genotype = c("A", "B", "C"), Sens = c(1, 1.1, NA)),
  ammi_genotype = data.frame(Genotype = c("A", "B", "C"), ASV = c(1, NA, NA)))
for (asv in list(c(1, NA, NA), c(NA, NA, NA), c(Inf, 2, NA), c(2, 2, NA), c(1, 2, 3))) {
  fake_result$ammi_genotype$ASV <- asv
  result <- build_met_integrated_ranking(trial, list(Yield = fake_result),
                                        c(mean = 1, fw = 0, asv = 0))
  stopifnot(nrow(result$ranking) == 3, all(is.finite(result$ranking$Integrated_MET_Index)),
            identical(result$ranking$Genotype, c("C", "B", "A")))
}
stopifnot(identical(met_standardize_component(c(1, NA, Inf)), c(0, NA_real_, NA_real_)),
          all(is.na(met_standardize_component(c(NA, Inf)))),
          identical(met_standardize_component(c(2, 2, NA)), c(0, 0, NA_real_)))
fake_result$ammi_genotype$ASV <- NA_real_
mean_only <- build_met_integrated_ranking(trial, list(Yield = fake_result), c(mean = 1, fw = 0, asv = 0))
missing_asv <- build_met_integrated_ranking(trial, list(Yield = fake_result), c(mean = 1, fw = 0, asv = 1))
stopifnot(identical(mean_only$ranking$Integrated_MET_Index, missing_asv$ranking$Integrated_MET_Index),
          grepl("ASV: 0/3 usable", missing_asv$trait_weights$Component_notes, fixed = TRUE))
cat("Sparse MET stability component checks passed.\n")
