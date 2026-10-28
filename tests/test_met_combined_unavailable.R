source("modules/module_5_met.R")
results <- list(met_trait_names = c("Failed", "Unavailable"),
  met_failed_traits = data.frame(Trait = "Failed", Error = "Model failed"))
combined <- met_combined_comparison_tables(results)
stopifnot(identical(combined$Means$Trait, c("Failed", "Unavailable")),
  combined$Means$Status[1] == "Model failed",
  combined$ANOVA$Status[1] == "Model failed",
  grepl("unavailable", combined$Means$Status[2]))
cat("Failed and unavailable traits remain in both combined tables.\n")
