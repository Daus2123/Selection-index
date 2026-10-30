source("modules/module_5_met.R")
results <- list(met_trait_names = c("Failed", "Unavailable"),
  met_failed_traits = data.frame(Trait = "Failed", Error = "Model failed"))
combined <- met_combined_comparison_tables(results)
stopifnot(all(c("Failed | Adjusted_mean_SE", "Unavailable | Adjusted_mean_SE") %in% names(combined$Means)),
  combined$Means[["Failed | Adjusted_mean_SE"]][1] == "Model failed",
  combined$ANOVA$Status[1] == "Model failed",
  grepl("unavailable", combined$Means[["Unavailable | Adjusted_mean_SE"]][1], ignore.case = TRUE))
cat("Failed and unavailable traits remain in both combined tables.\n")
