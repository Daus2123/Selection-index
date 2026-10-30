suppressMessages(invisible(source("app.R")))
set.seed(619)
d <- expand.grid(Genotype = paste0("G", 1:8), Environment = paste0("E", 1:4), Rep = factor(1:3))
missing <- expand.grid(Genotype = paste0("G", 4:8), Environment = c("E1", "E2"))
missing <- rbind(missing, data.frame(Genotype = paste0("G", 6:8), Environment = "E3"))
key <- function(x) paste(x$Genotype, x$Environment)
d <- d[!key(d) %in% key(missing), ]
d$Genotype <- factor(d$Genotype); d$Environment <- factor(d$Environment)
d$Weight <- 10 + as.integer(d$Genotype) + rnorm(nrow(d), sd = 0.5)
check_result <- function(result) {
  stopifnot(all(c("anova", "means", "pairwise", "simple_tests", "notes") %in% names(result)),
    any(!result$means$Estimable), any(is.finite(result$means$Adjusted_mean)),
    any(grepl("rank deficient", result$notes$Note)),
    setequal(result$means$View, unname(met_mean_comparison_views())))
  cells <- subset(result$means, View == "interaction" & Method == "tukey")
  stopifnot(nrow(cells) == 32, sum(!cells$Estimable) == nrow(missing))
  stopifnot(all(cells$Group[!cells$Estimable] == ""))
  stopifnot(all(result$means$Method == "tukey"), all(result$pairwise$Method == "tukey"))
  nonestimable_pairs <- result$pairwise[!is.finite(result$pairwise$p.value), ]
  stopifnot(all(is.na(nonestimable_pairs$Significant_0.05)))
}
previous_show_error <- getOption("show.error.messages")
blup <- met_build_mean_comparison(d, TRUE, FALSE)
stopifnot(all(blup$means$Estimable), sum(subset(blup$means, View == "interaction")$Evidence_status == "Untested - predicted") == nrow(missing))
stopifnot(identical(getOption("show.error.messages"), previous_show_error))
fixed <- lm(Weight ~ Environment + Environment:Rep + Genotype * Environment, data = d)
check_result(met_build_mean_comparison(d, TRUE, FALSE, fitted_model = fixed))
# Reproduce the exact matrix-inversion failure at the ANOVA boundary.
anova.met_singular_test <- function(object, ...) {
  stop("system is computationally singular: reciprocal condition number = 2.33266e-20")
}
class(fixed) <- c("met_singular_test", class(fixed))
recovered <- met_build_mean_comparison(d, TRUE, FALSE, fitted_model = fixed)
check_result(recovered)
stopifnot(any(grepl("computationally singular", recovered$notes$Note)),
  any(grepl("Not estimable", recovered$anova$Status)))
results <- list(met_trait_names = "Yield", met_by_trait = list(Yield = list(mean_comparison = recovered)))
wide <- met_combined_comparison_tables(results)$Means
stopifnot(nrow(wide) > 0, "Yield | Adjusted_mean_SE" %in% names(wide))
cat("Sparse comparisons preserve all cells and survive singular ANOVA tests.\n")
