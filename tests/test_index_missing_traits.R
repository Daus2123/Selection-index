suppressMessages(invisible(source("app.R")))
expect_missing <- function(expr, text = "missing required traits") {
  msg <- tryCatch({ force(expr); "NO ERROR" }, error = conditionMessage)
  stopifnot(grepl(text, msg, fixed = TRUE))
}
set.seed(84)
d <- expand.grid(Variety = LETTERS[1:5], Rep = 1:4, stringsAsFactors = FALSE)
d$Yield <- 10 + match(d$Variety, LETTERS) + rnorm(nrow(d))
d$Height <- 30 + 2 * match(d$Variety, LETTERS) + rnorm(nrow(d))
partial <- d; partial$Yield[1] <- NA_real_
si_require_index_coverage(partial, id_col, c("Yield", "Height"))
complete <- run_selection_pipeline(d, "A", "RCBD")
result <- run_selection_pipeline(partial, "A", "RCBD")
stopifnot(all(is.finite(result$index_ranking$Selection_Index)),
          nrow(result$index_ranking) == 5)
missing <- d; missing$Yield[missing$Variety == "B"] <- NA_real_
expect_missing(run_selection_pipeline(missing, "A", "RCBD"))
all_missing <- d; all_missing$Yield <- NA_real_
expect_missing(run_selection_pipeline(all_missing, "A", "RCBD"))
si_require_index_coverage(missing, id_col, c("Yield", "Height"), c(Yield = 0, Height = 1))
stopifnot(identical(standardize_trait(c(5, NA, 5)), c(0, NA_real_, 0)),
          all(is.na(standardize_trait(c(NA_real_, NA_real_)))))
met <- d; met$Environment <- "E1"
fake <- function(values) list(blups_main = data.frame(Genotype = LETTERS[1:5], BLUP_G = values),
  fw_results = data.frame(Genotype = LETTERS[1:5], Sens = 1), ammi_genotype = data.frame())
models <- list(Yield = fake(1:5), Height = fake(5:1))
ranked <- build_met_integrated_ranking(met, models)
stopifnot(nrow(ranked$ranking) == 5, all(ranked$ranking$Weight_coverage == 1))
models$Yield$blups_main$BLUP_G[2] <- NA_real_
expect_missing(build_met_integrated_ranking(met, models))
expect_missing(build_met_integrated_ranking(met, models["Height"]), "required MET trait models are missing")
cat("Missing-trait selection index checks passed.\n")
