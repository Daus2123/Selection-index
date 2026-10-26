suppressMessages(invisible(source("app.R")))
set.seed(92)
trial <- expand.grid(Variety = c("A", "B", "C", "D"), Rep = as.character(1:4),
  Location = c("E1", "E2", "E3"), stringsAsFactors = FALSE)
trial$Entry <- seq_len(nrow(trial))
trial$ID <- trial$Entry
trial$Type <- ifelse(trial$Variety == "A", "Check", "Candidate")
trial$Score <- rep(c(2, 3.5, 4, 4.5), length.out = nrow(trial))
trial$Yield <- 10 + match(trial$Variety, c("A", "B", "C", "D")) + rnorm(nrow(trial), sd = .5)
trial$Count <- rep(c(1, 2, 3, 4), length.out = nrow(trial))
meta <- trial[rep(1, 4), ]; meta[,] <- ""
meta$Variety <- c("DIRECTION", "WEIGHT", "TYPE", "UNIT")
meta$Score <- c("high", "2", "score", "score")
meta$Yield <- c("high", "3", "continuous", "kg")
meta$Count <- c("high", "1", "quantitative", "number")
book <- rbind(trial, meta)
prepared <- prepare_excel_input(book)
stopifnot(nrow(prepared$data) == nrow(trial),
  setequal(prepared$trait_cols, c("Score", "Yield", "Count")),
  identical(prepared$trait_metadata$Type, c("Score 1-5", "Quantitative", "Quantitative")),
  identical(prepared$trait_metadata$Unit, c("score", "kg", "number")),
  !any(c("Entry", "ID", "Location") %in% prepared$trait_cols))
met_prepared <- prepare_met_trait_settings(book)
stopifnot(nrow(met_prepared$data) == nrow(trial),
  identical(met_prepared$trait_metadata, prepared$trait_metadata),
  setequal(met_prepared$trait_cols, prepared$trait_cols))
stopifnot(identical(to_number(c("2500", "2.5", "2,5", "1,250", "1.250", "1,250.5")),
  c(2500, 2.5, 2.5, 1.25, 1.25, NA_real_)))
# Legacy UNIT-only files and type aliases remain supported.
unit_only <- book[book$Variety != "TYPE", ]
stopifnot(si_trait_metadata(unit_only, "Score")$Type == "Score 1-5")
unit_only$Yield[unit_only$Variety == "UNIT"] <- "continues"
stopifnot(si_trait_metadata(unit_only, "Yield")$Type == "Quantitative")
bad <- book; bad$Score[1] <- 6
review <- si_clean_screen(bad, "Score", id_col, c(id_col, rep_col), "iqr", 3)$review
stopifnot(grepl("Outside declared score range", review$Flag[1]),
  !any(si_clean_screen(book, "Score", id_col, c(id_col, rep_col), "iqr", 3)$review$Row > nrow(trial)))
duplicate <- rbind(book, meta[meta$Variety == "TYPE", ])
stopifnot(inherits(try(prepare_excel_input(duplicate), silent = TRUE), "try-error"))
conflict <- book; conflict$Score[conflict$Variety == "TYPE"] <- "quantitative"
stopifnot(inherits(try(prepare_excel_input(conflict), silent = TRUE), "try-error"))
single <- book[book$Location == "E1" | book$Variety %in% metadata_row_labels, ]
lpsi <- run_selection_pipeline(single, check_varieties = "A", model_type = "RCBD")
stopifnot(lpsi$trait_info$Type[lpsi$trait_info$Trait == "Score"] == "Score 1-5",
  lpsi$trait_info$Type[lpsi$trait_info$Trait == "Count"] == "Quantitative",
  lpsi$trait_info$Unit[lpsi$trait_info$Trait == "Yield"] == "kg",
  !any(metadata_row_labels %in% lpsi$cleaned_data$Original_ID))
stopifnot(any(grepl("trait_metadata", names(build_export_tables("LPSI", lpsi)))))
for (model_type in c("LMM", "ANOVA_RCBD", "SOMMER")) {
  result <- run_met_pipeline(book, "Yield", check_varieties = "A", replication_col = "Rep",
    min_envs_for_biplot = 4, model_type = model_type)
  stopifnot(nrow(result$met_data) == nrow(trial),
    !any(metadata_row_labels %in% as.character(result$met_data$Genotype)))
}
cat("TYPE/UNIT metadata, explicit types, decimal parsing, LPSI and all MET model checks passed.\n")
