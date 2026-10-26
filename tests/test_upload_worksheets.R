suppressMessages(invisible(source("app.R")))
trial <- expand.grid(Variety = c("A", "B", "C"), Rep = 1:3, stringsAsFactors = FALSE)
trial$Yield <- seq_len(nrow(trial)) + 0.5
legacy_rows <- trial[rep(1, 2), ]; legacy_rows[,] <- ""
legacy_rows$Variety <- c("DIRECTION", "WEIGHT"); legacy_rows$Yield <- c("high", "3")
legacy <- rbind(trial, legacy_rows)
updated <- trial; updated$ID <- seq_len(nrow(updated))
meta <- updated[rep(1, 4), ]; meta[,] <- ""
meta$Variety <- c("DIRECTION", "WEIGHT", "TYPE", "UNIT")
meta$Yield <- c("high", "3", "quantitative", "kg")
updated <- rbind(updated, meta)
path <- tempfile(fileext = ".xlsx")
writexl::write_xlsx(list(Legacy = legacy, Updated = updated), path)
first <- si_read_excel_upload(path, "synthetic.xlsx")
second <- si_read_excel_upload(path, "synthetic.xlsx", "Updated")
stopifnot(!"ID" %in% names(first), "ID" %in% names(second),
  nrow(second) == nrow(updated), "UNIT" %in% second$Variety,
  prepare_excel_input(second)$trait_metadata$Unit == "kg")
shiny::testServer(server, {
  session$setInputs(excel_file = list(name = "synthetic.xlsx", datapath = path, size = file.info(path)$size))
  stopifnot(selected_upload_sheet() == "Legacy", !"ID" %in% names(uploaded_data()),
    !grepl("UNIT rows detected", output$upload_message),
    !grepl("Worksheet:", output$upload_message))
  session$setInputs(excel_sheet = "Updated")
  stopifnot(selected_upload_sheet() == "Legacy", !"ID" %in% names(uploaded_data()),
    grepl("Scoring traits \\(1-5\\): 0", output$upload_message),
    grepl("Quantitative traits: 1", output$upload_message))
  cleaning$has_approved <- TRUE
  cleaning$settings <- data.frame(Trait = "Yield")
  next_path <- tempfile(fileext = ".xlsx")
  updated$Yield[updated$Variety == "TYPE"] <- "score"
  updated$Yield[updated$Variety == "UNIT"] <- "score"
  writexl::write_xlsx(list(Updated = updated, Legacy = legacy), next_path)
  session$setInputs(excel_file = list(name = "updated.xlsx", datapath = next_path, size = file.info(next_path)$size))
  stopifnot(!cleaning$has_approved, nrow(cleaning$settings) == 0)
  stopifnot(selected_upload_sheet() == "Updated", "ID" %in% names(uploaded_data()),
    grepl("Scoring traits \\(1-5\\): 1", output$upload_message),
    grepl("Quantitative traits: 0", output$upload_message))
  unlink(next_path)
})
unlink(path)
cat("Workbook import, worksheet selection, ID/metadata preservation and upload audit checks passed.\n")
