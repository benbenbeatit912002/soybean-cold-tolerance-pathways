# Preserve stored short-term scores exactly, with source-row identity checks.
source("R/core.R")

input <- data.frame(
  source_row = as.character(2:21),
  gene_id = sprintf("Glyma.01G%06d", pmax(0:19, 1)),
  raw_short = as.character(.02 + (0:19) * .035),
  raw_middle = as.character(.8 - (0:19) * .027),
  gene_length = as.character(c(100:118, 52247)),
  snp = as.character(c((0:18) %% 5 + 1, 0)),
  stringsAsFactors = FALSE)
reference <- input
reference$gene_id <- sub("Glyma.01G", "Glyma.01g", reference$gene_id, fixed = TRUE)
reference$reference_adjusted_short <- "0.9000000000000000"
reference$reference_score_short <- "0.045757490560675125409944193489"
reference$reference_snp_short <- "12"
reference$reference_log_row <- as.character(100:119)

expect_error <- function(expression) {
  failed <- tryCatch({ force(expression); FALSE }, error = function(e) TRUE)
  stopifnot(failed)
}

before <- input
matched <- match_short_reference(input, reference[20:1, , drop = FALSE])
stopifnot(nrow(matched) == 20,
          identical(as.character(matched$source_row), input$source_row),
          identical(matched$reference_log_row, reference$reference_log_row),
          identical(matched$reference_score_short, reference$reference_score_short),
          input$gene_id[1] == input$gene_id[2],
          matched$reference_log_row[1] != matched$reference_log_row[2],
          identical(input, before))
expect_error(match_short_reference(input, reference[-20, , drop = FALSE]))
extra <- reference[20, , drop = FALSE]
extra$source_row <- "999"
expect_error(match_short_reference(input, rbind(reference, extra)))

bad_identity <- list(source_row = c("999", "2", "21.5", "0", "NaN"),
                     gene_id = "Glyma.99G999999", raw_short = "0.999",
                     gene_length = "52248", snp = "1")
for (field in names(bad_identity)) for (value in bad_identity[[field]]) {
  changed <- reference
  changed[[field]][20] <- value
  expect_error(match_short_reference(input, changed))
}
for (value in c("2", "21.5", "0", "Inf")) {
  changed <- input
  changed$source_row[20] <- value
  expect_error(match_short_reference(changed, reference))
}
expect_error(match_short_reference(input[, names(input) != "source_row", drop = FALSE], reference))
expect_error(match_short_reference(input, reference[, names(reference) != "reference_score_short", drop = FALSE]))

invalid <- list(
  reference_adjusted_short = c("0", "-0.1", "1.1", "NaN", "Inf", ""),
  reference_score_short = c("-0.1", "NaN", "Inf", ""),
  reference_snp_short = c("0", "-1", "2.5", "NaN", "Inf", ""),
  reference_log_row = c("0", "-1", "2.5", "NaN", "Inf", ""))
for (field in names(invalid)) for (value in invalid[[field]]) {
  changed <- reference
  changed[[field]][1] <- value
  expect_error(match_short_reference(input, changed))
}

cat("Reference identity matching tests passed (verification only).\n")
