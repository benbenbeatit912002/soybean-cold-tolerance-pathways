# Parse archival scripts without executing their writes, installations, or analyses.
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
files <- list.files("legacy", full.names = TRUE)
if (!length(files)) stop("No legacy scripts found; check working directory and UTF-8 locale")
results <- lapply(files, function(f) {
  problem <- tryCatch({parse(f, encoding = "UTF-8"); ""}, error = function(e) conditionMessage(e))
  data.frame(file = basename(f), parses = identical(problem, ""), diagnostic = problem)
})
source("R/core.R")
dir.create("reports", showWarnings = FALSE)
write_tsv(do.call(rbind, results), "reports/legacy_parse.tsv")
print(do.call(rbind, results), row.names = FALSE)
