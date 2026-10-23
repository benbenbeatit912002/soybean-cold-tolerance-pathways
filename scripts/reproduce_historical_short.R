# Usage: Rscript scripts/reproduce_historical_short.R ORIGINAL_PROJECT OUTPUT_DIRECTORY
# Recalculate from four upstream files. Reference scores are not read.
# Exit status 2 means the archived procedure left nonfinite scores.
args <- commandArgs(trailingOnly = TRUE)
use_hypothesis <- "--zero-repair-hypothesis" %in% args
args <- args[args != "--zero-repair-hypothesis"]
if (length(args) != 2L) stop("Expected original project root and output directory")
script <- sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())])
source(file.path(dirname(normalizePath(script)), "../R/historical_short.R"), encoding = "UTF-8")
rebuilt <- reconstruct_historical_input(args[1L])
rows <- rebuilt$rows
imputation <- historical_snp_imputation(rows)
rows$historical_snp <- imputation$values
fit <- fit_historical_short(rows$raw_short, imputation$values)
counts <- imputation$values
adjusted <- fit$adjusted
if (use_hypothesis) {
  repair <- repair_zero_hypothesis(rows, counts, adjusted, args[1L])
  counts <- repair$counts
  adjusted <- repair$adjusted
}
export <- historical_score_export(adjusted)
rows$snp <- counts
rows$adjusted_short <- adjusted
rows$adjusted_short_saved_precision <- export$probability
rows$score_short_unrounded <- export$score
rows$score_short <- export$display
rows$calculation_status <- ifelse(is.finite(export$score), "finite", "unresolved_zero_probability")
output <- args[2L]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write_precise <- function(data, filename) {
  for (name in names(data)) {
    if (is.numeric(data[[name]])) data[[name]] <- sprintf("%.17g", data[[name]])
  }
  write.table(data, file.path(output, filename), sep = "\t", quote = FALSE,
              row.names = FALSE, fileEncoding = "UTF-8")
}
write_precise(rows, "recalculated_short.tsv")
write_precise(rebuilt$exclusions, "mapping_exclusions.tsv")
write_precise(imputation$audit, "historical_snp_imputation.tsv")
write_precise(fit$grid, "short_lambda_grid.tsv")
unresolved <- rows[!is.finite(export$score), , drop = FALSE]
write_precise(unresolved, "unresolved_rows.tsv")
if (use_hypothesis) {
  if (!is.null(repair$audit)) write_precise(repair$audit, "zero_repair_hypothesis.tsv")
  write_precise(repair$middle_grid, "middle_lambda_grid_for_hypothesis.tsv")
}
settings <- c(expression_records = rebuilt$expression_records, mapped_records = nrow(rows),
  unique_genes = length(unique(rows$gene_id)), lambda_short = fit$lambda,
  zero_markers_before = sum(rows$snp_markers == 0), zero_snp_after = sum(rows$snp == 0),
  nonfinite_scores = nrow(unresolved), arithmetic = "1-(1-p)^(SNP^lambda); ordinary floating-point subtraction",
  zero_policy = "original_bin_mean_including_zeros_no_fallback", reference_scores_read = "FALSE",
  export_rule = "15-significant-digit probability; zero score to 1e-16; inferred saved decimal format",
  zero_repair_hypothesis = use_hypothesis, repair_provenance_confirmed = "FALSE",
  R = R.version.string)
write_precise(data.frame(setting = names(settings), value = unname(settings)), "recalculation_settings.tsv")
write_precise(data.frame(source = names(rebuilt$paths), path = unname(rebuilt$paths),
  md5 = unname(tools::md5sum(rebuilt$paths))), "source_manifest.tsv")
cat("Recalculated", nrow(rows), "records; lambda=", fit$lambda, "; unresolved=", nrow(unresolved), "\n")
if (use_hypothesis) cat("HYPOTHESIS: matching saved values does not confirm the historical repair procedure.\n")
if (nrow(unresolved)) quit(status = 2L)
