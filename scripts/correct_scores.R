# Explicit correction entry point. Input: gene_id, raw_short, raw_middle, snp.
# Include gene_length (bp) when SNP counts contain zero. Use the same-group
# nonzero mean, or the global nonzero mean when the group has no donors.
# Raw inputs must be probabilities/pseudo-probabilities in (0,1], not logged scores.
# This command always recalculates scores. Use reproduce_historical_short.R
# for the archived arithmetic and input reconstruction.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) %% 2L != 0L) {
  stop("Usage: Rscript scripts/correct_scores.R input.tsv output_dir [--duplicate-policy error|preserve]")
}
duplicate_policy <- "error"
seen <- character()
if (length(args) > 2L) for (i in seq(3L, length(args), by = 2L)) {
  option <- args[i]
  if (option %in% seen) stop("Repeated option: ", option)
  seen <- c(seen, option)
  if (option == "--duplicate-policy" && args[i+1L] %in% c("error", "preserve")) {
    duplicate_policy <- args[i+1L]
  } else stop("Invalid option: ", option)
}
source("R/core.R")
d <- read_tsv(args[1])
required <- c("gene_id", "raw_short", "raw_middle", "snp")
if (!nrow(d) || !all(required %in% names(d))) stop("Required columns: ", paste(required, collapse = ", "))
d$gene_id <- canonical_id(d$gene_id)
if (any(is.na(d$gene_id) | d$gene_id == "")) stop("Missing gene IDs")
if (anyDuplicated(d$gene_id) && duplicate_policy == "error") {
  stop("Resolve duplicate gene IDs or explicitly use --duplicate-policy preserve")
}
d$snp_original <- d$snp
imputation <- impute_snp_counts(d$snp, d$gene_length)
d$snp <- imputation$values
d$snp_imputation_source <- "unchanged"
imputation$audit$gene_id <- d$gene_id[imputation$audit$row_number]
d$snp_imputation_source[imputation$audit$row_number] <- imputation$audit$source
out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_tsv(imputation$audit[, c("row_number", "gene_id", "bin_start", "bin_end_exclusive",
  "original_snp", "imputed_snp", "donor_count", "donor_mean", "source")],
  file.path(out, "snp_imputation.tsv"))
settings <- c(input = args[1], rows = nrow(d), unique_genes = length(unique(d$gene_id)),
  duplicate_policy = duplicate_policy, bin_width = 2000, bin_origin = 1,
  zero_policy = "group_nonzero_mean_then_global_nonzero_mean",
  donors = "positive_input_rows_before_any_imputation",
  rounding = "nearest_integer_ties_to_even", imputed_rows = nrow(imputation$audit),
  short_score_source = "recomputed_fosco")
write_tsv(data.frame(setting = names(settings), value = unname(settings)),
  file.path(out, "correction_settings.tsv"))
summary <- list()
for (phase in c("short", "middle")) {
  raw <- d[[paste0("raw_", phase)]]
  fit <- fit_fosco(raw, d$snp)
  d[[paste0("adjusted_", phase)]] <- fit$adjusted
  d[[paste0("score_", phase)]] <- -log10(fit$adjusted)
  score_source <- "recomputed_fosco"
  if (phase == "short") {
    d$short_score_source <- score_source
  }
  write_tsv(fit$grid, file.path(out, paste0(phase, "_lambda_grid.tsv")))
  summary[[phase]] <- data.frame(phase = phase, lambda = fit$lambda,
    score_source = score_source,
    raw_spearman = cor(d$snp, raw, method = "spearman"),
    adjusted_spearman = cor(d$snp, as.numeric(d[[paste0("adjusted_", phase)]]), method = "spearman"),
    recomputed_adjusted_spearman = cor(d$snp, fit$adjusted, method = "spearman"))
  png(file.path(out, paste0(phase, "_correction_diagnostics.png")), width = 1500, height = 1000, res = 150)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  point_size <- if (nrow(d) <= 1000) 0.8 else 0.35
  hist(-log10(raw), breaks = 50, col = "#337A8C", border = NA,
       main = paste(phase, "before correction"), xlab = "-log10(raw score)")
  hist(-log10(fit$adjusted), breaks = 50, col = "#CB7551", border = NA,
       main = paste(phase, "after correction"), xlab = "-log10(adjusted score)")
  plot(d$snp, -log10(raw), pch = 16, cex = point_size, col = "#337A8C90",
       xlab = "SNP count", ylab = "-log10(raw score)", main = "Raw score vs SNP count")
  plot(d$snp, -log10(fit$adjusted), pch = 16, cex = point_size, col = "#CB755190",
       xlab = "SNP count", ylab = "-log10(adjusted score)", main = paste("Fitted lambda =", fit$lambda))
  dev.off()
}
write_tsv(d, file.path(out, "corrected_scores.tsv"))
write_tsv(do.call(rbind, summary), file.path(out, "correction_summary.tsv"))
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
