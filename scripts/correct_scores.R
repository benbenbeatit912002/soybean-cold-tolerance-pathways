# Explicit correction entry point. Input: gene_id, raw_short, raw_middle, snp.
# Raw inputs must be probabilities/pseudo-probabilities in (0,1], not logged scores.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript scripts/correct_scores.R input.tsv new_output_directory")
source("R/core.R")
d <- read_tsv(args[1])
required <- c("gene_id", "raw_short", "raw_middle", "snp")
if (!all(required %in% names(d))) stop("Required columns: ", paste(required, collapse = ", "))
d$gene_id <- canonical_id(d$gene_id)
if (anyDuplicated(d$gene_id)) stop("Resolve duplicate gene mappings before correction")
if (any(!is.finite(d$snp) | d$snp <= 0)) stop("Resolve zero/missing SNP counts explicitly; no implicit imputation")
out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
summary <- list()
for (phase in c("short", "middle")) {
  raw <- d[[paste0("raw_", phase)]]
  fit <- fit_fosco(raw, d$snp)
  d[[paste0("adjusted_", phase)]] <- fit$adjusted
  d[[paste0("score_", phase)]] <- -log10(fit$adjusted)
  write_tsv(fit$grid, file.path(out, paste0(phase, "_lambda_grid.tsv")))
  summary[[phase]] <- data.frame(phase = phase, lambda = fit$lambda,
    raw_spearman = cor(d$snp, raw, method = "spearman"),
    adjusted_spearman = cor(d$snp, fit$adjusted, method = "spearman"))
  png(file.path(out, paste0(phase, "_correction_diagnostics.png")), width = 1500, height = 1000, res = 150)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
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
