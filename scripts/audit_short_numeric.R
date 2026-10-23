# Re-evaluate the archived R expression. Reference scores are not calculation inputs.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: audit_short_numeric.R INPUT.tsv OUTPUT_DIRECTORY")
data <- read.delim(args[1L], check.names = FALSE)
dir.create(args[2L], recursive = TRUE, showWarnings = FALSE)
p <- data$raw_short
snp <- data$snp
grid <- data.frame(lambda = (1:100)/100)
grid$absolute_spearman <- vapply(grid$lambda, function(lambda) {
  abs(cor(snp, 1-(1-p)^(snp^lambda), method = "spearman"))
}, 0.0)
lambda <- grid$lambda[which.min(grid$absolute_spearman)]
probability <- 1-(1-p)^(snp^lambda)
result <- data.frame(source_row = data$source_row,
  adjusted_short = sprintf("%.17g", probability),
  score_short = sprintf("%.17g", -log10(probability)))
write.table(result, file.path(args[2L], "r_naive_short.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(grid, file.path(args[2L], "r_naive_grid.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("Short lambda:", lambda, "\nSNP median:", median(snp),
    "\nPositive SNP mean:", mean(snp[snp > 0]), "\n")
