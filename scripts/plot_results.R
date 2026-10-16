# Plot available method results from one output directory.
# Usage: Rscript scripts/plot_results.R results/r_score5_full
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else "results/external_membership"
source("R/core.R", encoding = "UTF-8")
colors <- c(hypergeometric = "#337A8C", sumstat = "#527A4D", gsea = "#76558C", maxmean = "#B47043")
count <- 0L
for (phase in c("short", "middle")) {
 for (method in names(colors)) {
  input <- file.path(out, paste0(phase, "_", method, ".tsv"))
  if (!file.exists(input)) next
  d <- read_tsv(input)
  d <- d[d$ct_count > 0, , drop = FALSE]
  p_column <- if (method == "hypergeometric") "hypergeometric_p" else "p_value"
  q_column <- if (method == "hypergeometric") "hypergeometric_bh" else "p_bh"
  d <- head(d[order(d[[q_column]], d[[p_column]], d$go_id), ], 15L)
  if (!nrow(d)) next
  png(file.path(out, paste0(phase, "_", method, ".png")), width = 1500, height = 950, res = 150)
  par(mar = c(5, 19, 4, 2))
  description <- ifelse(nchar(d$description) > 42, paste0(substr(d$description, 1, 39), "..."), d$description)
  labels <- paste(d$go_id, description, sep = "  ")
  value <- -log10(pmax(d[[q_column]], .Machine$double.xmin))
  title <- paste(if (phase == "short") "Short-term" else "Mid-term", "|", toupper(method))
  if (method == "maxmean") title <- paste(title, "(experimental)")
  barplot(rev(value), names.arg = rev(labels), horiz = TRUE, las = 1, cex.names = 0.68,
          col = colors[method], border = NA, xlab = "-log10(BH-adjusted p-value)", main = title)
  mtext("Top 15 terms containing CTgenes; selection does not imply significance", side = 1, line = 3.7, cex = 0.7)
  dev.off()
  count <- count + 1L
 }
}
if (!count) stop("No results to plot: run the analysis first and choose its output folder")
message("Saved ", count, " plots to ", out)
