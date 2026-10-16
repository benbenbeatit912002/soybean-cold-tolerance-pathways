# Small deterministic reference: enumerate every 3-gene set, no RNG involved.
source("R/core.R")
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
scores <- c(5, 0, 2, 7, 1, 4, 8, 3, 6, 2)
hits <- c(1L, 4L, 7L)
null_sets <- combn(seq_along(scores), 3, simplify = FALSE)
z <- (scores-mean(scores))/sd(scores)
ranked <- sort(scores, decreasing = TRUE)
functions <- list(sumstat = function(i) sum(scores[i]),
                  maxmean = function(i) maxmean_stat(z[i]),
                  gsea = function(i) gsea_es(ranked, i))
result <- lapply(names(functions), function(method) {
  statistic <- functions[[method]]
  observed <- statistic(hits)
  null <- vapply(null_sets, statistic, numeric(1))
  data.frame(method = method, observed = observed, p_value = monte_carlo_p(observed, null))
})
write_tsv(do.call(rbind, result), file.path(out, "reference.tsv"))
raw <- c(.05, .15, .03, .55, .22, .81)
snp <- c(1, 3, 2, 4, 2, 5)
fit <- fit_fosco(raw, snp)
write_tsv(data.frame(raw = raw, snp = snp, adjusted = fosco(raw, snp, .47),
                     fitted = fit$adjusted, lambda = fit$lambda), file.path(out, "fosco.tsv"))
