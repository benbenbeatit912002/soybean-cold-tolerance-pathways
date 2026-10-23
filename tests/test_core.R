source("R/core.R")
check_error <- function(expr) {
  failed <- inherits(tryCatch(force(expr), error = identity), "error")
  stopifnot(failed)
}
stopifnot(identical(canonical_id(" Glyma.01g000100 "), "Glyma.01G000100"))
check_error(transform_score(c(0, 0.5), "pvalue"))
check_error(transform_score(c(1, Inf), "score"))
stopifnot(all.equal(transform_score(c(0.1, 1), "pvalue"), c(1, 0)))
stopifnot(all.equal(transform_score(c(2, 0.5), "fold_change"), rep(log(2), 2)))

# Direct combinatorial enumeration independently checks the upper-tail calculation.
universe <- LETTERS[1:10]
fg <- LETTERS[1:3]
sets <- data.frame(go_id = c("GO:1", "GO:2"), description = c("one", "two"),
                   genes = I(list(c("A", "B", "D", "X"), c("E", "F"))))
prepared <- prepare_membership(sets, universe, fg, 2, 4)
draws <- combn(universe, length(fg), simplify = FALSE)
exact <- mean(vapply(draws, function(g) sum(g %in% c("A", "B", "D")) >= 2, logical(1)))
stopifnot(all.equal(hypergeometric(prepared)[1], exact), hypergeometric(prepared)[2] == 1)
stopifnot(prepared$annotation_size[1] == 4, prepared$tested_size[1] == 3)
stopifnot(monte_carlo_p(5, 1:4) == 1/5, monte_carlo_p(1, rep(1, 4)) == 1)

# Stable correction agrees with algebra away from cancellation, retains tiny values.
stopifnot(all.equal(fosco(c(0.1, 0.5), c(4, 9), 0.5), 1-(1-c(0.1, 0.5))^c(2, 3)))
stopifnot(fosco(1e-30, 2, 0.5) > 0, fosco(1, 2, 0.5) == 1)
check_error(fosco(0.1, 0, 0.5))
check_error(fit_fosco(c(.1, .2), c(2, 2)))

# Independent expected means: group (2+5)/2 -> 4; global (2+5+20)/3 -> 9.
counts <- c(0, 2, 5, 0, 20, 0)
lengths <- c(1, 2, 2000, 2001, 4001, 120001)
imputed <- impute_snp_counts(counts, lengths)
stopifnot(identical(imputed$values, c(4, 2, 5, 9, 20, 9)),
          identical(imputed$audit$donor_count, c(2L, 3L, 3L)),
          identical(imputed$audit$source,
                    c("group_nonzero_mean", "global_nonzero_mean", "global_nonzero_mean")))
stopifnot(identical(rev(impute_snp_counts(rev(counts), rev(lengths))$values), imputed$values))
stopifnot(impute_snp_counts(c(0, 2, 3), c(1, 2, 3))$values[1] == 2)
stopifnot(identical(impute_snp_counts(c(1, 3))$values, c(1, 3)))
check_error(impute_snp_counts(c(0, 0), c(1, 2)))
check_error(impute_snp_counts(c(0, 1)))
check_error(impute_snp_counts(c(1, -1), c(1, 2)))
check_error(impute_snp_counts(c(1, NA), c(1, 2)))
check_error(impute_snp_counts(c(1, 2.5), c(1, 2)))
check_error(impute_snp_counts(c(0, 1), c(0, 2)))
check_error(impute_snp_counts(c(0, 1), 1))
check_error(impute_snp_counts(numeric(), numeric()))

# GSEA optimized statistic checked against a separate full running sum.
s <- c(5, 4, 3, 2, 1, 0)
for (hits in combn(seq_along(s), 2, simplify = FALSE)) {
  hit <- seq_along(s) %in% hits
  increments <- ifelse(hit, abs(s)/sum(abs(s[hit])), -1/sum(!hit))
  stopifnot(all.equal(gsea_es(s, hits), max(0, cumsum(increments))))
}
stopifnot(gsea_es(rep(0, 6), c(1, 2)) == 1)
check_error(gsea_es(s, seq_along(s)))
stopifnot(maxmean_stat(c(-2, 1, 0, 3)) == 1)
check_error(resolve_duplicate_scores(c("a", "a"), c(1, 3)))
stopifnot(unname(resolve_duplicate_scores(c("a", "a"), c(1, 3), "mean")) == 2)
scores <- setNames(1:10, universe)
for (method in c("sumstat", "gsea", "maxmean")) {
  x <- permutation_test(scores, prepared$genes, method, 39, 17)
  y <- permutation_test(scores, prepared$genes, method, 39, 17)
  stopifnot(identical(x, y), all(x$p_value >= 1/40), all(x$p_bh <= 1))
}
check_error(permutation_test(scores, prepared$genes, "sumstat", 0))
cat("PASS: exact hypergeometric, score domains, stable FOSCO, GSEA reference, ties, duplicates, seeded permutations.\n")
