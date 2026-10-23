# Validated building blocks. Base R only. See docs/METHODS.md for assumptions.
canonical_id <- function(x) {
  sub("^Glyma\\.([0-9]{2})[gG]", "Glyma.\\1G", trimws(as.character(x)))
}

read_tsv <- function(path) {
  read.delim(path, header = TRUE, sep = "\t", quote = "\"", comment.char = "",
             check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE, quote = TRUE, na = "NA",
              fileEncoding = "UTF-8")
}

read_gmt <- function(path) {
  lines <- strsplit(readLines(path, encoding = "UTF-8"), "\t", fixed = TRUE)
  if (!length(lines) || any(lengths(lines) < 3L)) stop("Invalid or empty GMT file")
  ids <- vapply(lines, `[`, "", 1)
  if (anyDuplicated(ids)) stop("Duplicate pathway IDs")
  data.frame(go_id = ids, description = vapply(lines, `[`, "", 2),
             genes = I(lapply(lines, function(x) unique(canonical_id(x[-c(1, 2)])))))
}

validate_values <- function(x, kind) {
  if (!is.numeric(x) || !length(x) || any(!is.finite(x))) stop("Scores must be finite numeric values")
  if (kind == "pvalue" && any(x <= 0 | x > 1)) stop("pvalue requires 0 < value <= 1; zeros need an explicit upstream policy")
  if (kind == "fold_change" && any(x <= 0)) stop("fold_change requires positive ratios")
  if (kind == "score" && any(x < 0)) stop("score requires nonnegative, larger-is-stronger values")
  if (!kind %in% c("pvalue", "fold_change", "score")) stop("Unknown score kind")
  invisible(TRUE)
}

transform_score <- function(x, kind) {
  validate_values(x, kind)
  switch(kind, pvalue = -log10(x), fold_change = abs(log(x)), score = x)
}

# Match saved short scores by source row; duplicate gene IDs remain distinct.
# Read the reference as character columns to preserve its printed score precision.
match_short_reference <- function(d, reference) {
  identity <- c("source_row", "gene_id", "raw_short", "gene_length", "snp")
  saved <- c("reference_adjusted_short", "reference_score_short",
             "reference_snp_short", "reference_log_row")
  if (!nrow(d) || !nrow(reference) || !all(identity %in% names(d)) ||
      !all(c(identity, saved) %in% names(reference))) {
    stop("Short reference and input require source-row identity and saved score columns")
  }
  positive_integer <- function(x) {
    x <- as.numeric(x)
    if (any(!is.finite(x) | x < 1 | x != floor(x))) {
      stop("Reference row numbers and SNP counts must be positive integers")
    }
    x
  }
  keys <- positive_integer(d$source_row)
  reference_keys <- positive_integer(reference$source_row)
  if (anyDuplicated(keys) || anyDuplicated(reference_keys) || !setequal(keys, reference_keys)) {
    stop("Short reference source rows must match input exactly, without duplicates")
  }
  matched <- reference[match(keys, reference_keys), , drop = FALSE]
  genes <- canonical_id(d$gene_id)
  old_genes <- canonical_id(matched$gene_id)
  if (anyNA(genes) || anyNA(old_genes) || any(!nzchar(genes)) || any(genes != old_genes)) {
    stop("Short reference gene identity differs from input")
  }
  for (field in c("raw_short", "gene_length", "snp")) {
    value <- as.numeric(d[[field]])
    original <- as.numeric(matched[[field]])
    if (any(!is.finite(value)) || any(!is.finite(original)) || any(value != original)) {
      stop("Short reference ", field, " differs from input")
    }
  }
  positive_integer(matched$reference_snp_short)
  positive_integer(matched$reference_log_row)
  validate_values(as.numeric(matched$reference_adjusted_short), "pvalue")
  validate_values(as.numeric(matched$reference_score_short), "score")
  matched
}

# Fill zeros using positive input counts in the same 2,000-bp group, else globally.
# Donors are frozen before filling; newly imputed values never become donors.
impute_snp_counts <- function(snp, gene_length = NULL, bin_width = 2000) {
  if (!is.numeric(snp) || !length(snp) ||
      any(!is.finite(snp) | snp < 0 | snp != floor(snp))) {
    stop("SNP counts must be finite nonnegative integers")
  }
  if (length(bin_width) != 1 || !is.finite(bin_width) ||
      bin_width < 1 || bin_width != floor(bin_width)) stop("Invalid bin width")
  positive <- snp > 0
  if (!any(positive)) stop("Cannot impute SNP counts: no positive input donors")
  audit <- data.frame(row_number = integer(), bin_start = numeric(),
    bin_end_exclusive = numeric(), original_snp = numeric(), imputed_snp = numeric(),
    donor_count = integer(), donor_mean = numeric(), source = character())
  if (is.null(gene_length)) {
    if (any(!positive)) stop("gene_length is required when SNP counts contain zero")
    return(list(values = snp, audit = audit))
  }
  if (!is.numeric(gene_length) || length(gene_length) != length(snp) ||
      any(!is.finite(gene_length) | gene_length < 1 | gene_length != floor(gene_length))) {
    stop("Gene lengths must be finite positive integers matching SNP counts")
  }
  bins <- floor((gene_length - 1) / bin_width)
  filled <- snp
  for (group in sort(unique(bins[!positive]))) {
    donors <- snp[positive & bins == group]
    method <- "group_nonzero_mean"
    if (!length(donors)) {
      donors <- snp[positive]
      method <- "global_nonzero_mean"
    }
    donor_mean <- mean(donors)
    replacement <- round(donor_mean)  # Nearest integer; ties to even.
    indices <- which(bins == group & !positive)
    filled[indices] <- replacement
    audit <- rbind(audit, data.frame(row_number = indices,
      bin_start = group * bin_width + 1, bin_end_exclusive = (group + 1) * bin_width + 1,
      original_snp = 0, imputed_snp = replacement, donor_count = length(donors),
      donor_mean = donor_mean, source = method))
  }
  list(values = filled, audit = audit[order(audit$row_number), , drop = FALSE])
}

# Stable evaluation of 1 - (1 - p)^(M^lambda).
# A fold-change-derived pseudo-probability is a score, not a calibrated p-value.
fosco <- function(p, snp, lambda) {
  validate_values(p, "pvalue")
  if (length(p) != length(snp) || any(!is.finite(snp) | snp <= 0)) stop("FOSCO needs a positive SNP count for every score")
  if (length(lambda) != 1L || !is.finite(lambda) || lambda < 0) stop("Invalid lambda")
  -expm1((snp^lambda) * log1p(-p))
}

fit_fosco <- function(p, snp, grid = seq(0.01, 1, 0.01)) {
  if (!length(grid) || any(!is.finite(grid) | grid < 0)) stop("Invalid lambda grid")
  objective <- vapply(grid, function(lam) {
    adjusted <- fosco(p, snp, lam)
    suppressWarnings(abs(cor(snp, adjusted, method = "spearman")))
  }, numeric(1))
  if (!any(is.finite(objective))) stop("Cannot fit FOSCO: constant scores or SNP counts")
  objective[!is.finite(objective)] <- Inf
  index <- which.min(objective)
  list(lambda = grid[index], adjusted = fosco(p, snp, grid[index]),
       grid = data.frame(lambda = grid, absolute_spearman = objective))
}

# Preserve the foreground membership independently of its prior ranking score.
prepare_membership <- function(pathways, universe, foreground, min_size = 6L, max_size = 1500L) {
  universe <- unique(canonical_id(universe))
  foreground <- intersect(unique(canonical_id(foreground)), universe)
  members <- lapply(pathways$genes, intersect, y = universe)
  n <- lengths(members)
  keep <- n >= min_size & n <= max_size
  result <- pathways[keep, , drop = FALSE]
  result$annotation_size <- lengths(pathways$genes)[keep]
  result$genes <- I(members[keep])
  result$tested_size <- n[keep]
  result$ct_count <- vapply(result$genes, function(g) length(intersect(g, foreground)), integer(1))
  result$ct_ratio <- result$ct_count / result$tested_size
  attr(result, "universe_size") <- length(universe)
  attr(result, "ct_in_universe") <- length(foreground)
  result
}

hypergeometric <- function(pathways) {
  N <- attr(pathways, "universe_size")
  K <- attr(pathways, "ct_in_universe")
  if (is.null(N) || is.null(K)) stop("Use prepare_membership first")
  phyper(pathways$ct_count - 1L, pathways$tested_size, N - pathways$tested_size,
         K, lower.tail = FALSE)
}

monte_carlo_p <- function(observed, null) {
  if (!is.finite(observed) || !length(null) || any(!is.finite(null))) stop("Invalid permutation statistics")
  # Include ties, with a small relative tolerance for floating-point sums.
  tolerance <- 1e-12 * max(1, abs(observed))
  (1 + sum(null >= observed - tolerance)) / (length(null) + 1)
}

maxmean_stat <- function(z) max(mean(pmax(z, 0)), mean(pmax(-z, 0)))

# Positive-tail weighted running sum; scores already sorted decreasingly.
# Evaluate hit positions only: between hits, the running sum can only decrease.
gsea_es <- function(ranked_scores, hits) {
  n <- length(ranked_scores)
  hits <- sort(unique(hits))
  if (!length(hits) || length(hits) >= n || any(hits < 1 | hits > n)) stop("GSEA needs 0 < set size < universe size")
  w <- abs(ranked_scores[hits])
  if (sum(w) == 0) w[] <- 1
  max(0, cumsum(w)/sum(w) - (hits-seq_along(hits))/(n-length(hits)))
}

# Compare each observed pathway with B random gene sets of the same size.
# Every draw contains distinct genes. The next draw starts from the whole universe.
permutation_test <- function(scores, members, method, B = 10000L, seed = 20261015L, verbose = FALSE) {
  if (!is.numeric(B) || length(B) != 1 || !is.finite(B) || B < 1 || B != as.integer(B)) stop("B must be a positive integer")
  if (!method %in% c("sumstat", "gsea", "maxmean")) stop("Unknown method")
  validate_values(unname(scores), "score")
  if (is.null(names(scores)) || anyDuplicated(names(scores))) stop("Need one named score per gene")
  set.seed(seed)
  if (method == "gsea") scores <- scores[order(-scores, names(scores))]
  if (method == "maxmean") {
    if (sd(scores) == 0) stop("MaxMean requires variable scores")
    scores <- (scores-mean(scores))/sd(scores)
  }
  indices <- lapply(members, function(g) match(g, names(scores)))
  if (anyNA(unlist(indices)) || any(lengths(indices) == 0)) stop("Pathway genes missing from score universe")
  statistic <- switch(method,
    sumstat = function(i) sum(scores[i]),
    maxmean = function(i) maxmean_stat(scores[i]),
    gsea = function(i) gsea_es(scores, i))
  observed <- vapply(indices, statistic, numeric(1))
  p <- numeric(length(indices))
  # Reuse null samples for equal-size pathways; each marginal null is unchanged.
  sizes <- sort(unique(lengths(indices)))
  for (size_index in seq_along(sizes)) {
    size <- sizes[size_index]
    # Hash sampling avoids allocating the whole genome for each small gene set.
    null <- replicate(B, statistic(sample.int(length(scores), size, replace = FALSE,
                                               useHash = size <= length(scores)/2)))
    selected <- which(lengths(indices) == size)
    p[selected] <- vapply(observed[selected], monte_carlo_p, numeric(1), null = null)
    if (verbose && (size_index %% 50L == 0L || size_index == length(sizes))) {
      message(method, ": completed ", size_index, "/", length(sizes), " distinct pathway sizes")
    }
  }
  data.frame(statistic = observed, p_value = p, p_bh = p.adjust(p, "BH"),
             null_exceedances = round(p * (B + 1) - 1), permutations = B, seed = seed)
}

resolve_duplicate_scores <- function(ids, values, policy = "error") {
  ids <- canonical_id(ids)
  if (!policy %in% c("error", "mean")) stop("Duplicate policy must be error or mean")
  if (anyDuplicated(ids)) {
    if (policy == "error") stop("Duplicate gene IDs: resolve provenance or explicitly select duplicate_policy=mean")
    values <- tapply(values, ids, mean)
    return(setNames(as.numeric(values), names(values)))
  }
  setNames(values, ids)
}
