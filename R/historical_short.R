# Reconstruct the original calculation. Saved scores are never inputs here.
historical_source_files <- c(
  expression = "\u7814\u7a76\u6cd5/\u5927\u8c46\u57fa\u56e0/\u6578\u64da/\u5916\u90e8\u8cc7\u6599\u532f\u5165/\u8010\u5bd2\u8cc7\u6599\u5eab0 1 24hr(\u8f49\u63db\u5f8c)\u57f9\u4fee\u7248.txt",
  correspondence = "\u7814\u7a76\u6cd5/\u5927\u8c46\u57fa\u56e0/\u5927\u8c46\u57fa\u56e0 \u5916\u90e8/Glyma_genes_correspondence_final.txt",
  gene_length = "\u7814\u7a76\u6cd5/\u5927\u8c46\u57fa\u56e0/\u6578\u64da/\u5916\u90e8\u57fa\u56e0\u6821\u6b63/Glyma2.0.txt",
  snp_markers = "\u7814\u7a76\u6cd5/\u5927\u8c46\u57fa\u56e0/\u6578\u64da/\u5916\u90e8\u57fa\u56e0\u6821\u6b63/snps_allgenes_20kb.set.txt"
)

resolve_historical_sources <- function(source_root) {
  # The portable manifest takes precedence; never silently use another copy.
  manifest <- file.path(source_root, "inputs.tsv")
  if (file.exists(manifest)) {
    entries <- read.delim(manifest, colClasses = "character", check.names = FALSE)
    if (nrow(entries) != length(historical_source_files) || anyDuplicated(entries$source) ||
        !setequal(entries$source, names(historical_source_files))) {
      stop("inputs.tsv must contain each of the four upstream sources once")
    }
    paths <- setNames(file.path(source_root, entries$file), entries$source)
    if (any(!file.exists(paths))) stop("Missing packaged upstream input")
    # Base R supplies MD5; the companion Python check also verifies SHA-256.
    if (any(unname(tools::md5sum(paths)) != entries$md5)) stop("Upstream checksum mismatch")
  } else {
    paths <- setNames(file.path(source_root, historical_source_files), names(historical_source_files))
    if (any(!file.exists(paths))) stop("Missing upstream inputs; use the included data/upstream directory")
  }
  paths
}

read_historical_table <- function(path) {
  read.delim(path, header = TRUE, check.names = FALSE, colClasses = "character",
             na.strings = NULL, comment.char = "", fileEncoding = "UTF-8")
}

count_historical_snp_markers <- function(path) {
  lines <- trimws(readLines(path, encoding = "UTF-8", warn = FALSE))
  lines <- lines[nzchar(lines)]
  starts <- which(substr(lines, 1L, 1L) == "G")
  ends <- which(lines == "END")
  if (length(starts) != length(ends) || any(ends <= starts) ||
      anyDuplicated(lines[starts]) || starts[1L] != 1L || tail(ends, 1L) != length(lines) ||
      (length(starts) > 1L && any(starts[-1L] != head(ends, -1L)+1L))) {
    stop("Malformed or duplicate SNP gene blocks")
  }
  marker_rows <- setdiff(seq_along(lines), c(starts, ends))
  if (any(substr(lines[marker_rows], 1L, 1L) != "A")) stop("Unexpected SNP marker line")
  setNames(ends-starts-1L, lines[starts])
}

reconstruct_historical_input <- function(source_root) {
  paths <- resolve_historical_sources(source_root)
  expression <- read_historical_table(paths[["expression"]])
  correspondence <- read_historical_table(paths[["correspondence"]])
  gene_length <- read_historical_table(paths[["gene_length"]])
  if (anyDuplicated(gene_length[[1L]])) stop("Duplicate gene-length key")
  lengths <- setNames(as.integer(gene_length[[5L]]), gene_length[[1L]])
  counts <- count_historical_snp_markers(paths[["snp_markers"]])
  # Deduplicate correspondence row numbers, not target genes: this preserves
  # the old rule that two matching correspondence rows are ambiguous.
  aliases <- split(rep(seq_len(nrow(correspondence)), 2L),
                   c(correspondence[[1L]], correspondence[[2L]]))
  target <- character(nrow(expression))
  reason <- character(nrow(expression))
  for (index in seq_len(nrow(expression))) {
    matches <- unique(aliases[[expression[[3L]][index]]])
    if (!length(matches)) reason[index] <- "no_correspondence"
    else if (length(matches) != 1L) reason[index] <- "ambiguous_correspondence_rows"
    else {
      target[index] <- correspondence[[3L]][matches]
      if (!nzchar(target[index])) reason[index] <- "blank_target"
    }
  }
  included <- which(reason == "")
  genes <- target[included]
  if (anyNA(lengths[genes]) || anyNA(counts[genes])) stop("Unmapped length or SNP gene")
  result <- data.frame(source_row = seq_along(included)+1L,
    expression_row = included+1L, original_gene_id = expression[[3L]][included],
    gene_id = sub("^Glyma\\.([0-9]{2})[gG]", "Glyma.\\1G", genes), raw_short = expression[[5L]][included],
    raw_middle = expression[[8L]][included],
    gene_length = as.integer(lengths[genes]), snp_markers = as.integer(counts[genes]))
  excluded <- which(reason != "")
  exclusions <- data.frame(expression_row = excluded+1L,
    original_gene_id = expression[[3L]][excluded], reason = reason[excluded])
  list(rows = result, exclusions = exclusions, paths = paths, expression_records = nrow(expression))
}

historical_snp_imputation <- function(rows) {
  values <- as.numeric(rows$snp_markers)
  audit <- list()
  # seq(1,100000,by=2000) ends at 98001. The archived loop's out-of-range
  # upper boundaries select no rows, so only these 49 bins change anything.
  for (start in seq(1, 96001, by = 2000)) {
    members <- which(rows$gene_length >= start & rows$gene_length < start+2000)
    zeros <- members[values[members] == 0]
    if (!length(zeros)) next
    group_mean <- mean(values[members])  # Includes zeros, as in the old code.
    replacement <- round(group_mean)
    audit[[length(audit)+1L]] <- data.frame(source_row = rows$source_row[zeros],
      gene_id = rows$gene_id[zeros], bin_start = start, bin_end_exclusive = start+2000,
      group_records = length(members), mean_including_zeros = group_mean,
      snp_before = 0, snp_after = replacement)
    values[zeros] <- replacement
  }
  audit <- if (length(audit)) do.call(rbind, audit) else data.frame(
    source_row = integer(), gene_id = character(), bin_start = numeric(),
    bin_end_exclusive = numeric(), group_records = integer(),
    mean_including_zeros = numeric(), snp_before = numeric(), snp_after = numeric())
  list(values = values, audit = audit)
}

historical_fosco <- function(probabilities, snp, exponent) {
  # Deliberately retain ordinary subtraction; log1p/expm1 changes old scores.
  1-(1-as.numeric(probabilities))^(as.numeric(snp)^exponent)
}

fit_historical_short <- function(probabilities, snp) {
  probabilities <- as.numeric(probabilities)
  snp <- as.numeric(snp)
  if (length(probabilities) != length(snp) || length(snp) < 2L ||
      any(!is.finite(probabilities) | probabilities <= 0 | probabilities > 1) ||
      any(!is.finite(snp) | snp < 0 | snp != round(snp))) stop("Invalid probability/SNP input")
  grid <- data.frame(lambda = (1:100)/100)
  grid$absolute_spearman <- vapply(grid$lambda, function(exponent) {
    abs(cor(snp, historical_fosco(probabilities, snp, exponent), method = "spearman"))
  }, 0.0)
  if (any(!is.finite(grid$absolute_spearman))) stop("Undefined Spearman correlation")
  best <- grid$lambda[which.min(grid$absolute_spearman)]
  adjusted <- historical_fosco(probabilities, snp, best)
  list(lambda = best, adjusted = adjusted, score = -log10(adjusted), grid = grid)
}

historical_score_export <- function(adjusted) {
  # These export rules are inferred from saved artifacts. No reference is read.
  probability <- as.numeric(sprintf("%.15g", adjusted))
  score <- -log10(probability)
  score[score == 0] <- 1e-16
  display <- vapply(score, function(value) {
    if (!is.finite(value)) return("Inf")
    if (value < 1e-4) return(sprintf("%.2E", value))
    text <- sprintf(if (value >= 10) "%.8f" else "%.9f", value)
    sub("\\.$", "", sub("0+$", "", text))
  }, "")
  list(probability = probability, score = score, display = display)
}

repair_zero_hypothesis <- function(rows, counts, adjusted, source_root) {
  # Explicitly unconfirmed: same whole-genome bin's positive mean, the last
  # phase's fitted lambda, then 7 significant digits (R console precision).
  last_fit <- fit_historical_short(rows$raw_middle, counts)
  paths <- resolve_historical_sources(source_root)
  lengths <- read_historical_table(paths[["gene_length"]])
  genome_counts <- count_historical_snp_markers(paths[["snp_markers"]])
  gene_lengths <- as.integer(lengths[[5L]])
  markers <- genome_counts[lengths[[1L]]]
  audit <- list()
  for (index in which(counts == 0)) {
    start <- 1+2000*((rows$gene_length[index]-1) %/% 2000)
    donors <- which(gene_lengths >= start & gene_lengths < start+2000 & !is.na(markers) & markers > 0)
    if (!length(donors)) next
    donor_mean <- mean(markers[donors])
    replacement <- round(donor_mean)
    probability <- historical_fosco(rows$raw_short[index], replacement, last_fit$lambda)
    rounded <- as.numeric(sprintf("%.7g", probability))
    counts[index] <- replacement
    adjusted[index] <- rounded
    audit[[length(audit)+1L]] <- data.frame(source_row = rows$source_row[index], gene_id = rows$gene_id[index],
      bin_start = start, bin_end_exclusive = start+2000,
      donor_genes = paste(lengths[[1L]][donors], collapse = ";"),
      donor_counts = paste(markers[donors], collapse = ";"), donor_mean = donor_mean,
      snp_after = replacement, lambda_from_middle_grid = last_fit$lambda,
      adjusted_before_rounding = probability, adjusted_7_significant_digits = rounded,
      evidence_status = "unconfirmed_reconstruction_hypothesis")
  }
  list(counts = counts, adjusted = adjusted, audit = do.call(rbind, audit), middle_grid = last_fit$grid)
}
