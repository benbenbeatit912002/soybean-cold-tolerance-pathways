# Exercise score-method and correction entry points on known synthetic inputs.
source("R/core.R")
source("config.R")
scratch <- "tests/tmp"
dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
ids <- sprintf("Glyma.01G%06d", 1:20)
input <- file.path(scratch, "input")
dir.create(input, recursive = TRUE, showWarnings = FALSE)
for (phase in c("short", "middle")) {
  write_tsv(data.frame(gene_id = ids[1:4], prior_score = 8:11), file.path(input, paste0("ct_", phase, ".tsv")))
  write_tsv(data.frame(gene_id = ids, source_col_2 = 1:20, source_col_3 = 1:20,
    source_col_4 = 1:20, source_col_5 = 1:20, source_col_6 = 1:20),
    file.path(input, paste0("external_", phase, ".tsv")))
}
writeLines(c(paste(c("GO:1", "first", ids[1:6]), collapse = "\t"),
             paste(c("GO:2", "second", ids[15:20]), collapse = "\t")), file.path(input, "pathways.gmt"))
config$input_dir <- input
config$output_dir <- file.path(scratch, "scores")
config$score_column <- 6L
config$score_kind <- "score"
config$permutations <- 39L
config$methods <- c("sumstat", "gsea", "maxmean")
config_file <- file.path(scratch, "config.R")
writeLines(c("config <-", capture.output(dput(config))), config_file)
rscript <- file.path(R.home("bin"), "Rscript.exe")
if (.Platform$OS.type != "windows") rscript <- file.path(R.home("bin"), "Rscript")
status <- system2(rscript, c("--vanilla", "scripts/run_analysis.R", config_file))
stopifnot(status == 0)
for (phase in c("short", "middle")) for (method in config$methods) {
  d <- read_tsv(file.path(config$output_dir, paste0(phase, "_", method, ".tsv")))
  stopifnot(nrow(d) == 2, all(d$p_value >= 1/40), all(is.finite(d$statistic)))
}
correction <- data.frame(gene_id = ids, raw_short = seq(.01, .9, length.out = 20),
  raw_middle = rev(seq(.02, .8, length.out = 20)), snp = rep(1:5, 4),
  gene_length = c(100:118, 52247))
correction$snp[c(1, 20)] <- 0
write_tsv(correction, file.path(scratch, "raw.tsv"))
status <- system2(rscript, c("--vanilla", "scripts/correct_scores.R",
                            file.path(scratch, "raw.tsv"), file.path(scratch, "correction")))
stopifnot(status == 0)
corrected <- read_tsv(file.path(scratch, "correction/corrected_scores.tsv"))
stopifnot(nrow(corrected) == 20, all(corrected$adjusted_short > 0),
          all(corrected$adjusted_middle <= 1), all(is.finite(corrected$score_short)))
audit <- read_tsv(file.path(scratch, "correction/snp_imputation.tsv"))
stopifnot(nrow(audit) == 2, all(audit$imputed_snp == 3),
          identical(audit$source, c("group_nonzero_mean", "global_nonzero_mean")),
          all(corrected$snp_original[c(1, 20)] == 0))
cat("PASS: both-phase score workflow and separate correction workflow on synthetic data.\n")
# Named recalculated scores must override the separate external input file.
fresh <- file.path(scratch,"recalculated.tsv")
write_tsv(data.frame(gene_id=ids,score_short=2*(1:20),score_middle=3*(1:20)),fresh)
config$score_input <- fresh; config$score_field <- "score_{phase}"
config$score_column <- NA_integer_; config$methods <- "sumstat"
config$output_dir <- file.path(scratch,"fresh")
writeLines(c("config <-",capture.output(dput(config))),config_file)
status <- system2(rscript,c("--vanilla","scripts/run_analysis.R",config_file))
stopifnot(status == 0)
for (phase in c("short","middle")) {
  result <- read_tsv(file.path(config$output_dir,paste0(phase,"_sumstat.tsv")))
  multiplier <- if (phase == "short") 2 else 3
  stopifnot(result$statistic[1] == multiplier*sum(1:6))
}
cat("PASS: both phases consume named recalculated scores.\n")
