# Full 10,000-draw reanalysis of score column 6; MaxMean remains experimental.
source("config.R")
config$score_column <- 6L
config$score_kind <- "score"
config$duplicate_policy <- "mean"
config$methods <- c("sumstat", "gsea", "maxmean")
config$permutations <- 10000L
config$output_dir <- "results/r_score6_full"
