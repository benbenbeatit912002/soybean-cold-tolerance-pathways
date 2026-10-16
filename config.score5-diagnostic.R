# Diagnostic sensitivity run, not a final inferential analysis.
source("config.R")
config$score_column <- 5L
config$score_kind <- "score"
config$duplicate_policy <- "mean"
config$permutations <- 199L
config$output_dir <- "results/score5_diagnostic"
