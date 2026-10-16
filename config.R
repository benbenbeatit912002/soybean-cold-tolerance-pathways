# Reanalysis settings, not an assertion of exact paper reproduction.
config <- list(
  input_dir = "data/input",
  phases = c("short", "middle"),
  min_size = 6L, max_size = 1500L,
  permutations = 10000L, seed = 20261015L,
  # Fill these only after confirming the headerless external file definitions.
  score_column = NA_integer_,       # Original external column number, e.g. 5 or 6
  score_kind = NA_character_,       # "score", "pvalue", or "fold_change"
  duplicate_policy = "error",     # "mean" must be chosen explicitly
  methods = c("sumstat", "gsea"),
  # MaxMean remains opt-in and experimental.
  output_dir = "results/external_membership"
)
