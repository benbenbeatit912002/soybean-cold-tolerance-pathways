config <- list(
  input_dir="data/input",
  score_input="results/recalculated_2026-10-21/r_scores/recalculated_scores.tsv",
  score_field="score_{phase}", score_column=NA_integer_, score_kind="score",
  phases=c("short","middle"), min_size=6L, max_size=1500L,
  permutations=10000L, seed=20261021L, duplicate_policy="mean",
  methods=c("sumstat","gsea","maxmean"),
  output_dir="results/recalculated_2026-10-21/r_full",
  score_provenance="Fresh calculation; explicit historical repair hypotheses enabled",
  analysis_protocol="Measured unique-gene universe, mean duplicate scores, BH over all measured-size-eligible pathways"
)
