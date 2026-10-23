# Run from the repository root. Each stage below writes a named, inspectable table.
# 1. Load a configuration and GO memberships.
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args)) args[1] else "config.R"
source("R/core.R", encoding = "UTF-8")
source(config_path, encoding = "UTF-8")
out <- config$output_dir
dir.create(out, recursive = TRUE, showWarnings = FALSE)
pathways <- read_gmt(file.path(config$input_dir, "pathways.gmt"))
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
dput(config, file = file.path(out, "settings.R"))
all_qc <- list()
for (phase in config$phases) {
  # 2. Load CT membership and the independent external scores for this phase.
  if (!phase %in% c("short", "middle")) stop("Only short and middle phases are supported")
  ct <- read_tsv(file.path(config$input_dir, paste0("ct_", phase, ".tsv")))
  score_input <- config$score_input
  ext <- read_tsv(if (is.null(score_input)) file.path(config$input_dir,paste0("external_",phase,".tsv")) else score_input)
  for (ids in list(ct$gene_id, ext$gene_id)) {
    if (!length(ids) || anyNA(ids) || any(!nzchar(ids))) stop("Missing gene identifiers")
  }
  if (anyDuplicated(ct$gene_id)) stop("Duplicate CTgene IDs")
  universe <- unique(ext$gene_id)
  missing <- ct[!ct$gene_id %in% universe, , drop = FALSE]
  write_tsv(missing, file.path(out, paste0(phase, "_ct_missing_external.tsv")))
  duplicated_ids <- unique(ext$gene_id[duplicated(ext$gene_id)])
  write_tsv(ext[ext$gene_id %in% duplicated_ids, , drop = FALSE],
            file.path(out, paste0(phase, "_duplicate_rows.tsv")))
  sets <- prepare_membership(pathways, universe, ct$gene_id, config$min_size, config$max_size)
  # 3. Test overlap using counts only; no numerical gene score is required here.
  result <- sets[, setdiff(names(sets), "genes"), drop = FALSE]
  result$hypergeometric_p <- hypergeometric(sets)
  # Include zero-overlap pathways in the testing family; filter for display afterwards.
  result$hypergeometric_bh <- p.adjust(result$hypergeometric_p, "BH")
  write_tsv(result, file.path(out, paste0(phase, "_hypergeometric.tsv")))
  all_qc[[phase]] <- data.frame(phase = phase, ct_input = nrow(ct),
    external_rows = nrow(ext), universe_size = length(universe),
    ct_in_universe = sum(ct$gene_id %in% universe),
    missing_ct = nrow(missing), duplicate_extra_rows = nrow(ext)-length(universe),
    tested_pathways = nrow(sets), pathways_with_ct = sum(sets$ct_count > 0))
  score_field <- if (is.null(config$score_field)) "" else gsub("{phase}",phase,config$score_field,fixed=TRUE)
  if ((!nzchar(score_field) && is.na(config$score_column)) || is.na(config$score_kind)) {
    message(phase, ": membership analysis complete; score analyses pending column definition.")
    next
  }
  column <- config$score_column
  if (!nzchar(score_field) && (length(column) != 1L || !column %in% 2:6)) stop("score_column must name an original external numeric column (2 through 6)")
  # Convert each source row, then aggregate only under an explicit duplicate policy.
  values <- transform_score(ext[[if (nzchar(score_field)) score_field else column]], config$score_kind)
  scores <- resolve_duplicate_scores(ext$gene_id, values, config$duplicate_policy)
  write_tsv(data.frame(gene_id = names(scores), analysis_score = unname(scores)),
            file.path(out, paste0(phase, "_analysis_scores.tsv")))
  for (method in config$methods) {
    # 4. Compare observed gene-set scores with matched-size random gene sets.
    message(phase, ": ", method, " with ", config$permutations, " random gene sets")
    if (method == "maxmean") warning("MaxMean is experimental; no paper-reproduction claim")
    test <- permutation_test(scores, sets$genes, method, config$permutations,
                             config$seed + match(phase, c("short", "middle")) - 1L,
                             verbose = TRUE)
    write_tsv(cbind(result[, 1:6], test), file.path(out, paste0(phase, "_", method, ".tsv")))
  }
}
write_tsv(do.call(rbind, all_qc), file.path(out, "input_summary.tsv"))
message("Outputs: ", normalizePath(out, winslash = "/"))
