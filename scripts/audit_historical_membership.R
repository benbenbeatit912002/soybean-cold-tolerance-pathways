# Recount CT membership and replay the archived full-annotation-size test.
# Only gene IDs are used; saved scores and pathway probabilities are not inputs.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
root <- dirname(dirname(normalizePath(script, winslash="/", mustWork=TRUE)))
source(file.path(root, "R/core.R"))
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Usage: Rscript audit_historical_membership.R GENE_TABLE OUTPUT")
out <- args[2]
if (file.exists(out) && (!dir.exists(out) || length(list.files(out, all.files=TRUE, no..=TRUE)))) {
  stop("Output directory must be new or empty")
}
rows <- read_tsv(args[1])
if (!nrow(rows) || !"gene_id" %in% names(rows) || anyNA(rows$gene_id) || any(!nzchar(trimws(rows$gene_id)))) {
  stop("A nonempty gene_id column is required")
}
universe <- unique(canonical_id(rows$gene_id))
pathway_path <- file.path(root, "data/input/pathways.gmt")
pathways <- read_gmt(pathway_path)
dir.create(out, recursive=TRUE, showWarnings=FALSE)
summaries <- list()
missing <- list()
sources <- c(args[1], pathway_path)
for (phase in c("short", "middle")) {
  ct_path <- file.path(root, paste0("data/input/ct_", phase, ".tsv"))
  sources <- c(sources, ct_path)
  primary <- unique(canonical_id(read_tsv(ct_path)$gene_id))
  measured <- intersect(primary, universe)
  excluded <- sort(setdiff(primary, measured))
  missing[[phase]] <- data.frame(phase=rep(phase, length(excluded)), gene_id=excluded)
  membership <- data.frame(phase=phase, go_id=pathways$go_id, description=pathways$description,
    annotation_size=lengths(pathways$genes),
    measured_size=vapply(pathways$genes, function(g) length(intersect(g, universe)), integer(1)),
    primary_ct_count=vapply(pathways$genes, function(g) length(intersect(g, primary)), integer(1)),
    measured_ct_count=vapply(pathways$genes, function(g) length(intersect(g, measured)), integer(1)),
    excluded_ct_ids=vapply(pathways$genes, function(g) paste(sort(intersect(g, excluded)), collapse=";"), ""))
  write_tsv(membership, file.path(out, paste0(phase, "_membership.tsv")))
  candidates <- membership[membership$annotation_size >= 5 & membership$annotation_size < 1500 &
                           membership$measured_ct_count > 0, , drop=FALSE]
  sizes <- candidates$annotation_size
  if (any(sizes > length(universe))) stop("Historical annotation size exceeds the supplied universe")
  p <- phyper(candidates$measured_ct_count-1, sizes, length(universe)-sizes, length(measured), lower.tail=FALSE)
  old <- 1-phyper(candidates$measured_ct_count-1, sizes, length(universe)-sizes, length(measured))
  calculated <- cbind(candidates, universe_size=length(universe), ct_population=length(measured),
    p_value=p, p_bh=p.adjust(p, "BH"), one_minus_cdf=old, one_minus_cdf_bh=p.adjust(old, "BH"))
  write_tsv(calculated, file.path(out, paste0(phase, "_historical_hypergeometric.tsv")))
  summaries[[phase]] <- data.frame(phase=phase, gene_records=nrow(rows), unique_genes=length(universe),
    primary_ct=length(primary), measured_ct=length(measured), excluded_ct=length(excluded),
    historical_candidates=nrow(candidates))
}
write_tsv(do.call(rbind, summaries), file.path(out, "input_summary.tsv"))
write_tsv(do.call(rbind, missing), file.path(out, "missing_ctgenes.tsv"))
write_tsv(data.frame(file=basename(sources), md5=unname(tools::md5sum(sources))), file.path(out, "input_checksums.tsv"))
writeLines(c("Historical membership audit; full annotation sizes are a legacy convention.",
             "Only gene_id is used from the supplied table. No expected results are calculation inputs."),
           file.path(out, "scope.txt"))
print(do.call(rbind, summaries))
