# Forensic historical replay; see the Python counterpart's module documentation.
# These intentionally labeled legacy count mistakes are NOT recommended tests.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript scripts/run_historical_diagnostics.R SCORES OUTPUT [B] [SEED]")
source("R/core.R", encoding = "UTF-8")
rows <- read_tsv(args[1]); out <- args[2]
B <- if (length(args) >= 3) as.integer(args[3]) else 10000L
seed_base <- if (length(args) >= 4) as.integer(args[4]) else 20261021L
if (is.na(B) || B < 1) stop("B must be positive")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
pathways <- read_gmt("data/input/pathways.gmt")
universe <- unique(rows$gene_id)
summaries <- list()
for (phase in c("short", "middle")) {
  ct <- read_tsv(paste0("data/input/ct_",phase,".tsv"))
  foreground <- unique(ct$gene_id)
  keep <- lengths(pathways$genes) >= 5 & lengths(pathways$genes) < 1500 &
    vapply(pathways$genes, function(g) any(g %in% foreground), logical(1))
  candidates <- pathways[keep,,drop=FALSE]
  annotation_union <- unique(unlist(candidates$genes))
  mixed_union <- union(annotation_union, as.character(lengths(candidates$genes)))
  values <- rows[[paste0("score_",phase)]]
  if (any(!is.finite(values))) stop("Finite recalculated scores required")
  indices <- lapply(candidates$genes, function(g) which(rows$gene_id %in% g))
  sizes <- lengths(indices)
  observed <- vapply(indices, function(i) sum(values[i]), numeric(1))
  counts <- integer(nrow(candidates))
  seed <- seed_base + match(phase,c("short","middle"))-1L
  set.seed(seed)
  distinct <- sort(unique(sizes))
  for (number in seq_along(distinct)) {
    size <- distinct[number]
    null <- replicate(B, sum(values[sample.int(length(values),size,replace=FALSE,useHash=TRUE)]))
    selected <- which(sizes == size)
    counts[selected] <- vapply(observed[selected],function(x) sum(null > x), integer(1))
    if (number %% 25L == 0L) message(phase,": historical SUMSTAT ",number,"/",length(distinct)," sizes")
  }
  base <- data.frame(go_id=candidates$go_id,description=candidates$description,
    annotation_size=lengths(candidates$genes),tested_size=sizes,
    ct_count=vapply(candidates$genes,function(g) sum(g %in% foreground),integer(1)))
  p <- counts/B
  write_tsv(cbind(base,statistic=observed,null_exceedances=counts,p_value=p,p_bh=p.adjust(p,"BH"),
    permutations=B,seed=seed),file.path(out,paste0(phase,"_sumstat.tsv")))
  profiles <- list(legacy_container_bug=c(length(mixed_union),2*length(foreground)),
    genome_with_primary_ct=c(length(universe),length(foreground)),
    genome_with_measured_ct=c(length(universe),length(intersect(universe,foreground))))
  variants <- lapply(names(profiles),function(name) {
    N <- profiles[[name]][1]; M <- profiles[[name]][2]
    p <- phyper(base$ct_count-1,base$annotation_size,N-base$annotation_size,M,lower.tail=FALSE)
    old <- 1-phyper(base$ct_count-1,base$annotation_size,N-base$annotation_size,M)
    cbind(base,profile=name,universe_size=N,ct_population=M,p_value=p,p_bh=p.adjust(p,"BH"),one_minus_cdf=old)
  })
  write_tsv(do.call(rbind,variants),file.path(out,paste0(phase,"_hypergeometric_variants.tsv")))
  summaries[[phase]] <- data.frame(phase=phase,records=nrow(rows),unique_genes=length(universe),
    candidates=nrow(candidates),distinct_record_sizes=length(distinct),annotation_union=length(annotation_union),
    mixed_union=length(mixed_union),ct_input=length(foreground),ct_measured=length(intersect(universe,foreground)))
}
write_tsv(do.call(rbind,summaries),file.path(out,"input_summary.tsv"))
dput(list(score_input=args[1],permutations=B,seed=seed_base,duplicate_policy="retain records",
  strict_exceedance=TRUE,no_exceedance_pseudocount=TRUE,paper_result_inputs=FALSE),file=file.path(out,"settings.R"))
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
message("Outputs: ",out)
