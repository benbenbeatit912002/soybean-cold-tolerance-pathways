# Offline pathway selection and overlap networks. Source R/core.R first.
phases <- c("short", "middle")
main_methods <- c("hypergeometric", "sumstat", "gsea")
audit_columns <- c("phase", "go_id", "description", "annotation_size", "tested_size",
  "ct_count", "ct_ratio", paste0(main_methods, "_bh"), "methods_below_alpha", "selected", "selection_basis")

checked_ids <- function(rows, field, label) {
  if (!field %in% names(rows)) stop("Missing ID column in ", label)
  ids <- canonical_id(rows[[field]])
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop("Missing or duplicate IDs in ", label)
  ids
}

select_pathways <- function(pathways, tables, universe, foreground, phase,
    alpha=.05, min_ct=1, min_ratio=0, method_rule="any", requested=NULL) {
  if (!is.finite(alpha) || alpha <= 0 || alpha > 1 || !is.finite(min_ct) ||
      min_ct < 0 || min_ct != floor(min_ct) || !is.finite(min_ratio) || min_ratio < 0 || min_ratio > 1) {
    stop("Invalid selection thresholds")
  }
  if (!method_rule %in% c("any", "all")) stop("method_rule must be any or all")
  checked_ids(pathways, "go_id", "pathways")
  ids <- lapply(main_methods, function(m) checked_ids(tables[[m]], "go_id", m))
  tested <- ids[[1]]
  if (!all(vapply(ids, setequal, logical(1), tested))) stop("Methods must contain the same complete tested GO family")
  if (!all(c(tested, requested) %in% pathways$go_id)) stop("Unknown GO ID in results or supplied selection")
  audit <- as.data.frame(setNames(replicate(length(audit_columns), character(), simplify=FALSE), audit_columns))
  for (go in sort(unique(c(tested, requested)))) {
    p <- pathways[match(go, pathways$go_id), , drop=FALSE]
    genes <- intersect(p$genes[[1]], universe)
    ct <- length(intersect(genes, foreground))
    ratio <- if (length(genes)) ct / length(genes) else 0
    values <- setNames(rep(NA_real_, 3), main_methods)
    for (method in main_methods) {
      table <- tables[[method]]
      index <- match(go, table$go_id)
      if (is.na(index)) next
      row <- table[index, , drop=FALSE]
      expected <- c(annotation_size=length(p$genes[[1]]), tested_size=length(genes), ct_count=ct)
      if (anyNA(row[names(expected)]) || any(as.numeric(row[names(expected)]) != expected)) {
        stop("Membership mismatch: ", phase, " ", method, " ", go)
      }
      if (!is.finite(row$ct_ratio) || abs(row$ct_ratio-ratio) > 1e-12) stop("CT ratio mismatch: ", go)
      value <- row[[if (method == "hypergeometric") "hypergeometric_bh" else "p_bh"]]
      if (!is.finite(value) || value < 0 || value > 1) stop("Invalid BH probability: ", method, " ", go)
      values[method] <- value
    }
    passing <- main_methods[!is.na(values) & values < alpha]
    accepted <- if (!is.null(requested)) go %in% requested else
      ct >= min_ct && ratio >= min_ratio && length(passing) >= if (method_rule == "any") 1 else 3
    row <- data.frame(phase=phase, go_id=go, description=p$description,
      annotation_size=length(p$genes[[1]]), tested_size=length(genes), ct_count=ct, ct_ratio=ratio,
      hypergeometric_bh=values[1], sumstat_bh=values[2], gsea_bh=values[3],
      methods_below_alpha=paste(passing, collapse=";"), selected=as.integer(accepted),
      selection_basis=if (is.null(requested)) "bh_reanalysis" else "provided_list", row.names=NULL)
    audit <- rbind(audit, row)
  }
  selected <- pathways[match(audit$go_id[audit$selected == 1], pathways$go_id), , drop=FALSE]
  list(audit=audit, selected=selected)
}

jaccard_edges <- function(memberships, keep_fraction) {
  if (!is.finite(keep_fraction) || keep_fraction < 0 || keep_fraction > 1) stop("Edge fraction must be between zero and one")
  if (any(lengths(memberships) == 0)) stop("A selected pathway has no genes under this membership policy")
  edges <- data.frame(source=character(), target=character(), intersection_size=integer(),
                      union_size=integer(), jaccard=numeric(), retained=integer())
  ids <- sort(names(memberships))
  if (length(ids) >= 2) {
    pairs <- combn(ids, 2)
    edges <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(i) {
      a <- pairs[1,i]; b <- pairs[2,i]
      common <- length(intersect(memberships[[a]], memberships[[b]]))
      total <- length(union(memberships[[a]], memberships[[b]]))
      data.frame(source=a, target=b, intersection_size=common, union_size=total,
                 jaccard=common/total, retained=0L)
    }))
  }
  positive <- sort(edges$jaccard[edges$intersection_size > 0], decreasing=TRUE)
  count <- ceiling(keep_fraction * length(positive))
  cutoff <- if (count) positive[count] else NA_real_
  if (count) edges$retained <- as.integer(edges$intersection_size > 0 & edges$jaccard >= cutoff)
  list(edges=edges, cutoff=cutoff)
}

# Empty values indicate unavailable statistics, never a fabricated zero p-value.
write_network_tsv <- function(x, path) {
  write.table(x, path, sep="\t", row.names=FALSE, quote=TRUE, na="", fileEncoding="UTF-8")
}

run_presoynet <- function(analysis, inputs, output, alpha=.05, min_ct=1, min_ratio=0,
    method_rule="any", membership="measured", keep_short=.8, keep_middle=.5,
    selected_pathways=NULL, selection_label="") {
  if (file.exists(output) && (!dir.exists(output) || length(list.files(output, all.files=TRUE, no..=TRUE)))) stop("Output must be new or empty")
  if (!membership %in% c("measured", "annotation")) stop("membership must be measured or annotation")
  requested <- NULL
  if (!is.null(selected_pathways)) {
    if (!nzchar(trimws(selection_label))) stop("A supplied selection requires --selection-label")
    provided <- read_tsv(selected_pathways)
    if (!"phase" %in% names(provided) || anyNA(provided$phase) || !all(provided$phase %in% phases)) stop("Selected pathway phase must be short or middle")
    requested <- setNames(lapply(phases, function(ph) checked_ids(provided[provided$phase == ph, , drop=FALSE], "go_id", "selection")), phases)
  }
  pathways <- read_gmt(file.path(inputs, "pathways.gmt"))
  data <- list(); summaries <- list()
  for (phase in phases) {
    scores <- read_tsv(file.path(analysis, paste0(phase, "_analysis_scores.tsv")))
    universe <- checked_ids(scores, "gene_id", "analysis scores")
    if (!length(universe) || any(!is.finite(scores$analysis_score) | scores$analysis_score < 0)) stop("Analysis scores must be nonempty and finite")
    ct <- checked_ids(read_tsv(file.path(inputs, paste0("ct_", phase, ".tsv"))), "gene_id", "CT list")
    tables <- setNames(lapply(main_methods, function(m) read_tsv(file.path(analysis, paste0(phase, "_", m, ".tsv")))), main_methods)
    selection <- select_pathways(pathways, tables, universe, ct, phase, alpha, min_ct, min_ratio, method_rule, requested[[phase]])
    members <- setNames(lapply(selection$selected$genes, function(g) if (membership == "measured") intersect(g, universe) else g), selection$selected$go_id)
    fraction <- if (phase == "short") keep_short else keep_middle
    network <- jaccard_edges(members, fraction)
    data[[phase]] <- c(selection, list(universe=universe, ct=ct, members=members, edges=network$edges))
    summaries[[phase]] <- data.frame(phase=phase, tested_pathways=nrow(tables$hypergeometric),
      selected_pathways=nrow(selection$selected), positive_edges=sum(network$edges$intersection_size > 0),
      retained_edges=sum(network$edges$retained), keep_fraction=fraction, jaccard_cutoff=network$cutoff)
  }
  dir.create(output, recursive=TRUE, showWarnings=FALSE)
  put <- function(x, name) write_network_tsv(x, file.path(output, name))
  membership_rows <- data.frame(phase=character(), go_id=character(), gene_id=character(), is_ct=integer(), is_measured=integer())
  for (phase in phases) {
    d <- data[[phase]]
    put(d$audit, paste0(phase, "_selection_audit.tsv"))
    put(d$audit[d$audit$selected == 1, , drop=FALSE], paste0(phase, "_selected_pathways.tsv"))
    put(d$edges, paste0(phase, "_all_pairs.tsv"))
    ids <- sort(as.character(names(d$members)))
    matrix <- matrix(0, length(ids), length(ids), dimnames=list(ids, ids)); diag(matrix) <- 1
    for (i in seq_len(nrow(d$edges))) {
      e <- d$edges[i,]; matrix[e$source,e$target] <- matrix[e$target,e$source] <- e$jaccard
    }
    put(data.frame(go_id=ids, matrix, check.names=FALSE), paste0(phase, "_jaccard.tsv"))
    for (go in ids) {
      genes <- sort(d$members[[go]])
      membership_rows <- rbind(membership_rows, data.frame(phase=phase, go_id=go, gene_id=genes,
        is_ct=as.integer(genes %in% d$ct), is_measured=as.integer(genes %in% d$universe)))
    }
  }
  put(membership_rows, "pathway_genes.tsv")
  node_ids <- sort(union(names(data$short$members), names(data$middle$members)))
  nodes <- data.frame(go_id=character(), description=character(), phase=character(),
                      network_gene_count=integer(), short_ct_count=integer(), middle_ct_count=integer())
  for (go in node_ids) {
    present <- phases[vapply(data, function(d) go %in% names(d$members), logical(1))]
    genes <- unique(unlist(lapply(data[present], function(d) d$members[[go]]), use.names=FALSE))
    nodes <- rbind(nodes, data.frame(go_id=go, description=pathways$description[match(go,pathways$go_id)],
      phase=if (length(present) == 2) "both" else present, network_gene_count=length(genes),
      short_ct_count=length(intersect(data$short$members[[go]], data$short$ct)),
      middle_ct_count=length(intersect(data$middle$members[[go]], data$middle$ct))))
  }
  put(nodes, "network_nodes.tsv")
  phase_edges <- lapply(phases, function(ph) {
    e <- data[[ph]]$edges; e <- e[e$retained == 1, c("source", "target", "jaccard"), drop=FALSE]
    names(e)[3] <- paste0(ph, "_jaccard"); e
  })
  edges <- merge(phase_edges[[1]], phase_edges[[2]], by=c("source", "target"), all=TRUE, sort=TRUE)
  put(edges, "network_edges.tsv")
  dir.create(file.path(output, "soynet_inputs"))
  shared <- intersect(names(data$short$members), names(data$middle$members))
  groups <- list(short=list(short=names(data$short$members)), middle=list(middle=names(data$middle$members)),
    shared_pathways=list(short=shared,middle=shared),
    combined=list(short=names(data$short$members),middle=names(data$middle$members)))
  exports <- list()
  for (group in names(groups)) {
    phase_ids <- groups[[group]]
    genes <- sort(unique(unlist(lapply(names(phase_ids), function(ph) unlist(data[[ph]]$members[phase_ids[[ph]]], use.names=FALSE)), use.names=FALSE)))
    ct <- sort(intersect(genes, unique(unlist(lapply(data[names(phase_ids)], `[[`, "ct"), use.names=FALSE))))
    writeLines(as.character(genes), file.path(output, "soynet_inputs", paste0(group, "_genes.txt")), useBytes=TRUE)
    writeLines(as.character(ct), file.path(output, "soynet_inputs", paste0(group, "_ctgenes.txt")), useBytes=TRUE)
    exports[[group]] <- data.frame(group=group,genes=length(genes),ctgenes=length(ct))
  }
  put(do.call(rbind, exports), "export_summary.tsv")
  summary <- do.call(rbind, summaries); put(summary, "network_summary.tsv")
  list(nodes=nodes, edges=edges, summary=summary)
}

plot_crosstalk <- function(nodes, edges, path) {
  colors <- c(short="#377EB8",both="#DDAA22",middle="#D86999")
  positions <- matrix(0,nrow(nodes),2,dimnames=list(nodes$go_id,c("x","y")))
  for (phase in names(colors)) {
    ids <- sort(nodes$go_id[nodes$phase == phase]); n <- length(ids)
    if (n) {
      angle <- 2*pi*(seq_len(n)-1)/n
      positions[ids,1] <- match(phase,names(colors))-1 + if (n > 2) .27*cos(angle) else 0
      positions[ids,2] <- if (n > 2) .5+.4*sin(angle) else 1-seq_len(n)/(n+1)
    }
  }
  png(path,width=1800,height=1050,res=150); on.exit(dev.off())
  par(mar=c(1,1,4,1))
  plot(NA,xlim=c(-.5,2.6),ylim=c(-.06,1.1),axes=FALSE,xlab="",ylab="",
       main="Pathway crosstalk from gene overlap")
  for (i in seq_len(nrow(edges))) {
    e <- edges[i,]; a <- positions[e$source,]; b <- positions[e$target,]
    weight <- max(e$short_jaccard,e$middle_jaccard,na.rm=TRUE)
    segments(a[1],a[2],b[1],b[2],col=adjustcolor("#777777",alpha.f=.3),lwd=.4+2*weight)
  }
  if (nrow(nodes)) {
    points(positions,pch=21,bg=colors[nodes$phase],col=colors[nodes$phase],cex=sqrt((40+2*sqrt(nodes$network_gene_count))/40))
    text(positions[,1]+.04,positions[,2],labels=nodes$go_id,adj=c(0,.5),cex=.65)
  } else text(1,.5,"No pathways pass the selected rule")
  text(0:2,1.03,c("Short","Both","Middle"),col=colors,font=2)
}
