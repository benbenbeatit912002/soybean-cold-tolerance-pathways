# Complete candidate tables and explicit rules; no expected GO IDs or counts.
select_candidates <- function(rows, config) {
  methods <- c("hypergeometric","sumstat","gsea")
  field <- config$probability_field; alpha <- config$alpha
  if (!field %in% c("p_value","p_bh") || !is.finite(alpha) || alpha <= 0 || alpha > 1) stop("Invalid probability field or threshold")
  if (!config$method_rule %in% c("any","all") || !config$count_method %in% methods) stop("Invalid method rule or count source")
  if (!nrow(rows) || anyNA(rows[c("phase","go_id","method")]) || any(!nzchar(trimws(rows$go_id))) ||
      !all(rows$phase %in% c("short","middle")) || !all(rows$method %in% methods)) stop("Invalid phase, method or GO ID")
  keys <- paste(rows$phase,rows$go_id,rows$method,sep="/")
  if (anyDuplicated(keys)) stop("Duplicate phase/GO/method row")
  audit <- list()
  for (phase in c("short","middle")) {
    rule <- config$phases[[phase]]
    if (!is.finite(rule$min_ct) || rule$min_ct < 0 || rule$min_ct != floor(rule$min_ct) ||
        !is.finite(rule$min_ratio) || rule$min_ratio < 0 || rule$min_ratio > 1 ||
        !rule$biological_rule %in% c("count","ratio","and","or")) stop("Invalid biological rule")
    ids <- sort(unique(rows$go_id[rows$phase==phase]))
    if (!length(ids)) stop("Candidate input must contain both phases")
    for (go in ids) {
      group <- rows[rows$phase==phase & rows$go_id==go,,drop=FALSE]
      if (!setequal(group$method,methods)) stop("Incomplete method family for ",phase," ",go)
      group <- group[match(methods,group$method),,drop=FALSE]
      sizes <- group$annotation_size; counts <- group$ct_count; values <- group[[field]]
      if (any(!is.finite(c(sizes,counts,values))) || any(sizes < 1 | sizes != floor(sizes)) ||
          any(counts < 0 | counts > sizes | counts != floor(counts)) || any(values < 0 | values > 1)) stop("Invalid candidate metadata or probability")
      if (length(unique(sizes)) != 1) stop("Annotation sizes differ across methods")
      source <- match(config$count_method,methods)
      count <- counts[source]; size <- sizes[source]; ratio <- count/size
      count_pass <- count >= rule$min_ct; ratio_pass <- ratio >= rule$min_ratio
      biological <- switch(rule$biological_rule,count=count_pass,ratio=ratio_pass,
        and=count_pass && ratio_pass,or=count_pass || ratio_pass)
      passing <- methods[values < alpha]
      statistic <- length(passing) >= if (config$method_rule=="any") 1 else 3
      audit[[length(audit)+1]] <- data.frame(phase=phase,go_id=go,description=group$description[source],
        annotation_size=size,ct_count=count,ct_ratio=ratio,count_method=config$count_method,
        method_ct_counts_agree=as.integer(length(unique(counts))==1), probability_field=field,
        hypergeometric_probability=values[1],sumstat_probability=values[2],gsea_probability=values[3],
        count_pass=as.integer(count_pass),ratio_pass=as.integer(ratio_pass),biological_pass=as.integer(biological),
        methods_passing=paste(passing,collapse=";"),statistic_pass=as.integer(statistic),
        selected=as.integer(biological && statistic),selection_basis=config$label)
    }
  }
  do.call(rbind,audit)
}
