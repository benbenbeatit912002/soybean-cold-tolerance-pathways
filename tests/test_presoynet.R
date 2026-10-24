source("R/core.R")
source("R/presoynet.R")
members <- list(a=c("x","y"), b=c("y","z"), c=c("x","z"), d="q")
network <- jaccard_edges(members,.5)
stopifnot(nrow(network$edges)==6, network$cutoff==1/3, sum(network$edges$retained)==3,
          !any(network$edges$retained[network$edges$target=="d"]),
          sum(jaccard_edges(members,0)$edges$retained)==0,
          sum(jaccard_edges(members,1)$edges$retained)==3,
          nrow(jaccard_edges(list(),.8)$edges)==0,
          nrow(jaccard_edges(list(a="x"),.8)$edges)==0)
fails <- function(expr) inherits(tryCatch({force(expr); NULL},error=identity),"error")
stopifnot(fails(jaccard_edges(list(a=character()),.8)), fails(jaccard_edges(list(),NaN)))
paths <- data.frame(go_id=c("A","B"),description=c("A","B"),genes=I(list(c("x","y"),"z")))
metadata <- data.frame(go_id="A",annotation_size=2,tested_size=2,ct_count=1,ct_ratio=.5)
tables <- list(hypergeometric=cbind(metadata,hypergeometric_bh=.01),
               sumstat=cbind(metadata,p_bh=.8),gsea=cbind(metadata,p_bh=.8))
selected <- select_pathways(paths,tables,c("x","y","z"),"x","short",requested="B")
stopifnot(identical(selected$selected$go_id,"B"),is.na(selected$audit$sumstat_bh[2]),
          selected$audit$selection_basis[2]=="provided_list",
          nrow(select_pathways(paths,tables,c("x","y","z"),"x","short",method_rule="all")$selected)==0)
tables$sumstat$ct_count <- 2
stopifnot(fails(select_pathways(paths,tables,c("x","y","z"),"x","short")))
tables$sumstat$ct_count <- 1; tables$sumstat$p_bh <- NaN
stopifnot(fails(select_pathways(paths,tables,c("x","y","z"),"x","short")))
cat("Crosstalk known-set and selection tests passed.\n")
