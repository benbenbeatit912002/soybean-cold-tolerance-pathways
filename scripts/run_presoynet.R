# Independent base-R entry point for offline pathway crosstalk.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
root <- dirname(dirname(normalizePath(script,winslash="/",mustWork=TRUE)))
source(file.path(root,"R/core.R"))
source(file.path(root,"R/presoynet.R"))
args <- commandArgs(trailingOnly=TRUE)
usage <- paste("Usage: Rscript --vanilla scripts/run_presoynet.R ANALYSIS OUTPUT",
  "[--input-dir DIR] [--selected-pathways TSV --selection-label LABEL]",
  "[--membership measured|annotation] [--keep-short 0.8] [--keep-middle 0.5]",
  "[--alpha 0.05] [--min-ct 1] [--min-ratio 0] [--method-rule any|all]",sep="\n")
if ("--help" %in% args) { cat(usage,"\n"); quit(status=0L) }
if (length(args)<2) stop(usage)
settings <- list(analysis=args[1],output=args[2],inputs=file.path(root,"data/input"),
  alpha=.05,min_ct=1,min_ratio=0,method_rule="any",membership="measured",keep_short=.8,
  keep_middle=.5,selected_pathways=NULL,selection_label="")
numeric_args <- c("alpha","min_ct","min_ratio","keep_short","keep_middle")
i <- 3L
while (i <= length(args)) {
  key <- gsub("-","_",sub("^--","",args[i]),fixed=TRUE)
  if (key=="input_dir") key <- "inputs"
  if (!startsWith(args[i],"--") || !key %in% names(settings) || key %in% c("analysis","output") || i==length(args)) stop("Unknown or incomplete argument: ",args[i])
  settings[[key]] <- if (key %in% numeric_args) as.numeric(args[i+1]) else args[i+1]
  i <- i+2L
}
files <- c(file.path(settings$inputs,"pathways.gmt"),settings$selected_pathways)
for (phase in phases) files <- c(files,file.path(settings$inputs,paste0("ct_",phase,".tsv")),
  file.path(settings$analysis,paste0(phase,"_analysis_scores.tsv")),
  file.path(settings$analysis,paste0(phase,"_",main_methods,".tsv")))
hashes <- tools::md5sum(files)
if (anyNA(hashes)) stop("Missing input file")
result <- do.call(run_presoynet,settings)
plot_crosstalk(result$nodes,result$edges,file.path(settings$output,"crosstalk.png"))
write_network_tsv(data.frame(file=normalizePath(files,winslash="/"),md5=unname(hashes)),file.path(settings$output,"input_checksums.tsv"))
values <- c(status="completed",language="R",selection_mode=if (is.null(settings$selected_pathways)) "bh_reanalysis" else "provided_list",
  unlist(settings),paper_reproduction_verified=FALSE,soynet_contacted=FALSE,
  edge_rule="ceil(fraction * positive undirected pairs), include all cutoff ties")
write_network_tsv(data.frame(setting=names(values),value=unname(values)),file.path(settings$output,"manifest.tsv"))
cat("Completed offline crosstalk:",nrow(result$nodes),"nodes,",nrow(result$edges),"edges. Gene lists:",file.path(settings$output,"soynet_inputs"),"\n")
