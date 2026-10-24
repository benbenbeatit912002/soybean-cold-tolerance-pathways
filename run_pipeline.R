# Complete base-R workflow. Locate this file rather than the caller's directory.
# Usage: Rscript --vanilla run_pipeline.R --mode demo
#        Rscript --vanilla run_pipeline.R --mode full --use-reconstruction-hypotheses
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE))
args <- commandArgs(trailingOnly = TRUE)
usage <- paste(
  "Usage: Rscript --vanilla run_pipeline.R [--mode demo|quick|full]",
  "       [--output NEW_DIRECTORY] [--use-reconstruction-hypotheses] [--include-maxmean]",
  "       [--crosstalk reference|reanalysis]",
  "Demo uses invented data. Quick/full recalculate the included real inputs.",
  "Quick uses 199 draws; full uses 10000. Real-data repairs require the explicit flag.",
  "Output paths are relative to your current directory. Existing results are not overwritten.", sep = "\n")
if ("--help" %in% args) { cat(usage, "\n"); quit(status = 0L) }
mode <- "demo"
output <- NULL
hypotheses <- FALSE
maxmean <- FALSE
crosstalk <- NULL
i <- 1L
while (i <= length(args)) {
  arg <- args[i]
  if (arg %in% c("--mode", "--output", "--crosstalk")) {
    if (i == length(args)) stop("Missing value for ", arg)
    if (arg == "--mode") mode <- args[i+1L] else if (arg == "--output") output <- args[i+1L] else crosstalk <- args[i+1L]
    i <- i+2L
  } else if (arg == "--use-reconstruction-hypotheses") { hypotheses <- TRUE; i <- i+1L
  } else if (arg == "--include-maxmean") { maxmean <- TRUE; i <- i+1L
  } else stop("Unknown argument: ", arg, "\n", usage)
}
if (!mode %in% c("demo", "quick", "full")) stop("Mode must be demo, quick, or full")
if (mode == "demo" && hypotheses) stop("Demo uses invented scores and does not accept reconstruction hypotheses")
if (!is.null(crosstalk) && !crosstalk %in% c("reference", "reanalysis")) stop("Crosstalk must be reference or reanalysis")
if (mode == "demo" && identical(crosstalk, "reference")) stop("The historical soybean selection cannot be applied to invented demo pathways")
if (is.null(output)) output <- file.path(root, "results", paste0("r_", mode))
if (file.exists(output) && (!dir.exists(output) || length(list.files(output, all.files=TRUE, no..=TRUE)))) {
  stop("Output must be new or empty. Use --output with a different directory.")
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output, winslash = "/", mustWork = TRUE)
setwd(root)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
steps <- data.frame(script=character(), exit_code=integer(), elapsed_seconds=numeric())
status <- "running"
save_manifest <- function() {
  write.table(steps, file.path(output, "execution_steps.tsv"), sep="\t", row.names=FALSE, quote=TRUE)
  values <- c(language="R", mode=mode, status=status, hypotheses_enabled=hypotheses,
              historical_repairs_confirmed=FALSE, paper_reproduction_verified=FALSE,
              permutations=if (mode == "full") 10000L else 199L, include_maxmean=maxmean,
              crosstalk=if (is.null(crosstalk)) "none" else crosstalk)
  write.table(data.frame(setting=names(values), value=unname(values)),
              file.path(output, "run_manifest.tsv"), sep="\t", row.names=FALSE, quote=TRUE)
}
run_step <- function(script, arguments=character()) {
  cat("\nRunning", script, "\n")
  start <- proc.time()[["elapsed"]]
  result <- system2(rscript, shQuote(c("--vanilla", file.path(root, script), arguments)))
  steps <<- rbind(steps, data.frame(script=script, exit_code=result,
                                   elapsed_seconds=proc.time()[["elapsed"]]-start))
  save_manifest()
  if (result != 0L) {
    status <<- if (result == 2L && grepl("recalculate", script)) "unresolved_scores" else "failed"
    save_manifest()
    message("Stopped. Inspect ", output, "; no downstream success is being reported.")
    quit(status=as.integer(result))
  }
}

# 1. Rebuild scores independently in R; no Python runtime is needed.
save_manifest()
if (mode == "demo") {
  inputs <- file.path(root, "examples/data")
  scores <- file.path(inputs, "scores.tsv")
} else {
  inputs <- file.path(root, "data/input")
  flags <- if (hypotheses) "--use-reconstruction-hypotheses" else character()
  run_step("scripts/recalculate_phase_scores.R", c(file.path(root, "data/upstream"), file.path(output, "scores"), flags))
  scores <- file.path(output, "scores/recalculated_scores.tsv")
}

# 2. Save every setting, including the exact named score columns and null size.
config <- list(input_dir=inputs, score_input=scores, score_field="score_{phase}",
  score_column=NA_integer_, score_kind="score", phases=c("short", "middle"),
  min_size=if (mode == "demo") 2L else 6L, max_size=if (mode == "demo") 18L else 1500L,
  permutations=if (mode == "full") 10000L else 199L, seed=20261021L,
  duplicate_policy="mean", methods=c("sumstat", "gsea", if (maxmean) "maxmean"),
  output_dir=file.path(output, "analysis"),
  score_provenance=if (mode == "demo") "Invented teaching data" else
    "Recalculated historical p-value branch with explicit repair hypotheses",
  paper_reproduction_verified=FALSE)
config_path <- file.path(output, "analysis_config.R")
writeLines(c("config <-", capture.output(dput(config))), config_path)

# 3. Exact overlap probabilities, score simulations, and labelled figures.
run_step("scripts/run_analysis.R", config_path)
run_step("scripts/plot_results.R", file.path(output, "analysis"))
if (!is.null(crosstalk)) {
  flags <- if (crosstalk == "reference") c("--selected-pathways",file.path(root,"data/reference/crosstalk/selected_pathways.tsv"),
    "--selection-label","Historical Excel selection","--membership","annotation","--edge-rule","historical-code") else character()
  run_step("scripts/run_presoynet.R",c(file.path(output,"analysis"),file.path(output,"crosstalk"),"--input-dir",inputs,flags))
}
status <- "completed"
save_manifest()
cat("\nCompleted", mode, ". Tables and plots:", file.path(output, "analysis"), "\n")
if (mode != "full") cat("This is an execution/teaching run; use full mode for the 10,000-draw analysis.\n")
