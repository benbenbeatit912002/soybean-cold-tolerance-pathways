# Recalculate both phases. Expected/result score files are never read.
args <- commandArgs(trailingOnly = TRUE)
hypotheses <- "--use-reconstruction-hypotheses" %in% args
args <- args[args != "--use-reconstruction-hypotheses"]
if (length(args) != 2L) stop("Usage: recalculate_phase_scores.R SOURCE_ROOT OUTPUT [--use-reconstruction-hypotheses]")
source("R/historical_short.R")
rebuilt <- reconstruct_historical_input(args[1L])
rows <- rebuilt$rows
imputation <- historical_snp_imputation(rows)
counts <- imputation$values
rows$historical_snp <- counts
fits <- lapply(c("short", "middle"), function(phase) fit_historical_short(rows[[paste0("raw_",phase)]], counts))
names(fits) <- c("short", "middle")
final_counts <- counts
short_adjusted <- fits$short$adjusted
if (hypotheses) {
  repair <- repair_zero_hypothesis(rows, counts, short_adjusted, args[1L])
  final_counts <- repair$counts
  short_adjusted <- repair$adjusted
}
output <- args[2L]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write_precise <- function(data, name) {
  for (field in names(data)) if (is.numeric(data[[field]])) data[[field]] <- sprintf("%.17g", data[[field]])
  write.table(data, file.path(output,name), sep="\t", quote=FALSE, row.names=FALSE, fileEncoding="UTF-8")
}
exceptions <- list()
for (phase in names(fits)) {
  fit <- fits[[phase]]
  adjusted <- if (phase == "short") short_adjusted else fit$adjusted
  repaired <- which(counts == 0 & final_counts > 0)
  if (phase == "middle" && hypotheses && length(repaired)) {
    adjusted[repaired] <- as.numeric(sprintf("%.7g", historical_fosco(rows$raw_middle[repaired], final_counts[repaired], fit$lambda)))
  }
  export <- historical_score_export(adjusted)
  finite_max <- max(export$score[is.finite(export$score)])
  status <- rep("recalculated", nrow(rows))
  status[repaired] <- "zero_snp_repair_hypothesis"
  overflow <- which(!is.finite(export$score) & final_counts > 0)
  if (hypotheses && length(overflow)) {
    export$display[overflow] <- sprintf("%.7g", finite_max)
    export$score[overflow] <- as.numeric(export$display[overflow])
    status[overflow] <- "maximum_finite_score_repair_hypothesis"
  }
  status[!is.finite(export$score)] <- "unresolved_nonfinite"
  rows$snp <- final_counts
  rows[[paste0("adjusted_",phase)]] <- adjusted
  rows[[paste0("adjusted_",phase,"_serialized")]] <- export$probability
  rows[[paste0("score_",phase,"_unrounded")]] <- export$score
  rows[[paste0("score_",phase)]] <- export$display
  rows[[paste0("status_",phase)]] <- status
  affected <- which(status != "recalculated")
  if (length(affected)) exceptions[[phase]] <- data.frame(phase=phase, source_row=rows$source_row[affected],
    gene_id=rows$gene_id[affected], status=status[affected], historical_snp=counts[affected],
    final_snp=final_counts[affected], adjusted_probability=adjusted[affected], score=export$display[affected],
    calculated_finite_max=finite_max)
  write_precise(fit$grid,paste0(phase,"_lambda_grid.tsv"))
}
write_precise(rows,"recalculated_scores.tsv")
write_precise(do.call(rbind,exceptions),"exception_audit.tsv")
write_precise(imputation$audit,"historical_snp_imputation.tsv")
if (hypotheses) write_precise(repair$audit,"snp_repair_donors.tsv")
settings <- c(hypotheses_enabled=hypotheses, historical_repairs_confirmed=FALSE, reference_scores_read=FALSE,
              lambda_short=fits$short$lambda, lambda_middle=fits$middle$lambda)
write_precise(data.frame(setting=names(settings),value=unname(settings)),"score_settings.tsv")
write_precise(data.frame(source=names(rebuilt$paths),path=unname(rebuilt$paths),md5=unname(tools::md5sum(rebuilt$paths))),"source_manifest.tsv")
cat("Both phases recalculated:",nrow(rows),"records; lambdas",fits$short$lambda,fits$middle$lambda,"\n")
if (any(rows$status_short == "unresolved_nonfinite" | rows$status_middle == "unresolved_nonfinite")) quit(status=2L)
