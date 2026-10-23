# Read saved R workspaces as data, without running archived functions.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Usage: inspect_saved_r_state.R OUTPUT_DIRECTORY WORKSPACE...")
dir.create(args[1L], recursive = TRUE, showWarnings = FALSE)
for (file_index in seq_along(args[-1L])) {
  file <- args[-1L][file_index]
  objects <- new.env(parent = emptyenv())
  load(file, envir = objects)
  inventory <- lapply(ls(objects, all.names = TRUE), function(name) {
    value <- get(name, envir = objects)
    data.frame(file = file, object = name, type = typeof(value),
               class = paste(class(value), collapse = ","),
               length = length(value), dimensions = paste(dim(value), collapse = "x"),
               columns = paste(head(colnames(value), 15L), collapse = ","))
  })
  inventory <- do.call(rbind, inventory)
  stem <- paste0("workspace_", file_index)
  write.table(inventory, file.path(args[1L], paste0(stem, "_inventory.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Workspace", file_index, ":", nrow(inventory), "objects\n")
  print(head(inventory[, c("object", "type", "length", "dimensions", "columns")], 80L), row.names = FALSE)
  for (name in ls(objects, all.names = TRUE)) {
    value <- get(name, envir = objects)
    # Wide GO membership tables are not correction-session numeric tables.
    if ((is.data.frame(value) || is.matrix(value)) && ncol(value) <= 30L) {
      hits <- which(apply(value, 1L, function(row) any(grepl("14[gG]129500", as.character(row)))))
      if (length(hits)) {
        cat("Target records in", file, name, "\n")
        print(value[hits, , drop = FALSE])
        write.table(value[hits, , drop = FALSE],
                    file.path(args[1L], paste0(stem, "_", name, "_target.tsv")),
                    sep = "\t", row.names = TRUE, quote = TRUE)
      }
    }
  }
}
