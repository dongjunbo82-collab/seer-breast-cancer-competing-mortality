# Run from repository root. Check dependencies without installing software.
needed <- c("data.table", "splines", "nnet", "MASS", "survival", "ggplot2", "patchwork", "svglite", "ragg", "posterior")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("Install required R packages: ", paste(missing, collapse=", "))
for (p in c("work", "work/models", "work/lifetables", "outputs/tables", "outputs/figures")) dir.create(p, recursive=TRUE, showWarnings=FALSE)
cat("Dependencies available; local runtime directories prepared.\n")
