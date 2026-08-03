if (!exists("PROJECT_ROOT")) source("scripts/00_setup.R")
if (!exists("pbmc")) pbmc <- load_checkpoint("08_annotation_validated.rds")

if (!RUN_COVERAGE_PLOTS) {
  message("Coverage plots skipped: RUN_COVERAGE_PLOTS = FALSE")
} else if (!fragment_files_available()) {
  warning("Coverage plots skipped: fragment file/index unavailable.")
} else {
  expression_assay <- if ("SCT" %in% Seurat::Assays(pbmc)) "SCT" else "RNA"
  for (gene in COVERAGE_GENES) {
    tryCatch({
      p <- Signac::CoveragePlot(
        object = pbmc, region = gene, features = gene,
        assay = "ATAC", expression.assay = expression_assay,
        peaks = TRUE, extend.upstream = 10000, extend.downstream = 10000
      )
      save_plot(p, file.path(COVERAGE_DIR, paste0("coverage_", gene, ".png")), 12, 7)
    }, error = function(e) {
      msg <- paste0("CoveragePlot failed for ", gene, ": ", conditionMessage(e))
      warning(msg)
      writeLines(msg, file.path(VALIDATION_DIR, paste0("coverage_", gene, "_error.txt")))
    })
  }
}
