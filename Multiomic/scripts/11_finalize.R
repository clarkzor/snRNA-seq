if (!exists("PROJECT_ROOT")) source("scripts/00_setup.R")
if (!exists("pbmc")) {
  candidate <- if (file.exists(file.path(OBJECT_DIR, "10_peak_gene_linked.rds"))) {
    "10_peak_gene_linked.rds"
  } else {
    "08_annotation_validated.rds"
  }
  pbmc <- load_checkpoint(candidate)
}

save_checkpoint(pbmc, "St11_multiome_final.rds")
metadata <- pbmc[[]]
metadata$cell_barcode <- rownames(metadata)
write_table(metadata, file.path(ANNOTATION_TABLE_DIR, "final_cell_metadata.csv"))

packages <- c("R", "Seurat", "Signac", "SeuratObject", "ggplot2", "patchwork",
              "Matrix", "GenomeInfoDb", "GenomicRanges", "ensembldb", ENSDB_PACKAGE)
versions <- vapply(packages, function(pkg) {
  if (pkg == "R") return(as.character(getRversion()))
  if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_
}, character(1))
write_table(data.frame(package = packages, version = versions),
            file.path(VALIDATION_DIR, "package_versions.csv"))
writeLines(capture.output(sessionInfo()), file.path(VALIDATION_DIR, "sessionInfo.txt"))

run_summary <- c(
  paste("Final cells:", ncol(pbmc)),
  paste("Final RNA features:", nrow(pbmc[["RNA"]])),
  paste("Final ATAC peaks:", nrow(pbmc[["ATAC"]])),
  paste("Final WNN resolution:", FINAL_WNN_RESOLUTION),
  paste("Unresolved cells:", if ("final_annotation" %in% colnames(pbmc[[]])) sum(as.character(pbmc$final_annotation) == UNASSIGNED_LABEL) else NA),
  paste("Annotation conflicts:", if ("expression_annotation_conflict" %in% colnames(pbmc[[]])) sum(pbmc$expression_annotation_conflict) else NA)
)
writeLines(run_summary, file.path(VALIDATION_DIR, "run_summary.txt"))
message("Pipeline finalized. Final object: ", file.path(OBJECT_DIR, "St11_multiome_final.rds"))
