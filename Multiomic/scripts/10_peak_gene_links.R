if (!exists("PROJECT_ROOT")) source("scripts/00_setup.R")
if (!exists("pbmc")) pbmc <- load_checkpoint("08_annotation_validated.rds")

if (!RUN_PEAK_GENE_LINKS) {
  message("Peak-gene linking skipped: RUN_PEAK_GENE_LINKS = FALSE")
} else if (!fragment_files_available()) {
  warning("Peak-gene linking skipped: fragment file/index unavailable.")
} else if (!requireNamespace(BSGENOME_PACKAGE, quietly = TRUE)) {
  warning("Peak-gene linking skipped: BSgenome package not installed: ", BSGENOME_PACKAGE)
} else {
  suppressPackageStartupMessages(library(BSGENOME_PACKAGE, character.only = TRUE))
  genome_object <- get(BSGENOME_PACKAGE)
  DefaultAssay(pbmc) <- "ATAC"
  pbmc <- Signac::RegionStats(pbmc, genome = genome_object)
  expression_assay <- if ("SCT" %in% Seurat::Assays(pbmc)) "SCT" else "RNA"
  pbmc <- Signac::LinkPeaks(
    object = pbmc, peak.assay = "ATAC",
    expression.assay = expression_assay,
    genes.use = LINKPEAKS_GENES
  )
  link_ranges <- Signac::Links(pbmc[["ATAC"]])
  write_table(as.data.frame(link_ranges), file.path(PEAK_LINK_DIR, "peak_gene_links.csv"))
  save_checkpoint(pbmc, "10_peak_gene_linked.rds")
}
