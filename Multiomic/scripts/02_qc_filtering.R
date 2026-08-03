if (!exists("PROJECT_ROOT")) source("scripts/00_setup.R")
if (!exists("pbmc")) pbmc <- load_checkpoint("01_raw_multiome.rds")

DefaultAssay(pbmc) <- "ATAC"
has_fragments <- fragment_files_available()

if (RUN_FRAGMENT_QC && has_fragments) {
  # Deprecated in recent Signac but retained as a Windows-compatible fallback
  # while ATACqc requires fragtk.
  pbmc <- suppressWarnings(Signac::NucleosomeSignal(pbmc))
  pbmc <- suppressWarnings(Signac::TSSEnrichment(pbmc, fast = FALSE))
}

qc_features <- c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "nFeature_ATAC")
if ("nucleosome_signal" %in% colnames(pbmc[[]])) qc_features <- c(qc_features, "nucleosome_signal")
if ("TSS.enrichment" %in% colnames(pbmc[[]])) qc_features <- c(qc_features, "TSS.enrichment")

plots <- Seurat::VlnPlot(pbmc, features = qc_features, pt.size = 0.1, combine = FALSE)
plots <- lapply(plots, function(p) {
  p + ggplot2::theme_classic() + ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )
})
qc_panel <- patchwork::wrap_plots(plots, ncol = 3)
save_plot(qc_panel, file.path(QC_FIG_DIR, "pre_filter_qc_violin.png"), 12, 8)

metadata_before <- pbmc[[]]
metadata_before$cell_barcode <- rownames(metadata_before)
write_table(metadata_before, file.path(QC_TABLE_DIR, "cell_metadata_before_qc.csv"))

keep_cells <- pbmc$nCount_RNA > MIN_RNA_COUNTS &
  pbmc$nCount_RNA < MAX_RNA_COUNTS &
  pbmc$nCount_ATAC > MIN_ATAC_COUNTS &
  pbmc$nCount_ATAC < MAX_ATAC_COUNTS
if ("nucleosome_signal" %in% colnames(pbmc[[]])) {
  keep_cells <- keep_cells & pbmc$nucleosome_signal < MAX_NUCLEOSOME_SIGNAL
}
if ("TSS.enrichment" %in% colnames(pbmc[[]])) {
  keep_cells <- keep_cells & pbmc$TSS.enrichment > MIN_TSS_ENRICHMENT
}

qc_decisions <- data.frame(cell_barcode = colnames(pbmc), pass_qc = as.logical(keep_cells))
write_table(qc_decisions, file.path(QC_TABLE_DIR, "cell_qc_decisions.csv"))

n_before <- ncol(pbmc)
pbmc <- subset(pbmc, cells = colnames(pbmc)[keep_cells])
n_after <- ncol(pbmc)
write_table(
  data.frame(metric = c("cells_before_qc", "cells_after_qc", "cells_removed"),
             value = c(n_before, n_after, n_before - n_after)),
  file.path(QC_TABLE_DIR, "qc_summary.csv")
)

metadata_after <- pbmc[[]]
metadata_after$cell_barcode <- rownames(metadata_after)
write_table(metadata_after, file.path(QC_TABLE_DIR, "cell_metadata_after_qc.csv"))
save_checkpoint(pbmc, "02_qc_filtered.rds")
