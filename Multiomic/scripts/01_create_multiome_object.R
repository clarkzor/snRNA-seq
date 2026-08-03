if (!exists("PROJECT_ROOT")) source("scripts/00_setup.R")

assert_file(FEATURE_MATRIX_H5, "Filtered feature-barcode H5 matrix")
assert_file(RULES_FILE, "Expression annotation rules")

counts <- Seurat::Read10X_h5(FEATURE_MATRIX_H5)
if (!is.list(counts)) stop("Read10X_h5 did not return a modality list.")
required_modalities <- c("Gene Expression", "Peaks")
missing_modalities <- setdiff(required_modalities, names(counts))
if (length(missing_modalities) > 0) {
  stop("H5 is missing modalities: ", paste(missing_modalities, collapse = ", "))
}

rna_counts <- counts[["Gene Expression"]]
atac_counts <- counts[["Peaks"]]
common_cells <- intersect(colnames(rna_counts), colnames(atac_counts))
if (length(common_cells) == 0) stop("RNA and ATAC matrices share no barcodes.")
rna_counts <- rna_counts[, common_cells, drop = FALSE]
atac_counts <- atac_counts[, common_cells, drop = FALSE]

ensdb_object <- get(ENSDB_PACKAGE)
annotation <- Signac::GetGRangesFromEnsDb(ensdb = ensdb_object)
annotation <- harmonize_annotation_seqlevels(annotation, rownames(atac_counts))
has_fragments <- fragment_files_available()

if (!has_fragments) {
  warning(
    "Fragment file/index not found. Matrix-based RNA/ATAC processing can run, ",
    "but fragment QC, coverage plots, tile plots, and LinkPeaks will be skipped."
  )
}

pbmc <- Seurat::CreateSeuratObject(
  counts = rna_counts, assay = "RNA", project = "Xtrop_St11_Multiome",
  min.cells = 0, min.features = 0
)
pbmc[["ATAC"]] <- Signac::CreateChromatinAssay(
  counts = atac_counts, sep = c(":", "-"),
  fragments = if (has_fragments) FRAGMENTS_FILE else NULL,
  annotation = annotation, min.cells = 0, min.features = 0,
  validate.fragments = has_fragments
)

input_summary <- data.frame(
  metric = c("cells", "RNA_features", "ATAC_peak_features",
             "fragment_file_available", "fragment_index_available"),
  value = c(ncol(pbmc), nrow(pbmc[["RNA"]]), nrow(pbmc[["ATAC"]]),
            file.exists(FRAGMENTS_FILE), file.exists(FRAGMENTS_INDEX))
)
write_table(input_summary, file.path(VALIDATION_DIR, "input_summary.csv"))
save_checkpoint(pbmc, "01_raw_multiome.rds")
