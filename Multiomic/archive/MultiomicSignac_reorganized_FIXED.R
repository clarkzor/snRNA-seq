#!/usr/bin/env Rscript

# ==============================================================================
# Xenopus tropicalis 10x Multiome: organized Signac + Seurat workflow
# ==============================================================================
# This script:
#   1. reads the paired RNA/ATAC Cell Ranger ARC matrix,
#   2. creates one Seurat object with RNA and ATAC assays,
#   3. performs QC before dimensional reduction and clustering,
#   4. builds RNA PCA, ATAC LSI, WNN neighbors, WNN UMAP, and Leiden clusters,
#   5. assigns auditable cell-level annotations from normalized RNA expression,
#   6. saves plots, metadata, rule summaries, and a final RDS object.
#
# Important:
# - Package installation is intentionally NOT performed inside this analysis.
# - Annotation thresholds are starting points copied/adapted from the Stage 11
#   Python notebook. Review marker plots and edit the CSV rule file as needed.
# - The original manual cluster map is preserved in config/ for comparison only.
#   Do not apply it blindly to newly recalculated clusters.
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(1234)

# ==============================================================================
# 0. User configuration
# ==============================================================================

PROJECT_DIR <- Sys.getenv(
  "MULTIOME_PROJECT_DIR",
  unset = "C:/Users/coron/OneDrive/Desktop/10xGenomics"
)

INPUT_DIR <- PROJECT_DIR
OUTPUT_DIR <- file.path(PROJECT_DIR, "results", "St11_multiome")

FEATURE_MATRIX_H5 <- file.path(INPUT_DIR, "filtered_feature_bc_matrix.h5")
FRAGMENTS_FILE <- file.path(INPUT_DIR, "atac_fragments.tsv.gz")
FRAGMENTS_INDEX <- paste0(FRAGMENTS_FILE, ".tbi")

ENSDB_PACKAGE <- "EnsDb.Xtropicalis.v111"

# RNA and ATAC QC thresholds from the original script.
MIN_RNA_COUNTS <- 1000
MAX_RNA_COUNTS <- 25000
MIN_ATAC_COUNTS <- 1000
MAX_ATAC_COUNTS <- 100000
MAX_NUCLEOSOME_SIGNAL <- 2
MIN_TSS_ENRICHMENT <- 1

USE_SCTRANSFORM <- TRUE
RUN_FRAGMENT_QC <- TRUE
RUN_RNA_ONLY_CLUSTERING <- TRUE

RNA_PCS <- 50
ATAC_LSI_COMPONENTS <- 50
RNA_DIMS_WNN <- 1:30
ATAC_DIMS_WNN <- 2:40

# Several resolutions are retained as separate metadata columns.
WNN_LEIDEN_RESOLUTIONS <- c(0.5, 1.0, 2.0, 5.0, 13.0)
FINAL_WNN_RESOLUTION <- 1.0

# Optional propagation of expression labels to currently unresolved cells using
# consensus within each WNN Leiden cluster. The direct expression seeds are
# always preserved in expression_annotation_seed.
USE_CLUSTER_CONSENSUS_FILL <- FALSE
CLUSTER_FILL_MIN_LABELED_CELLS <- 10
CLUSTER_FILL_MIN_FRACTION <- 0.70

# ==============================================================================
# 1. Locate script resources and create output folders
# ==============================================================================

get_script_directory <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    source_path <- tryCatch(
      rstudioapi::getSourceEditorContext()$path,
      error = function(e) ""
    )
    if (nzchar(source_path)) {
      return(dirname(normalizePath(source_path)))
    }
  }

  normalizePath(getwd())
}

SCRIPT_DIR <- get_script_directory()
CONFIG_DIR <- file.path(SCRIPT_DIR, "config")
RULES_FILE <- file.path(CONFIG_DIR, "stage11_expression_annotation_rules.csv")
MANUAL_MAP_FILE <- file.path(CONFIG_DIR, "original_manual_cluster_map.csv")

FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
OBJECT_DIR <- file.path(OUTPUT_DIR, "objects")
QC_DIR <- file.path(OUTPUT_DIR, "qc")
VALIDATION_DIR <- file.path(OUTPUT_DIR, "validation")

for (folder in c(
  OUTPUT_DIR, FIGURE_DIR, TABLE_DIR, OBJECT_DIR, QC_DIR, VALIDATION_DIR
)) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
}

# ==============================================================================
# 2. Package checks
# ==============================================================================

required_packages <- c(
  "Seurat",
  "Signac",
  "SeuratObject",
  "ggplot2",
  "patchwork",
  "Matrix",
  "GenomeInfoDb",
  "GenomicRanges",
  "ensembldb",
  ENSDB_PACKAGE
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script. Package installation is kept ",
    "separate from the analysis for reproducibility."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(SeuratObject)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(GenomeInfoDb)
  library(GenomicRanges)
  library(ensembldb)
})

suppressPackageStartupMessages(
  library(ENSDB_PACKAGE, character.only = TRUE)
)

# ==============================================================================
# 3. General helper functions
# ==============================================================================

assert_file <- function(path, description) {
  if (!file.exists(path)) {
    stop(description, " not found:\n", path)
  }
}

save_plot <- function(plot_object, filename, width, height) {
  ggplot2::ggsave(
    filename = file.path(FIGURE_DIR, filename),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
}

resolution_slug <- function(x) {
  gsub("\\.", "_", format(x, scientific = FALSE, trim = TRUE))
}

split_markers <- function(x) {
  if (is.na(x) || !nzchar(trimws(x))) {
    return(character(0))
  }
  trimws(strsplit(x, ";", fixed = TRUE)[[1]])
}

append_token <- function(existing, new_value) {
  ifelse(
    is.na(existing) | existing == "",
    new_value,
    paste(existing, new_value, sep = ";")
  )
}

count_unique_tokens <- function(x) {
  if (is.na(x) || x == "") {
    return(0L)
  }
  tokens <- strsplit(x, ";", fixed = TRUE)[[1]]
  length(unique(tokens[tokens != ""]))
}

harmonize_annotation_seqlevels <- function(annotation, peak_names) {
  peak_chromosomes <- unique(sub(":.*$", "", peak_names))
  annotation_chromosomes <- GenomeInfoDb::seqlevels(annotation)

  common <- intersect(annotation_chromosomes, peak_chromosomes)

  if (length(common) == 0) {
    stripped <- sub("^chr", "", annotation_chromosomes)
    if (length(intersect(stripped, peak_chromosomes)) > 0) {
      mapping <- stats::setNames(stripped, annotation_chromosomes)
      annotation <- GenomeInfoDb::renameSeqlevels(annotation, mapping)
      annotation_chromosomes <- GenomeInfoDb::seqlevels(annotation)
      common <- intersect(annotation_chromosomes, peak_chromosomes)
    }
  }

  if (length(common) == 0) {
    stop(
      "No chromosome names overlap between the EnsDb annotation and peak ",
      "matrix. Peak examples: ",
      paste(head(peak_chromosomes), collapse = ", "),
      "; annotation examples: ",
      paste(head(annotation_chromosomes), collapse = ", ")
    )
  }

  GenomeInfoDb::keepSeqlevels(
    annotation,
    value = common,
    pruning.mode = "coarse"
  )
}

# ==============================================================================
# 4. Read Cell Ranger ARC matrix and create the multimodal Seurat object
# ==============================================================================

assert_file(FEATURE_MATRIX_H5, "Filtered feature-barcode H5 matrix")
assert_file(RULES_FILE, "Expression-annotation rule table")

counts <- Seurat::Read10X_h5(FEATURE_MATRIX_H5)

if (!is.list(counts)) {
  stop("Read10X_h5 did not return a modality list.")
}

required_modalities <- c("Gene Expression", "Peaks")
missing_modalities <- setdiff(required_modalities, names(counts))

if (length(missing_modalities) > 0) {
  stop(
    "The H5 file is missing modalities: ",
    paste(missing_modalities, collapse = ", ")
  )
}

rna_counts <- counts[["Gene Expression"]]
atac_counts <- counts[["Peaks"]]

common_cells <- intersect(colnames(rna_counts), colnames(atac_counts))

if (length(common_cells) == 0) {
  stop("RNA and ATAC matrices do not share any cell barcodes.")
}

rna_counts <- rna_counts[, common_cells, drop = FALSE]
atac_counts <- atac_counts[, common_cells, drop = FALSE]

ensdb_object <- get(ENSDB_PACKAGE)
annotation <- Signac::GetGRangesFromEnsDb(ensdb = ensdb_object)
annotation <- harmonize_annotation_seqlevels(annotation, rownames(atac_counts))

has_fragments <- file.exists(FRAGMENTS_FILE) && file.exists(FRAGMENTS_INDEX)

if (!has_fragments) {
  warning(
    "The fragment file and/or its tabix index were not found. RNA-based ",
    "annotation and matrix-based RNA/ATAC processing can run, but fragment-",
    "dependent QC, CoveragePlot, TilePlot, and LinkPeaks cannot run.\n",
    "Expected files:\n", FRAGMENTS_FILE, "\n", FRAGMENTS_INDEX
  )
}

pbmc <- Seurat::CreateSeuratObject(
  counts = rna_counts,
  assay = "RNA",
  project = "Xtrop_St11_Multiome",
  min.cells = 0,
  min.features = 0
)

pbmc[["ATAC"]] <- Signac::CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  fragments = if (has_fragments) FRAGMENTS_FILE else NULL,
  annotation = annotation,
  min.cells = 0,
  min.features = 0,
  validate.fragments = has_fragments
)

input_summary <- data.frame(
  metric = c(
    "cells",
    "RNA_features",
    "ATAC_peak_features",
    "fragment_file_available",
    "fragment_index_available"
  ),
  value = c(
    ncol(pbmc),
    nrow(pbmc[["RNA"]]),
    nrow(pbmc[["ATAC"]]),
    file.exists(FRAGMENTS_FILE),
    file.exists(FRAGMENTS_INDEX)
  )
)

write.csv(
  input_summary,
  file.path(VALIDATION_DIR, "input_summary.csv"),
  row.names = FALSE
)

saveRDS(
  pbmc,
  file.path(OBJECT_DIR, "01_multimodal_object_unfiltered.rds")
)

# ==============================================================================
# 5. QC before clustering
# ==============================================================================

DefaultAssay(pbmc) <- "ATAC"

if (RUN_FRAGMENT_QC && has_fragments) {
  # These functions are deprecated in recent Signac releases but remain useful
  # on native Windows when fragtk/ATACqc is unavailable.
  pbmc <- suppressWarnings(Signac::NucleosomeSignal(pbmc))
  pbmc <- suppressWarnings(Signac::TSSEnrichment(pbmc, fast = FALSE))
}

qc_features <- c(
  "nCount_RNA",
  "nFeature_RNA",
  "nCount_ATAC",
  "nFeature_ATAC"
)

if ("nucleosome_signal" %in% colnames(pbmc[[]])) {
  qc_features <- c(qc_features, "nucleosome_signal")
}
if ("TSS.enrichment" %in% colnames(pbmc[[]])) {
  qc_features <- c(qc_features, "TSS.enrichment")
}

qc_plots <- Seurat::VlnPlot(
  object = pbmc,
  features = qc_features,
  pt.size = 0.1,
  combine = FALSE,
  ncol = 3
)

qc_plots <- lapply(
  qc_plots,
  function(plot_object) {
    plot_object +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
  }
)

qc_panel <- patchwork::wrap_plots(qc_plots, ncol = 3)

ggplot2::ggsave(
  filename = file.path(QC_DIR, "pre_filter_qc_violin.png"),
  plot = qc_panel,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

metadata_before_qc <- pbmc[[]]
metadata_before_qc$cell_barcode <- rownames(metadata_before_qc)
write.csv(
  metadata_before_qc,
  file.path(QC_DIR, "cell_metadata_before_qc.csv"),
  row.names = FALSE
)

keep_cells <- (
  pbmc$nCount_RNA > MIN_RNA_COUNTS &
  pbmc$nCount_RNA < MAX_RNA_COUNTS &
  pbmc$nCount_ATAC > MIN_ATAC_COUNTS &
  pbmc$nCount_ATAC < MAX_ATAC_COUNTS
)

if ("nucleosome_signal" %in% colnames(pbmc[[]])) {
  keep_cells <- keep_cells &
    pbmc$nucleosome_signal < MAX_NUCLEOSOME_SIGNAL
}

if ("TSS.enrichment" %in% colnames(pbmc[[]])) {
  keep_cells <- keep_cells &
    pbmc$TSS.enrichment > MIN_TSS_ENRICHMENT
}

qc_decisions <- data.frame(
  cell_barcode = colnames(pbmc),
  pass_qc = as.logical(keep_cells)
)
write.csv(
  qc_decisions,
  file.path(QC_DIR, "cell_qc_decisions.csv"),
  row.names = FALSE
)

pbmc <- subset(pbmc, cells = colnames(pbmc)[keep_cells])

metadata_after_qc <- pbmc[[]]
metadata_after_qc$cell_barcode <- rownames(metadata_after_qc)
write.csv(
  metadata_after_qc,
  file.path(QC_DIR, "cell_metadata_after_qc.csv"),
  row.names = FALSE
)

saveRDS(
  pbmc,
  file.path(OBJECT_DIR, "02_multimodal_object_qc_filtered.rds")
)

# ==============================================================================
# 6. RNA processing
# ==============================================================================

DefaultAssay(pbmc) <- "RNA"

pbmc <- Seurat::NormalizeData(
  pbmc,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

pbmc <- Seurat::FindVariableFeatures(
  pbmc,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

if (USE_SCTRANSFORM) {
  pbmc <- Seurat::SCTransform(
    pbmc,
    assay = "RNA",
    new.assay.name = "SCT",
    vst.flavor = "v2",
    verbose = FALSE
  )

  rna_reduction_assay <- "SCT"
} else {
  pbmc <- Seurat::ScaleData(
    pbmc,
    assay = "RNA",
    features = Seurat::VariableFeatures(pbmc[["RNA"]]),
    verbose = FALSE
  )
  rna_reduction_assay <- "RNA"
}

pbmc <- Seurat::RunPCA(
  pbmc,
  assay = rna_reduction_assay,
  npcs = RNA_PCS,
  reduction.name = "rna.pca",
  reduction.key = "rnaPC_",
  verbose = FALSE
)

if (RUN_RNA_ONLY_CLUSTERING) {
  pbmc <- Seurat::FindNeighbors(
    pbmc,
    reduction = "rna.pca",
    dims = RNA_DIMS_WNN,
    graph.name = c("rna_nn", "rna_snn"),
    verbose = FALSE
  )

  # This recreates the type of high-resolution RNA-only clustering used in the
  # original script. algorithm = 1 is Louvain, not Leiden.
  pbmc <- Seurat::FindClusters(
    pbmc,
    graph.name = "rna_snn",
    resolution = 13,
    algorithm = 1,
    cluster.name = "rna_louvain_res_13",
    random.seed = 1234,
    verbose = FALSE
  )
}

# ==============================================================================
# 7. ATAC processing using the Cell Ranger ARC peak matrix already in the H5
# ==============================================================================

DefaultAssay(pbmc) <- "ATAC"

pbmc <- Signac::FindTopFeatures(
  pbmc,
  min.cutoff = 5
)

pbmc <- Signac::RunTFIDF(pbmc)

pbmc <- Signac::RunSVD(
  pbmc,
  n = ATAC_LSI_COMPONENTS,
  reduction.name = "atac.lsi",
  reduction.key = "atacLSI_"
)

# ==============================================================================
# 8. WNN integration, WNN UMAP, and true Leiden clustering
# ==============================================================================

pbmc <- Seurat::FindMultiModalNeighbors(
  object = pbmc,
  reduction.list = list("rna.pca", "atac.lsi"),
  dims.list = list(RNA_DIMS_WNN, ATAC_DIMS_WNN),
  knn.graph.name = "wknn",
  snn.graph.name = "wsnn",
  weighted.nn.name = "weighted.nn",
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

pbmc <- Seurat::RunUMAP(
  object = pbmc,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_",
  seed.use = 1234,
  verbose = TRUE
)

for (resolution in WNN_LEIDEN_RESOLUTIONS) {
  cluster_column <- paste0(
    "wsnn_leiden_res_",
    resolution_slug(resolution)
  )

  pbmc <- Seurat::FindClusters(
    object = pbmc,
    graph.name = "wsnn",
    resolution = resolution,
    algorithm = 4,
    leiden_method = "igraph",
    leiden_objective_function = "modularity",
    cluster.name = cluster_column,
    random.seed = 1234,
    verbose = FALSE
  )
}

FINAL_CLUSTER_COLUMN <- paste0(
  "wsnn_leiden_res_",
  resolution_slug(FINAL_WNN_RESOLUTION)
)

if (!FINAL_CLUSTER_COLUMN %in% colnames(pbmc[[]])) {
  stop("Final cluster column was not generated: ", FINAL_CLUSTER_COLUMN)
}

Idents(pbmc) <- pbmc[[FINAL_CLUSTER_COLUMN, drop = TRUE]]

cluster_plot <- Seurat::DimPlot(
  pbmc,
  reduction = "wnn.umap",
  group.by = FINAL_CLUSTER_COLUMN,
  label = TRUE,
  repel = TRUE,
  pt.size = 1
) + ggplot2::theme_classic()

save_plot(
  cluster_plot,
  paste0(FINAL_CLUSTER_COLUMN, "_wnn_umap.png"),
  width = 9,
  height = 7
)

# ==============================================================================
# 9. Expression-based annotation rule engine
# ==============================================================================

apply_expression_rules <- function(
    object,
    rules_file,
    assay = "RNA",
    layer = "data",
    unassigned_label = "Other") {

  rules <- read.csv(
    rules_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  rules <- rules[order(rules$priority), , drop = FALSE]

  all_rule_genes <- unique(unlist(lapply(
    seq_len(nrow(rules)),
    function(i) {
      c(
        split_markers(rules$positive_all[i]),
        split_markers(rules$positive_any[i]),
        split_markers(rules$negative_all[i])
      )
    }
  )))

  available_genes <- intersect(all_rule_genes, rownames(object[[assay]]))
  missing_genes <- setdiff(all_rule_genes, available_genes)

  marker_inventory <- data.frame(
    gene = all_rule_genes,
    present = all_rule_genes %in% available_genes
  )

  if (length(available_genes) == 0) {
    stop("None of the annotation marker genes are present in the RNA assay.")
  }

  # Pull expression directly from the requested assay/layer.
  # This avoids ambiguity when the same genes are present in both RNA and SCT.
  expression_matrix <- SeuratObject::LayerData(
    object = object[[assay]],
    layer = layer,
    features = available_genes
  )

  # LayerData returns genes x cells; annotation rules use cells x genes.
  expression <- as.data.frame(
    as.matrix(t(expression_matrix)),
    check.names = FALSE
  )

  # Force exact Seurat cell order.
  expression <- expression[colnames(object), , drop = FALSE]

  labels <- rep(unassigned_label, nrow(expression))
  names(labels) <- rownames(expression)

  assignment_rule <- rep(NA_character_, nrow(expression))
  names(assignment_rule) <- rownames(expression)

  candidate_rules <- rep("", nrow(expression))
  candidate_labels <- rep("", nrow(expression))

  rule_summaries <- vector("list", nrow(rules))

  for (i in seq_len(nrow(rules))) {
    positive_all <- split_markers(rules$positive_all[i])
    positive_any <- split_markers(rules$positive_any[i])
    negative_all <- split_markers(rules$negative_all[i])

    required_genes <- unique(c(positive_all, negative_all))
    missing_required <- setdiff(required_genes, available_genes)
    available_any <- intersect(positive_any, available_genes)

    if (length(missing_required) > 0 ||
        (length(positive_any) > 0 && length(available_any) == 0)) {
      rule_summaries[[i]] <- data.frame(
        priority = rules$priority[i],
        rule_name = rules$rule_name[i],
        label = rules$label[i],
        status = "skipped_missing_markers",
        missing_markers = paste(
          unique(c(missing_required, setdiff(positive_any, available_any))),
          collapse = ";"
        ),
        matched_cells = 0,
        newly_assigned_cells = 0
      )
      next
    }

    positive_threshold <- as.numeric(rules$positive_threshold[i])
    negative_threshold <- as.numeric(rules$negative_threshold[i])

    mask <- rep(TRUE, nrow(expression))

    if (length(positive_all) > 0) {
      positive_all_matrix <- as.matrix(
        expression[, positive_all, drop = FALSE]
      )
      mask <- mask &
        rowSums(positive_all_matrix > positive_threshold) ==
        length(positive_all)
    }

    if (length(available_any) > 0) {
      positive_any_matrix <- as.matrix(
        expression[, available_any, drop = FALSE]
      )
      mask <- mask &
        rowSums(positive_any_matrix > positive_threshold) >= 1
    }

    if (length(negative_all) > 0) {
      negative_matrix <- as.matrix(
        expression[, negative_all, drop = FALSE]
      )
      mask <- mask &
        rowSums(negative_matrix > negative_threshold) == 0
    }

    candidate_rules[mask] <- append_token(
      candidate_rules[mask],
      rules$rule_name[i]
    )
    candidate_labels[mask] <- append_token(
      candidate_labels[mask],
      rules$label[i]
    )

    newly_assigned <- mask & labels == unassigned_label
    labels[newly_assigned] <- rules$label[i]
    assignment_rule[newly_assigned] <- rules$rule_name[i]

    rule_summaries[[i]] <- data.frame(
      priority = rules$priority[i],
      rule_name = rules$rule_name[i],
      label = rules$label[i],
      status = "evaluated",
      missing_markers = "",
      matched_cells = sum(mask),
      newly_assigned_cells = sum(newly_assigned)
    )
  }

  candidate_label_count <- vapply(
    candidate_labels,
    count_unique_tokens,
    integer(1)
  )

  # Add metadata with explicit cell-barcode row names so Seurat aligns rows safely.
  annotation_metadata <- data.frame(
    expression_annotation_seed = labels,
    expression_annotation_rule = assignment_rule,
    expression_annotation_candidates = candidate_labels,
    expression_annotation_candidate_rules = candidate_rules,
    expression_annotation_conflict = candidate_label_count > 1,
    expression_annotation_n_candidate_labels = candidate_label_count,
    stringsAsFactors = FALSE,
    row.names = colnames(object)
  )

  annotation_metadata$expression_annotation_seed <- factor(
    annotation_metadata$expression_annotation_seed
  )

  object <- SeuratObject::AddMetaData(
    object = object,
    metadata = annotation_metadata
  )

  list(
    object = object,
    rules = rules,
    rule_summary = do.call(rbind, rule_summaries),
    marker_inventory = marker_inventory,
    missing_genes = missing_genes
  )
}

annotation_result <- apply_expression_rules(
  pbmc,
  rules_file = RULES_FILE,
  assay = "RNA",
  layer = "data"
)

pbmc <- annotation_result$object

write.csv(
  annotation_result$rule_summary,
  file.path(TABLE_DIR, "expression_annotation_rule_summary.csv"),
  row.names = FALSE
)

write.csv(
  annotation_result$marker_inventory,
  file.path(VALIDATION_DIR, "expression_annotation_marker_inventory.csv"),
  row.names = FALSE
)

# ==============================================================================
# 10. Optional cluster-consensus fill for cells remaining Other
# ==============================================================================

fill_unassigned_by_cluster <- function(
    object,
    cluster_col,
    seed_col,
    output_col,
    unassigned_label = "Other",
    min_labeled_cells = 10,
    min_fraction = 0.70) {

  clusters <- as.character(object[[cluster_col, drop = TRUE]])
  seed_labels <- as.character(object[[seed_col, drop = TRUE]])
  filled_labels <- seed_labels

  cluster_summaries <- lapply(
    sort(unique(clusters)),
    function(cluster_id) {
      in_cluster <- clusters == cluster_id
      labeled <- seed_labels[in_cluster]
      labeled <- labeled[labeled != unassigned_label]

      if (length(labeled) == 0) {
        return(data.frame(
          cluster = cluster_id,
          labeled_cells = 0,
          consensus_label = NA_character_,
          consensus_fraction = NA_real_,
          filled_cells = 0,
          status = "no_expression_seeds"
        ))
      }

      counts <- sort(table(labeled), decreasing = TRUE)
      consensus_label <- names(counts)[1]
      consensus_fraction <- as.numeric(counts[1]) / sum(counts)

      can_fill <- (
        length(labeled) >= min_labeled_cells &&
        consensus_fraction >= min_fraction
      )

      fill_mask <- (
        in_cluster &
        seed_labels == unassigned_label &
        can_fill
      )

      filled_labels[fill_mask] <<- consensus_label

      data.frame(
        cluster = cluster_id,
        labeled_cells = length(labeled),
        consensus_label = consensus_label,
        consensus_fraction = consensus_fraction,
        filled_cells = sum(fill_mask),
        status = if (can_fill) "eligible" else "insufficient_consensus"
      )
    }
  )

  object[[output_col]] <- factor(filled_labels)

  list(
    object = object,
    summary = do.call(rbind, cluster_summaries)
  )
}

cluster_fill_result <- fill_unassigned_by_cluster(
  pbmc,
  cluster_col = FINAL_CLUSTER_COLUMN,
  seed_col = "expression_annotation_seed",
  output_col = "expression_annotation_cluster_filled",
  min_labeled_cells = CLUSTER_FILL_MIN_LABELED_CELLS,
  min_fraction = CLUSTER_FILL_MIN_FRACTION
)

pbmc <- cluster_fill_result$object

write.csv(
  cluster_fill_result$summary,
  file.path(TABLE_DIR, "cluster_consensus_fill_summary.csv"),
  row.names = FALSE
)

FINAL_ANNOTATION_COLUMN <- if (USE_CLUSTER_CONSENSUS_FILL) {
  "expression_annotation_cluster_filled"
} else {
  "expression_annotation_seed"
}

pbmc$final_annotation <- pbmc[[FINAL_ANNOTATION_COLUMN, drop = TRUE]]
Idents(pbmc) <- "final_annotation"

# ==============================================================================
# 11. Annotation validation plots and tables
# ==============================================================================

annotation_plot <- Seurat::DimPlot(
  pbmc,
  reduction = "wnn.umap",
  group.by = "final_annotation",
  label = TRUE,
  repel = TRUE,
  pt.size = 1
) + ggplot2::theme_classic()

save_plot(
  annotation_plot,
  "final_expression_annotation_wnn_umap.png",
  width = 11,
  height = 8
)

marker_panel <- c(
  "sox2", "sox3", "pax3",
  "krt7", "grhl1", "tp63",
  "tbxt", "fgf8", "chrd.1", "wnt8a",
  "darmin", "gata6", "slc5a8.1",
  "agr2", "ag1",
  "mcidas", "foxj1", "dnah5",
  "tp73", "dscaml1"
)

available_marker_panel <- intersect(marker_panel, rownames(pbmc[["RNA"]]))

# FeaturePlot does not use an assay argument reliably when genes are duplicated
# across RNA/SCT assays. Explicitly make RNA the active assay first.
DefaultAssay(pbmc) <- "RNA"

feature_plots <- Seurat::FeaturePlot(
  pbmc,
  features = available_marker_panel,
  reduction = "wnn.umap",
  order = TRUE,
  pt.size = 0.8,
  combine = FALSE
)

feature_plots <- lapply(
  feature_plots,
  function(plot_object) {
    plot_object +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.position = "right")
  }
)

feature_panel <- patchwork::wrap_plots(feature_plots, ncol = 4)

save_plot(
  feature_panel,
  "annotation_marker_featureplot_panel.png",
  width = 16,
  height = 20
)

dot_plot <- Seurat::DotPlot(
  pbmc,
  features = available_marker_panel,
  assay = "RNA",
  group.by = "final_annotation"
) +
  Seurat::RotatedAxis() +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 9)
  )

save_plot(
  dot_plot,
  "annotation_marker_dotplot.png",
  width = 14,
  height = 8
)

annotation_counts <- as.data.frame(table(pbmc$final_annotation))
colnames(annotation_counts) <- c("final_annotation", "n_cells")
annotation_counts$percentage <- 100 * annotation_counts$n_cells / ncol(pbmc)

write.csv(
  annotation_counts,
  file.path(TABLE_DIR, "final_annotation_counts.csv"),
  row.names = FALSE
)

cluster_annotation_table <- as.data.frame.matrix(table(
  pbmc[[FINAL_CLUSTER_COLUMN, drop = TRUE]],
  pbmc$final_annotation
))
cluster_annotation_table[[FINAL_CLUSTER_COLUMN]] <- rownames(cluster_annotation_table)
cluster_annotation_table <- cluster_annotation_table[
  c(FINAL_CLUSTER_COLUMN, setdiff(colnames(cluster_annotation_table), FINAL_CLUSTER_COLUMN))
]

write.csv(
  cluster_annotation_table,
  file.path(TABLE_DIR, "wnn_leiden_by_final_annotation_counts.csv"),
  row.names = FALSE
)

annotation_metadata <- pbmc[[]]
annotation_metadata$cell_barcode <- rownames(annotation_metadata)
write.csv(
  annotation_metadata,
  file.path(TABLE_DIR, "final_cell_metadata.csv"),
  row.names = FALSE
)

# Preserve the original manual mapping file as a reference, but do not apply it
# automatically because QC/reclustering can change numeric cluster identities.
if (file.exists(MANUAL_MAP_FILE)) {
  file.copy(
    MANUAL_MAP_FILE,
    file.path(VALIDATION_DIR, "original_manual_cluster_map_reference.csv"),
    overwrite = TRUE
  )
}

# ==============================================================================
# 12. Optional fragment-dependent downstream analysis examples
# ==============================================================================

if (has_fragments) {
  # Example coverage plot. This uses the existing Cell Ranger ARC ATAC assay;
  # no second peak assay is needed unless you deliberately call a new peak set.
  coverage_sox2 <- Signac::CoveragePlot(
    object = pbmc,
    region = "sox2",
    features = "sox2",
    assay = "ATAC",
    expression.assay = if ("SCT" %in% Seurat::Assays(pbmc)) "SCT" else "RNA",
    peaks = TRUE,
    extend.upstream = 10000,
    extend.downstream = 10000
  )

  save_plot(
    coverage_sox2,
    "coverage_sox2.png",
    width = 12,
    height = 7
  )

  # LinkPeaks additionally requires RegionStats to have been run with a matching
  # Xenopus BSgenome package. Keep this disabled until the genome package is
  # installed and chromosome names have been checked.
  #
  # library(BSgenome.Xtropicalis.NCBI.UCBXtro10.0)
  # DefaultAssay(pbmc) <- "ATAC"
  # pbmc <- RegionStats(
  #   pbmc,
  #   genome = BSgenome.Xtropicalis.NCBI.UCBXtro10.0
  # )
  # pbmc <- LinkPeaks(
  #   object = pbmc,
  #   peak.assay = "ATAC",
  #   expression.assay = if ("SCT" %in% Seurat::Assays(pbmc)) "SCT" else "RNA",
  #   genes.use = c("sox17a", "tbxt")
  # )
}

# ==============================================================================
# 13. Save the final object and reproducibility information
# ==============================================================================

saveRDS(
  pbmc,
  file.path(OBJECT_DIR, "St11_multiome_expression_annotated.rds")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(VALIDATION_DIR, "sessionInfo.txt")
)

writeLines(
  c(
    paste("Final cells:", ncol(pbmc)),
    paste("Final RNA features:", nrow(pbmc[["RNA"]])),
    paste("Final ATAC peaks:", nrow(pbmc[["ATAC"]])),
    paste("Final cluster column:", FINAL_CLUSTER_COLUMN),
    paste("Final annotation source:", FINAL_ANNOTATION_COLUMN),
    paste("Unresolved cells:", sum(pbmc$final_annotation == "Other")),
    paste("Annotation conflicts:", sum(pbmc$expression_annotation_conflict))
  ),
  file.path(VALIDATION_DIR, "run_summary.txt")
)

message("Pipeline complete.")
message("Final object: ", file.path(OBJECT_DIR, "St11_multiome_expression_annotated.rds"))
message("Results directory: ", OUTPUT_DIR)
