if (!exists("PROJECT_ROOT")) source("scripts/00_setup.R")
if (!exists("pbmc")) pbmc <- load_checkpoint("02_qc_filtered.rds")

DefaultAssay(pbmc) <- "RNA"
pbmc <- Seurat::NormalizeData(
  pbmc, assay = "RNA", normalization.method = "LogNormalize",
  scale.factor = 10000, verbose = FALSE
)
pbmc <- Seurat::FindVariableFeatures(
  pbmc, assay = "RNA", selection.method = "vst",
  nfeatures = RNA_VARIABLE_FEATURES, verbose = FALSE
)

if (USE_SCTRANSFORM) {
  pbmc <- Seurat::SCTransform(
    pbmc, assay = "RNA", new.assay.name = "SCT",
    vst.flavor = "v2", verbose = FALSE
  )
  rna_reduction_assay <- "SCT"
} else {
  pbmc <- Seurat::ScaleData(
    pbmc, assay = "RNA", features = Seurat::VariableFeatures(pbmc[["RNA"]]),
    verbose = FALSE
  )
  rna_reduction_assay <- "RNA"
}

pbmc <- Seurat::RunPCA(
  pbmc, assay = rna_reduction_assay, npcs = RNA_PCS,
  reduction.name = "rna.pca", reduction.key = "rnaPC_", verbose = FALSE
)

if (RUN_RNA_ONLY_CLUSTERING) {
  pbmc <- Seurat::FindNeighbors(
    pbmc, reduction = "rna.pca", dims = RNA_DIMS_WNN,
    graph.name = c("rna_nn", "rna_snn"), verbose = FALSE
  )
  pbmc <- Seurat::FindClusters(
    pbmc, graph.name = "rna_snn", resolution = RNA_LOUVAIN_RESOLUTION,
    algorithm = 1, cluster.name = "rna_louvain_res_13",
    random.seed = SEED, verbose = FALSE
  )
}

pca_plot <- Seurat::DimPlot(pbmc, reduction = "rna.pca", group.by = "orig.ident") + ggplot2::theme_classic()
save_plot(pca_plot, file.path(RNA_FIG_DIR, "rna_pca.png"), 8, 6)
save_checkpoint(pbmc, "03_rna_processed.rds")
