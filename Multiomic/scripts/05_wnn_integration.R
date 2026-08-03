if (!exists("PROJECT_ROOT")) source("scripts/00_setup.R")
if (!exists("pbmc")) pbmc <- load_checkpoint("04_atac_processed.rds")

pbmc <- Seurat::FindMultiModalNeighbors(
  object = pbmc,
  reduction.list = list("rna.pca", "atac.lsi"),
  dims.list = list(RNA_DIMS_WNN, ATAC_DIMS_WNN),
  knn.graph.name = "wknn", snn.graph.name = "wsnn",
  weighted.nn.name = "weighted.nn", modality.weight.name = "RNA.weight",
  verbose = TRUE
)
pbmc <- Seurat::RunUMAP(
  object = pbmc, nn.name = "weighted.nn",
  reduction.name = "wnn.umap", reduction.key = "wnnUMAP_",
  seed.use = SEED, verbose = TRUE
)

wnn_plot <- Seurat::DimPlot(pbmc, reduction = "wnn.umap", pt.size = 1) + ggplot2::theme_classic()
save_plot(wnn_plot, file.path(INTEGRATION_FIG_DIR, "wnn_umap_unclustered.png"), 8, 6)
save_checkpoint(pbmc, "05_wnn_integrated.rds")
