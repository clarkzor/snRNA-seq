if (!exists("PROJECT_ROOT")) source("scripts/00_setup.R")
if (!exists("pbmc")) pbmc <- load_checkpoint("05_wnn_integrated.rds")

for (resolution in WNN_LEIDEN_RESOLUTIONS) {
  cluster_column <- paste0("wsnn_leiden_res_", resolution_slug(resolution))
  pbmc <- Seurat::FindClusters(
    object = pbmc, graph.name = "wsnn", resolution = resolution,
    algorithm = 4, leiden_method = "igraph",
    leiden_objective_function = "modularity",
    cluster.name = cluster_column, random.seed = SEED, verbose = FALSE
  )
}

FINAL_CLUSTER_COLUMN <- paste0("wsnn_leiden_res_", resolution_slug(FINAL_WNN_RESOLUTION))
if (!FINAL_CLUSTER_COLUMN %in% colnames(pbmc[[]])) stop("Missing final cluster column: ", FINAL_CLUSTER_COLUMN)
Idents(pbmc) <- pbmc[[FINAL_CLUSTER_COLUMN, drop = TRUE]]

cluster_plot <- Seurat::DimPlot(
  pbmc, reduction = "wnn.umap", group.by = FINAL_CLUSTER_COLUMN,
  label = TRUE, repel = TRUE, pt.size = 1
) + ggplot2::theme_classic()
save_plot(
  cluster_plot,
  file.path(INTEGRATION_FIG_DIR, paste0(FINAL_CLUSTER_COLUMN, "_wnn_umap.png")),
  9, 7
)

cluster_counts <- as.data.frame(table(pbmc[[FINAL_CLUSTER_COLUMN, drop = TRUE]]))
colnames(cluster_counts) <- c("cluster", "n_cells")
write_table(cluster_counts, file.path(INTEGRATION_TABLE_DIR, "final_cluster_counts.csv"))
save_checkpoint(pbmc, "06_wnn_clustered.rds")
