if (!exists("PROJECT_ROOT")) source("scripts/00_setup.R")
if (!exists("pbmc")) pbmc <- load_checkpoint("07_expression_annotated.rds")
FINAL_CLUSTER_COLUMN <- paste0("wsnn_leiden_res_", resolution_slug(FINAL_WNN_RESOLUTION))

annotation_plot <- Seurat::DimPlot(
  pbmc, reduction = "wnn.umap", group.by = "final_annotation",
  label = TRUE, repel = TRUE, pt.size = 1
) + ggplot2::theme_classic()
save_plot(annotation_plot,
          file.path(ANNOTATION_FIG_DIR, "final_expression_annotation_wnn_umap.png"), 11, 8)

conflict_plot <- Seurat::DimPlot(
  pbmc, reduction = "wnn.umap", group.by = "expression_annotation_conflict",
  pt.size = 1
) + ggplot2::theme_classic()
save_plot(conflict_plot,
          file.path(ANNOTATION_FIG_DIR, "annotation_conflicts_wnn_umap.png"), 8, 6)

marker_table <- read.csv(MARKER_PANEL_FILE, stringsAsFactors = FALSE)
marker_panel <- unique(marker_table$gene)
available_marker_panel <- intersect(marker_panel, rownames(pbmc[["RNA"]]))
missing_marker_panel <- setdiff(marker_panel, available_marker_panel)
write_table(
  data.frame(gene = marker_panel, present = marker_panel %in% available_marker_panel),
  file.path(VALIDATION_DIR, "marker_panel_inventory.csv")
)

# FeaturePlot does not take an assay argument; set the default assay explicitly.
DefaultAssay(pbmc) <- "RNA"
feature_plots <- Seurat::FeaturePlot(
  pbmc, features = available_marker_panel, reduction = "wnn.umap",
  order = TRUE, pt.size = 0.8, combine = FALSE
)
feature_plots <- lapply(feature_plots, function(p) {
  p + ggplot2::theme_classic() + ggplot2::theme(legend.position = "right")
})
feature_panel <- patchwork::wrap_plots(feature_plots, ncol = 4)
save_plot(feature_panel,
          file.path(ANNOTATION_FIG_DIR, "annotation_marker_featureplot_panel.png"), 16, 20)

dot_plot <- Seurat::DotPlot(
  pbmc, features = available_marker_panel, group.by = "final_annotation"
) + Seurat::RotatedAxis() + ggplot2::theme_classic() +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9))
save_plot(dot_plot,
          file.path(ANNOTATION_FIG_DIR, "annotation_marker_dotplot.png"), 14, 8)

cluster_annotation <- as.data.frame.matrix(table(
  pbmc[[FINAL_CLUSTER_COLUMN, drop = TRUE]], pbmc$final_annotation
))
cluster_annotation[[FINAL_CLUSTER_COLUMN]] <- rownames(cluster_annotation)
cluster_annotation <- cluster_annotation[
  c(FINAL_CLUSTER_COLUMN, setdiff(colnames(cluster_annotation), FINAL_CLUSTER_COLUMN))
]
write_table(cluster_annotation,
            file.path(ANNOTATION_TABLE_DIR, "wnn_leiden_by_final_annotation_counts.csv"))

if (file.exists(MANUAL_MAP_FILE)) {
  file.copy(MANUAL_MAP_FILE,
            file.path(VALIDATION_DIR, "original_manual_cluster_map_reference.csv"),
            overwrite = TRUE)
}

validation_summary <- data.frame(
  metric = c("cells", "unresolved_other", "annotation_conflicts", "marker_panel_missing"),
  value = c(
    ncol(pbmc),
    sum(as.character(pbmc$final_annotation) == UNASSIGNED_LABEL),
    sum(pbmc$expression_annotation_conflict),
    length(missing_marker_panel)
  )
)
write_table(validation_summary,
            file.path(VALIDATION_DIR, "annotation_validation_summary.csv"))
save_checkpoint(pbmc, "08_annotation_validated.rds")
