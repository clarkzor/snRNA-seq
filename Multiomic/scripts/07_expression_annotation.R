if (!exists("PROJECT_ROOT")) source("scripts/00_setup.R")
if (!exists("pbmc")) pbmc <- load_checkpoint("06_wnn_clustered.rds")

FINAL_CLUSTER_COLUMN <- paste0("wsnn_leiden_res_", resolution_slug(FINAL_WNN_RESOLUTION))
annotation_result <- apply_expression_rules(
  pbmc, rules_file = RULES_FILE, assay = "RNA", layer = "data",
  unassigned_label = UNASSIGNED_LABEL
)
pbmc <- annotation_result$object

write_table(annotation_result$rule_summary,
            file.path(ANNOTATION_TABLE_DIR, "expression_annotation_rule_summary.csv"))
write_table(annotation_result$marker_inventory,
            file.path(VALIDATION_DIR, "expression_annotation_marker_inventory.csv"))

cluster_fill_result <- fill_unassigned_by_cluster(
  pbmc, cluster_col = FINAL_CLUSTER_COLUMN,
  seed_col = "expression_annotation_seed",
  output_col = "expression_annotation_cluster_filled",
  unassigned_label = UNASSIGNED_LABEL,
  min_labeled_cells = CLUSTER_FILL_MIN_LABELED_CELLS,
  min_fraction = CLUSTER_FILL_MIN_FRACTION
)
pbmc <- cluster_fill_result$object
write_table(cluster_fill_result$summary,
            file.path(ANNOTATION_TABLE_DIR, "cluster_consensus_fill_summary.csv"))

FINAL_ANNOTATION_COLUMN <- if (USE_CLUSTER_CONSENSUS_FILL) {
  "expression_annotation_cluster_filled"
} else {
  "expression_annotation_seed"
}

final_metadata <- data.frame(
  final_annotation = pbmc[[FINAL_ANNOTATION_COLUMN, drop = TRUE]],
  row.names = colnames(pbmc)
)
pbmc <- SeuratObject::AddMetaData(pbmc, metadata = final_metadata)
Idents(pbmc) <- "final_annotation"

annotation_counts <- as.data.frame(table(pbmc$final_annotation))
colnames(annotation_counts) <- c("final_annotation", "n_cells")
annotation_counts$percentage <- 100 * annotation_counts$n_cells / ncol(pbmc)
write_table(annotation_counts, file.path(ANNOTATION_TABLE_DIR, "final_annotation_counts.csv"))
save_checkpoint(pbmc, "07_expression_annotated.rds")
