# Pipeline

## Core workflow

1. Read paired Cell Ranger ARC Gene Expression and Peaks matrices.
2. Create a Seurat object with RNA and ATAC assays and Xenopus gene annotations.
3. Calculate RNA/ATAC QC metrics and apply configured thresholds.
4. Normalize RNA, identify variable features, optionally SCTransform, and run RNA PCA.
5. Run ATAC feature selection, TF-IDF, and LSI.
6. Construct weighted nearest neighbors (WNN) from RNA PCA and ATAC LSI.
7. Compute WNN UMAP and Leiden clusters on the `wsnn` graph.
8. Apply prioritized cell-level RNA expression annotation rules.
9. Validate annotations with marker UMAPs, dot plots, conflicts, and cluster cross-tabs.
10. Optionally generate coverage plots and peak-to-gene links.
11. Save a final Seurat RDS, metadata, package versions, and session information.

## Checkpoints

The scripts save numbered RDS checkpoints under `results/objects/`. These are ignored by Git.
Each numbered script can reload the preceding checkpoint when run independently.
