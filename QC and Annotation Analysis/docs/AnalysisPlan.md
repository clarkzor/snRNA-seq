# Analysis Plan

## Purpose

This project processes and annotates a developmental single-nucleus RNA-seq
time course from _Xenopus tropicalis_ gastrulation. The completed notebooks are
organized around stage-specific analyses rather than a single integrated
object.

## Analysis units

Seven independent datasets are analyzed:

- Stage 10/10.5 WT
- Stage 11.25 WT
- Stage 11.25 Foxi2 morpholino
- Stage 12.0 WT
- Stage 12.0 Foxi2 morpholino
- Stage 12.5 WT
- Stage 13 WT

WT and morpholino data are processed in separate notebooks at Stage 11.25 and
Stage 12.0. The current notebook series does not perform a formal WT-versus-MO
differential-expression test; Section 5 focuses on annotation-level contrasts
within each dataset.

## Objectives

1. Remove low-complexity nuclei and nuclei with outlying mitochondrial or
   ribosomal count fractions.
2. predict and remove likely doublets.
3. construct PCA, t-SNE, neighbor-graph, and UMAP representations.
4. identify broad and fine cell populations using Leiden clustering and marker
   expression.
5. refine manual annotations using local nearest-neighbor information.
6. quantify annotation abundance and annotation-to-cluster correspondence.
7. identify marker genes for each final annotation and for selected pairwise
   tissue-layer comparisons.
8. save reproducible figures, tables, metadata, and processed AnnData objects.

## Completed notebook sections

### Section 1 — preprocessing and QC

- load an unfiltered stage-specific `.h5ad`;
- retain cells with at least 1,000 detected genes;
- annotate mitochondrial and ribosomal genes;
- compute per-cell QC metrics;
- retain `pct_counts_mt < 20` and `pct_counts_ribo < 3.5`;
- run Scrublet and store doublet scores and predictions;
- normalize to 10,000 counts per cell and apply `log1p`;
- select highly variable genes.

### Section 2 — dimensionality reduction

- calculate 100 principal components with the ARPACK solver;
- inspect PCA projections, variance ratios, and loadings;
- generate exploratory t-SNE embeddings;
- compare multiple neighbor-graph parameterizations;
- generate UMAP embeddings.

### Section 3 — marker survey and doublet validation

- visualize canonical ectoderm, mesoderm, and endoderm marker panels;
- overlay predicted doublets and total counts;
- remove Scrublet-predicted doublets after marker-based review.

### Section 4 — clustering and annotation

- compare Leiden resolutions;
- subset broad germ-layer regions;
- assign labels using marker-expression thresholds;
- preserve or overwrite prior labels according to explicit rules;
- fill remaining `Other` cells using nearby labeled cells;
- smooth selected boundaries while retaining pre-refinement labels;
- transfer subobject labels back to the full AnnData object.

### Section 5 — differential expression

The cleaned notebooks use Wilcoxon tests and save all statistical outputs.

- all sufficiently populated annotations versus the rest;
- reciprocal pairwise tests for selected layer or region pairs;
- complete, significant, and top-ranked result tables;
- group-size and parameter tables;
- ranked-gene plots, heatmaps, dotplots, and matrix plots.

### Section 6 — final naming and outputs

- rename temporary layer labels to final biological names where specified;
- save final annotation summaries and cluster cross-tabulations;
- save the annotated `.h5ad` and complete `obs` metadata table.

## Primary outputs

- one final annotated AnnData object per dataset;
- one full cell-metadata CSV per dataset;
- annotation counts and percentages;
- cluster-to-annotation count and row-percentage tables;
- label-transition tables for KNN refinement or smoothing;
- all-group and pairwise DEG tables;
- section-specific figures.

## Acceptance criteria

A dataset analysis is considered complete when:

- all code cells run from a fresh kernel without manual path edits;
- no final cell is labeled `Other` or `Unknown`, unless explicitly retained and
  documented;
- every annotation contains enough cells for its intended statistical use;
- every `rank_genes_groups` key used for plotting was generated on the same
  AnnData object;
- complete DEG results and test parameters are saved;
- final annotation counts sum to `adata.n_obs`;
- processed data and metadata are written to the configured results directory.
