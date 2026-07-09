# Analysis Plan

## Project Title

Single-Cell RNA-seq Preprocessing and Cell Type Annotation of Early Gastrula *Xenopus tropicalis* Embryos

---

# Objective

The objective of this analysis is to process raw single-cell RNA sequencing data from early gastrula stage *Xenopus tropicalis* embryos (stages 10.5, 11.25, 12.0, 12.5, 13.0) to identify high-quality cells, characterize transcriptionally distinct cell populations, assign biological cell type identities, and identify marker genes associated with each population.

---

# Biological Question

What transcriptionally distinct cell populations are present during the early gastrula stage of *Xenopus tropicalis* embryonic development, and what genes define each cell population?

---

# Input Data

Input data consist of reformatted PARSE output matrices containing:

* Gene expression count matrix
* Cell barcode information

Associated sample metadata are incorporated during preprocessing in R. The packages "Matrix", "Seurat" and "SeuratDisk" are used to combine the raw and unfiltered "count_matrix.mtx", "all_genes.csv" and "cell_metadata.csv" into a SeuratObject, saved as a .h5Seurat file and converted into an .h5ad file. (R code can be found in /notebooks/PARSEconversion.R)

---

# Expected Outputs

The pipeline generates:

* Filtered AnnData (.h5ad) object
* Quality control metrics
* PCA embeddings
* UMAP visualization
* Leiden clusters
* Cell type annotations
* Marker gene tables
* Differential expression results

---

# Software

Primary analysis software:

* Python
* Scanpy
* AnnData
* NumPy
* Pandas
* SciPy
* Matplotlib
* Seaborn
* Scrublet

---

# Analysis Workflow

1. Import Cell Ranger count matrices.
2. Create AnnData object.
3. Calculate quality control metrics.
4. Remove low-quality cells.
5. Detect and remove predicted doublets.
6. Normalize gene expression.
7. Log-transform normalized counts.
8. Identify highly variable genes.
9. Scale expression matrix.
10. Perform Principal Component Analysis (PCA).
11. Construct nearest-neighbor graph.
12. Generate UMAP embedding.
13. Perform Leiden graph-based clustering.
14. Identify cluster marker genes.
15. Assign biological cell identities using known marker genes.
16. Perform differential expression analyses between selected cell populations.
17. Export processed datasets, figures, and marker gene tables.

---

# Quality Control Strategy

Quality control is performed prior to downstream analyses to minimize technical artifacts.

QC includes assessment of:

* Number of genes detected per cell
* Total UMI counts
* Percentage of mitochondrial transcripts
* Detection of likely doublets using Scrublet

Cells failing quality thresholds are excluded from downstream analyses.

---

# Normalization

Gene expression counts are normalized to account for differences in sequencing depth between cells.

Normalized counts are log-transformed prior to downstream analyses.

---

# Feature Selection

Highly variable genes are identified and retained for dimensionality reduction and clustering.

Low-information genes are excluded from downstream analyses.

---

# Dimensionality Reduction

Principal Component Analysis (PCA) is performed using highly variable genes.

PCA coordinates are used to construct the neighborhood graph for downstream analyses.

UMAP is used to generate a two-dimensional visualization of transcriptional similarity between cells.

---

# Clustering

Graph-based clustering is performed using the Leiden community detection algorithm.

Clusters represent transcriptionally distinct cellular populations.

---

# Cell Type Annotation

Clusters are manually annotated through comparison of cluster-specific marker genes with established developmental marker genes reported in the literature.

---

# Differential Gene Expression

Differential expression analyses are performed using Scanpy's `rank_genes_groups()` function to identify genes enriched within individual cell populations and between selected biological groups.

Results include:

* Ranked marker genes
* Log fold changes
* P-values
* Adjusted p-values

---

# Statistical Methods

The analysis includes:

* Highly variable gene selection based on gene dispersion
* Principal Component Analysis (PCA)
* k-nearest neighbor graph construction
* Leiden graph clustering
* UMAP dimensionality reduction
* Differential expression analysis using Student's t-test (`rank_genes_groups`)
* Multiple-testing correction performed by Scanpy
* Doublet detection using Scrublet

---

# Assumptions

* PARSE sequenced count matrices accurately represent transcript abundance.
* Quality control thresholds remove technical artifacts while preserving biological variability.
* Highly variable genes adequately capture biological signal.
* Leiden clustering identifies biologically meaningful cellular populations.
* Published marker genes accurately define embryonic cell identities.

---

# Limitations

* Cell annotations depend on currently available marker genes.
* Clustering results may vary with parameter selection.
* Differential expression assumes sufficient cell numbers within each cluster.
* Lowly expressed genes may not be detected reliably.
* UMAP is intended for visualization and should not be interpreted as preserving true biological distances.

---

# Deliverables

Final outputs include:

* Processed AnnData object
* QC summaries
* Cluster annotations
* Marker gene tables
* Differential expression tables
* UMAP visualizations
* Heatmaps
* Dot plots
* Documentation describing the complete analysis workflow
