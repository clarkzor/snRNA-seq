# Pipeline

## Project

Single-nucleus RNA-seq preprocessing, clustering, and annotation of *Xenopus tropicalis* gastrula-stage embryos.

## Overview

This pipeline processes single-nucleus RNA-seq data from *Xenopus tropicalis* embryos and performs quality control, normalization, dimensionality reduction, clustering, marker-gene visualization, doublet detection, cell type annotation, and differential gene expression analysis.

The Stage 10.5 notebook contains the full preprocessing workflow, while the Stage 13 notebook uses an existing AnnData object and focuses on cluster annotation and refinement.

---

# Input Data

## Stage 10.5 Early Gastrula

Input file:

```text
PARSE-WTSt105-Unfiltered_FromSeurat.h5ad
```

The dataset is loaded into Scanpy as an AnnData object.

Initial dataset size:

```text
215,209 cells × 21,821 genes
```

The notebook metadata describes this as Stage 10.5 single-nucleus RNA-seq data generated using Parse Biosciences technology.

## Stage 13 Late Gastrula

Input file:

```text
PARSE-Unfiltered_St13_Compatible.h5ad
```

The Stage 13 notebook begins from a pre-existing AnnData object and performs additional clustering, subsetting, marker-based annotation, annotation refinement, and export.

---

# Software

The analysis uses the following Python packages:

```text
scanpy
anndata
numpy
pandas
scipy
matplotlib
scrublet
seaborn
stream2
```

Primary analysis environment:

```text
Python / Jupyter Notebook
```

---

# Workflow Summary

```text
Load AnnData object
        ↓
Initial cell filtering
        ↓
Mitochondrial and ribosomal QC metric calculation
        ↓
Doublet prediction with Scrublet
        ↓
Normalization and log transformation
        ↓
Highly variable gene selection
        ↓
PCA
        ↓
Neighborhood graph construction
        ↓
UMAP / tSNE visualization
        ↓
Marker-gene visualization
        ↓
Doublet removal
        ↓
Recompute neighborhood graph and UMAP
        ↓
Leiden clustering
        ↓
Marker-based cell type annotation
        ↓
Nearest-neighbor label refinement
        ↓
Differential gene expression analysis
        ↓
Export annotated AnnData object and figures
```

---

# Stage 10.5 Processing Workflow

## 1. Load Data

The Stage 10.5 AnnData object is loaded using:

```python
adataR = sc.read_h5ad(...)
```

Initial object:

```text
215,209 cells × 21,821 genes
```

---

## 2. Initial Cell Filtering

Cells expressing fewer than 1,000 genes are removed.

```python
sc.pp.filter_cells(adataR, min_genes=1000)
```

After filtering:

```text
15,920 cells × 21,821 genes
```

Genes are not filtered at this stage. The notebook notes that lowly expressed genes were intentionally retained until after clustering to avoid removing genes that may define small cell populations.

---

## 3. Mitochondrial and Ribosomal QC Metrics

Mitochondrial genes are identified using gene-name prefixes such as:

```text
cox, COX, ND, rnr1, mt-, atp6, atp8, timm
```

Ribosomal genes are identified using:

```text
rps, rpl
```

QC metrics are calculated using:

```python
sc.pp.calculate_qc_metrics(
    adataR,
    qc_vars=['mt', 'ribo'],
    percent_top=None,
    log1p=False,
    inplace=True
)
```

The following QC metrics are generated:

```text
n_genes_by_counts
total_counts
pct_counts_mt
pct_counts_ribo
```

Mitochondrial and ribosomal thresholds are discussed in the notebook, but the actual filtering commands are commented out. Therefore, mitochondrial and ribosomal content were evaluated but not used to remove cells in the uploaded Stage 10.5 workflow.

---

## 4. Doublet Detection

Doublets are predicted using Scrublet.

```python
scrub = scr.Scrublet(adataR.raw.X)
adataR.obs['doublet_scores'], adataR.obs['predicted_doublets'] = scrub.scrub_doublets()
```

Scrublet automatically selected a doublet score threshold of:

```text
0.20
```

Scrublet output:

```text
Detected doublet rate: 11.4%
Expected doublet rate: 10.0%
Estimated overall doublet rate: 16.0%
Predicted doublets: 1,820 cells
```

Predicted doublets are stored in:

```text
adataR.obs['predicted_doublets']
adataR.obs['doublet_scores']
adataR.obs['doublet_info']
```

Doublet predictions are visualized using violin plots and UMAP plots.

---

## 5. Normalization and Log Transformation

Expression values are normalized to 10,000 counts per cell.

```python
sc.pp.normalize_per_cell(adataR, counts_per_cell_after=1e4)
```

The normalized matrix is log-transformed.

```python
sc.pp.log1p(adataR)
```

The normalized/log-transformed object is stored in:

```python
adataR.raw = adataR
```

---

## 6. Highly Variable Gene Selection

Highly variable genes are identified using:

```python
sc.pp.highly_variable_genes(
    adataR,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.05
)
```

Number of highly variable genes identified:

```text
6,439
```

The AnnData object is then subset to highly variable genes.

```python
adataR = adataR[:, adataR.var['highly_variable']]
```

---

## 7. PCA

Principal Component Analysis is performed using:

```python
sc.tl.pca(adataR, svd_solver='arpack', n_comps=100)
```

PCA variance ratio and PCA loadings are visualized to assess major axes of variation.

---

## 8. tSNE

tSNE is calculated using different values of `n_pcs`, including:

```text
30 PCs
99 PCs
```

Example:

```python
sc.tl.tsne(adataR, n_pcs=30)
```

tSNE is used as an exploratory visualization.

---

## 9. Neighborhood Graph and UMAP

A nearest-neighbor graph is computed from PCA coordinates.

Different `n_pcs` values are tested:

```text
10 PCs
20 PCs
50 PCs
99 PCs
```

A commonly used setting in the notebook is:

```python
sc.pp.neighbors(adataR, n_pcs=99, n_neighbors=50)
sc.tl.umap(adataR, n_components=2)
```

UMAP is used for visualization of cellular transcriptional relationships.

---

## 10. Marker-Gene Visualization

Known marker genes are visualized on UMAP to assess biological structure and doublet predictions.

Marker categories include:

### Epidermal markers

```text
krt7
grhl1
```

### Neuroectoderm markers

```text
sox2
sox3
pax3
```

### Mesodermal markers

```text
tbxt
fgf8
```

### Endodermal markers

```text
darmin
gata6
```

---

## 11. Doublet Removal

Predicted doublets are removed using:

```python
adataR = adataR[adataR.obs['doublet_info'] == 'False', :]
```

Post-doublet-removal dataset size:

```text
14,100 cells × 21,821 genes
```

UMAP is recalculated after doublet removal.

---

## 12. Leiden Clustering

Leiden clustering is performed at multiple resolutions.

```python
sc.tl.leiden(adataR, resolution=0.1, key_added="leiden_0.1")
sc.tl.leiden(adataR, resolution=0.3, key_added="leiden_0.3")
sc.tl.leiden(adataR, resolution=0.8, key_added="leiden_0.8")
sc.tl.leiden(adataR, resolution=1.2, key_added="leiden_1.2")
sc.tl.leiden(adataR, resolution=10.0, key_added="leiden_10.0")
```

The notebook primarily uses low-resolution Leiden clusters for broad germ-layer annotation.

---

## 13. Stage 10.5 Cell Type Annotation

Cells are annotated using marker-gene thresholding.

For broad annotation, the notebook uses marker combinations such as:

| Category      | Marker genes      |
| ------------- | ----------------- |
| Mesoderm      | `tbxt`, `fgf8`    |
| Epidermis     | `krt7`, `grhl1`   |
| Neuroectoderm | `sox2`, `sox3`    |
| Endoderm      | `darmin`, `gata6` |

Cells are assigned to categories if marker expression exceeds a threshold of approximately:

```text
0.5
```

Unassigned cells are initially labeled:

```text
Unknown
```

Unknown labels are then refined using nearest-neighbor propagation from nearby labeled cells.

---

## 14. Differential Gene Expression Analysis

Differential expression is performed using Scanpy's `rank_genes_groups()` function.

Primary method used:

```python
method='t-test'
```

Examples include:

```python
sc.tl.rank_genes_groups(adataR, 'category', method='t-test')
```

Pairwise comparisons include:

```text
Mesoderm1 vs Mesoderm2
Epidermis1 vs Epidermis2
```

Outputs include ranked marker genes, test statistics, p-values, adjusted p-values, and log fold changes.

---

# Stage 13 Annotation Workflow

The Stage 13 notebook starts from a pre-existing AnnData object:

```python
corr_data = sc.read_h5ad(...)
```

The uploaded Stage 13 notebook does not contain the full preprocessing/QC workflow. Instead, it performs cluster annotation and annotation refinement.

---

## 1. Leiden Clustering

Leiden clustering is run at multiple resolutions:

```python
sc.tl.leiden(corr_data, resolution=0.2, key_added="leiden_0.2")
sc.tl.leiden(corr_data, resolution=0.3, key_added="leiden_0.3")
sc.tl.leiden(corr_data, resolution=0.4, key_added="leiden_0.4")
sc.tl.leiden(corr_data, resolution=1.2, key_added="leiden_1.2")
sc.tl.leiden(corr_data, resolution=1.4, key_added="leiden_1.4")
```

---

## 2. Inner Ectoderm Subsetting and Annotation

The following Leiden 0.3 clusters are subset for inner ectoderm analysis:

```text
3, 7, 2, 8, 9, 4
```

Marker-based annotations include:

| Annotation                        | Marker genes      |
| --------------------------------- | ----------------- |
| Anterior Neural Plate             | `rax`, `pax6`     |
| Neural Plate                      | `sox2`, `sox3`    |
| Neural Plate Border               | `pax3`            |
| Neural Crest                      | `snai2`, `sox9`   |
| six1/eya1 PPR                     | `six1`, `eya1`    |
| pax8 PPR                          | `pax8`            |
| pitx1 PPR                         | `pitx1`           |
| Anterior Neural Ridge             | `foxg1`, `emx1l`  |
| Inner NNE Prog.                   | `tp63`            |
| Notoplate                         | `shh`, `ptch2`    |
| Early Neuron                      | `ebf2`            |
| Neural Plate HB                   | `pax2`            |
| Basal Cell                        | `foxa1`, `pou2f3` |
| Lateral Mesoderm                  | `tbxt`            |
| Cranial NCC                       | `rpe65`, `mafb`   |
| Ciliated L/R Organizer Primordium | `foxj1`           |

Nearest-neighbor label smoothing is applied using the KNN graph.

---

## 3. Outer Ectoderm and Endoderm Subsetting and Annotation

The following Leiden 0.3 clusters are subset:

```text
0, 10, 11, 14, 15, 13
```

Marker-based annotations include:

| Annotation                     | Marker genes     |
| ------------------------------ | ---------------- |
| Cement Gland                   | `agr2`, `pitx1`  |
| Goblet Cell                    | `itln1`, `itln2` |
| Outer NNE                      | `krt70`          |
| Ciliated Epidermal Prog.       | `foxj1`, `tp73`  |
| Ionocyte                       | `foxi1`          |
| Foregut Primordium             | `slc5a8`         |
| Midgut Primordium              | `ndrg1`          |
| Hindgut Primordium             | `gatm`           |
| Ventral Involuted Mesoderm     | `cdx4`, `tbxt`   |
| Outer Neural Plate Border      | `pax3`           |
| Liver Primordium               | `hhex`           |
| Anterior Pharyngeal Primordium | `gsc`            |
| Pharyngeal Primordium          | `nkx2-6`         |
| Pancreatic Primordium          | `pdx1`           |

Additional Leiden clustering is performed within this subset, and clusters 15 and 16 are removed from the outer ectoderm subset.

---

## 4. Mesoderm Subsetting and Annotation

The following Leiden 0.3 clusters are subset for mesoderm analysis:

```text
5, 6, 1, 12
```

Marker-based annotations include:

| Annotation                  | Marker genes    |
| --------------------------- | --------------- |
| Pronephric Primordium       | `pax8`, `lhx1`  |
| Anterior Paraxial Mesoderm  | `ank1`, `actc1` |
| Posterior Paraxial Mesoderm | `myf5`, `hoxd3` |
| Pharyngeal Mesoderm         | `tbx1`          |
| Cardiac Mesoderm            | `nkx2-5`        |
| Ventral Blood Island        | `runx1`         |
| Ventral Mesoderm            | `gata6`         |
| Notochord                   | `not`           |
| Lateral Mesoderm            | `foxc1`, `tbxt` |
| L/R Organizer Primordium    | `cirop`         |

Nearest-neighbor label refinement is also applied to mesoderm annotations.

---

## 5. Combining Annotations

Annotations from separate subset analyses are transferred back into the full `corr_data` object.

Intermediate annotation columns include:

```text
region_annotation_meso
region_annotation_outerecto
region_annotation_innerecto
```

These are combined into a final annotation column:

```text
refined_region_annotation
```

Cells without a refined label are assigned:

```text
Other
```

---

## 6. Export

The final annotated Stage 13 object is saved as:

```text
Stage13_WT_Annotated.h5ad
```

---

# Main Outputs

Expected outputs from the complete workflow include:

| Output                          | Description                                             |
| ------------------------------- | ------------------------------------------------------- |
| `adata_processed.h5ad`          | Processed Stage 10.5 AnnData object                     |
| `Stage13_WT_Annotated.h5ad`     | Annotated Stage 13 AnnData object                       |
| QC violin plots                 | Gene counts, total counts, mitochondrial %, ribosomal % |
| Scrublet histogram              | Doublet score distribution                              |
| UMAP plots                      | Visualization of clusters, markers, and annotations     |
| Leiden cluster labels           | Cluster assignments across multiple resolutions         |
| Marker gene plots               | UMAPs showing known lineage markers                     |
| Differential expression results | Ranked genes by cell category or pairwise comparison    |
| Annotation metadata             | Final cell type labels stored in `.obs`                 |

---

# Notes and Recommended Improvements

The uploaded notebooks perform the major steps expected in a single-cell analysis. To make the pipeline more reproducible and industry-ready, the following improvements are recommended:

1. Save all key parameters in a `config.yaml` file.
2. Save all run-specific values in `results/run_config.yaml`.
3. Save QC summaries as CSV files.
4. Re-run normalization, highly variable gene selection, PCA, neighbors, UMAP, and clustering after doublet removal.
5. Use fixed random seeds for PCA, UMAP, tSNE, and Leiden clustering where possible.
6. Export marker-gene and differential-expression tables instead of only plotting them.
7. Replace manual notebook paths with relative paths or config-based paths.
8. Record package versions automatically.
9. Consider using `method='wilcoxon'` for marker-gene discovery in addition to, or instead of, `method='t-test'`.
10. Convert the notebooks into modular scripts or a Snakemake/Nextflow workflow for better reproducibility.
