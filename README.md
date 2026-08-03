# Xenopus tropicalis Single-Nucleus(Cell) and Multiomic Analysis

This repository contains analysis workflows for **gastrula-stage *Xenopus tropicalis*** single-nuclei and multiomic sequencing data. The repository is organized into two complementary analysis branches:

1. **Multiomic RNA + ATAC analysis** using Seurat and Signac.
2. **snRNA-seq quality control and cell-type annotation** using Scanpy.

The goal of this repository is to provide a reproducible framework for processing, quality controlling, integrating, visualizing, and annotating early embryonic nuclei populations across transcriptomic and chromatin-accessibility measurements.

---

## Repository Structure

```text
.
├── Multiomic/
│   ├── scripts/
│   ├── config/
│   ├── results/
│   ├── archive/
|   ├── scripts/
│   └── docs/
│
├── QC_Annotation/
|   ├── config/
|   ├── tables/
│   ├── notebooks/
│   ├── results/
│   ├── reproducibility/
│   └── docs/
│
└── README.md
```

> Folder names may differ slightly depending on the local clone. Each analysis folder contains its own documentation describing the corresponding workflow in more detail.

---

# 1. Multiomic RNA + ATAC Analysis

The `Multiomic/` folder contains the paired RNA-seq and ATAC-seq analysis workflow for 10x Genomics Multiome data.

The workflow is implemented primarily in **R** using:

* Seurat
* Signac
* SeuratObject
* EnsDb / ensembldb
* ggplot2
* patchwork

The analysis begins with the filtered Cell Ranger ARC feature matrix and ATAC fragment file and constructs a Seurat object containing both RNA and chromatin-accessibility assays.

## Multiomic workflow

```text
10x Genomics Multiome output
        │
        ├── Gene Expression matrix
        └── ATAC peak matrix + fragments
                ↓
        Create Seurat object
                ↓
        RNA and ATAC quality control
                ↓
        Cell filtering
                ↓
        RNA normalization and PCA
                ↓
        ATAC TF-IDF and LSI
                ↓
        Weighted Nearest Neighbor integration
                ↓
        Joint UMAP
                ↓
        Leiden clustering
                ↓
        RNA expression-based cell annotation
                ↓
        Annotation validation
                ↓
        ATAC coverage and regulatory analyses
```

## RNA processing

RNA measurements are processed using standard Seurat workflows, including normalization, variable-feature selection, dimensionality reduction, and SCTransform where appropriate.

RNA expression is also used to assign biological cell identities following marker-based annotation rules.

## ATAC processing

Chromatin accessibility is analyzed using Signac. Major steps include:

* ATAC fragment-based quality control
* nucleosome signal calculation
* transcription start site enrichment
* peak accessibility matrix processing
* TF-IDF normalization
* latent semantic indexing
* coverage visualization
* peak-associated regulatory analysis

The Cell Ranger ARC peak matrix is used as the primary accessibility feature set unless a new peak set is deliberately generated.

## Joint RNA + ATAC integration

RNA PCA and ATAC LSI representations are integrated using Seurat's **Weighted Nearest Neighbor (WNN)** framework.

The resulting weighted neighbor graph is used for:

* joint UMAP visualization
* multimodal clustering
* comparison of transcriptional and chromatin states

## Expression-based annotation

Cell identities are assigned using normalized RNA expression and a configurable set of marker-gene rules.

Rather than relying only on manual renaming of numerical clusters, the annotation workflow records:

* marker-based cell assignments
* the rule responsible for each assignment
* cells satisfying multiple annotation rules
* unresolved cells
* cluster-to-annotation relationships

The annotation rules are stored separately from the analysis code so that marker combinations and thresholds can be inspected and modified without rewriting the main workflow.

Example lineage markers are used to identify populations associated with:

* Neural Plate
* Neural Plate Border
* Non-neural Ectoderm
* Ciliated Epidermal Progenitors
* Cement Gland
* Dorsal and Ventral Mesoderm
* Dorsal and Ventral Endoderm

## Regulatory analysis

The multiomic workflow also supports downstream chromatin analyses such as:

* locus-specific accessibility visualization with `CoveragePlot`
* read-density visualization
* GC-content calculation with `RegionStats`
* peak-to-gene association with `LinkPeaks`

These analyses require appropriate *X. tropicalis* genomic annotation and genome-reference resources.

---

# 2. snRNA-seq Quality Control and Annotation

The `QC_Annotation/` folder contains the single-nucleus RNA-seq preprocessing, quality-control, clustering, visualization, and annotation workflow.

This analysis is implemented primarily in **Python** using Scanpy and AnnData.

## snRNA-seq workflow

```text
Raw / filtered expression matrix
        ↓
AnnData construction
        ↓
Quality-control metrics
        ↓
Cell and gene filtering
        ↓
Doublet detection / removal
        ↓
Normalization and log transformation
        ↓
Highly variable gene selection
        ↓
PCA
        ↓
Neighbor graph
        ↓
UMAP
        ↓
Leiden clustering
        ↓
Marker-gene inspection
        ↓
Cell-type annotation
        ↓
Annotation validation
```

## Quality control

QC analysis includes assessment of features such as:

* reads or counts per cell
* genes detected per cell
* doublet identification
* expression distributions before and after filtering
* cell retention following QC

QC figures and summary tables are saved separately from downstream biological analyses.

## Cell-type annotation

Cell identities are assigned using lineage-specific marker expression and cluster structure.

The annotation workflow is designed to preserve both:

1. the underlying numerical cluster identity, and
2. the biologically interpreted cell-type annotation.

Marker expression is evaluated through UMAP visualization, expression plots, and targeted marker panels.

---

# Relationship Between the Two Analysis Branches

The two folders address complementary aspects of the same biological system.

The snRNA-seq workflow provides a transcriptome-focused framework for:

* quality control
* cluster identification
* marker discovery
* cell-type annotation

The multiomic workflow extends this framework by jointly examining:

* gene expression
* chromatin accessibility

Together, these analyses allow transcriptional cell states to be compared with their underlying regulatory landscapes.

---

# Data

Large sequencing files are intentionally excluded from GitHub.

Typical required input files include:

```text
filtered_feature_bc_matrix.h5
atac_fragments.tsv.gz
atac_fragments.tsv.gz.tbi
```

Additional reference resources may include:

```text
Xenopus tropicalis genome annotation
EnsDb annotation package
BSgenome reference package
```
---

# Reproducibility

The repository is organized so that:

* analysis code is separated from raw data
* large intermediate objects are not committed to GitHub
* marker-based annotation rules are stored explicitly
* generated figures and tables are separated by analysis stage
* validation outputs document whether scripts and notebooks are structurally valid
* package versions can be recorded for reproducibility

For the R multiomic workflow, `renv` can be used to preserve package versions.

For the Python snRNA-seq workflow, package versions should be recorded using the active conda or Python environment.

---

# Software

## R / Multiomic analysis

Core packages include:

```text
Seurat
Signac
SeuratObject
ggplot2
patchwork
ensembldb
GenomicRanges
```

## Python / snRNA-seq analysis

Core packages include:

```text
scanpy
anndata
numpy
pandas
matplotlib
scipy
scrublet
```

Exact package versions should be recorded with each analysis release.

---

# Project Status

This repository is under active development.

The current workflows support:

* snRNA-seq preprocessing and QC
* doublet filtering
* dimensionality reduction
* Leiden clustering
* marker-based cell annotation
* paired RNA + ATAC processing
* WNN multimodal integration
* multimodal UMAP visualization
* expression-driven annotation
* ATAC locus visualization
* peak-to-gene analysis
---

# Organism

**Species:** *Xenopus tropicalis*
**Developmental context:** Early gastrula-stage embryos
**Modalities:** snRNA-seq and paired RNA + ATAC multiome sequencing

