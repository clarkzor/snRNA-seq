## Overview

Quality control (QC) and cell-type annotation analysis for single cell/nucleus RNA-sequencing datasets. Here, unfiltered [Gene X Expression] matrices produced from alignment are passed through different quality control metrics to remove unwanted cells and utilize gene expression signatures for cell-type annotation.

## Dataset

snRNAseq from X. tropicalis gastrulation, stages 10-13
- 86,350 cells across 5 temporal samples
- Samples and libraries prepared via PARSE Evercode WT-mini, Sequenced on Illumina NovaSeq 6000

## Pipeline overview

Pre-processing and quality control
- Filtering genes/cell
- Filtering cells with abberant mitochondrial and ribosome gene expression
- Doublet identification
- Normalization (including batch effects)
- High Variable Gene calling
- Regression analysis

Linear and non-linear dimension reduction analysis techniques
- Principle Component Analysis (PCA)
- t-Stochastic Neighbor Embedding (tSNE)
- Uniform Manifold Approximation Projection (UMAP)

Marker gene expression and doublet identification
- Cell-type marker gene expression visualization
- Doublet validation via marker gene expression

Clustering algorithms and Cell-type annoation
- Leiden clustering
- Gene expression based annotation

Differential gene expression analysis
- Statistical test utilization
- Plotting Differential gene expression visualization


## Requirements

- python       3.11.5,
- scanpy       1.9.6,
- pandas       2.1.4,
- scipy        1.11.4,
- anndata      0.10.3,
- numpy        1.26.2,
- matplotlib   3.8.2,
- scrublet     0.2.3,
- seaborn      0.12.2,

# Outputs

The pipeline generates:
- Processed AnnData object (.h5ad) with field standard genes/cell and high variable gene call thresholds
- Quality control attributes including mitochondrial, ribosomal expressing cells and doublet detection
- Both linear and non-linear dimension reduced (UMAP) visualizations of gene expression
- Clustering algorithm cell-type predictions
- Marker gene cell-type annotations
- Visualizations of differential gene expression analysis

# Quality Control 

Quality control includes filtering cells based on:
- Total numbers of genes/cell
- Proportion of mitochondrial and ribosomal reads
- Doublet detection

# Validation

 Pipeline validation was performed using the publicly available bechmark dataset found at GSEX. Reproducibility of QC, marker gene detection, clustering consistency and cell-type annotation. 

 See **docs/Validation_Report.md**
