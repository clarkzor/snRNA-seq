## Overview

Quality Control and Cell-Type Annotation analysis for single cell/nucleus RNA-sequencing datasets. Here, unfiltered [Gene X Expression] matrices produced from alignment are passed through different quality control metrics to remove unwanted cells and utilize gene expression signatures for cell-type annotation.

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
