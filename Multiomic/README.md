# Xenopus tropicalis Stage 11 RNA+ATAC multiome

Reproducible Seurat/Signac workflow for paired 10x Genomics RNA + ATAC multiome data from
**Xenopus tropicalis Stage 11**. The repository separates raw inputs, configurable analysis
decisions, modular processing scripts, expression-based cell-type annotation, outputs, and
validation.

## Analysis overview

```text
Cell Ranger ARC RNA + ATAC matrices
            |
            v
   Seurat + ChromatinAssay
            |
            v
       RNA / ATAC QC
            |
      +-----+-----+
      |           |
      v           v
 RNA PCA      ATAC TF-IDF/LSI
      |           |
      +-----+-----+
            |
            v
      WNN integration
            |
            v
       WNN UMAP + Leiden
            |
            v
 RNA expression-rule annotation
            |
            v
 annotation validation / conflicts
            |
            v
 coverage + peak-gene links (optional)
```

## Repository structure

```text
config/      Analysis parameters, local path template, annotation rules, marker panels
scripts/     Numbered modular pipeline plus shared utility functions
notebooks/   Narrative RMarkdown entry point and original PBMC demo reference
results/     Generated QC, RNA, ATAC, integration, annotation, coverage, links, checkpoints
validation/  Static repository checks and runtime reproducibility summaries
docs/        Pipeline, annotation, QC, ATAC, and reproducibility documentation
data/        Input-data documentation; large data are ignored by Git
archive/     Original exploratory and single-script reorganized workflows for provenance
```

## Required input

Core matrix processing requires:

```text
filtered_feature_bc_matrix.h5
```

Fragment-level QC/coverage additionally require:

```text
atac_fragments.tsv.gz
atac_fragments.tsv.gz.tbi
```

Large data files are intentionally excluded from Git.

## Local paths

The downloadable bundle includes an ignored `config/paths.local.R` configured for:

```text
C:/Users/coron/OneDrive/Desktop/10xGenomics
```

For another machine, copy `config/paths_example.R` to `config/paths.local.R` and edit it,
or set the `MULTIOME_DATA_DIR` environment variable.

## Run the pipeline

From the repository root:

```bash
Rscript scripts/run_all.R
```

Or run the numbered scripts one at a time. Each stage writes an RDS checkpoint to
`results/objects/`, allowing later scripts to reload the previous stage.

## Annotation strategy

The original numeric-cluster manual map is preserved in
`config/original_manual_cluster_map.csv`. The primary workflow instead performs prioritized
cell-level marker-expression annotation using
`config/stage11_expression_annotation_rules.csv`.

The rule engine records both the selected label and every candidate label so ambiguous cells
can be reviewed rather than hidden. Cluster-consensus filling of unresolved cells is optional
and disabled by default.

## Validation

Run static repository validation with:

```bash
Rscript scripts/validate_repository.R
```

Runtime analysis also writes input summaries, annotation validation summaries, package
versions, `sessionInfo()`, and a final run summary under `validation/`.

## Reproducibility

Package installation is deliberately separated from the analysis. See
`docs/Reproducibility.md`. Once the working environment is confirmed, create and commit an
`renv.lock` using `renv::snapshot()`.

## Important optional steps

Coverage plots and peak-to-gene links are disabled by default. They require fragment files;
peak-to-gene links also require the matching Xenopus BSgenome package. Keeping them optional
prevents environment-specific plotting or genome-package issues from stopping the core RNA,
ATAC, WNN, clustering, and annotation workflow.
