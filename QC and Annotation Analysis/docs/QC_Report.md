# Quality-Control Report

## Scope

This report summarizes the thresholds implemented in the uploaded notebooks
and the values recorded in their existing execution outputs. Because several
notebooks were edited after execution, all values should be regenerated from a
clean top-to-bottom run before manuscript submission or archival release.

## Implemented thresholds

| Step | Implemented criterion |
|---|---|
| Cell complexity | at least 1,000 detected genes |
| Mitochondrial fraction | `pct_counts_mt < 20` |
| Ribosomal fraction | `pct_counts_ribo < 3.5` |
| Library-size normalization | 10,000 counts per cell |
| Transformation | natural-log `log1p` |
| Doublet detection | Scrublet automatic threshold |
| All-group DEG inclusion | at least 7 cells per annotation |
| Pairwise DEG inclusion | at least 2 cells in each group |

## Notebook-derived retention

| Dataset       |   Input barcodes |   Genes |   ≥1,000 genes |   MT/ribo filters |   Predicted doublets |   Retained singlets |   Final annotated |
|:--------------|-----------------:|--------:|---------------:|------------------:|---------------------:|--------------------:|------------------:|
| Stage10 WT    |           215209 |   21821 |          15920 |             15920 |                 1820 |               14100 |             14100 |
| Stage11.25 WT |           110592 |   21821 |           7971 |              7971 |                 1446 |                6525 |              6516 |
| Stage11.25 MO |           110590 |   21821 |           6533 |              6533 |                  503 |                6030 |              5995 |
| Stage12.0 WT  |           110592 |   21821 |          14754 |             14754 |                    0 |               14754 |             10789 |
| Stage12.0 MO  |           110590 |   21821 |          13967 |             13967 |                 1369 |               12598 |              8598 |
| Stage12.5 WT  |           110592 |   29231 |          11908 |             11908 |                 1608 |               10300 |             10237 |
| Stage13 WT    |           110584 |   29231 |          10074 |             10074 |                 1166 |                8908 |              8837 |

The seven notebook outputs sum to:

- **878,749** input barcodes;
- **81,127** cells after the implemented complexity and
  mitochondrial/ribosomal filters;
- **7,912** Scrublet-predicted doublets;
- **73,215** retained singlets immediately after doublet removal;
- **65,072** cells represented in the final annotation-count
  outputs parsed from the notebooks.

The difference between retained singlets and final annotated cells reflects
later dataset-specific filtering, removal of `Unknown`/`Other`, or subcluster
selection.

## Dataset-specific Scrublet results

| Dataset       |   Automatic threshold |   Detected doublet rate (%) |   Predicted doublets |
|:--------------|----------------------:|----------------------------:|---------------------:|
| Stage10 WT    |                  0.2  |                        11.4 |                 1820 |
| Stage11.25 WT |                  0.15 |                        18.1 |                 1446 |
| Stage11.25 MO |                  0.24 |                         7.7 |                  503 |
| Stage12.0 WT  |                  0.76 |                         0   |                    0 |
| Stage12.0 MO  |                  0.22 |                         9.8 |                 1369 |
| Stage12.5 WT  |                  0.17 |                        13.5 |                 1608 |
| Stage13 WT    |                  0.16 |                        11.6 |                 1166 |

## QC observations requiring attention

1. **Stage 12.0 WT reports zero predicted doublets.** Scrublet selected a
   threshold of 0.76 and returned a 0% detected doublet rate. This differs
   strongly from the other datasets and should be inspected using the score
   histogram, marker co-expression, and expected multiplet rate.
2. **The mitochondrial and ribosomal filters did not reduce the post-minimum-
   gene cell counts in the stored outputs.** This may be correct, but should be
   verified after rerun.
3. **The ribosomal narrative and code need consistent wording.** The notebooks
   retain cells below 3.5% ribosomal counts. Documentation should not describe
   the criterion as removal of low-ribosomal cells unless the code changes.
4. **Doublets are removed after exploratory normalization and embeddings.**
   Some notebooks recompute neighbors and UMAP afterward, while others retain
   earlier representations. A stricter production pipeline would remove
   doublets before the final HVG/PCA/neighbor computation.
5. **Executed outputs contain older absolute Windows paths.** The clean
   notebooks remove outputs so that GitHub does not display local user paths.

## QC outputs written by each notebook

All Section 1 figures now write directly to the stage-specific `qc/` directory.
The notebooks explicitly save:

- combined violin plots for detected genes, total counts, mitochondrial counts,
  and ribosomal counts;
- top expressed, mitochondrial, and ribosomal gene bar graphs;
- mitochondrial and ribosomal percentage distributions;
- the Scrublet doublet-score histogram;
- the doublet/singlet detected-gene violin plot;
- the highly variable gene plot;
- a post-doublet histogram with detected reads per cell on the x-axis and number
  of retained cells on the y-axis;
- a post-doublet mean-versus-median reads-per-cell bar graph.

The post-doublet cells also save `post_doublet_qc_summary.csv` and
`post_doublet_reads_per_cell_summary.csv` in the same `qc/` directory.

## QC tables generated by the cleaned notebooks

The following tables are recommended and implemented where possible:

- final annotation count and percentage table;
- Leiden-cluster by annotation count table;
- Leiden-cluster by annotation row-percentage table;
- before/after label-refinement transition table;
- DEG group-count table;
- list of groups excluded from DEG because of small sample size;
- pairwise-comparison inventory with group sizes and completion status.

## Rerun checklist

- restart the kernel;
- run every cell in order;
- confirm input shape;
- confirm QC counts after each filter;
- inspect the Scrublet histogram and marker overlays;
- verify final annotation counts sum to `adata.n_obs`;
- verify every selected DEG group has the intended number of cells;
- archive the regenerated CSV summaries alongside the notebook version.
