# Notebook Audit

## Files reviewed

Seven notebooks were reviewed, covering Stage 10/10.5, Stage 11.25 WT/MO,
Stage 12.0 WT/MO, Stage 12.5 WT, and Stage 13 WT.

The machine-readable inventory is available in:

- `tables/notebook_inventory.csv`
- `tables/notebook_qc_summary.csv`
- `tables/final_annotation_counts_from_notebook_outputs.csv`
- `tables/figure_filename_changes.csv`

## Main findings

### 1. Figure filenames depended on notebook cell numbers

The notebooks used suffixes such as:

```text
_cell088_rank_genes_groups_1.png
_cell109_umap_1.png
```

Cell numbers are unstable and do not explain the biological content. The
cleaned copies replace **293 non-Section-5 save suffixes**
with section-based descriptive names. Section 5 now saves exact filenames
through `scripts/output_utils.py`.

### 2. Differential-expression tables were not exported

The original notebooks stored `rank_genes_groups` results only in memory and
saved plots. The replacement Section 5 exports:

- all genes;
- significant genes at FDR ≤ 0.05 and |log2FC| ≥ 1;
- top 25 genes per group;
- test parameters;
- group counts;
- comparison completion status.

### 3. Some DEG plots used the wrong key

Several heatmap, dotplot, and matrix-plot cells requested `key="wilcoxon"` on a
pairwise subset even though that key had not been generated on the subset.
Other notebooks generated reciprocal comparisons under separate keys and then
attempted to plot both groups from a single unrelated key.

The replacement Section 5 runs both groups under one key on a two-group subset,
so every plot and export uses a result generated on the same object.

### 4. Section headings described t-tests while code used Wilcoxon

The original text called the analysis a t-test, but active code specified
`method="wilcoxon"`. The cleaned notebooks now describe the Wilcoxon test.

### 5. Pairwise subsection labels were duplicated or inaccurate

Later-stage notebooks used a second heading reading “Neuroectoderm 1 and
Neuroectoderm 2” for the Neural Plate Border comparison. The replacement code
uses explicit comparison names and an inventory table.

### 6. Path templates were inconsistent

- six source notebooks omitted `Stage12` from the stage-folder mapping;
- Stage 10 loaded from a hard-coded `St10` subfolder rather than
  `STAGE_RAW_DATA_DIR`;
- executed outputs still displayed older absolute Windows paths and older stage
  names.

The cleaned copies use one consistent project-relative path template.

### 7. Final output names were inconsistent

Examples included `Stage10_5`, `Stage12`, and `Stage12_0` naming differences.
The cleaned notebooks define one `DATASET_SLUG` and use it for figures, tables,
metadata, and AnnData outputs.

### 8. The README cell total was not supported by notebook outputs

The prior README stated 86,350 cells. Existing notebook outputs record
73,215 retained singlets and 65,072 cells in the
parsed final annotation tables. This discrepancy is documented rather than
silently resolved.

### 9. Regression and batch correction were not active

The notebooks contain commented regression, MNN, and batch-aware HVG examples.
They should not be listed as completed pipeline steps unless activated and
validated.

### 10. Stored outputs are stale relative to source code

Some output text reports stage/path values that differ from the current source
cells. Clean notebooks have outputs removed and must be rerun.

## Changes made in the GitHub-ready copies

- standardized stage and condition paths;
- standardized dataset slugs;
- corrected the 1,000-gene filtering comment;
- renamed figure-save suffixes outside Section 5;
- replaced Section 5 with consistent Wilcoxon analysis and CSV exports;
- added final annotation, cluster cross-tab, and label-transition tables;
- standardized final `.h5ad` and metadata filenames;
- removed execution outputs and counters;
- preserved all annotation logic outside the replaced Section 5.

## Items intentionally not changed

The audit did not reinterpret marker thresholds or change the biological
annotation rules. Scientific changes to thresholds, protected labels, or
lineage definitions should be made deliberately and documented separately.
