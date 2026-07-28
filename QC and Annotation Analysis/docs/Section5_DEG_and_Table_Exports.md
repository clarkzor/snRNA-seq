# Section 5: DEG and Table Exports

## Why Section 5 was replaced

The original Section 5 saved visualizations but did not save the statistical
results. Some pairwise plots also referenced a DEG key that did not exist on
the pairwise AnnData object.

The GitHub-ready notebooks use one consistent workflow implemented with
`scripts/output_utils.py`.

## All-annotation output files

For a dataset such as `Stage12_0-WT`, the all-group test creates:

```text
differential_expression/
├── stage12_0-wt_all_annotations_vs_rest_wilcoxon_all_results.csv
├── stage12_0-wt_all_annotations_vs_rest_wilcoxon_significant_fdr0p05_abslog2fc1p0.csv
├── stage12_0-wt_all_annotations_vs_rest_wilcoxon_top25_per_group.csv
├── stage12_0-wt_all_annotations_vs_rest_wilcoxon_group_counts.csv
└── stage12_0-wt_all_annotations_vs_rest_wilcoxon_parameters.json
```

The exact prefix is slugified to remain filesystem safe.

## Pairwise output files

Each configured comparison creates the same five files, using a comparison
slug such as:

```text
stage12_0-wt_neural_plate_1_vs_2_wilcoxon_*
```

The test is run after subsetting to exactly the two specified groups. Both
groups are included under one key with `reference="rest"`.

## Statistical columns

The full result table includes available Scanpy fields:

- `group`
- `names`
- `scores`
- `logfoldchanges`
- `pvals`
- `pvals_adj`
- `pct_nz_group`
- `pct_nz_reference`

The percentage columns require `pts=True`, which is enabled in the revised
tests.

## Filtering used for the significant table

The helper applies:

```text
adjusted p-value ≤ 0.05
absolute log2 fold change ≥ 1.0
```

The full unfiltered result is always saved, so thresholds can be changed later
without rerunning the statistical test.

## Summary tables

Section 5 also saves:

- final annotation counts and percentages;
- annotations excluded because they contain fewer than seven cells;
- pairwise comparison group sizes;
- whether each pair completed or was skipped.

Section 6 saves:

- final annotation counts after renaming;
- every retained Leiden-to-annotation count crosstab;
- row-normalized Leiden-to-annotation percentages;
- before/after transition tables for label-refinement backups.

## Highest-priority tables for publication or supplement

1. final cell counts and percentages by stage and annotation;
2. complete all-annotation marker table;
3. top 25 markers per annotation;
4. pairwise layer-comparison DEG tables;
5. marker-expression summary by annotation;
6. label-refinement transition table;
7. cluster-to-annotation correspondence table.

## Adding a new comparison

Add a tuple to `PAIRWISE_COMPARISONS`:

```python
PAIRWISE_COMPARISONS = [
    ("Epidermis 1", "Epidermis 2", "epidermis_1_vs_2"),
    ("Neural Plate 1", "Neural Plate 2", "neural_plate_1_vs_2"),
]
```

The first two values must exactly match `region_annotation`. The third value is
the output filename slug.
