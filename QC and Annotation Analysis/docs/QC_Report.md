# Quality Control Report

## Project

Single-nucleus RNA-seq preprocessing and annotation of *Xenopus tropicalis* gastrula-stage embryos.

## Dataset

This QC report describes the quality control workflow applied to the Stage 10.5 early gastrula single-nucleus RNA-seq dataset.

Input file:

```text
PARSE-WTSt105-Unfiltered_FromSeurat.h5ad
```

Initial dataset size:

```text
215,209 cells × 21,821 genes
```

---

# QC Objectives

The objectives of quality control were to:

1. Remove low-quality barcodes with insufficient gene detection.
2. Evaluate mitochondrial and ribosomal transcript content.
3. Predict and remove likely doublets.
4. Confirm doublet predictions using marker-gene visualization.
5. Generate a filtered dataset suitable for normalization, dimensionality reduction, clustering, and annotation.

---

# QC Workflow Summary

```text
Load raw AnnData object
        ↓
Filter low-gene-count cells
        ↓
Calculate mitochondrial and ribosomal QC metrics
        ↓
Visualize QC metrics
        ↓
Predict doublets with Scrublet
        ↓
Visualize doublet scores
        ↓
Validate doublet predictions using lineage marker genes
        ↓
Remove predicted doublets
        ↓
Proceed to downstream analysis
```

---

# 1. Initial Cell Filtering

Cells were filtered based on the number of detected genes.

Filtering command:

```python
sc.pp.filter_cells(adataR, min_genes=1000)
```

Filtering criterion:

| Metric                  | Threshold | Action      |
| ----------------------- | --------: | ----------- |
| Genes detected per cell |  `< 1000` | Remove cell |

Result:

| Metric         |   Value |
| -------------- | ------: |
| Initial cells  | 215,209 |
| Cells removed  | 199,289 |
| Cells retained |  15,920 |
| Genes retained |  21,821 |

This filtering step removed likely background barcodes and low-quality nuclei with insufficient detected genes.

Genes were not filtered at this stage. The notebook notes that lowly expressed genes were intentionally retained until after clustering to avoid removing genes that may define rare cell populations.

---

# 2. Mitochondrial and Ribosomal QC Metrics

Mitochondrial and ribosomal genes were labeled in the gene metadata.

Mitochondrial genes were identified using gene-name prefixes including:

```text
cox, COX, ND, rnr1, mt-, atp6, atp8, timm
```

Ribosomal genes were identified using:

```text
rps, rpl
```

QC metrics were calculated using:

```python
sc.pp.calculate_qc_metrics(
    adataR,
    qc_vars=['mt', 'ribo'],
    percent_top=None,
    log1p=False,
    inplace=True
)
```

The following metrics were calculated:

| Metric              | Description                                |
| ------------------- | ------------------------------------------ |
| `n_genes_by_counts` | Number of genes detected per cell          |
| `total_counts`      | Total counts per cell                      |
| `pct_counts_mt`     | Percent of counts from mitochondrial genes |
| `pct_counts_ribo`   | Percent of counts from ribosomal genes     |

These metrics were visualized using violin plots.

---

# 3. Mitochondrial and Ribosomal Filtering

The notebook discusses common single-cell RNA-seq thresholds:

| Metric              |                Potential threshold discussed |
| ------------------- | -------------------------------------------: |
| Mitochondrial reads | Remove cells with `>20%` mitochondrial reads |
| Ribosomal reads     |    Remove cells with `<3.5%` ribosomal reads |

However, these filtering commands were commented out in the uploaded notebook.

Therefore:

```text
No cells were removed based on mitochondrial percentage.
No cells were removed based on ribosomal percentage.
```

The rationale stated in the notebook is that this dataset is single-nucleus RNA-seq, where appreciable mitochondrial and ribosomal recovery is not necessarily expected in the same way as whole-cell scRNA-seq.

Cells retained after mitochondrial/ribosomal evaluation:

```text
15,920 cells
```

---

# 4. Doublet Detection

Doublets were predicted using Scrublet.

Command:

```python
scrub = scr.Scrublet(adataR.raw.X)
adataR.obs['doublet_scores'], adataR.obs['predicted_doublets'] = scrub.scrub_doublets()
```

Scrublet results:

| Metric                                         |       Value |
| ---------------------------------------------- | ----------: |
| Automatically selected doublet score threshold |        0.20 |
| Detected doublet rate                          |       11.4% |
| Expected doublet rate                          |       10.0% |
| Estimated detectable doublet fraction          |       71.6% |
| Estimated overall doublet rate                 |       16.0% |
| Predicted doublets                             | 1,820 cells |

Doublet annotations were stored in:

```text
doublet_scores
predicted_doublets
doublet_info
```

Doublet score distributions were visualized using a Scrublet histogram.

---

# 5. Doublet Validation

Before removing doublets, predicted doublets were evaluated visually.

The notebook uses known germ-layer marker genes to determine whether predicted doublets show mixed lineage marker expression.

Marker genes visualized include:

| Lineage / category | Marker genes           |
| ------------------ | ---------------------- |
| Epidermal          | `krt7`, `grhl1`        |
| Neuroectoderm      | `sox2`, `sox3`, `pax3` |
| Mesodermal         | `tbxt`, `fgf8`         |
| Endodermal         | `darmin`, `gata6`      |

Predicted doublets were visualized on UMAP together with total counts.

This provides a biological sanity check for Scrublet predictions.

---

# 6. Doublet Removal

Predicted doublets were removed using:

```python
adataR = adataR[adataR.obs['doublet_info'] == 'False', :]
```

Result after doublet removal:

| Metric                               |  Value |
| ------------------------------------ | -----: |
| Cells before doublet removal         | 15,920 |
| Predicted doublets removed           |  1,820 |
| Cells retained after doublet removal | 14,100 |
| Genes retained                       | 21,821 |

Final post-QC dataset:

```text
14,100 cells × 21,821 genes
```

---

# 7. QC Summary

| QC Step                  | Filtering applied? | Result                                    |
| ------------------------ | ------------------ | ----------------------------------------- |
| Minimum detected genes   | Yes                | Removed cells with fewer than 1,000 genes |
| Gene filtering           | No                 | All 21,821 genes retained                 |
| Mitochondrial percentage | No                 | Metrics calculated, no filtering applied  |
| Ribosomal percentage     | No                 | Metrics calculated, no filtering applied  |
| Doublet prediction       | Yes                | Scrublet used                             |
| Doublet removal          | Yes                | 1,820 predicted doublets removed          |
| Final retained cells     | Yes                | 14,100 cells retained                     |

---

# 8. Downstream Readiness

After QC, the dataset was used for:

* normalization
* log transformation
* highly variable gene selection
* PCA
* UMAP
* Leiden clustering
* marker-gene visualization
* manual cell type annotation
* differential gene expression analysis

---

# 9. QC Limitations

Several limitations should be noted:

1. Mitochondrial and ribosomal filtering thresholds were discussed but not applied.
2. QC thresholds were manually selected rather than determined from sample-specific distributions.
3. Scrublet predictions were used for doublet removal, but no orthogonal doublet detection method was used.
4. The notebook relies heavily on visual validation of marker genes.
5. Full QC summaries were not exported as structured CSV files in the uploaded notebook.
6. Reproducibility would be improved by saving all QC parameters and cell counts automatically.

---

# 10. Recommended QC Improvements

For future versions of this pipeline, the following improvements are recommended:

1. Export a `qc_summary.csv` file containing the number of cells removed at each step.
2. Export per-cell QC metrics to `qc_metrics.csv`.
3. Save doublet scores to `doublet_scores.csv`.
4. Save all QC plots as PDF or PNG files.
5. Record all QC thresholds in `config.yaml`.
6. Record the exact QC results in `results/run_config.yaml`.
7. Re-run normalization, highly variable gene selection, PCA, neighbors, UMAP, and clustering after doublet removal.
8. Include final pass/fail acceptance criteria for each QC step.
9. Add a small test dataset to confirm that the QC pipeline runs successfully.
10. Add package version capture for reproducibility.
