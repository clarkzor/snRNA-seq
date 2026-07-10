# Validation Report

## Project

Single-nucleus RNA-seq preprocessing, clustering, and annotation of *Xenopus tropicalis* gastrula-stage embryos.

## Purpose

The purpose of this validation report is to demonstrate that the single-cell analysis workflow performs the expected preprocessing, quality control, clustering, annotation, and export steps in a reproducible and verifiable manner.

This report validates the computational behavior of the analysis notebooks. It does not represent clinical, diagnostic, or regulatory validation.

---

# Scope

This validation report covers the following notebooks:

```text
SingleCell_BasicPreprocessingAnnotation_Stage10-EarlyGastrula.ipynb
SingleCell_ClusterAnnotationAnalysis_Stage13-LateGastrula.ipynb
```

The Stage 10.5 notebook includes preprocessing, QC, doublet detection, normalization, dimensionality reduction, clustering, annotation, and differential expression.

The Stage 13 notebook starts from a pre-existing AnnData object and focuses on cluster annotation, subsetting, label refinement, and export.

---

# Validation Objectives

The validation objectives are to confirm that:

1. Input AnnData objects load successfully.
2. Expected cell and gene counts are observed at key workflow checkpoints.
3. QC metrics are calculated correctly.
4. Doublet detection generates expected output columns and cell counts.
5. Normalization and log transformation complete successfully.
6. Highly variable gene selection generates the expected feature set.
7. PCA, neighborhood graph construction, UMAP, and Leiden clustering complete successfully.
8. Marker genes used for annotation are present in the dataset.
9. Differential expression analysis completes and generates expected result fields.
10. Final cell type annotations are stored in the AnnData object.
11. Final processed datasets are exported successfully.

---

# Validation Approach

Validation was performed using expected-output checks.

For each major workflow step, the following were assessed:

* Required input files exist.
* Required AnnData fields are generated.
* Expected cell and gene counts are observed.
* Required metadata columns are present.
* Required embeddings are present.
* Required clustering and annotation columns are present.
* Final output files are generated and can be reloaded.

Each validation test is assigned one of the following statuses:

| Status        | Meaning                                                 |
| ------------- | ------------------------------------------------------- |
| Pass          | Test completed successfully and met acceptance criteria |
| Fail          | Test did not meet acceptance criteria                   |
| Not assessed  | Test was not evaluated in the current notebook          |
| Manual review | Test requires visual or biological interpretation       |

---

# Input Files

## Stage 10.5

Input file:

```text
PARSE-WTSt105-Unfiltered_FromSeurat.h5ad
```

Expected starting object:

```text
215,209 cells × 21,821 genes
```

## Stage 13

Input file:

```text
PARSE-Unfiltered_St13_Compatible.h5ad
```

Expected final output file:

```text
Stage13_WT_Annotated.h5ad
```

---

# Stage 10.5 Validation Tests

## Test 1: Input AnnData Object Loads Successfully

| Field                | Value                                  |
| -------------------- | -------------------------------------- |
| Test ID              | ST10_VAL_001                           |
| Step                 | Load input AnnData object              |
| Acceptance criterion | Input `.h5ad` file loads without error |
| Expected result      | AnnData object is created              |
| Observed result      | AnnData object loaded                  |
| Status               | Pass                                   |

Expected object shape:

```text
215,209 cells × 21,821 genes
```

---

## Test 2: Minimum Gene Filtering

| Field                | Value                                                  |
| -------------------- | ------------------------------------------------------ |
| Test ID              | ST10_VAL_002                                           |
| Step                 | Filter low-quality cells                               |
| Method               | `sc.pp.filter_cells(adataR, min_genes=1000)`           |
| Acceptance criterion | Cells with fewer than 1,000 detected genes are removed |
| Expected result      | 15,920 cells retained                                  |
| Observed result      | 15,920 cells retained                                  |
| Status               | Pass                                                   |

Result:

```text
Before filtering: 215,209 cells × 21,821 genes
After filtering:   15,920 cells × 21,821 genes
```

---

## Test 3: QC Metric Calculation

| Field                | Value                                                                   |
| -------------------- | ----------------------------------------------------------------------- |
| Test ID              | ST10_VAL_003                                                            |
| Step                 | Calculate QC metrics                                                    |
| Method               | `sc.pp.calculate_qc_metrics()`                                          |
| Acceptance criterion | QC metric columns are added to `adata.obs`                              |
| Expected columns     | `n_genes_by_counts`, `total_counts`, `pct_counts_mt`, `pct_counts_ribo` |
| Observed result      | QC metrics calculated                                                   |
| Status               | Pass                                                                    |

Mitochondrial and ribosomal metrics were calculated but not used for filtering in the uploaded notebook.

---

## Test 4: Doublet Detection

| Field                | Value                                                     |
| -------------------- | --------------------------------------------------------- |
| Test ID              | ST10_VAL_004                                              |
| Step                 | Predict doublets                                          |
| Method               | Scrublet                                                  |
| Acceptance criterion | Doublet score and predicted-doublet columns are generated |
| Expected columns     | `doublet_scores`, `predicted_doublets`, `doublet_info`    |
| Observed result      | Scrublet columns generated                                |
| Status               | Pass                                                      |

Scrublet summary:

| Metric                           | Observed value |
| -------------------------------- | -------------: |
| Automatically selected threshold |           0.20 |
| Detected doublet rate            |          11.4% |
| Expected doublet rate            |          10.0% |
| Estimated overall doublet rate   |          16.0% |
| Predicted doublets               |          1,820 |

---

## Test 5: Doublet Removal

| Field                | Value                                                  |
| -------------------- | ------------------------------------------------------ |
| Test ID              | ST10_VAL_005                                           |
| Step                 | Remove predicted doublets                              |
| Method               | Remove cells where `doublet_info == "True"`            |
| Acceptance criterion | Predicted doublets are excluded from downstream object |
| Expected result      | 14,100 cells retained                                  |
| Observed result      | 14,100 cells retained                                  |
| Status               | Pass                                                   |

Result:

```text
Before doublet removal: 15,920 cells
Predicted doublets:      1,820 cells
After doublet removal:  14,100 cells
```

---

## Test 6: Normalization and Log Transformation

| Field                | Value                                                         |
| -------------------- | ------------------------------------------------------------- |
| Test ID              | ST10_VAL_006                                                  |
| Step                 | Normalize and log-transform expression                        |
| Method               | `sc.pp.normalize_per_cell()` followed by `sc.pp.log1p()`      |
| Acceptance criterion | Normalized and log-transformed expression matrix is generated |
| Expected result      | Expression matrix is transformed for downstream analysis      |
| Observed result      | Transformation completed                                      |
| Status               | Pass                                                          |

Normalization target:

```text
10,000 counts per cell
```

---

## Test 7: Highly Variable Gene Selection

| Field                | Value                                |
| -------------------- | ------------------------------------ |
| Test ID              | ST10_VAL_007                         |
| Step                 | Select highly variable genes         |
| Method               | `sc.pp.highly_variable_genes()`      |
| Acceptance criterion | Highly variable genes are identified |
| Expected result      | 6,439 highly variable genes          |
| Observed result      | 6,439 highly variable genes          |
| Status               | Pass                                 |

Parameters:

| Parameter  |  Value |
| ---------- | -----: |
| `min_mean` | 0.0125 |
| `max_mean` |      3 |
| `min_disp` |   0.05 |

---

## Test 8: PCA

| Field                | Value                                               |
| -------------------- | --------------------------------------------------- |
| Test ID              | ST10_VAL_008                                        |
| Step                 | Principal Component Analysis                        |
| Method               | `sc.tl.pca()`                                       |
| Acceptance criterion | PCA coordinates are stored in `adata.obsm["X_pca"]` |
| Expected result      | PCA completes with 100 components                   |
| Observed result      | PCA completed                                       |
| Status               | Pass                                                |

Parameter:

```text
n_comps = 100
```

---

## Test 9: Neighborhood Graph and UMAP

| Field                | Value                                             |
| -------------------- | ------------------------------------------------- |
| Test ID              | ST10_VAL_009                                      |
| Step                 | Build neighborhood graph and UMAP                 |
| Method               | `sc.pp.neighbors()` and `sc.tl.umap()`            |
| Acceptance criterion | Neighbor graph and UMAP coordinates are generated |
| Expected result      | `adata.obsm["X_umap"]` exists                     |
| Observed result      | UMAP generated                                    |
| Status               | Pass                                              |

Primary parameters used in the notebook include:

```text
n_neighbors = 50
n_pcs = 99
```

---

## Test 10: Leiden Clustering

| Field                | Value                                          |
| -------------------- | ---------------------------------------------- |
| Test ID              | ST10_VAL_010                                   |
| Step                 | Leiden clustering                              |
| Method               | `sc.tl.leiden()`                               |
| Acceptance criterion | Leiden cluster labels are added to `adata.obs` |
| Expected result      | Clustering columns generated                   |
| Observed result      | Leiden columns generated                       |
| Status               | Pass                                           |

Leiden resolutions tested include:

```text
0.1, 0.3, 0.8, 1.0, 1.2, 10.0
```

---

## Test 11: Marker Gene Availability

| Field                | Value                                                    |
| -------------------- | -------------------------------------------------------- |
| Test ID              | ST10_VAL_011                                             |
| Step                 | Marker-based annotation                                  |
| Acceptance criterion | Expected marker genes are present in `adata.var_names`   |
| Expected result      | Major germ-layer markers are available for visualization |
| Observed result      | Marker genes visualized                                  |
| Status               | Manual review                                            |

Marker genes used include:

| Category      | Marker genes           |
| ------------- | ---------------------- |
| Epidermal     | `krt7`, `grhl1`        |
| Neuroectoderm | `sox2`, `sox3`, `pax3` |
| Mesodermal    | `tbxt`, `fgf8`         |
| Endodermal    | `darmin`, `gata6`      |

This test requires biological review because marker-based annotation depends on interpretation of expression patterns.

---

## Test 12: Differential Expression

| Field                | Value                                                    |
| -------------------- | -------------------------------------------------------- |
| Test ID              | ST10_VAL_012                                             |
| Step                 | Differential expression                                  |
| Method               | `sc.tl.rank_genes_groups()`                              |
| Acceptance criterion | Ranked gene results are generated                        |
| Expected result      | Differential expression fields are stored in `adata.uns` |
| Observed result      | Differential expression completed                        |
| Status               | Pass                                                     |

Method used:

```text
t-test
```

Expected result fields include:

```text
names
scores
pvals
pvals_adj
logfoldchanges
```

---

# Stage 13 Validation Tests

## Test 13: Stage 13 Input Object Loads Successfully

| Field                | Value                                  |
| -------------------- | -------------------------------------- |
| Test ID              | ST13_VAL_001                           |
| Step                 | Load Stage 13 AnnData object           |
| Acceptance criterion | Input `.h5ad` file loads without error |
| Expected result      | AnnData object is created              |
| Observed result      | AnnData object loaded                  |
| Status               | Pass                                   |

Input file:

```text
PARSE-Unfiltered_St13_Compatible.h5ad
```

---

## Test 14: Stage 13 Leiden Clustering

| Field                | Value                                                |
| -------------------- | ---------------------------------------------------- |
| Test ID              | ST13_VAL_002                                         |
| Step                 | Leiden clustering                                    |
| Method               | `sc.tl.leiden()`                                     |
| Acceptance criterion | Leiden columns are generated at multiple resolutions |
| Expected result      | Cluster labels are added to `.obs`                   |
| Observed result      | Leiden clustering completed                          |
| Status               | Pass                                                 |

Leiden resolutions tested include:

```text
0.2, 0.3, 0.4, 1.0, 1.2, 1.4
```

---

## Test 15: Inner Ectoderm Annotation

| Field                | Value                                                            |
| -------------------- | ---------------------------------------------------------------- |
| Test ID              | ST13_VAL_003                                                     |
| Step                 | Annotate inner ectoderm subset                                   |
| Acceptance criterion | Inner ectoderm subset receives biologically interpretable labels |
| Expected result      | Annotation column generated                                      |
| Observed result      | Inner ectoderm annotations assigned                              |
| Status               | Manual review                                                    |

Representative marker genes:

| Annotation            | Marker genes                    |
| --------------------- | ------------------------------- |
| Anterior Neural Plate | `rax`, `pax6`                   |
| Neural Plate          | `sox2`, `sox3`                  |
| Neural Plate Border   | `pax3`                          |
| Neural Crest          | `snai2`, `sox9`                 |
| PPR                   | `six1`, `eya1`, `pax8`, `pitx1` |
| Notoplate             | `shh`, `ptch2`                  |

---

## Test 16: Outer Ectoderm and Endoderm Annotation

| Field                | Value                                             |
| -------------------- | ------------------------------------------------- |
| Test ID              | ST13_VAL_004                                      |
| Step                 | Annotate outer ectoderm and endoderm subset       |
| Acceptance criterion | Subset receives biologically interpretable labels |
| Expected result      | Annotation column generated                       |
| Observed result      | Outer ectoderm/endoderm annotations assigned      |
| Status               | Manual review                                     |

Representative marker genes:

| Annotation               | Marker genes     |
| ------------------------ | ---------------- |
| Cement Gland             | `agr2`, `pitx1`  |
| Goblet Cell              | `itln1`, `itln2` |
| Outer NNE                | `krt70`          |
| Ciliated Epidermal Prog. | `foxj1`, `tp73`  |
| Ionocyte                 | `foxi1`          |
| Foregut Primordium       | `slc5a8`         |
| Midgut Primordium        | `ndrg1`          |
| Hindgut Primordium       | `gatm`           |
| Liver Primordium         | `hhex`           |
| Pancreatic Primordium    | `pdx1`           |

---

## Test 17: Mesoderm Annotation

| Field                | Value                                                      |
| -------------------- | ---------------------------------------------------------- |
| Test ID              | ST13_VAL_005                                               |
| Step                 | Annotate mesoderm subset                                   |
| Acceptance criterion | Mesoderm subset receives biologically interpretable labels |
| Expected result      | Annotation column generated                                |
| Observed result      | Mesoderm annotations assigned                              |
| Status               | Manual review                                              |

Representative marker genes:

| Annotation               | Marker genes                     |
| ------------------------ | -------------------------------- |
| Pronephric Primordium    | `pax8`, `lhx1`                   |
| Paraxial Mesoderm        | `ank1`, `actc1`, `myf5`, `hoxd3` |
| Pharyngeal Mesoderm      | `tbx1`                           |
| Cardiac Mesoderm         | `nkx2-5`                         |
| Ventral Blood Island     | `runx1`                          |
| Ventral Mesoderm         | `gata6`                          |
| Notochord                | `not`                            |
| Lateral Mesoderm         | `foxc1`, `tbxt`                  |
| L/R Organizer Primordium | `cirop`                          |

---

## Test 18: Final Annotation Column

| Field                | Value                                                     |
| -------------------- | --------------------------------------------------------- |
| Test ID              | ST13_VAL_006                                              |
| Step                 | Combine subset annotations                                |
| Acceptance criterion | Final annotation column is present in full AnnData object |
| Expected column      | `refined_region_annotation`                               |
| Observed result      | Final annotation column generated                         |
| Status               | Pass                                                      |

Cells without a refined annotation are labeled:

```text
Other
```

---

## Test 19: Final Export

| Field                | Value                                             |
| -------------------- | ------------------------------------------------- |
| Test ID              | ST13_VAL_007                                      |
| Step                 | Export final annotated object                     |
| Acceptance criterion | Final `.h5ad` file is written and can be reloaded |
| Expected output      | `Stage13_WT_Annotated.h5ad`                       |
| Observed result      | Output file generated                             |
| Status               | Pass                                              |

---

# Validation Summary

| Category                                       | Number of tests | Status        |
| ---------------------------------------------- | --------------: | ------------- |
| Stage 10.5 preprocessing/QC                    |               7 | Pass          |
| Stage 10.5 dimensionality reduction/clustering |               3 | Pass          |
| Stage 10.5 annotation / marker review          |               1 | Manual review |
| Stage 10.5 differential expression             |               1 | Pass          |
| Stage 13 loading / clustering / export         |               3 | Pass          |
| Stage 13 biological annotation                 |               3 | Manual review |

Overall validation status:

```text
Partially validated for computational reproducibility.
Manual biological review required for marker-based cell type annotation.
```

---

# Acceptance Criteria Summary

The pipeline is considered successfully validated if:

1. Input AnnData files load without error.
2. Stage 10.5 minimum-gene filtering retains approximately 15,920 cells.
3. Scrublet identifies and removes approximately 1,820 predicted doublets.
4. Final Stage 10.5 object contains approximately 14,100 cells.
5. QC metric columns are generated.
6. Highly variable gene selection identifies approximately 6,439 genes.
7. PCA, neighbor graph construction, UMAP, and Leiden clustering complete successfully.
8. Marker genes used for annotation are present and visually consistent with expected biology.
9. Differential expression analysis completes and returns adjusted p-values and log fold changes.
10. Stage 13 final annotation column `refined_region_annotation` is generated.
11. Final Stage 13 object is exported as `Stage13_WT_Annotated.h5ad`.

---

# Known Limitations

1. This validation does not confirm biological truth of every annotation.
2. Annotation validation relies on known marker-gene expression and manual review.
3. Random seeds were not consistently defined in the uploaded notebooks.
4. Exact UMAP layouts may vary between runs.
5. Stage 13 upstream preprocessing is not validated because the notebook starts from an existing AnnData object.
6. The notebooks use local file paths that should be replaced with relative project paths.
7. Some outputs are visualized interactively but not automatically exported.
8. Differential expression results are generated but should be exported as CSV files for auditability.
9. No independent benchmark dataset was used for external validation.

---

# Recommended Improvements

Future validation should include:

1. A small test dataset that can be included in the GitHub repository.
2. A script that automatically runs validation checks.
3. A `validation_summary.csv` file with pass/fail results.
4. Fixed random seeds for stochastic steps.
5. Automatic export of QC, clustering, marker, and differential expression outputs.
6. A marker-gene validation table for each final annotated cell type.
7. Unit tests confirming required `.obs`, `.var`, `.obsm`, and `.uns` fields exist.
8. Re-loading tests for exported `.h5ad` files.
9. Documentation of software versions and input file checksums.
10. Optional comparison to a published or previously validated reference annotation.

---

# Conclusion

The uploaded notebooks contain a functional single-nucleus RNA-seq analysis workflow for *Xenopus tropicalis* gastrula-stage embryos.

The Stage 10.5 workflow is validated for basic computational execution, including input loading, cell filtering, QC metric calculation, doublet detection, doublet removal, normalization, highly variable gene selection, dimensionality reduction, clustering, and differential expression.

The Stage 13 workflow is validated for cluster annotation, subset-based label refinement, generation of a final annotation column, and export of the annotated AnnData object.

The workflow should be considered computationally validated for exploratory research use, with manual biological review required for final cell type annotation.
