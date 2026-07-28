# Pipeline

## 1. Project paths and output routing

Each notebook defines `STAGE`, `CONDITION`, and `DATASET_SLUG`. Inputs are read
from `data/raw/<stage-folder>/`. Outputs are written beneath
`results/<stage-folder>/`, with WT/MO subdirectories for Stage 11.25 and Stage
12.0.

Figures are routed by notebook section. The cleaned notebook copies use
descriptive suffixes such as:

```text
S04_12_umap_region_annotation.png
S05_03_Stage12_0-WT_neural_plate_1_vs_2_heatmap.png
```

This replaces filenames based on mutable notebook cell numbers.

## 2. Input filtering and QC

The notebooks retain cells with at least 1,000 detected genes:

```python
sc.pp.filter_cells(adata, min_genes=1000)
```

Mitochondrial and ribosomal genes are identified and included in Scanpy QC
metrics. Cells are retained with:

```python
adata = adata[adata.obs["pct_counts_mt"] < 20, :]
adata = adata[adata.obs["pct_counts_ribo"] < 3.5, :]
```

These are the thresholds implemented by the code. The biological rationale and
the code comments should be kept synchronized if the ribosomal criterion is
changed.

## 3. Doublet detection

Scrublet is run on the matrix stored in `adata.raw.X`. The notebooks save:

- `doublet_scores`
- `predicted_doublets`
- string-form `doublet_info`

Predicted doublets are visualized against marker expression and total counts,
then excluded:

```python
adata = adata[adata.obs["doublet_info"] == "False", :]
```

## 4. Normalization and highly variable genes

The active normalization code uses 10,000 counts per cell followed by log
transformation:

```python
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
```

Most notebooks call highly variable genes with:

```python
sc.pp.highly_variable_genes(
    adata,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.05,
)
```

The Stage 13 notebook uses broader thresholds:

```python
sc.pp.highly_variable_genes(
    adata,
    min_mean=0.005,
    max_mean=5,
    min_disp=0.005,
)
```

Regression, scaling, MNN correction, and batch-aware HVG selection appear as
commented examples. They are not active in the current completed analyses.

## 5. Dimensionality reduction

The common PCA call is:

```python
sc.tl.pca(adata, svd_solver="arpack", n_comps=100)
```

The notebooks compare t-SNE embeddings based on 30 and 99 PCs. Neighbor graphs
and UMAPs are explored over several parameter combinations, including 10, 20,
25, 50, 70, or 99 PCs and approximately 25–75 neighbors depending on the
dataset.

Because multiple exploratory graphs are generated, the graph used for final
clustering should be clearly identified in the notebook narrative.

## 6. Clustering and annotation

Leiden clustering is evaluated at multiple resolutions. Common resolutions
include 0.05, 0.1, 0.2, 0.3, 0.8, 1.0, 1.2, and 10.0, with stage-specific
variation.

Annotation then proceeds through germ-layer subsets and marker-expression
rules. The notebooks include labels for:

- ectodermal and neural populations;
- epidermal, ionocyte, goblet, ciliated, placodal, neural-crest, and neural
  plate/border populations;
- mesodermal subdivisions, including paraxial, pharyngeal, pronephric,
  ventral-blood-island, notochord, and organizer-related populations;
- endodermal and organ-primordium populations.

Marker thresholds are intentionally stage- and lineage-specific. Assignment
rules explicitly determine whether a new label:

- overwrites every existing label;
- applies only to `Other`;
- applies only to a named set of existing labels; or
- protects selected labels.

## 7. KNN label refinement

Several subobjects use existing UMAP coordinates to fill `Other` labels from
nearby annotated cells. The implemented variants include:

- five nearest neighbors;
- inverse-distance weighting;
- iterative assignment;
- confidence thresholds;
- optional minimum labeled-neighbor requirements;
- preservation of pre-refinement labels.

Some notebooks also perform multi-iteration boundary smoothing across all
labels. Transition tables should be reviewed because aggressive smoothing can
erase biologically meaningful boundaries.

## 8. Differential expression

The GitHub-ready notebooks replace the prior Section 5 with a consistent
implementation.

### All annotations versus rest

Annotations with at least seven cells are tested with:

```python
sc.tl.rank_genes_groups(
    adata_deg,
    groupby="region_annotation",
    groups=valid_groups,
    reference="rest",
    method="wilcoxon",
    pts=True,
    key_added="all_annotations_vs_rest_wilcoxon",
)
```

### Pairwise comparisons

Each pair is subset to exactly two groups and tested under one key with both
groups included. In a two-group subset, each group versus `rest` is the
reciprocal comparison.

### Exported statistics

For every DEG key, the helper module saves:

- gene name;
- tested group;
- score;
- log fold change;
- raw p-value;
- adjusted p-value;
- percentage expressing in the group and reference when available.

## 9. Final outputs

The final notebooks save:

```text
<DATASET_SLUG>_processed_annotated.h5ad
<DATASET_SLUG>_cell_metadata.csv
<DATASET_SLUG>_final_annotation_counts.csv
```

They also save crosstabs for every retained Leiden column and transition tables
for any `region_annotation_before*` backup columns.
