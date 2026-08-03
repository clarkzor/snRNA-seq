# Validation Report

## Validation strategy implemented in the notebooks

Cell identities are supported through several complementary checks rather than
a single automated classifier.

### Marker-expression validation

The notebooks visualize broad germ-layer and lineage markers on UMAPs before
and during annotation. Fine annotations are assigned with explicit
gene-expression thresholds, including logical combinations of markers and
rules that protect previously assigned identities.

### Doublet validation

Scrublet predictions are visualized with total counts and marker panels.
Predicted doublets are removed after reviewing whether cells co-express
incompatible lineage markers.

### Clustering validation

Multiple Leiden resolutions are compared. Broad clusters are subset into
lineage-specific AnnData objects, where marker thresholds and local structure
are used for finer annotation.

### Local-neighborhood validation

Unassigned cells are filled from nearby labeled cells using KNN procedures.
Several notebooks preserve the original labels and record which cells changed,
allowing transition-table review.

### Differential-expression validation

Section 5 evaluates:

- all final annotations versus the remaining cells;
- layer-specific or region-specific pairs;
- ranked genes, heatmaps, dotplots, and matrix plots.

The cleaned notebooks save the underlying statistical tables rather than
retaining only visual summaries.

### Literature-based interpretation

Section 6 interprets inner/outer tissue-layer distinctions using differential
expression and published _Xenopus_ developmental literature. Temporary labels
such as `Epidermis 1/2`, `Neural Plate 1/2`, and `Neural Plate Border 1/2` are
renamed to outer/inner identities in notebooks where this interpretation is
applied.

## Recommended validation tables

### Essential

1. **Final annotation abundance**  
   Cell count and percent for every final annotation.

2. **Leiden-to-annotation cross-tabulation**  
   Demonstrates whether labels correspond to coherent clusters and reveals
   annotations distributed across multiple clusters.

3. **Marker summary by annotation**  
   For each marker and annotation, save mean expression, median expression,
   percent expressing, and cell count.

4. **Complete DEG table**  
   Save every tested gene, score, log fold change, raw p-value, adjusted
   p-value, and fraction expressing.

5. **Label-transition table**  
   Crosstab of labels before and after KNN filling or boundary smoothing.

### Useful additional tables

- annotation counts before and after removal of `Unknown` or `Other`;
- marker-threshold assignment counts;
- protected-label conflict counts;
- pairwise comparison sample sizes;
- top 25 markers per annotation;
- genes significant in both reciprocal comparison directions;
- overlap of top markers across stages.

## Repository-level notebook validation

Notebook-file validation is performed independently of biological validation:

```bash
python scripts/validate_notebooks.py
```

The command validates every file in `notebooks/` with `nbformat.validate` and
performs static Python syntax compilation of each code cell after applying
IPython input transformations. Reports are written to the top-level
`validation/` directory as CSV, JSON, and readable text files.

A `PASS` result confirms that a notebook has valid Jupyter structure and that
its code cells are syntactically compilable. It does **not** confirm that the
notebook executes successfully, that required `.h5ad` files are available, or
that the resulting biological analysis is correct. Those checks require a full
fresh-kernel execution with the intended data and software environment.

## Interpretation cautions

1. **Threshold-based labels are parameter dependent.** Thresholds differ across
   genes and stages and should be recorded in a machine-readable table.
2. **UMAP-space KNN is a visualization-space heuristic.** It can improve label
   continuity but may blur biological boundaries. Results should be compared
   with PCA- or graph-based neighbors.
3. **Boundary smoothing can propagate common labels.** The stored transition
   tables should be inspected for loss of rare populations.
4. **Small annotations produce unstable DEG estimates.** The revised workflow
   excludes groups with fewer than seven cells from the all-group test and
   records exclusions.
5. **One-versus-rest is not always a pure pairwise comparison.** The revised
   pairwise workflow subsets to exactly two groups before testing.
6. **Final renaming is notebook-specific.** Analyses performed before Section 6
   use temporary labels; result filenames and tables should state whether they
   use pre-renaming or final labels.
7. **Stage and condition comparisons are not yet formalized.** WT and MO
   notebooks are currently analyzed independently. Cross-condition inference
   requires a separate model that accounts for biological replication and
   batch structure.

## Validation completion criteria

- canonical markers are enriched in the expected annotation;
- no final annotation is defined by only one weakly expressed gene without
  supporting context;
- cluster-to-annotation tables show interpretable correspondence;
- smoothing transitions do not eliminate rare biologically plausible groups;
- DEG tables contain adjusted p-values and sample-size metadata;
- final labels, figures, and saved tables use consistent nomenclature.
