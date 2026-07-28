# Xenopus tropicalis gastrula snRNA-seq annotation

Quality control, dimensionality reduction, clustering, marker-guided cell-type
annotation, label refinement, and differential gene-expression analysis for a
single-nucleus RNA-sequencing time course spanning _Xenopus tropicalis_
gastrulation.

## Dataset

The repository contains seven stage-condition analysis notebooks representing
five developmental time points:

| Dataset       |   Input barcodes |   Genes |   ≥1,000 genes |   MT/ribo filters |   Predicted doublets |   Retained singlets |   Final annotated |
|:--------------|-----------------:|--------:|---------------:|------------------:|---------------------:|--------------------:|------------------:|
| Stage10 WT    |           215209 |   21821 |          15920 |             15920 |                 1820 |               14100 |             14100 |
| Stage11.25 WT |           110592 |   21821 |           7971 |              7971 |                 1446 |                6525 |              6516 |
| Stage11.25 MO |           110590 |   21821 |           6533 |              6533 |                  503 |                6030 |              5995 |
| Stage12.0 WT  |           110592 |   21821 |          14754 |             14754 |                    0 |               14754 |             10789 |
| Stage12.0 MO  |           110590 |   21821 |          13967 |             13967 |                 1369 |               12598 |              8598 |
| Stage12.5 WT  |           110592 |   29231 |          11908 |             11908 |                 1608 |               10300 |             10237 |
| Stage13 WT    |           110584 |   29231 |          10074 |             10074 |                 1166 |                8908 |              8837 |

The executed notebook outputs currently report **73,215 retained
singlet nuclei** after Scrublet filtering and **65,072 nuclei in the
final annotation tables**.

Libraries were prepared using PARSE Evercode WT-mini and sequenced on an
Illumina NovaSeq 6000, according to the project metadata supplied with this
repository.

## Workflow

1. **Preprocessing and quality control**
   - retain nuclei expressing at least 1,000 genes;
   - calculate mitochondrial and ribosomal count fractions;
   - retain cells with `pct_counts_mt < 20` and `pct_counts_ribo < 3.5`;
   - normalize to 10,000 counts per cell and log-transform;
   - identify highly variable genes;
   - predict and remove doublets with Scrublet.

2. **Dimensionality reduction**
   - principal-component analysis;
   - exploratory t-SNE projections;
   - neighbor-graph construction and UMAP visualization.

3. **Clustering and annotation**
   - Leiden clustering across multiple resolutions;
   - marker-panel visualization;
   - germ-layer subsetting;
   - threshold-based assignment of region and cell-type labels;
   - nearest-neighbor filling of unassigned cells and boundary refinement;
   - transfer of subcluster annotations back to the full AnnData object.

4. **Differential expression**
   - Wilcoxon tests for each final annotation versus the remaining cells;
   - selected reciprocal two-group comparisons;
   - ranked-gene plots, heatmaps, dotplots, and matrix plots;
   - CSV export of complete, significant, and top-ranked DEG results.

Regression and batch-correction examples remain in the notebooks as commented
exploratory code; they are **not active steps in the completed workflow**.

## Repository layout

```text
.
├── README.md
├── config/
│   └── run_config.yaml
├── data/
│   ├── raw/                 # local input matrices; not committed
│   └── metadata/
├── docs/
│   ├── AnalysisPlan.md
│   ├── Notebook_Audit.md
│   ├── Pipeline.md
│   ├── QC_Report.md
│   ├── Reproducibility.md
│   ├── Section5_DEG_and_Table_Exports.md
│   └── Validation_Report.md
├── notebooks/               # cleaned, GitHub-ready notebooks
├── results/                 # generated outputs; not committed by default
├── scripts/
│   └── output_utils.py
├── tables/                  # audit tables derived from uploaded outputs
├── environment.yml
└── requirements.txt
```

## Installation

```bash
conda env create -f environment.yml
conda activate xtrop-gastrula-snrnaseq
jupyter lab
```

Place the unfiltered `.h5ad` inputs in the stage-specific directories described
in `data/README.md`, then run each notebook from top to bottom.

## Generated outputs

Each dataset writes to a stage- and condition-specific results directory:

```text
results/<stage>[/<condition>]/
├── qc/
├── figures/
│   ├── Section1-Preprocessing/
│   ├── Section2-DimensionReduction/
│   ├── Section3-MarkerGeneExpression/
│   ├── Section4-ClusteringAnnotation/
│   ├── Section5-DifferentialExpression/
│   └── Section6-FinalAnnotation/
├── processed_data/
├── differential_expression/
├── validation/
└── tables/
```

Important generated files include:

- final annotated AnnData objects;
- complete cell metadata;
- annotation count and percentage tables;
- Leiden-to-annotation cross-tabulations;
- label-refinement transition tables;
- complete and filtered DEG tables;
- DEG comparison inventories and parameter records;
- descriptively named figures without notebook cell numbers.

## Documentation

- [`docs/AnalysisPlan.md`](docs/AnalysisPlan.md): completed analysis aims,
  units, comparisons, and deliverables.
- [`docs/Pipeline.md`](docs/Pipeline.md): implemented workflow and parameters.
- [`docs/QC_Report.md`](docs/QC_Report.md): QC thresholds and notebook-derived
  retention statistics.
- [`docs/Validation_Report.md`](docs/Validation_Report.md): annotation evidence,
  validation strategy, limitations, and recommended checks.
- [`docs/Reproducibility.md`](docs/Reproducibility.md): environment, paths,
  execution, and version-control guidance.
- [`docs/Notebook_Audit.md`](docs/Notebook_Audit.md): issues found and changes
  made to the GitHub-ready notebook copies.
- [`docs/Section5_DEG_and_Table_Exports.md`](docs/Section5_DEG_and_Table_Exports.md):
  exact DEG and table outputs produced by the replacement Section 5.

## Data and version-control policy

Raw matrices and generated `.h5ad` files are intentionally excluded from Git.
The repository tracks code, configuration, documentation, and lightweight audit
tables. Selected publication-quality figures can be copied to a dedicated
tracked directory if needed.

## Citation

Add the associated manuscript citation, data accession, and a repository DOI
after archival release. The Stage 10/10.5 notebook also references the published
dataset DOI `10.1371/journal.pbio.3003476`.
