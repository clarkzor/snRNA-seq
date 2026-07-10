# Reproducibility

## Project

Single-nucleus RNA-seq preprocessing, clustering, and annotation of *Xenopus tropicalis* gastrula-stage embryos.

This document describes how to reproduce the analyses performed in the Stage 10.5 early gastrula and Stage 13 late gastrula notebooks.

---

# Reproducibility Status

The current workflow is partially reproducible.

The computational preprocessing, quality control, normalization, dimensionality reduction, clustering, and differential expression steps can be reproduced if the required input `.h5ad` files and Python environment are available.

However, several aspects of the current workflow require manual review or future improvement:

* Input files are referenced using local absolute Windows paths.
* Some cell type annotations are assigned manually based on marker-gene inspection.
* Random seeds are not explicitly set in the uploaded notebooks.
* Most figures are displayed interactively rather than saved automatically.
* Some intermediate QC and differential expression tables are not exported as structured files.
* The Stage 13 notebook begins from a pre-existing AnnData object and does not show the full upstream preprocessing workflow.

Future versions of this workflow should use relative paths, a configuration file, fixed random seeds, automatic output export, and a documented environment file.

---

# Input Data

## Stage 10.5 Early Gastrula

Input file used in the notebook:

```text
PARSE-WTSt105-Unfiltered_FromSeurat.h5ad
```

Original local path in notebook:

```text
C:/Users/coron/OneDrive/Desktop/snRNAseq/St10/PARSE-WTSt105-Unfiltered_FromSeurat.h5ad
```

Initial dataset size:

```text
215,209 cells × 21,821 genes
```

After filtering cells with fewer than 1,000 detected genes:

```text
15,920 cells × 21,821 genes
```

After Scrublet doublet removal:

```text
14,100 cells × 21,821 genes
```

---

## Stage 13 Late Gastrula

Input file used in the notebook:

```text
PARSE-Unfiltered_St13_Compatible.h5ad
```

Original local path in notebook:

```text
C:/Users/coron/OneDrive/Desktop/snRNAseq/St13/X__tropicalis_st13_NUCLEI/DGE_unfiltered/PARSE-Unfiltered_St13_Compatible.h5ad
```

The Stage 13 notebook starts from this existing AnnData object and performs cluster annotation, subsetting, label refinement, and export.

Final output file:

```text
Stage13_WT_Annotated.h5ad
```

---

# Recommended Repository Paths

For reproducible reruns, local absolute paths should be replaced with project-relative paths.

Recommended structure:

```text
project/
├── notebooks/
│   ├── SingleCell_BasicPreprocessingAnnotation_Stage10-EarlyGastrula.ipynb
│   └── SingleCell_ClusterAnnotationAnalysis_Stage13-LateGastrula.ipynb
├── data/
│   ├── raw/
│   └── processed/
├── results/
│   ├── qc/
│   ├── figures/
│   ├── differential_expression/
│   └── processed_data/
├── docs/
│   ├── AnalysisPlan.md
│   ├── Pipeline.md
│   ├── QC_Report.md
│   └── Reproducibility.md
├── config/
│   └── sample_config.yaml
└── environment.yml
```

Example path replacement:

```python
from pathlib import Path

PROJECT_DIR = Path.cwd().parent
DATA_DIR = PROJECT_DIR / "data"
RESULTS_DIR = PROJECT_DIR / "results"

adata = sc.read_h5ad(DATA_DIR / "raw" / "PARSE-WTSt105-Unfiltered_FromSeurat.h5ad")
```

---

# Software Environment

The analysis was performed using Python and Scanpy-based single-cell analysis packages.

Required packages include:

```text
python
scanpy
anndata
numpy
pandas
scipy
matplotlib
scrublet
seaborn
stream2
```

Exact versions should be captured from the active environment using the following code:

```python
import sys
import platform
from importlib.metadata import version, PackageNotFoundError

packages = [
    "scanpy",
    "anndata",
    "numpy",
    "pandas",
    "scipy",
    "matplotlib",
    "scrublet",
    "seaborn",
    "stream2",
]

print("Python:", sys.version)
print("Platform:", platform.platform())

for pkg in packages:
    try:
        print(f"{pkg}: {version(pkg)}")
    except PackageNotFoundError:
        print(f"{pkg}: not installed")
```

These versions should be recorded in either:

```text
environment.yml
requirements.txt
results/run_config.yaml
```

---

# Environment Re-Creation

A Conda environment file should be provided to allow others to recreate the analysis environment.

Example:

```yaml
name: xenopus_scrnaseq
channels:
  - conda-forge
  - bioconda
  - defaults

dependencies:
  - python=3.11
  - scanpy
  - anndata
  - numpy
  - pandas
  - scipy
  - matplotlib
  - seaborn
  - pip
  - pip:
      - scrublet
      - stream2
```

To recreate the environment:

```bash
conda env create -f environment.yml
conda activate xenopus_scrnaseq
```

---

# Random Seeds

The uploaded notebooks do not explicitly define random seeds.

For improved reproducibility, future versions should set seeds before PCA, tSNE, UMAP, Leiden clustering, and any other stochastic steps.

Recommended seed setup:

```python
import random
import numpy as np
import scanpy as sc

RANDOM_SEED = 42

random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)
sc.settings.seed = RANDOM_SEED
```

Where supported, pass the seed directly:

```python
sc.tl.umap(adata, random_state=RANDOM_SEED)
sc.tl.leiden(adata, random_state=RANDOM_SEED)
sc.tl.tsne(adata, random_state=RANDOM_SEED)
```

---

# Key Analysis Parameters

## Stage 10.5

| Step                    |              Parameter | Value                         |
| ----------------------- | ---------------------: | ----------------------------- |
| Cell filtering          | Minimum genes per cell | 1,000                         |
| Normalization           |        Counts per cell | 10,000                        |
| Highly variable genes   |             `min_mean` | 0.0125                        |
| Highly variable genes   |             `max_mean` | 3                             |
| Highly variable genes   |             `min_disp` | 0.05                          |
| PCA                     |   Number of components | 100                           |
| Neighbors               |          `n_neighbors` | 50                            |
| Neighbors               |         `n_pcs` tested | 10, 20, 50, 99                |
| Main neighbors setting  |                `n_pcs` | 99                            |
| tSNE                    |         `n_pcs` tested | 30, 99                        |
| Leiden clustering       |     Resolutions tested | 0.1, 0.3, 0.8, 1.0, 1.2, 10.0 |
| Differential expression |          Scanpy method | `t-test`                      |

---

## Scrublet Doublet Detection

Scrublet was used to identify predicted doublets.

| Scrublet metric                                | Value |
| ---------------------------------------------- | ----: |
| Automatically selected doublet score threshold |  0.20 |
| Detected doublet rate                          | 11.4% |
| Expected doublet rate                          | 10.0% |
| Estimated detectable doublet fraction          | 71.6% |
| Estimated overall doublet rate                 | 16.0% |
| Predicted doublets removed                     | 1,820 |

Doublets were removed using:

```python
adata = adata[adata.obs["doublet_info"] == "False", :]
```

---

## Stage 13

The Stage 13 notebook performs cluster annotation and refinement using an existing AnnData object.

| Step                                     | Parameter / setting                        |
| ---------------------------------------- | ------------------------------------------ |
| Leiden clustering resolutions            | 0.2, 0.3, 0.4, 1.0, 1.2, 1.4               |
| Main broad-cluster column                | `leiden_0.3`                               |
| Inner ectoderm KNN refinement            | `n_neighbors=20`, `use_rep="X_pca"`        |
| Outer ectoderm / endoderm KNN refinement | `n_neighbors=5`, `use_rep="X_pca"`         |
| Mesoderm KNN refinement                  | `n_neighbors=5` or `10`, `use_rep="X_pca"` |
| Final annotation column                  | `refined_region_annotation`                |
| Final exported object                    | `Stage13_WT_Annotated.h5ad`                |

---

# Expected Outputs

A successful rerun should generate the following outputs.

## Stage 10.5

```text
results/processed_data/
    Stage10_5_processed.h5ad

results/qc/
    qc_summary.csv
    qc_metrics.csv
    doublet_scores.csv

results/figures/
    qc_violin_plots.pdf
    scrublet_doublet_histogram.pdf
    umap_clusters.pdf
    umap_celltype_annotations.pdf
    marker_gene_umaps.pdf

results/differential_expression/
    category_rank_genes_groups.csv
    Mesoderm1_vs_Mesoderm2.csv
    Epidermis1_vs_Epidermis2.csv
```

## Stage 13

```text
results/processed_data/
    Stage13_WT_Annotated.h5ad

results/figures/
    stage13_refined_region_annotation_umap.pdf
    stage13_inner_ectoderm_umap.pdf
    stage13_outer_ectoderm_endoderm_umap.pdf
    stage13_mesoderm_umap.pdf
```

---

# Suggested Output Validation

After rerunning the workflow, confirm that the following checks pass.

## Stage 10.5

| Check                             | Expected result |
| --------------------------------- | --------------: |
| Raw object loads successfully     |            Pass |
| Initial cell count                |         215,209 |
| Initial gene count                |          21,821 |
| Cells after min-gene filtering    |          15,920 |
| Highly variable genes             |           6,439 |
| Predicted doublets removed        |           1,820 |
| Final cells after doublet removal |          14,100 |
| Final genes after doublet removal |          21,821 |
| PCA completes                     |            Pass |
| UMAP generates                    |            Pass |
| Leiden labels generated           |            Pass |
| Differential expression completes |            Pass |

## Stage 13

| Check                                              | Expected result |
| -------------------------------------------------- | --------------- |
| Input AnnData object loads successfully            | Pass            |
| Leiden clustering completes                        | Pass            |
| Inner ectoderm subset annotated                    | Pass            |
| Outer ectoderm/endoderm subset annotated           | Pass            |
| Mesoderm subset annotated                          | Pass            |
| Final `refined_region_annotation` column generated | Pass            |
| Final `.h5ad` object exported                      | Pass            |

---

# File Integrity

For full reproducibility, checksums should be recorded for all input files.

Example checksum command:

```bash
sha256sum data/raw/PARSE-WTSt105-Unfiltered_FromSeurat.h5ad
sha256sum data/raw/PARSE-Unfiltered_St13_Compatible.h5ad
```

On Windows PowerShell:

```powershell
Get-FileHash data/raw/PARSE-WTSt105-Unfiltered_FromSeurat.h5ad -Algorithm SHA256
Get-FileHash data/raw/PARSE-Unfiltered_St13_Compatible.h5ad -Algorithm SHA256
```

The resulting hashes should be recorded here:

| File                                       | SHA256      |
| ------------------------------------------ | ----------- |
| `PARSE-WTSt105-Unfiltered_FromSeurat.h5ad` | To be added |
| `PARSE-Unfiltered_St13_Compatible.h5ad`    | To be added |

---

# Running the Notebooks

The notebooks can be run interactively using JupyterLab:

```bash
jupyter lab
```

Recommended execution order:

```text
1. SingleCell_BasicPreprocessingAnnotation_Stage10-EarlyGastrula.ipynb
2. SingleCell_ClusterAnnotationAnalysis_Stage13-LateGastrula.ipynb
```

If notebooks are converted into scripts, they should be executable from the command line.

Example:

```bash
python scripts/preprocess_stage10_5.py
python scripts/annotate_stage13.py
```

---

# Reproducibility Limitations

The following limitations apply to the uploaded notebook versions:

1. The notebooks use local absolute Windows paths.
2. The full data files are not stored in the GitHub repository due to file size.
3. Random seeds are not consistently defined.
4. Several annotations depend on manual marker-gene interpretation.
5. Some plots are displayed but not saved automatically.
6. Some intermediate result tables are not exported.
7. Stage 13 preprocessing steps prior to the loaded AnnData object are not documented in the uploaded notebook.
8. The current notebooks are exploratory rather than fully parameterized pipeline scripts.

---

# Recommended Improvements

To improve reproducibility, future versions should:

1. Replace absolute file paths with relative project paths.
2. Add `environment.yml` or `requirements.txt`.
3. Add `config/sample_config.yaml`.
4. Automatically generate `results/run_config.yaml`.
5. Automatically export `qc_summary.csv`.
6. Automatically export differential expression tables.
7. Save all final figures to `results/figures/`.
8. Set explicit random seeds.
9. Record input file checksums.
10. Convert repeated code into reusable functions or scripts.
11. Add a small example dataset for testing.
12. Add a `tests/` folder with simple validation checks.





# What you should still add to your notebook at the beginning

import sys
import platform
import random
import numpy as np
import scanpy as sc
from pathlib import Path
from importlib.metadata import version, PackageNotFoundError

RANDOM_SEED = 42

random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)
sc.settings.seed = RANDOM_SEED

PROJECT_DIR = Path.cwd().parent
DATA_DIR = PROJECT_DIR / "data"
RESULTS_DIR = PROJECT_DIR / "results"

for folder in [
    RESULTS_DIR / "qc",
    RESULTS_DIR / "figures",
    RESULTS_DIR / "processed_data",
    RESULTS_DIR / "differential_expression",
]:
    folder.mkdir(parents=True, exist_ok=True)

packages = [
    "scanpy",
    "anndata",
    "numpy",
    "pandas",
    "scipy",
    "matplotlib",
    "scrublet",
    "seaborn",
    "stream2",
]

print("Python:", sys.version)
print("Platform:", platform.platform())

for pkg in packages:
    try:
        print(f"{pkg}: {version(pkg)}")
    except PackageNotFoundError:
        print(f"{pkg}: not installed")



# What you should still add to your notebook at the end

import yaml
from datetime import datetime

run_config = {
    "analysis_date": datetime.now().isoformat(),
    "random_seed": RANDOM_SEED,
    "input_file": "PARSE-WTSt105-Unfiltered_FromSeurat.h5ad",
    "parameters": {
        "min_genes": 1000,
        "normalization_counts_per_cell": 10000,
        "hvg_min_mean": 0.0125,
        "hvg_max_mean": 3,
        "hvg_min_disp": 0.05,
        "pca_n_comps": 100,
        "neighbors_n_pcs": 99,
        "neighbors_n_neighbors": 50,
        "differential_expression_method": "t-test",
    },
}

with open(RESULTS_DIR / "run_config.yaml", "w") as f:
    yaml.dump(run_config, f, sort_keys=False)

