# Reproducibility

## Software environment

The project records the following environment:

| Package | Version |
|---|---:|
| Python | 3.11.5 |
| Scanpy | 1.9.6 |
| pandas | 2.1.4 |
| SciPy | 1.11.4 |
| AnnData | 0.10.3 |
| NumPy | 1.26.2 |
| Matplotlib | 3.8.2 |
| Scrublet | 0.2.3 |
| seaborn | 0.12.2 |

Create the environment with:

```bash
conda env create -f environment.yml
conda activate xtrop-gastrula-snrnaseq
```

## Input placement

Raw matrices are not committed. Place each input in its stage directory:

```text
data/raw/Stage10/PARSE-WTSt105-Unfiltered_FromSeurat.h5ad
data/raw/Stage11_25/PARSE-WTSt1125-Unfiltered_FromSeurat.h5ad
data/raw/Stage11_25/PARSE-MOSt1125-Unfiltered_FromSeurat.h5ad
data/raw/Stage12_0/PARSE-WTSt12-Unfiltered_FromSeurat.h5ad
data/raw/Stage12_0/PARSE-MOSt12-Unfiltered_FromSeurat.h5ad
data/raw/Stage12_5/PARSE-WTSt12_5-Unfiltered_FromSeurat.h5ad
data/raw/Stage13/PARSE-WTSt13-Unfiltered_FromSeurat.h5ad
```

## Notebook execution

1. clone the repository;
2. create and activate the Conda environment;
3. place the input files;
4. start Jupyter from inside the repository;
5. open one cleaned notebook;
6. restart the kernel and run all cells in order;
7. review QC and annotation summaries;
8. commit code and lightweight tables, not large generated matrices.

The notebooks locate the repository root using project-specific markers:
`README.md`, `notebooks/`, `scripts/output_utils.py`, `data/`, and `results/`.
This prevents the `notebooks/` directory from being mistaken for the project
root. No user-specific absolute path is required.

Before running a notebook, the repository layout can be checked with:

```bash
python scripts/check_project_layout.py
```

Data and generated outputs must exist only at the repository level:

```text
data/
results/
```

The following locations are invalid and should not exist:

```text
notebooks/data/
notebooks/results/
```

## Determinism

The uploaded notebooks do not consistently specify random seeds for every
stochastic operation. For stricter reproducibility, define a project seed and
pass it to UMAP, Leiden, and any other function exposing `random_state`.

Recommended convention:

```python
SEED = 0
np.random.seed(SEED)
```

Record parameter values and package versions for every release.

## Notebook output policy

The GitHub-ready notebook copies have cell outputs and execution counters
removed. This:

- avoids committing local Windows paths;
- reduces repository size;
- prevents stale output from appearing current;
- makes code review easier.

For releases, executed HTML reports can be generated separately and attached as
release assets.

## File naming

Generated figures use stable section-based names rather than notebook cell
numbers. Cell numbers change whenever a cell is inserted, while names such as
`S05_03_neural_plate_1_vs_2_heatmap.png` remain interpretable.

## Data integrity

For archival releases, record:

- input file names;
- file sizes;
- SHA-256 checksums;
- genome and annotation versions;
- alignment/quantification software and parameters;
- sample metadata and biological replicate identifiers.

## Version-control recommendations

- keep raw data and `.h5ad` outputs outside normal Git tracking;
- use Git LFS only if large versioned files are essential;
- tag analysis releases;
- archive a release with Zenodo or a comparable service;
- add a `CITATION.cff` after authorship and publication metadata are finalized.

## Reproducibility limitations found during audit

- stored notebook outputs were generated with older path configurations;
- several source notebooks were edited after execution;
- some final counts differ from the immediately retained singlet counts because
  of later filtering;
- regression and batch correction are commented, not executed;
- the original Section 5 reused incompatible DEG keys for some plots.

The cleaned notebook copies correct the path templates and Section 5 output
handling, but the analyses must be rerun to produce authoritative final values.

## QC-directory revision

Section 1 preprocessing figures are routed directly to
`results/<stage>[/<condition>]/qc/`. Post-doublet read-count distribution and
mean/median plots are generated after singlet filtering. The general `figures/`
directory now contains Sections 2-6 only.


## Static notebook validation

Run the repository-level notebook validator after editing notebook files:

```bash
python scripts/validate_notebooks.py
```

Review the generated reports in `validation/`. Static validation should be
followed by `Restart Kernel and Run All Cells` for each notebook when the
required input data are available.
