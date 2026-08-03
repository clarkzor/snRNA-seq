# Results

This directory is populated by the notebooks and is ignored by Git by default.

Each dataset generates:

```text
results/<stage>[/<condition>]/
├── qc/                      # Section 1 and post-doublet QC plots/tables
├── figures/                 # Sections 2-6
├── processed_data/
├── differential_expression/
├── validation/
└── tables/
```

The `qc/` directory contains the preprocessing violin plots, top-gene plots,
mitochondrial/ribosomal distributions, Scrublet score histogram, post-doublet
reads-per-cell histogram, mean/median reads-per-cell bar graph, and associated
CSV summaries.

Keep large `.h5ad` files and complete figure collections outside normal Git
history. Copy selected final figures to a separate tracked manuscript directory
when needed.
