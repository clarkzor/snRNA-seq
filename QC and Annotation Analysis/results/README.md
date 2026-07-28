# Results

This directory is populated by the notebooks and is ignored by Git by default.

Each dataset generates:

```text
results/<stage>[/<condition>]/
├── qc/
├── figures/
├── processed_data/
├── differential_expression/
├── validation/
└── tables/
```

Keep large `.h5ad` files and complete figure collections outside normal Git
history. Copy selected final figures to a separate tracked manuscript directory
when needed.
