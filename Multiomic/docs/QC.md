# Quality control

Configured thresholds are stored in `config/analysis_config.R`.

Current defaults:

- RNA counts: >1,000 and <25,000
- ATAC counts: >1,000 and <100,000
- nucleosome signal: <2 when available
- TSS enrichment: >1 when available

`NucleosomeSignal()` and `TSSEnrichment()` are retained as a Windows-compatible fallback
while the newer `ATACqc()` route requires an accessible `fragtk` executable.

QC decisions and pre/post metadata are written under `results/qc/tables/`.
