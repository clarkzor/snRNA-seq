# Data

Large sequencing files are intentionally excluded from Git.

Place local Cell Ranger ARC outputs in `data/raw/`, or edit the ignored
`config/paths.local.R` file to point to an existing data directory.

Required for the core RNA+ATAC matrix workflow:

- `filtered_feature_bc_matrix.h5`

Required for fragment-level QC and chromatin visualization:

- `atac_fragments.tsv.gz`
- `atac_fragments.tsv.gz.tbi`

The current local paths template points to the existing Windows data directory.
