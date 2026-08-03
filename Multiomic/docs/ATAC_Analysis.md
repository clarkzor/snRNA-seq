# ATAC analysis

The Cell Ranger ARC peak matrix already present in `filtered_feature_bc_matrix.h5` is used
as the ATAC assay. The core workflow does not recreate a second peak assay from the same
Cell Ranger peak set.

ATAC dimensional reduction uses:

1. `FindTopFeatures()`
2. `RunTFIDF()`
3. `RunSVD()`

The resulting `atac.lsi` reduction is combined with `rna.pca` using WNN.

## Fragment-dependent steps

Coverage plots require `atac_fragments.tsv.gz` plus its `.tbi` index. They are disabled by
default in `config/analysis_config.R` so plotting-library compatibility issues cannot stop
the core pipeline.

Peak-to-gene linking additionally requires the matching Xenopus BSgenome package and
`RegionStats()` before `LinkPeaks()`.
