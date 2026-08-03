# Run manually in R. Package installation is deliberately separated from analysis.
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
install.packages(c("Seurat", "Signac", "ggplot2", "patchwork", "remotes", "renv"))
BiocManager::install(c("GenomeInfoDb", "GenomicRanges", "ensembldb", "BSgenome", "biovizBase"))

message("Custom Xenopus packages are not installed automatically here.")
message("Install EnsDb.Xtropicalis.v111 before running the core pipeline.")
message("Install BSgenome.Xtropicalis.NCBI.UCBXtro10.0 before enabling LinkPeaks.")
