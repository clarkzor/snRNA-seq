if (!exists("PROJECT_ROOT")) source("scripts/00_setup.R")
if (!exists("pbmc")) pbmc <- load_checkpoint("03_rna_processed.rds")

DefaultAssay(pbmc) <- "ATAC"
pbmc <- Signac::FindTopFeatures(pbmc, min.cutoff = 5)
pbmc <- Signac::RunTFIDF(pbmc)
pbmc <- Signac::RunSVD(
  pbmc, n = ATAC_LSI_COMPONENTS,
  reduction.name = "atac.lsi", reduction.key = "atacLSI_"
)
save_checkpoint(pbmc, "04_atac_processed.rds")
