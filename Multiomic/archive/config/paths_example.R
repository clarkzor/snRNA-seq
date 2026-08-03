# Portable path template. Copy to config/paths.local.R and edit for your machine.
# config/paths.local.R is ignored by Git.

MULTIOME_DATA_DIR <- Sys.getenv(
  "MULTIOME_DATA_DIR",
  unset = file.path(PROJECT_ROOT, "data", "raw")
)

FEATURE_MATRIX_H5 <- file.path(MULTIOME_DATA_DIR, "filtered_feature_bc_matrix.h5")
FRAGMENTS_FILE <- file.path(MULTIOME_DATA_DIR, "atac_fragments.tsv.gz")
FRAGMENTS_INDEX <- paste0(FRAGMENTS_FILE, ".tbi")

# Custom Xenopus annotation packages installed in the active R library.
ENSDB_PACKAGE <- Sys.getenv("XTROP_ENSDB_PACKAGE", unset = "EnsDb.Xtropicalis.v111")
BSGENOME_PACKAGE <- Sys.getenv(
  "XTROP_BSGENOME_PACKAGE",
  unset = "BSgenome.Xtropicalis.NCBI.UCBXtro10.0"
)
