# Shared setup. Run all analysis scripts from the repository root.
options(stringsAsFactors = FALSE)

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    markers <- c(
      file.exists(file.path(current, "README.md")),
      dir.exists(file.path(current, "scripts")),
      dir.exists(file.path(current, "config")),
      dir.exists(file.path(current, "results"))
    )
    if (all(markers)) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  stop("Could not locate repository root. Run from the cloned repository.")
}

PROJECT_ROOT <- find_project_root()
setwd(PROJECT_ROOT)

source(file.path(PROJECT_ROOT, "config", "analysis_config.R"))
paths_file <- file.path(PROJECT_ROOT, "config", "paths.local.R")
if (!file.exists(paths_file)) {
  paths_file <- file.path(PROJECT_ROOT, "config", "paths_example.R")
}
source(paths_file)
set.seed(SEED)

RESULTS_DIR <- file.path(PROJECT_ROOT, "results")
QC_FIG_DIR <- file.path(RESULTS_DIR, "qc", "figures")
QC_TABLE_DIR <- file.path(RESULTS_DIR, "qc", "tables")
RNA_FIG_DIR <- file.path(RESULTS_DIR, "rna", "figures")
RNA_TABLE_DIR <- file.path(RESULTS_DIR, "rna", "tables")
ATAC_FIG_DIR <- file.path(RESULTS_DIR, "atac", "figures")
ATAC_TABLE_DIR <- file.path(RESULTS_DIR, "atac", "tables")
INTEGRATION_FIG_DIR <- file.path(RESULTS_DIR, "integration", "figures")
INTEGRATION_TABLE_DIR <- file.path(RESULTS_DIR, "integration", "tables")
ANNOTATION_FIG_DIR <- file.path(RESULTS_DIR, "annotation", "figures")
ANNOTATION_TABLE_DIR <- file.path(RESULTS_DIR, "annotation", "tables")
COVERAGE_DIR <- file.path(RESULTS_DIR, "coverage")
PEAK_LINK_DIR <- file.path(RESULTS_DIR, "peak_gene_links")
OBJECT_DIR <- file.path(RESULTS_DIR, "objects")
VALIDATION_DIR <- file.path(PROJECT_ROOT, "validation")

for (d in c(QC_FIG_DIR, QC_TABLE_DIR, RNA_FIG_DIR, RNA_TABLE_DIR,
            ATAC_FIG_DIR, ATAC_TABLE_DIR, INTEGRATION_FIG_DIR,
            INTEGRATION_TABLE_DIR, ANNOTATION_FIG_DIR, ANNOTATION_TABLE_DIR,
            COVERAGE_DIR, PEAK_LINK_DIR, OBJECT_DIR, VALIDATION_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

RULES_FILE <- file.path(PROJECT_ROOT, "config", "stage11_expression_annotation_rules.csv")
MANUAL_MAP_FILE <- file.path(PROJECT_ROOT, "config", "original_manual_cluster_map.csv")
MARKER_PANEL_FILE <- file.path(PROJECT_ROOT, "config", "marker_panels.csv")

required_packages <- c(
  "Seurat", "Signac", "SeuratObject", "ggplot2", "patchwork", "Matrix",
  "GenomeInfoDb", "GenomicRanges", "ensembldb", ENSDB_PACKAGE
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ", paste(missing_packages, collapse = ", "),
    "\nSee docs/Reproducibility.md and scripts/install_dependencies.R.", call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(SeuratObject)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(GenomeInfoDb)
  library(GenomicRanges)
  library(ensembldb)
})
suppressPackageStartupMessages(library(ENSDB_PACKAGE, character.only = TRUE))

source(file.path(PROJECT_ROOT, "scripts", "utils", "io_utils.R"))
source(file.path(PROJECT_ROOT, "scripts", "utils", "plotting_utils.R"))
source(file.path(PROJECT_ROOT, "scripts", "utils", "genomics_utils.R"))
source(file.path(PROJECT_ROOT, "scripts", "utils", "annotation_utils.R"))

message("Project root: ", PROJECT_ROOT)
message("Data directory: ", MULTIOME_DATA_DIR)
