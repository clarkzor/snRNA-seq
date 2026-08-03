# Run from the repository root:
#   Rscript scripts/run_all.R

source("scripts/00_setup.R")
pipeline_scripts <- c(
  "01_create_multiome_object.R",
  "02_qc_filtering.R",
  "03_rna_processing.R",
  "04_atac_processing.R",
  "05_wnn_integration.R",
  "06_clustering.R",
  "07_expression_annotation.R",
  "08_annotation_validation.R",
  "09_coverage_plots.R",
  "10_peak_gene_links.R",
  "11_finalize.R"
)
for (script_name in pipeline_scripts) {
  message("\n===== Running ", script_name, " =====")
  source(file.path("scripts", script_name), local = .GlobalEnv)
}
message("\nAll enabled pipeline stages completed.")
