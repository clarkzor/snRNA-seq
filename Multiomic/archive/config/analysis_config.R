# Analysis parameters for Stage 11 Xenopus tropicalis RNA+ATAC multiome.
# Edit deliberately and commit changes so analysis decisions remain auditable.

SEED <- 1234

# QC thresholds inherited from the original Signac workflow.
MIN_RNA_COUNTS <- 1000
MAX_RNA_COUNTS <- 25000
MIN_ATAC_COUNTS <- 1000
MAX_ATAC_COUNTS <- 100000
MAX_NUCLEOSOME_SIGNAL <- 2
MIN_TSS_ENRICHMENT <- 1

# Processing choices
USE_SCTRANSFORM <- TRUE
RUN_FRAGMENT_QC <- TRUE
RUN_RNA_ONLY_CLUSTERING <- TRUE

RNA_VARIABLE_FEATURES <- 3000
RNA_PCS <- 50
ATAC_LSI_COMPONENTS <- 50
RNA_DIMS_WNN <- 1:30
ATAC_DIMS_WNN <- 2:40

# Clustering
RNA_LOUVAIN_RESOLUTION <- 13
WNN_LEIDEN_RESOLUTIONS <- c(0.5, 1.0, 2.0, 5.0, 13.0)
FINAL_WNN_RESOLUTION <- 1.0

# Expression annotation
UNASSIGNED_LABEL <- "Other"
USE_CLUSTER_CONSENSUS_FILL <- FALSE
CLUSTER_FILL_MIN_LABELED_CELLS <- 10
CLUSTER_FILL_MIN_FRACTION <- 0.70

# Optional fragment/genome-dependent analyses.
# Disabled by default so plotting or BSgenome issues do not stop the core pipeline.
RUN_COVERAGE_PLOTS <- FALSE
RUN_PEAK_GENE_LINKS <- FALSE
COVERAGE_GENES <- c("sox2", "pax3", "tp63", "tbxt")
LINKPEAKS_GENES <- c("sox17a", "tbxt")
