install.packages("seuratdisk")

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")


library(Seurat)
library(Matrix)
library(SeuratDisk)
library(reticulate)
mat = t(readMM("count_matrix.mtx"))
genes = read.csv("all_genes.csv")
meta = read.csv("cell_metadata.csv")


# Set gene names as rownames — ensuring uniqueness
gene_names <- make.unique(as.character(genes$gene_name))
rownames(mat) <- gene_names

# Assign cell barcodes
rownames(meta) <- meta$bc_wells
colnames(mat) <- meta$bc_wells

# Create Seurat object
obj <- CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0, meta.data = meta)


library(Seurat)

# Extract data from the v5 Assay5 object
counts <- GetAssayData(obj, slot = "counts")
data <- GetAssayData(obj, slot = "data")
scale.data <- GetAssayData(obj, slot = "scale.data")

# Create a classic Assay object
old_assay <- CreateAssayObject(counts = counts)
old_assay@data <- data
#old_assay@scale.data <- scale.data


# Replace the RNA assay in your Seurat object with the classic one
obj@assays$RNA <- old_assay
DefaultAssay(obj) <- "RNA"

DefaultAssay(obj) <- "RNA"
obj[["RNA"]] <- CreateAssayObject(counts = obj@assays$RNA@counts)


# Normalize data
#obj <- NormalizeData(obj)

# Save again
SaveH5Seurat(obj, filename = "PARSE-Unfiltered_St10-5_Compatible.h5Seurat", overwrite = TRUE)

# Convert to h5ad
Convert("PARSE-Unfiltered_St1-5_Compatible.h5Seurat", dest = "h5ad")

str(obj@meta.data)








######## NEW MO BETTA ##############
######## NEW MO BETTA ##############
######## NEW MO BETTA ##############
######## NEW MO BETTA ##############
######## NEW MO BETTA ##############
######## NEW MO BETTA ##############


meta <- read.csv("cell_metadata.csv", stringsAsFactors = FALSE)

# Make sure all columns are atomic (not lists or weird types)
meta[] <- lapply(meta, function(x) {
  if (is.logical(x) || is.numeric(x) || is.character(x) || is.factor(x)) return(x)
  return(as.character(x))  # Fallback for unknown types like Date
})

rownames(meta) <- meta$bc_wells


library(Matrix)
library(Seurat)
library(SeuratDisk)

# Read matrix and metadata
mat <- t(readMM("count_matrix.mtx"))
genes <- read.csv("all_genes.csv", stringsAsFactors = FALSE)
meta <- read.csv("cell_metadata.csv", stringsAsFactors = FALSE)

# Sanitize gene names
gene_names <- make.unique(as.character(genes$gene_name))
rownames(mat) <- gene_names

# Assign barcodes to matrix columns (cells)
colnames(mat) <- meta$bc_wells
rownames(meta) <- meta$bc_wells


obj <- CreateSeuratObject(counts = mat, meta.data = meta, project = "parse_data")

SaveH5Seurat(obj, filename = "PARSE-Unfiltered_St10-5-WT_Compatible.h5Seurat", overwrite = TRUE)

Convert("PARSE-Unfiltered_St10-5-WT_Compatible.h5Seurat", dest = "h5ad", overwrite = TRUE)
