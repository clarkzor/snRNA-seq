setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("BSgenomeForge")
BiocManager::install("Seurat")

install.packages("BSgenomeForge")
library(BSgenomeForge)
install.packages("Signac")
library("Signac")
library("ggplot2")
library("Seurat")
install.packages("biovizBase")
install.packages("BSgenome")
library(BSgenome)
library(EnsDb.Xtropicalis.v111)
library(ensembldb)
install.packages("biovizBase")
library(biovizBase)
library(devtools)
install("C:/Users/coron/OneDrive/Desktop/10xGenomics/EnsDb.Xtropicalis.v111")

BiocManager::install("devtools")
library(MACSr)

install.packages("installr")
library(installr)
updateR()


install.packages("hdf5r")
library(hdf5r)


# Install and load BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
library(BiocManager)

# Install EnsDb package (replace 'Xtropicalis' with the appropriate species name)
BiocManager::install("ensembldb")


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biovizBase")

# Load the installed EnsDb package
library(ensembldb)
library(biovizBase)

gtf <- "/Users/coron/OneDrive/Desktop/10xGenomics/Xenopus_tropicalis.UCB_Xtro_10.0.111.gtf"


### BUILDING THE FIRST REFERENCE FROM GTF ANNOTATION ###

#gtf
#ensDbFromGtf(gtf, , , )

#ensdb <- "/Users/coron/OneDrive/Desktop/10xGenomics/Xenopus_tropicalis.UCB_Xtro_10.0.111.sqlite"

#makeEnsembldbPackage(ensdb, version="0.0.1", maintainer='zor <zor@zor.com>', author="Clickityclack", destDir="/Users/coron/OneDrive/Desktop/10xGenomics/")

install.packages("/Users/coron/OneDrive/Desktop/10xGenomics/EnsDb.Xtropicalis.v111", repos = NULL, type = "source")


.libPaths()



######### UPLOADING DATA AND FORMING A SERUAT OBJECT###########
library(Signac)
library(Seurat)
library(EnsDb.Xtropicalis.v111)
# load the RNA and ATAC data
counts <- Read10X_h5("/Users/coron/OneDrive/Desktop/10xGenomics/filtered_feature_bc_matrix.h5")

fragpath <- "/Users/coron/OneDrive/Desktop/10xGenomics/atac_fragments.tsv.gz"

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb =EnsDb.Xtropicalis.v111)
#seqlevels(annotation) <- paste0('chr', seqlevels(annotation))




# create a Seurat object containing the RNA adata
pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)


##ANNOTATING THE CELL TYPES VIA RNAseq DATA#####


DefaultAssay(pbmc) <- "RNA"
DefaultAssay(pbmc)

# perform visualization and clustering steps
pbmc <- NormalizeData(pbmc, assay="RNA")
pbmc <- FindVariableFeatures(pbmc, assay="RNA")
pbmc <- ScaleData(pbmc, assay = "RNA")
pbmc <- RunPCA(pbmc, verbose = FALSE, assay ="RNA")
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 13, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:30)
DimPlot(pbmc, label = TRUE)




# Visualize gene expression on UMAP for the chosen assay
FeaturePlot(pbmc2, features = ("incenp"), pt.size = 1.25, reduction = "umap") + NoLegend()

# Remove the cells belonging to the cluster "Unknown"
pbmc <- subset(pbmc, idents = "Unknown")

genes_of_interest <- c("sox2", "foxi1", "sox17a",'grhl1', 'grhl3', 'tbxt',
                       'tp63', 'tp73', 'darmin', 'epas1', 'sox3', 'pax3', 'krt7',
                       'krt70', 'ventx2.1', 'chrd.1', 'wnt8a', 'mcidas', 'dnah5',
                       'foxj1', 'lrp2', 'dscaml1', 'slc5a8.1'
                       )
genes_of_interest <- c("COX3", 'ND5', 'COX1', 'ND3')



genes_of_interest <- c("agr2", 'itln1')

FeaturePlot(pbmc, features = genes_of_interest, pt.size = 1, reduction = "umap") + NoLegend()


DimPlot(new_obj, label = TRUE, repel = TRUE, group.by = "seurat_clusters", reduction = "umap", pt.size=1)


### ANNOTATING CLUSTER IDs ###


new.cluster.ids <- c("Novel Cluster (tp73/dscaml1)", "Non-neural Ectoderm (tp63+)", "Non-neural Ectoderm (tp63+)", "Neural Plate", "Non-neural Ectoderm (outer)", "Non-neural Ectoderm (tp63+)", "Non-neural Ectoderm (outer)", "Cement Gland", "Non-neural Ectoderm (outer)", "Unknown", "Non-neural Ectoderm (tp63+)",
                        "Neural Plate", "Dorsal Mesoderm", "Ventral Mesoderm", "Ventral Endoderm", "Dorsal Mesoderm", "Neural Plate (outer)", "Dorsal Mesoderm", "Neural Plate", "Neural Plate Border", "Unknown",
                     "Non-neural Ectoderm (outer)", "Ventral Mesoderm", "Unknown", "Dorsal Endoderm", "Neural Plate", "Involuting Endoderm", "Neural Plate Border", "Ventral Mesoderm", "Non-neural Ectoderm (tp63+)", "Neural Plate",
                        "Unknown", "Cement Gland", "Neural Plate Border", "Neural Plate", "Neural Plate (outer)", "Neural Plate", "Neural Plate", "Non-neural Ectoderm (outer)", "Dorsal Mesoderm", "Dorsal Endoderm",
                     "Non-neural Ectoderm (outer)", "Neural Plate Border", "Unknown", "Neural Plate", "Ciliated Epidermal Prog. (mcidas)", "Non-neural Ectoderm (tp63+)", "Involuting Ventral Mesoderm", "Dorsal Endoderm", "Non-neural Ectoderm (outer)", "Non-neural Ectoderm (tp63+)", "Involuting Dorsal Mesoderm",
                     "Dorsal Mesoderm", "Dorsal Mesoderm", "Non-neural Ectoderm (tp63+)", "Neural Plate", "Non-neural Ectoderm (outer)", "Non-neural Ectoderm (outer)")


names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
current_names <- levels(pbmc2)

# Order the cluster names alphabetically
desired_order <- sort(current_names)

# Create a named vector with the desired order
new_names <- setNames(desired_order, desired_order)

# Rename the clusters with the sorted order
pbmc2 <- RenameIdents(pbmc2, new_names)

plot <- DimPlot(pbmc2, reduction = "umap", label = FALSE, pt.size = 2.5)
plot + theme(text = element_text(size = 19))


# Assuming your Seurat object is called pbmc

# Rename the cluster
pbmc2 <- RenameIdents(pbmc2, "Non-neural Ectoderm (tp63+)" = "Non-neural Ectoderm (inner)")
pbmc2 <- RenameIdents(pbmc2, "Novel Cluster (tp73/dscaml1)" = "Novel Ecto. (tp73/dscaml1)")
pbmc2 <- RenameIdents(pbmc2, "Ciliated Epidermal Prog. (mcidas)" = "Ciliated Epidmeral Prog.")
pbmc2 <- RenameIdents(pbmc2, "Ciliated Epidmeral Prog." = "Ciliated Epidermal Prog.")


pbmc <- subset(pbmc, idents = c("Novel Ecto. (tp73/dscaml1)", "Non-neural Ectoderm (tp63+)", "Neural Plate", "Non-neural Ectoderm (outer)", "Cement Gland", "Dorsal Mesoderm", "Ventral Mesoderm", "Ventral Endoderm", "Neural Plate (outer)", "Neural Plate Border", "Dorsal Endoderm", "Involuting Endoderm", "Ciliated Epidermal Prog. (mcidas)", "Involuting Ventral Mesoderm", "Involuting Dorsal Mesoderm"
))




common_barcodes <- intersect(colnames(pbmc), colnames(pbmc2))
new_obj <- AddMetaData(pbmc, metadata=pbmc2@meta.data[common_barcodes, c('seurat_clusters'), drop=FALSE])




pbmc3 <- subset(pbmc2, idents = c("Novel Ecto. (tp73/dscaml1)", "Non-neural Ectoderm (inner)", "Neural Plate", "Non-neural Ectoderm (outer)", "Cement Gland", "Neural Plate (outer)", "Neural Plate Border", "Ciliated Epidermal Prog."
))


# find all markers of cluster 2
Unknown.markers <- FindMarkers(pbmc2, ident.1 = "Unknown")
head(Unknown.markers, n = 5)



###PRELIM ATAC CALCULATIONS###

DefaultAssay(pbmc) <- "ATAC"

pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)


VlnPlot(
  object = pbmc,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)


# filter out low quality cells
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)


###PEAK CALLING#####

# call peaks using MACS2
peaks <- CallPeaks(pbmc)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
#peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)



# In lieu of MACS2, we're going to use the already called peaks from the cellranger arc output
# To do this we'll need to upload the .bed file and then convert it into a granges object using rtracklayer
if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  install.packages("rtracklayer")
}
library(rtracklayer)

peaks <- import.bed("/Users/coron/OneDrive/Desktop/10xGenomics/atac_peaks.narrowPeak.bed")



# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(pbmc),
  features = peaks,
  cells = colnames(pbmc)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
pbmc[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)


##### GENE EXPRESSION DATA PROCESSING #####
BiocManager::install('glmGamPoi')

DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc)
pbmc <- RunPCA(pbmc)



##### DNA ACCESSIBILITY DATA PROCESSING #####

DefaultAssay(pbmc) <- "peaks"
pbmc <- FindTopFeatures(pbmc, min.cutoff = 5)
pbmc <- RunTFIDF(pbmc)

#Installing matrix and irlba from the source fixed the error
#install.packages("Matrix", type = "source")
#install.packages("irlba", type = "source")

pbmc <- RunSVD(pbmc)



##### JOINT UMAP VISUALIZATION ######

# build a joint neighbor graph using both assays
pbmc <- FindMultiModalNeighbors(
  object = pbmc,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)


# build a joint UMAP visualization
pbmc <- RunUMAP(
  object = pbmc,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

DimPlot(pbmc2, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend()


# Assuming your Seurat object is called pbmc

# Save the Seurat object to a .RDS file
saveRDS(pbmc, file = "St12_Multiomic_seurat_object.rds")

# Load the Seurat object back into R
pbmc2 <- readRDS("St12_Multiomic_seurat_object.rds")


# Run Leiden clustering on UMAP coordinates
#First check which graphs you have calculated
names(pbmc@graphs)
#Can use wknn or wsnn based on above result
pbmc <- FindClusters(pbmc, graph.name="wsnn", resolution = 13, algorithm=1)

DimPlot(pbmc, label = TRUE, repel = TRUE, reduction = "umap", pt.size = 2.5) + NoLegend()

DimPlot(pbmc, label = FALSE, repel = TRUE, reduction = "umap", pt.size = 2)


#### VISUALIZING PARTICULAR GENES #####


CoveragePlot(pbmc, region = 'sox2', features = 'sox2', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE)
CoveragePlot(pbmc3, region = 5:104010000, features = 'tp63', assay = 'ATAC',
             expression.assay = 'SCT', peaks = TRUE,
             extend.upstream = 10000,
             extend.downstream = 10000)

zeb2
tp63
krt7



# Assuming your Seurat object is called pbmc

# Define the desired order of cluster names
desired_order <- c("Neural Plate", "Neural Plate Border", "Neural Plate (outer)", 
                   "Novel Ecto. (tp73/dscaml1)", "Ciliated Epidermal Prog.", 
                   "Non-neural Ectoderm (inner)", "Non-neural Ectoderm (outer)", 
                   "Cement Gland")

# Get the current cluster names
current_names <- levels(pbmc3)

# Create a named vector with the desired order
new_names <- setNames(desired_order, desired_order)

# Rename the clusters with the specified order
pbmc3 <- RenameIdents(pbmc3, new_names)



##VISUALIZING READ DENSITY AT A GENOMIC LOCI#####

tile_plot <- TilePlot(
  object = pbmc,
  region = "chr2-1000-10000",
  idents = c("Ventral Mesoderm", "Dorsal Mesoderm")
)
tile_plot




##### LINKING GENES TO PEAKS #####

BiocManager::install("BSgenomeForge")
library(BSgenomeForge)
forgeBSgenomeDataPkgFromNCBI("GCF_000004195.4", pkg_maintainer="zor <clarklh@uci.edu>", destdir="/Users/coron/OneDrive/Desktop/10xGenomics/")
install.packages("/Users/coron/OneDrive/Desktop/10xGenomics/BSgenome.Xtropicalis.NCBI.UCBXtro10.0", repos = NULL, type = "source")
library(BSgenome.Xtropicalis.NCBI.UCBXtro10.0)

DefaultAssay(pbmc) <- "peaks"

# first compute the GC content for each peak
pbmc <- RegionStats(pbmc, genome = BSgenome.Xtropicalis.NCBI.UCBXtro10.0)

install.packages("devtools")
library(devtools)
devtools::install_github("cysouw/qlcMatrix")



# link peaks to genes
pbmc <- LinkPeaks(
  object = pbmc,
  peak.assay = "peaks",
  expression.assay = "SCT",
  genes.use = c("sox17a", "tbxt",
                gene.id=TRUE)
)


pbmc <- LinkPeaks(
  object = pbmc,
  peak.assay = "peaks",
  expression.assay = "SCT",
)



CoveragePlot(
  object = pbmc,
  region = "sox2",
  features = "sox2",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 500,
  extend.downstream = 10000
)
