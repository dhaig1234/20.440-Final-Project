# Pset 4 - Dana Haig
library(Seurat)

# Install and load packages quietly (not getting a bunch of messages spit out in console)
library(Seurat)
suppressMessages(require(Seurat))
suppressMessages(require(dplyr))
suppressMessages(require(ggplot2))
suppressMessages(require(patchwork))
suppressMessages(require(Matrix))
install.packages("Matrix", repos = "http://cran.r-project.org") # Matrix was causign errors, so download twice..?
install.packages("remotes")
library(remotes)
install.packages("fields")
library(fields)
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", upgrade = FALSE,
                        dependencies = FALSE)
suppressMessages(require(DoubletFinder))

# specify path containing the 10x data
data_path <- "/Users/danahaig/Desktop/440_Final_Project_GSE166326"
# rename downloaded data to remove GEO Index (GSE166326)
data <- Read10X("/Users/danahaig/Desktop/440_Final_Project_GSE166326", strip.suffix = TRUE)

# create seurat object, include gene expression 
c_seurat_object <- CreateSeuratObject(counts = data$`Gene Expression`)
Assays(c_seurat_object)
dim(c_seurat_object)
# include antibody capture
adt_assay <- CreateAssayObject(counts = data$`Antibody Capture`)
c_seurat_object[["ADT"]] <- adt_assay
dim(adt_assay)
# include custom assay 
custom_assay <- CreateAssayObject(counts = data$`Custom`)
c_seurat_object[["Custom"]] <- custom_assay
dim(custom_assay)

# name columns after percent mitochondria, ribosome, hb, and plat
c_seurat_object <- PercentageFeatureSet(c_seurat_object, "^MT-", col.name = "percent_mito")
c_seurat_object <- PercentageFeatureSet(c_seurat_object, "^RP[SL]", col.name = "percent_ribo")
c_seurat_object <- PercentageFeatureSet(c_seurat_object, "^HB[^(P)]", col.name = "percent_hb")
c_seurat_object <- PercentageFeatureSet(c_seurat_object, "PECAM1|PF4", col.name = "percent_plat")

# add features 
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
VlnPlot(c_seurat_object, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) +
  NoLegend()
Assays(c_seurat_object)
DefaultAssay(c_seurat_object) <- "RNA"
DefaultAssay(c_seurat_object)

# normalize, scale, create a PCE plot, cluster, and run UMAP
c_seurat_object <- NormalizeData(c_seurat_object)
c_seurat_object <- FindVariableFeatures(c_seurat_object)
c_seurat_object <- ScaleData(c_seurat_object)
c_seurat_object <- RunPCA(c_seurat_object, verbose = FALSE)
c_seurat_object <- FindNeighbors(c_seurat_object, dims = 1:30)
c_seurat_object <- FindClusters(c_seurat_object, resolution = 0.8, verbose = FALSE)
c_seurat_object <- RunUMAP(c_seurat_object, dims = 1:30)

# display plot! 
DimPlot(c_seurat_object, label = TRUE)