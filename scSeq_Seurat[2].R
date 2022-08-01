library(dplyr)
library(patchwork)
library(NMF)
library(RColorBrewer)
library(statmod)
library(ggplot2)
library(cowplot)
library(Seurat)
library(Matrix)
library(stringr)
library(plyr)
library(ggsignif)
library(reshape2)
library(quantreg)
library(R.utils)

### Load all the required functions for this analysis
source("C:/Users/LehtinenLab/Dropbox/code/single_cell_airway/Fxns.R")

setwd("E:/sc data/")

#load and format metadata
header <- scan("scEmbryo_cell_metadata.txt", nlines = 1, what = character())
metadata <- read.delim("scEmbryo_cell_metadata.txt", skip = 1, row.names = 1, sep = "\t")
names(metadata) <- header[(-1)]

#load preprocessed, batch corrected data and make seurat object (with metadata)
expression_matrix <- ReadMtx(mtx = "sc_rna_data.mtx", cells = "sc_rna_barcodes.tsv", features = "sc_rna_genes.tsv", feature.column = 1)
seurat_object <- CreateSeuratObject(counts=expression_matrix, meta.data = metadata)

#use Montero et al. script to get variable genes by batch
thresh = -0.15
v <- get.variable.genes.umis(expression_matrix, residual.threshold = thresh, batch = seurat_object$donor_id, ret.plot = T)
outliers = subset(v$fit.data, residual < thresh)

#use only genes that are in all three batches
g = outliers$gene
g_dup = g[duplicated(g)]
g_trip = g_dup[duplicated(g_dup)]
VariableFeatures(seurat_object) <- g_trip

# A<-as.numeric(v$var.genes)
# B<-v$fit.data$gene[A]?
# VariableFeatures(seurat_object) <- B

#scale the data
seurat_object <- ScaleData(seurat_object)

VlnPlot(seurat_object, features = "Sod3", group.by = "cell_type__ontology_label", log = TRUE)

#Run PCA
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(seurat_object))

#Run tSNE embedding
seurat_object <- RunTSNE(seurat_object, ndims = 1:5)

#plot tSNE
DimPlot(seurat_object, reduction = "tsne", group.by = "cell_type__ontology_label")
DimPlot(seurat_object, reduction = "tsne", group.by = "donor_id")
DimPlot(seurat_object, reduction = "tsne", group.by = "ventricle")

# clustering (skip this step) ----
# sig_PCs <- Reductions(seurat_object, slot = "pca")
# sig_PCs <- sig_PCs[,1:10]
seurat_object <- FindNeighbors(seurat_object, reduction = "pca", dims = 1:10, return.neighbor = FALSE)
seurat_object <- FindClusters(seurat_object, method = "igraph", resolution = 0.1)

#Run tSNE embedding
seurat_object <- RunTSNE(seurat_object, ndims = 1:5)
DimPlot(seurat_object, reduction = 'tsne')

cluster4_markers <- FindMarkers(seurat_object, ident.1 = 4, ident.2 = 0, min.pct = 0.25)
head(cluster4_markers, n=10)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

seurat_object.markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster_top <- seurat_object.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)












