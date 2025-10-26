#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 09/30/2025
# Description: Reference-based annotation for mouse tumors
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

library(Seurat)
library(SeuratDisk)
library(googlesheets4)
library(tidyverse)
library(plyr)
library(ggrepel)
library(patchwork)
library(anndata)
library(biomaRt)

setwd("/home/hnatri/SPP1_mouse_scRNAseq/")
set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)
reduction <- "integratedSCTumap"

source("/home/hnatri/SPP1_mouse_scRNAseq/code/CART_plot_functions.R")
source("/home/hnatri/SPP1_mouse_scRNAseq/code/colors_themes.R")

#==============================================================================#
# Import data
#==============================================================================#

seurat_data <- readRDS("/tgen_labs/banovich/BCTCSF/SPP1_mouse_scRNAseq/scRNAseq_Seurat_dim8_nCount1k_nFeature500_dblrate15.rds_scGSEAmicroglia_reclustered.rds")

ref_data <- anndata::read_h5ad("/tgen_labs/banovich/BCTCSF/SPP1_mouse_scRNAseq/Reference_data/86afbf8c-e49c-4e71-9ba3-e07b571f1acf.h5ad")
ref_data_seurat <- CreateSeuratObject(counts = t(as.matrix(ref_data$raw)),
                                      meta.data = ref_data$obs)

#saveRDS(ref_data_seurat, "/tgen_labs/banovich/BCTCSF/SPP1_mouse_scRNAseq/Reference_data/86afbf8c-e49c-4e71-9ba3-e07b571f1acf.rds")
#ref_data_seurat <- readRDS("/tgen_labs/banovich/BCTCSF/SPP1_mouse_scRNAseq/Reference_data/86afbf8c-e49c-4e71-9ba3-e07b571f1acf.rds")

head(ref_data_seurat@meta.data)

#==============================================================================#
# Transfering labels
#==============================================================================#

# Gene names
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

gene_names <- getBM(attributes = c("mgi_symbol", "ensembl_gene_id"), 
                    filters = "ensembl_gene_id", 
                    values = rownames(ref_data_seurat),
                    mart = ensembl)

# Removing problematic gene names
gene_names$mgi_symbol[gene_names$mgi_symbol==""]

gene_names <- gene_names[-which(gene_names$mgi_symbol == ""),]
gene_names <- gene_names[-which(duplicated(gene_names$mgi_symbol)),]

counts <- LayerData(ref_data_seurat,
                    assay = "RNA",
                    layer = "counts")
counts <- as.matrix(counts)
counts <- counts[gene_names$ensembl_gene_id,]

identical(rownames(counts), gene_names$ensembl_gene_id)

rownames(counts) <- gene_names$mgi_symbol
ref_data_seurat[["new_RNA"]] <- CreateAssay5Object(counts = counts)

DefaultAssay(ref_data_seurat) <- "new_RNA"
ref_data_seurat <- SCTransform(ref_data_seurat,
                               assay = "new_RNA")
ref_data_seurat <- RunPCA(ref_data_seurat,
                          assay = "SCT")

anchors <- FindTransferAnchors(reference = ref_data_seurat,
                               query = seurat_data,
                               normalization.method = "SCT",
                               reference.assay = "SCT",
                               query.assay = "SCT",
                               dims = 1:30,
                               reference.reduction = "pca",
                               features = VariableFeatures(seurat_data))
predictions <- TransferData(anchorset = anchors,
                            refdata = ref_data_seurat$cell_type,
                            dims = 1:30)
seurat_data <- AddMetaData(seurat_data, metadata = predictions)

head(seurat_data@meta.data)

saveRDS(seurat_data, "/tgen_labs/banovich/BCTCSF/SPP1_mouse_scRNAseq/scRNAseq_Seurat_dim8_nCount1k_nFeature500_dblrate15.rds_scGSEAmicroglia_reclustered_refAnnot.rds")