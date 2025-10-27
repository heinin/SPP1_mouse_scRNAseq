#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 04/10/2025
# Description: Gene Set Enrichment Analysis for CAR+anti-SPP1 mouse tumors
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(patchwork)
library(tidyr)
library(tidyverse)
library(googlesheets4)
library(scGSVA)

#==============================================================================#
# Helper functions and eEnvironment variables
#==============================================================================#

setwd("/home/hnatri/SPP1_mouse_scRNAseq/")
set.seed(1234)
options(future.globals.maxSize = 30000 * 1024^2)
reduction <- "umap"

source("/home/hnatri/SPP1_mouse_scRNAseq/code/CART_plot_functions.R")
source("/home/hnatri/SPP1_mouse_scRNAseq/code/colors_themes.R")

#==============================================================================#
# Import data
#==============================================================================#

seurat_data <- readRDS("/tgen_labs/banovich/BCTCSF/SPP1_mouse_scRNAseq/scRNAseq_Seurat_dim8_nCount1k_nFeature500_dblrate15_scGSEAmicroglia_reclustered_refAnnot.rds")
DefaultAssay(seurat_data) <- "RNA"

head(seurat_data@meta.data)

DimPlot(seurat_data,
        group.by = "integratedSCTsnn_res.0.8")

# Gene sets
gs4_deauth()
genesets  <- gs4_get("https://docs.google.com/spreadsheets/d/127J6C4KF7uBGKUnrPuC1mcsb_wNCN6k1zXKSCbJ6q0M/edit?usp=sharing")
genesets <- read_sheet(genesets, sheet = "Myeloid modules")

#==============================================================================#
# Run scGSVA
#==============================================================================#

unique(genesets$gene)

setdiff(genesets$gene, rownames(seurat_data))

grep("H2", rownames(seurat_data), value = T)

# scGSVA analysis
module_annot <- buildAnnot(species="mouse", keytype="SYMBOL", anntype="GO")
module_annot@anntype <- "custom"
module_annot@keytype

names(module_annot)
typeof(module_annot@annot)
head(module_annot@annot)

module_annot@annot <- genesets

module_annot@annot$Dummy <- module_annot@annot$module
module_annot@annot <- as.data.frame(module_annot@annot)
colnames(module_annot@annot) <- c("GeneID", "Dummy", "Annot")

DefaultAssay(seurat_data) <- "RNA"
res <- scgsva(seurat_data,
              annot = module_annot,
              verbose = T)

res_df <- as.data.frame(res)

rownames(res_df)

identical(rownames(res_df), colnames(seurat_data))

seurat_data$DAM <- res_df$DAM
seurat_data$Homeostatic <- res_df$Homeostatic
seurat_data$MHC_I_machinery <- res_df$MHC_I_machinery
seurat_data$MHC_II_machinery <- res_df$MHC_II_machinery
seurat_data$Costimulation_APCB <- res_df$Costimulation_ligands_APCB
seurat_data$Adhesion_immunesynapse <- res_df$Adhesion_immunesynapse

unique(genesets$module)

saveRDS(seurat_data, "/tgen_labs/banovich/BCTCSF/SPP1_mouse_scRNAseq/scRNAseq_Seurat_dim8_nCount1k_nFeature500_dblrate15_scGSEAmodules_reclustered_refAnnot.rds")

