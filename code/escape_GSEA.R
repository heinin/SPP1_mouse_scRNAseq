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

seurat_data <- readRDS("/tgen_labs/banovich/BCTCSF/SPP1_mouse_scRNAseq/scRNAseq_Seurat_dim8_nCount1k_nFeature500_dblrate15.rds")
DefaultAssay(seurat_data) <- "RNA"

head(seurat_data@meta.data)

DimPlot(seurat_data,
        group.by = "integratedSCTsnn_res.0.8")

#==============================================================================#
# Run GSEA
#==============================================================================#


