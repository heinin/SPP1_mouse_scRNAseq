#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 04/10/2025
# Description: CellChat analysis on CAR+anti-SPP1 treated Kluc tumors
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(CellChat)

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
# Importing data and creating inputs
#==============================================================================#

seurat_object <- readRDS("/tgen_labs/banovich/BCTCSF/SPP1_mouse_scRNAseq/scRNAseq_Seurat_dim8_granular_annot.rds")

# annot_granular
seurat_object$celltype <- seurat_object$annot_granular
seurat_object$celltype <- factor(seurat_object$celltype,
                                 levels = sort(unique(seurat_object$celltype)))

# Cluster annotations
#gs4_deauth()
#cluster_annot  <- gs4_get("https://docs.google.com/spreadsheets/d/127J6C4KF7uBGKUnrPuC1mcsb_wNCN6k1zXKSCbJ6q0M/edit?usp=sharing")
#cluster_annot <- read_sheet(cluster_annot, sheet = "Cluster annotation")
#
## Updating annotations
#seurat_object$annot <- mapvalues(x = seurat_object$sub.cluster,
#                                 from = cluster_annot$sub.cluster,
#                                 to = cluster_annot$annot)
#
#seurat_object$celltype <- factor(seurat_object$annot,
#                                 levels = sort(unique(seurat_object$annot)))
#
Idents(seurat_object) <- seurat_object$celltype
seurat_object$Sample <- seurat_object$orig.ident
seurat_object$Sample <- gsub("SPP1\\+CAR", "CAR_SPP1", seurat_object$Sample)
seurat_object$Sample <- gsub("TUMOR", "CTRL", seurat_object$Sample)

table(seurat_object$celltype, seurat_object$Sample)

#DefaultAssay(seurat_object) <- "RNA_human"
#seurat_object <- NormalizeData(seurat_object)
#seurat_object <- ScaleData(seurat_object)

DefaultAssay(seurat_object) <- "RNA"
seurat_object <- NormalizeData(seurat_object)
seurat_object <- ScaleData(seurat_object)

# Inputs for CellChat:
# Log-normalized counts and cell labels (cluster/celltype)
# Creating CC objects for each group, and a merged object for comparative analysis
all_counts <- GetAssayData(seurat_object, assay = "RNA", slot = "data")
all_CC <- createCellChat(object = all_counts,
                         meta = seurat_object@meta.data,
                         group.by = "celltype")

## Control
CTRL_seurat <- subset(seurat_object, subset = Sample == "CTRL")
CTRL_counts <- GetAssayData(CTRL_seurat, assay = "RNA",
                            slot = "data")
CTRL_CC <- createCellChat(object = CTRL_counts,
                          meta = CTRL_seurat@meta.data,
                          group.by = "celltype")

## CAR only
CART_seurat <- subset(seurat_object, subset = Sample == "CAR")
CART_counts <- GetAssayData(CART_seurat, assay = "RNA",
                              slot = "data")
CART_CC <- createCellChat(object = CART_counts,
                          meta = CART_seurat@meta.data,
                          group.by = "celltype")

## anti-SPP1 only
SPP1_seurat <- subset(seurat_object, subset = Sample == "SPP1")
SPP1_counts <- GetAssayData(SPP1_seurat, assay = "RNA",
                             slot = "data")
SPP1_CC <- createCellChat(object = SPP1_counts,
                          meta = SPP1_seurat@meta.data,
                          group.by = "celltype")

## CAR + anti-SPP1
CART_SPP1_seurat <- subset(seurat_object, subset = Sample == "CAR_SPP1")
CAR_SPP1_counts <- GetAssayData(CART_SPP1_seurat, assay = "RNA",
                          slot = "data")
CAR_SPP1_CC <- createCellChat(object = CAR_SPP1_counts,
                              meta = CART_SPP1_seurat@meta.data,
                              group.by = "celltype")


# Adding metadata
CTRL_CC <- addMeta(CTRL_CC, meta = CTRL_seurat@meta.data)
CTRL_CC <- setIdent(CTRL_CC, ident.use = "celltype")
# Dropping an unused ident
CTRL_CC@idents <- factor(CTRL_CC@idents, levels = unique(CTRL_CC@idents))

CART_CC <- addMeta(CART_CC, meta = CART_seurat@meta.data)
CART_CC <- setIdent(CART_CC, ident.use = "celltype")
CART_CC@idents <- factor(CART_CC@idents, levels = unique(CART_CC@idents))

SPP1_CC <- addMeta(SPP1_CC, meta = SPP1_seurat@meta.data)
SPP1_CC <- setIdent(SPP1_CC, ident.use = "celltype")
SPP1_CC@idents <- factor(SPP1_CC@idents, levels = unique(SPP1_CC@idents))

CART_SPP1_CC <- addMeta(CAR_SPP1_CC, meta = CART_SPP1_seurat@meta.data)
CART_SPP1_CC <- setIdent(CART_SPP1_CC, ident.use = "celltype")
CART_SPP1_CC@idents <- factor(CART_SPP1_CC@idents, levels = unique(CART_SPP1_CC@idents))

all_CC <- addMeta(all_CC, meta = seurat_object@meta.data)
all_CC <- setIdent(all_CC, ident.use = "celltype")

setdiff(CART_SPP1_CC@idents, CART_CC@idents)
setdiff(CART_CC@idents, CART_SPP1_CC@idents)

#==============================================================================#
# CellChat analysis
#==============================================================================#

# Setting the ligand-receptor database
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

CC_objects <- list("all_CC" = all_CC,
                   "CART_CC" = CART_CC,
                   "CART_SPP1_CC" = CART_SPP1_CC)

CC_objects <- lapply(CC_objects, function(xx){
  # Set the used database in the object
  xx@DB <- CellChatDB.use
  
  # Preprocessing the expression data for cell-cell communication analysis
  # Subset the expression data of signaling genes for saving computation cost
  xx <- subsetData(xx) # This step is necessary even if using the whole database
  
  xx <- identifyOverExpressedGenes(xx)
  xx <- identifyOverExpressedInteractions(xx)
  
  # Project gene expression data onto PPI (Optional: when running it, USER should
  # set `raw.use = FALSE` in the function `computeCommunProb()` in order to use
  # the projected data)
  xx <- projectData(xx, PPI.human)
  
  ## Inference of cell-cell communication network
  # Compute the communication probability and infer cellular communication network
  xx <- computeCommunProb(xx) # raw.use = FALSE
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  xx <- filterCommunication(xx, min.cells = 10)
  
  # Infer the cell-cell communication at a signaling pathway level
  unique(xx@idents)
  
  xx <- computeCommunProbPathway(xx)
  
  # Calculate the aggregated cell-cell communication network
  xx <- aggregateNet(xx)
  
  xx
})

# Merging objects for comparative analysis
CC_objects_compare <- CC_objects[c("CART_SPP1_CC", "CART_CC")]
CC_merged_object <- mergeCellChat(CC_objects_compare, add.names = names(CC_objects_compare))

saveRDS(CC_objects_compare, "/scratch/hnatri/CART/CART_SPP1_Kluc_Cellchat_compare.rsd")
saveRDS(CC_merged_object, "/scratch/hnatri/CART/CART_SPP1_Kluc_Cellchat_merged.rds")

