#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 09/06/2025
# Description: Saving celltype top markers
#==============================================================================#

#==============================================================================#
# Loading libraries
#==============================================================================#

library(tidyverse)
library(Seurat)

#==============================================================================#
# Environment variables and helper functions
#==============================================================================#

setwd("/home/hnatri/SPP1_mouse_scRNAseq/")
set.seed(1234)

#==============================================================================#
# Import data and save top markers
#==============================================================================#

seurat_object <- readRDS("/tgen_labs/banovich/BCTCSF/SPP1_mouse_scRNAseq/scRNAseq_Seurat_dim8_granular_annot.rds")

markers <- presto::wilcoxauc(seurat_object,
                             group_by = "annot_granular",
                             assay = "data",
                             seurat_assay = "RNA_human")

output_cluster_markers <- markers %>%
  arrange(dplyr::desc(logFC)) %>%
  group_by(group) %>%
  dplyr::slice(1:30)

write.table(output_cluster_markers, "/home/hnatri/SPP1_mouse_scRNAseq/annot_granular_top30_markers.tsv",
            quote = F, row.names = F, sep = "\t")

