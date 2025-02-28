#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 01/04/2025
# Description: CAR+SPP1 mouse tumor colors and themes
#==============================================================================#

library(ggplot2)
library(RColorBrewer)
library(plyr)
library(circlize)
library(googlesheets4)

#==============================================================================
# Environment variables
#==============================================================================

set.seed(1234)

#==============================================================================
# Colors and themes
#==============================================================================

# Colors for plotting
# Define colors for each level of categorical variables

# Cluster colors
data_clusters <- as.factor(c(0, seq(1:16)))

carspp1_cluster_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(data_clusters))
names(carspp1_cluster_col) <- levels(data_clusters)

# Cell type colors
# Cluster annotations
gs4_deauth()
cluster_annot  <- gs4_get("https://docs.google.com/spreadsheets/d/127J6C4KF7uBGKUnrPuC1mcsb_wNCN6k1zXKSCbJ6q0M/edit?usp=sharing")
cluster_annot <- read_sheet(cluster_annot, sheet = "Cluster annotation")

carspp1_celltypes <- cluster_annot$annot
carspp1_celltype_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(carspp1_celltypes))
names(carspp1_celltype_col) <- carspp1_celltypes

myeloid_clusters <- as.factor(c(0, seq(1:16)))
myeloid_cluster_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(myeloid_clusters))
names(myeloid_cluster_col) <- myeloid_clusters

myeloid_celltypes <- c("M1_suppressive", "M6", "M9", "M2_suppressive", "M4",
                       "M5_suppressive", "M7", "M3_suppressive",
                       "M8_suppressive_G2MS", "M10", "M11")
#myeloid_colors <- tinter("aquamarine3", steps = 8, crop = 3,
#                         direction = "both", adjust = 0)
myeloid_colors <- colorRampPalette(c("aquamarine1", "darkgreen"))(nb.cols <- length(myeloid_celltypes))
names(myeloid_colors) <- myeloid_celltypes
#myeloid_colors <- myeloid_colors[myeloid_celltypes]

lymphoid_clusters <- as.factor(c(0, seq(1:10)))
lymphoid_cluster_col <- colorRampPalette(brewer.pal(10, "Paired"))(nb.cols <- length(lymphoid_clusters))
names(lymphoid_cluster_col) <- lymphoid_clusters

# Sample type
#sample_type_col <- c("Mock" = "azure3",
#                     "Pre_CAR" = "lightskyblue3",
#                     "5050_CAR" = "darkslategray2",
#                     "CAR" = "turquoise2")


