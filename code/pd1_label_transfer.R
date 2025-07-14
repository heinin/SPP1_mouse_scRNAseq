## Label transfer from the SPP1 data to PD1 and vice versa

library(Seurat)

pd1_seurat <- readRDS("/tgen_labs/banovich/BCTCSF/PD1_mm_Seurat/PD1_mm_Seurat_merged.Rds")
spp1_seurat <- readRDS("/tgen_labs/banovich/BCTCSF/SPP1_mouse_scRNAseq/scRNAseq_Seurat_dim8_annotated.rds")

DefaultAssay(pd1_seurat) <- "RNA"
DefaultAssay(spp1_seurat) <- "RNA"

pd1_seurat <- NormalizeData(pd1_seurat)
pd1_seurat <- FindVariableFeatures(pd1_seurat)
pd1_seurat <- ScaleData(pd1_seurat)
pd1_seurat <- RunPCA(pd1_seurat)

spp1_seurat <- NormalizeData(spp1_seurat)
spp1_seurat <- FindVariableFeatures(spp1_seurat)
spp1_seurat <- ScaleData(spp1_seurat)
spp1_seurat <- RunPCA(spp1_seurat)
#spp1_seurat <- RunUMAP()

# SPP1 to PD1
anchors <- FindTransferAnchors(reference = spp1_seurat,
                               query = pd1_seurat,
                               dims = 1:30,
                               reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors,
                            refdata = spp1_seurat$annot_granular,
                            dims = 1:30)
pd1_seurat <- AddMetaData(pd1_seurat,
                          metadata = predictions)

#spp1_seurat <- RunUMAP(spp1_seurat, dims = 1:30, reduction = "integrated.cca", return.model = TRUE)
pd1_seurat <- MapQuery(anchorset = anchors,
                       reference = spp1_seurat,
                       query = pd1_seurat,
                       refdata = list(annot_granular = "annot_granular"),
                       reference.reduction = "pca",
                       reduction.model = "integratedSCTumap")

DimPlot(pd1_seurat,
        group.by = "predicted.annot_granular",
        reduction = "umap",
        label = T) +
  coord_fixed() +
  theme_classic() +
  NoLegend()

# PD1 to SPP1
anchors <- FindTransferAnchors(reference = pd1_seurat,
                               query = spp1_seurat,
                               dims = 1:30,
                               reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors,
                            refdata = pd1_seurat$celltype,
                            dims = 1:30)
spp1_seurat <- AddMetaData(spp1_seurat,
                           metadata = predictions)

#pd1_seurat <- RunUMAP(pd1_seurat, dims = 1:30, reduction = "integrated.cca", return.model = TRUE)
spp1_seurat <- MapQuery(anchorset = anchors,
                        reference = pd1_seurat,
                        query = spp1_seurat,
                        refdata = list(celltype = "celltype"),
                        reference.reduction = "pca",
                        reduction.model = "umap")

head(spp1_seurat)

#saveRDS(spp1_seurat, "/tgen_labs/banovich/BCTCSF/SPP1_mouse_scRNAseq/scRNAseq_Seurat_dim8_annotated.rds")
spp1_seurat <- readRDS("/tgen_labs/banovich/BCTCSF/SPP1_mouse_scRNAseq/scRNAseq_Seurat_dim8_annotated.rds")

DimPlot(spp1_seurat,
        group.by = "predicted.celltype",
        reduction = "integratedSCTumap",
        label = T) +
  coord_fixed() +
  theme_classic() +
  NoLegend()

