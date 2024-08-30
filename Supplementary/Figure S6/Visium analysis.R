library(Seurat)
library(ggplot2)
library(dplyr)
library(viridis)
library(UCell)
library(patchwork)

#loading and processing lung cancer Visium HD dataset (Fig S7A-B)####

localdir <- "/Visium_HD_Human_Lung_Cancer_binned_outputs/binned_outputs/square_008um" 
object <- Load10X_Spatial(data.dir = localdir)
object <- NormalizeData(object)

#canonical markers
SpatialFeaturePlot(object, features = c("CD19", "MS4A1", "SELL", "CD27",
                                             "TNFRSF13B", "MT1X", "JCHAIN", "BCL6"),
                        pt.size.factor = 3, ncol = 4, image.alpha = 0.5) *
  scale_fill_viridis(option = "magma")

obj_b = subset(object, subset = CD19 > 0 | MS4A1 > 0)
obj_p = subset(object, subset = JCHAIN > 0)

lung_sig = read.csv("/LC_top10_0.2_RNA_rPCA.csv")
atlas_sig = read.csv("/top10_DEGs.csv")

b.gene.sets.l = list(naiveBcell_sig = c("CD19", "MS4A1", "SELL", atlas_sig[atlas_sig$cluster == "Naive B cells",]$gene, lung_sig[lung_sig$cluster == "Naive B cells",]$gene),
                     memoryBcell_sig = c("CD19", "MS4A1", "CD27", atlas_sig[atlas_sig$cluster == "Activated B cells",]$gene, lung_sig[lung_sig$cluster == "Activated B cells",]$gene),
                     atypicalBcell_sig = c("CD19", "MS4A1", "TNFRSF13B", atlas_sig[atlas_sig$cluster == "Atypical Memory B cells",]$gene, lung_sig[lung_sig$cluster == "Atypical Memory B cells",]$gene),
                     GC_sig = c("BCL6", atlas_sig[atlas_sig$cluster == "GC B cells",]$gene, lung_sig[lung_sig$cluster == "GC B cells",]$gene))
b.gene.sets.l <- lapply(b.gene.sets.l, function(genes) {
  genes[!grepl("^HLA-|ACT|MYO", genes)]
})


p.gene.sets.l = list(conv_plasma_sig = c("JCHAIN", atlas_sig[atlas_sig$cluster == "Conventional Plasma cells",]$gene, lung_sig[lung_sig$cluster == "Conventional Plasma cells",]$gene),
                     mt1x_plasma_sig = c("JCHAIN", "MT1X", atlas_sig[atlas_sig$cluster == "MT1X-high Plasma/Plasmablasts",]$gene, lung_sig[lung_sig$cluster == "MT1X-high Plasma/Plasmablasts",]$gene))
p.gene.sets.l <- lapply(p.gene.sets.l, function(genes) {
  genes[!grepl("^IGK|HB|IGHG1|APOE", genes)]
})

scores.l <- ScoreSignatures_UCell(obj_b@assays$Spatial$counts, features = b.gene.sets.l)
obj_b@meta.data$naiveBcell_sig_UCell <- scores.l[,"naiveBcell_sig_UCell"]
obj_b@meta.data$memoryBcell_sig_UCell <- scores.l[,"memoryBcell_sig_UCell"]
obj_b@meta.data$atypicalBcell_sig_UCell <- scores.l[,"atypicalBcell_sig_UCell"]
obj_b@meta.data$GC_sig_UCell <- scores.l[,"GC_sig_UCell"]

scores.l <- ScoreSignatures_UCell(obj_p@assays$Spatial$counts, features = p.gene.sets.l)
obj_p@meta.data$conv_plasma_sig_UCell <- scores.l[,"conv_plasma_sig_UCell"]
obj_p@meta.data$mt1x_plasma_sig_UCell <- scores.l[,"mt1x_plasma_sig_UCell"]

#atlas-derived + lung cancer-specific signatures
SpatialFeaturePlot(obj_b, features = c("naiveBcell_sig_UCell",
                                       "memoryBcell_sig_UCell",
                                       "atypicalBcell_sig_UCell",
                                       "GC_sig_UCell"),
                   pt.size.factor = 3, ncol = 2, image.alpha = 0.5, keep.scale = "all") *
  scale_fill_viridis(option = "magma", limits = c(0, 0.5))

SpatialFeaturePlot(obj_p, features = c("conv_plasma_sig_UCell",
                                       "mt1x_plasma_sig_UCell"),
                   pt.size.factor = 3, ncol = 2, image.alpha = 0.5, keep.scale = "all") *
  scale_fill_viridis(option = "magma", limits = c(0, 0.5))




#loading and processing breast cancer dataset (Fig S7C-D)####
data.dir <- "Visium_Breast"
img.dir <- "/Visium_Breast/spatial"
s13_img = Read10X_Image(img.dir, image.name = "tissue_hires_image.png")
s13 <- Load10X_Spatial(data.dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", image = s13_img)
s13@images[["slice1"]]@scale.factors[["lowres"]] = s13@images[["slice1"]]@scale.factors[["hires"]]
s13 <- SCTransform(s13, assay = "Spatial", verbose = FALSE)

s13_b = subset(s13, subset = CD19 > 0 | MS4A1 > 0)
s13_p = subset(s13, subset = JCHAIN > 0)

breast_sig = read.csv("/BC_0.2_RNA_rPCA.csv")
atlas_sig = read.csv("/top10_DEGs.csv")

b.gene.sets.4 = list(naiveBcell_sig = c("CD19", "MS4A1", "SELL", atlas_sig[atlas_sig$cluster == "Naive B cells",]$gene, breast_sig[breast_sig$cluster == "Naive B cells",]$gene),
                     memoryBcell_sig = c("CD19", "MS4A1", "CD27", atlas_sig[atlas_sig$cluster == "Activated Memory B cells",]$gene, breast_sig[breast_sig$cluster == "Activated memory B cells",]$gene[1:10]),
                     atypicalBcell_sig = c("CD19", "MS4A1", "TNFRSF13B", atlas_sig[atlas_sig$cluster == "Atypical Memory B cells",]$gene, breast_sig[breast_sig$cluster == "Atypical memory B cells",]$gene[1:10]),
                     GC_sig = c("BCL6", atlas_sig[atlas_sig$cluster == "Transitional/Proliferative 2",]$gene, breast_sig[breast_sig$cluster == "Transitional/proliferative 2",]$gene[1:10]))

b.gene.sets.4 <- lapply(b.gene.sets.4, function(genes) {
  genes[!grepl("^HLA-|ACT", genes)]
})


p.gene.sets.1 = list(conv_plasma_sig = c("JCHAIN", atlas_sig[atlas_sig$cluster == "Conventional Plasma cells",]$gene, breast_sig[breast_sig$cluster == "Conventional plasma cells",]$gene[1:10]),
                     mt1x_plasma_sig = c("JCHAIN", "MT1X", atlas_sig[atlas_sig$cluster == "MT1X-high Plasma/Plasmablasts",]$gene, breast_sig[breast_sig$cluster == "MT1X-high plasma/plasmablasts",]$gene[1:10]))
p.gene.sets.1 <- lapply(p.gene.sets.1, function(genes) {
  genes[!grepl("^IGK", genes)]
})


scores4 <- ScoreSignatures_UCell(s13_b@assays$SCT@counts, features = b.gene.sets.4)
s13_b@meta.data$naiveBcell_sig_UCell <- scores4[,"naiveBcell_sig_UCell"]
s13_b@meta.data$memoryBcell_sig_UCell <- scores4[,"memoryBcell_sig_UCell"]
s13_b@meta.data$atypicalBcell_sig_UCell <- scores4[,"atypicalBcell_sig_UCell"]
s13_b@meta.data$GC_sig_UCell <- scores4[,"GC_sig_UCell"]

scores4 <- ScoreSignatures_UCell(s13_p@assays$SCT@counts, features = p.gene.sets.1)
s13_p@meta.data$conv_plasma_sig_UCell <- scores4[,"conv_plasma_sig_UCell"]
s13_p@meta.data$mt1x_plasma_sig_UCell <- scores4[,"mt1x_plasma_sig_UCell"]

#canonical markers
SpatialFeaturePlot(s13, features = c("CD19", "MS4A1", "SELL", "CD27",
                                     "TNFRSF13B", "MT1X", "JCHAIN", "BCL6"),
                   pt.size.factor = 2.5, ncol = 4, image.alpha = 0.5) *
  scale_fill_viridis(option = "magma")

#atlas-derived signatures
SpatialFeaturePlot(s13_b, features = c("naiveBcell_sig_UCell",
                                       "memoryBcell_sig_UCell",
                                       "atypicalBcell_sig_UCell",
                                       "GC_sig_UCell"),
                   pt.size.factor = 2.5, image.alpha = 0.5, keep.scale = "all") *
  scale_fill_viridis(option = "magma", limits = c(0, 1))

SpatialFeaturePlot(s13_p, features = c("conv_plasma_sig_UCell",
                                       "mt1x_plasma_sig_UCell"),
                   pt.size.factor = 2.5, image.alpha = 0.5, keep.scale = "all") *
  scale_fill_viridis(option = "magma", limits = c(0, 1))
