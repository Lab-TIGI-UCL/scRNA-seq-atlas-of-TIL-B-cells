#Differential Gene Expression


library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(dplyr)

Bcell_integrated<-readRDS("Tumourlowres_subsetted_integrated.rds")

Idents(object = Bcell_integrated) <- "integrated_snn_res.0.2"
DefaultAssay(Bcell_integrated) <- "RNA"
markers <- FindAllMarkers(Bcell_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers,file="tumoursubset__0.2_RNA_rPCA.csv")
top10<-markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
write.csv(top10,file="tumourtop10subset_0.2_RNA_rPCA.csv")