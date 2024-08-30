#QC Metrics

library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(viridisLite)
library(viridis)

Bcell_integrated<-readRDS("Tumourlowres_subsetted_integrated.rds")
DefaultAssay(Bcell_integrated) <- "SCT"
Idents(object = Bcell_integrated) <- "integrated_snn_res.0.2"
Bcell_integrated <- RenameIdents(object = Bcell_integrated, '3' = "Naive B cells", '1' = "Activated Memory B cells", '0' = "Resting Memory B cells", '6' = "IGKC-high Plasma/Plasmablasts", '2' = "Conventional Plasma cells", '4' = "Stressed Plasma cells", '5' = "Atypical Memory B cells", '7' = "Transitional/Proliferative 1", '8' = "Transitional/Proliferative 2", '9' = "MT1X-high Plasma/Plasmablasts")


VlnPlot(Bcell_integrated, features = c("nFeature_RNA", "percent.mt", "nCount_RNA"), pt.size = 0) + theme(text = element_text(size = 15)) 



### Elbow plot to choose UMAP resolution

library(ggplot2)
a <- SSE_results
plot(a$nclust, a$var, xlab = "Number of clusters", ylab = "var")
lines(a$nclust[order(a$nclust)], a$var[order(a$nclust)], xlim=range(a$nclust), ylim=range(a$var), pch=16)


ggplot(a, aes(x = nclust, y = var)) +
  geom_line() + geom_point()+
  scale_x_continuous(breaks = 1:10)
