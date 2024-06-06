#Merge Clusters

library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(viridisLite)
library(viridis)

Bcell_integrated<-readRDS("Tumourlowres_subsetted_integrated.rds")
DefaultAssay(Bcell_integrated) <- "SCT"
Idents(object = Bcell_integrated) <- "integrated_snn_res.0.2"

Bcell_integrated <- RenameIdents(object = Bcell_integrated, '3' = "Naive B cells", '1' = "Memory B cells", '0' = "Memory B cells", '6' = "Plasma cells", '2' = "Plasma cells", '4' = "Plasma cells", '5' = "Memory B cells", '7' = "Transitional/Proliferative cells", '8' = "Transitional/Proliferative cells", '9' = "Plasma cells")

umap.df <- FetchData(Bcell_integrated, vars = c("UMAP_1", "UMAP_2", "integrated_snn_res.0.2", "ident"))

ggplot(umap.df, aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point(aes(fill = ident), size = 4.5, shape = 21) + 
  geom_label_repel(umap.df %>% group_by(ident) %>% summarise_at(c("UMAP_1", "UMAP_2"), mean), mapping = aes(x = UMAP_1, y = UMAP_2, label = ident)) + 
  scale_fill_manual(values = c("Naive B cells" = "#FF977E","Memory B cells" = "#D64550",
                               "Transitional/Proliferative cells" = "#9A64A0",
                               "Plasma cells" = "#5ECBC8")) + 
  guides(colour = guide_legend(override.aes = list(size=2))) + 
  theme_classic() + NoLegend() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + theme(text = element_text(size = 15))


ggsave(
  "UMAP_paper_merged.pdf",
  plot = last_plot(),
  width = 250,
  height = 200,
  units = "mm",
  dpi = 300) 


levels(Bcell_integrated) <- c("Naive B cells", "Memory B cells", "Transitional/Proliferative cells", "Plasma cells")
features <- c("CD19", "MS4A1", "TCL1A", "FCER2", "FCMR", "SELL", "CD69", "CD27", "CCR7", "CD83", "NR4A2", "CHMP1B", "NFKBID", "DUSP4", "CD82", "CD86", "ISG15", "IFITM1", "DRAP1", "TNFRSF13B", "LGALS9", "HLA-DPB1", "HLA-DRA", "HLA-DPA1", "BCL6", "RGS13", "LCK", "MYO1E", "IRAG2", "MARCKSL1", "SERPINA9", "CD22", "AURKB", "UBE2C", "MKI67", "TUBA1B", "MT1X","MT2A", "IGKC", "IGHA1", "IGHG1", "IGHG4", "DERL3", "FKBP11", "MZB1", "XBP1", "JCHAIN", "JUN", "CD38", "SDC1", "HSPA1B", "HSPA1A")
DEG <- c("CD19", "TCL1A", "FCER2", "SELL",  "CNPY3", "SNX9", "ST13", "FAM177A1", "UCP2", "EIF2S3", "CDC37", "CCT7", "TMEM230", "CCT4", "CD27", "CD69", "CD82", "CCR7", "REL", "RILPL2", "SLC2A3", "BCL2A1", "MS4A1", "CAPG", "SPIB", "H4C3", "VPREB3", "STMN1", "EIF6", "HMGB2", "PSMA7", "FABP5", "CRIP1", "PTMA", "TUBB", "IGLV3-1", "IGHG1", "IGHG3", "IGLC2", "IGKC", "IGHA1", "IGLC3", "IGHG4", "XBP1", "IGLL5", "SDC1", "JCHAIN")
DotPlot(Bcell_integrated, features = DEG) + theme(axis.text.y = element_text(size=10)) + theme(axis.text.x = element_text(angle= 90, vjust = 0.5, hjust = 1, size=8)) + scale_colour_viridis(option="magma")

