#UMAP

library(Seurat)
library(ggplot2)
library(dplyr)

Bcell_integrated<-readRDS("Tumourlowres_subsetted_integrated.rds")
DefaultAssay(Bcell_integrated) <- "SCT"
Idents(object = Bcell_integrated) <- "integrated_snn_res.0.2"
Bcell_integrated <- RenameIdents(object = Bcell_integrated, '3' = "Naive B cells", '1' = "Activated Memory B cells", '0' = "Resting Memory B cells",
                                 '6' = "IGKC-high Plasma/Plasmablasts", '2' = "Conventional Plasma cells", '4' = "Stressed Plasma cells", 
                                 '5' = "Atypical Memory B cells", '7' = "Transitional/Proliferative 1", '8' = "Transitional/Proliferative 2", 
                                 '9' = "MT1X-high Plasma/Plasmablasts")


umap.df <- FetchData(Bcell_integrated, vars = c("UMAP_1", "UMAP_2", "integrated_snn_res.0.2", "ident"))

ggplot(umap.df, aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point(aes(fill = ident), size = 4.5, shape = 21) + 
  geom_label_repel(umap.df %>% group_by(ident) %>% summarise_at(c("UMAP_1", "UMAP_2"), mean), mapping = aes(x = UMAP_1, y = UMAP_2, label = ident)) + 
  scale_fill_manual(values = c("Naive B cells" = "#E66C37","Resting Memory B cells" = "#D64550",
                               "Activated Memory B cells" = "#FF977E","Atypical Memory B cells" = "#A43B76",
                               "Transitional/Proliferative 1" = "#9A64A0","Transitional/Proliferative 2"="#750985",
                               "Conventional Plasma cells" = "#5ECBC8", "Stressed Plasma cells" = "#28788D", "IGKC-high Plasma/Plasmablasts" = "#37A794", 
                               "MT1X-high Plasma/Plasmablasts" = "#4C5D8A")) + 
  guides(colour = guide_legend(override.aes = list(size=2))) + 
  theme_classic() + NoLegend() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + theme(text = element_text(size = 15))