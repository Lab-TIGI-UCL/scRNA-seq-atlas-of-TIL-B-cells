#Excluded clusters

Bcell_integrated<-readRDS("Tumour0.2_B_integrated.rds")

Idents(object = Bcell_integrated) <- "integrated_snn_res.0.4"
Bcell_integrated <- RenameIdents(object = Bcell_integrated, '3' = "Naive B cells", '1' = "Activated Memory B cells", '0' = "Resting Memory B cells", '6' = "Atypical Memory B cells", '2' = "Conventional Plasma cells", '4' = "Stressed Plasma cells", '5' = "IGKC-high Plasma/Plasmablasts" , '7' = "T/B cell doublets", '8' = "Transitional/Proliferative 1", '9' = "Transitional/Proliferative 2", '11' = "MNP/B cell doublets", '13' = "Artifacts", '12' = 'Technical', '10' = "MT1X-high Plasma/Plasmablasts")#Old version

umap.df <- FetchData(Bcell_integrated, vars = c("UMAP_1", "UMAP_2", "integrated_snn_res.0.4", "ident"))

ggplot(umap.df, aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point(aes(fill = ident), size = 3, shape = 21) + 
  geom_label_repel(umap.df %>% group_by(ident) %>% summarise_at(c("UMAP_1", "UMAP_2"), mean), mapping = aes(x = UMAP_1, y = UMAP_2, label = ident)) + 
  scale_fill_manual(values = c("Naive B cells" = "#E66C37","Resting Memory B cells" = "#D64550",
                               "Activated Memory B cells" = "#FF977E","Atypical Memory B cells" = "#A43B76",
                               "Transitional/Proliferative 1" = "#9A64A0","Transitional/Proliferative 2"="#750985",
                               "Conventional Plasma cells" = "#5ECBC8", "Stressed Plasma cells" = "#28788D", "IGKC-high Plasma/Plasmablasts" = "#37A794", 
                               "MT1X-high Plasma/Plasmablasts" = "#4C5D8A", "T/B cell doublets" = "#ff93ac", "MNP/B cell doublets" = "#b3d4cf", "Artifacts" = "#819090", "Technical" = "#f4d4d4")) + 
  guides(colour = guide_legend(override.aes = list(size=2))) + 
  theme_classic() + NoLegend() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + theme(text = element_text(size = 15))

ggsave(
  "UMAP_withdoublets.pdf",
  plot = last_plot(),
  width = 250,
  height = 200,
  units = "mm",
  dpi = 300)


Idents(object = Bcell_integrated) <- "integrated_snn_res.0.4"
Bcell_integrated <- RenameIdents(object = Bcell_integrated, '3' = "Naive B cells", '1' = "Activated Memory B cells", '0' = "Resting Memory B cells", '6' = "Atypical Memory B cells", '2' = "Conventional Plasma cells", '4' = "Stressed Plasma cells", '5' = "IGKC-high Plasma/Plasmablasts" , '7' = "T cell doublets/misannotation", '8' = "Transitional/Proliferative 1", '9' = "Transitional/Proliferative 2", '11' = "MNP cell doublets/misannotation", '13' = "Artifacts", '12' = 'Technical', '10' = "MT1X-high Plasma/Plasmablasts")#Old version
levels(Bcell_integrated) <- c("Naive B cells", "Resting Memory B cells", "Activated Memory B cells", "Atypical Memory B cells", "Transitional/Proliferative 2", "Transitional/Proliferative 1", "MT1X-high Plasma/Plasmablasts", "IGKC-high Plasma/Plasmablasts", "Conventional Plasma cells", "Stressed Plasma cells", "T cell doublets/misannotation", "MNP cell doublets/misannotation", "Artifacts", "Technical")
features <- c("CD19", "MS4A1", "SDC1", "IGHG4","CD2", "CD3D", "CD3E", "CD3G", "LYZ", "FCER1G", "AIF1", "S100A9", "CLEC9A", "ARRB2", "GZMB", "CLIC3", "MPEG1", "ACTN4")
DotPlot(Bcell_integrated, features = features) + theme(axis.text.y = element_text(size=10)) + theme(axis.text.x = element_text(angle= 90, vjust = 0.5, hjust = 1, size=8)) + scale_colour_viridis(option="magma")



cluster_donors <- read.csv("~/Documents/cluster_donors.csv")
ident <- cluster_donors$Ident
cells <- cluster_donors$Number.of.Cells
donors <- cluster_donors$Number.of.Donors
data <- data.frame(ident, cells, donors)

melt <- melt(data, id.var = "ident")

ggplot(melt, aes(fill=variable, y=value, x=ident)) + 
  geom_bar(position="dodge", stat="identity")

ggplot(data=data, mapping=aes(x=cells, y=donors)) +
  geom_point() + 
  theme_bw() +
  geom_label_repel(aes(x=cells, y=donors, label = ident))

