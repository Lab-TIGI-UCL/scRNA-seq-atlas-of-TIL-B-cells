#Dotplot

Bcell_integrated <- RenameIdents(object = Bcell_integrated, '3' = "Naive B cells", '1' = "Activated Memory B cells", '0' = "Resting Memory B cells", 
                                 '6' = "IGKC-high Plasma/Plasmablasts", '2' = "Conventional Plasma cells", '4' = "Stressed Plasma cells", 
                                 '5' = "Atypical Memory B cells", '7' = "Transitional/Proliferative 1", '8' = "Transitional/Proliferative 2", 
                                 '9' = "MT1X-high Plasma/Plasmablasts")

levels(Bcell_integrated) <- c("Naive B cells", "Resting Memory B cells", "Activated Memory B cells", "Atypical Memory B cells", "Transitional/Proliferative 2", 
                              "Transitional/Proliferative 1", "MT1X-high Plasma/Plasmablasts", "IGKC-high Plasma/Plasmablasts", "Conventional Plasma cells", 
                              "Stressed Plasma cells")


features <- c("CD19", "MS4A1", "IL4R", "TCL1A", "FCER2", "FCMR", "SELL", "LY86", "BLK", "LTB", "BANK1", "CD69", "CCR7", "CD83", "NR4A2", "CHMP1B", "NFKBID", 
              "DUSP4", "CD82", "CD86", "ISG15", "CAPG", "ACP5", "COTL1", "IFITM1", "DRAP1", "TNFRSF13B", "LGALS9", "SHFL", "PLSCR1", "MX2", "IFI6", "IFI44L", 
              "HLA-DPB1", "HLA-DRA", "HLA-DPA1", "SUGCT", "BCL6", "SSBP2", "STAG3", "RGS13", "LCK", "MYO1E", "IRAG2", "HLA_DMB", "LAT2", "RRAS2", "MARCKSL1", 
              "SERPINA9", "PRPSAP2", "LMO2", "CD22","H4C3", "VPREB3", "STMN1", "HMGB2", "EIF6", "PSMA7", "CRIP1", "PTMA", "TYMS", "TUBB", "AURKB", "UBE2C", 
              "MKI67", "TUBA1B", "MT1X","MT2A", "MT1E", "MT1F", "IGKC", "IGHM", "IGHD", "IGHA1","IGHA2", "IGHG1", "IGHG2", "IGHG4", "IGHGP", "IGLC2", "IGLC3", 
              "DERL3", "FKBP11", "PRDX4", "TMED2", "DPP7", "FIS1", "ISG20", "EZR", "PSME2", "PRDM1", "MCL1", "MZB1", "XBP1", "VIM", "SSR4", "JCHAIN", "JUN", 
              "CD38", "SDC1", "CD27", "FOS", "DNAJB1","DUSP1", "ZFAND2A", "HSPA1B", "HSPA1A", "HSPB1", "HSPA6")


#Load Module IDs
Module_IDs <- read.csv("~/Module_IDs.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
GeneList <- vector(mode = "list", length = length(unique(Module_IDs$module)))

for (i in 1:length(unique(Module_IDs$module)))
{
  names(GeneList)[[i]] <- unique(Module_IDs$module)[i]
  GeneList[[i]] <- Module_IDs[Module_IDs$module == unique(Module_IDs$module)[i],"genes"]
}

GeneList

DotPlot(Bcell_integrated, features = GeneList) + theme(axis.text.y = element_text(size=10)) + 
  theme(axis.text.x = element_text(angle= 90, vjust = 0.5, hjust = 1, size=8)) + 
  scale_colour_viridis(option="magma")

