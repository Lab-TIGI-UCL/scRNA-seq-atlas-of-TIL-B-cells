#Mapping signatures back to each study and cancer type

library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

Bcell_integrated<-readRDS("Tumourlowres_subsetted_integrated.rds")
DefaultAssay(Bcell_integrated) <- "SCT"
Idents(object = Bcell_integrated) <- "integrated_snn_res.0.2"
Bcell_integrated <- RenameIdents(object = Bcell_integrated, '3' = "Naive B cells", '1' = "Activated Memory B cells", '0' = "Resting Memory B cells", '6' = "IGKC-high Plasma/Plasmablasts", '2' = "Conventional Plasma cells", '4' = "Stressed Plasma cells", '5' = "Atypical Memory B cells", '7' = "Transitional/Proliferative 1", '8' = "Transitional/Proliferative 2", '9' = "MT1X-high Plasma/Plasmablasts")

#Naive
sig = list(c('IL4R', 'TCL1A', 'FCER2', 'FCMR', 'SELL'))
Bcell_integrated <- AddModuleScore(object = Bcell_integrated, features = sig, name = "naive_sig")
FeaturePlot(object = Bcell_integrated, 
            features = "naive_sig1", 
            raster = "true", 
            split.by = "cancer_type",
            order = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(object = Bcell_integrated, 
            features = "naive_sig1", 
            raster = "true", 
            split.by = "study",
            order = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


#Resting Memory
sig = list(c('LY86', 'BLK', 'LTB', 'MARCHF1'))
Bcell_integrated <- AddModuleScore(object = Bcell_integrated, features = sig, name = "restingmem_sig")
FeaturePlot(object = Bcell_integrated, 
            features = "restingmem_sig1", 
            raster = "true", 
            split.by = "cancer_type",
            order = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(object = Bcell_integrated, 
            features = "restingmem_sig1", 
            raster = "true", 
            split.by = "study",
            order = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


#Activated Memory
sig = list(c('CD69', 'CCR7', 'CD83', 'SLC2A3', 'REL', 'NR4A2', 'CHMP1B', 'NFKBID'))
Bcell_integrated <- AddModuleScore(object = Bcell_integrated, features = sig, name = "activatedmem_sig")
FeaturePlot(object = Bcell_integrated, 
            features = "activatedmem_sig1", 
            raster = "true", 
            split.by = "cancer_type",
            order = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(object = Bcell_integrated, 
            features = "activatedmem_sig1", 
            raster = "true", 
            split.by = "study",
            order = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


#Atypical Memory
sig = list(c('DUSP4', 'CD82', 'CD86', 'ISG15', 'CAPG', 'ACP5', 'COTL1', 'IFITM1', 'DRAP1', 'TNFRSF13B', 'LGALS9', 'SHFL', 'PLSCR1', 'MX2', 'IFI6', 'IFI44L'))
Bcell_integrated <- AddModuleScore(object = Bcell_integrated, features = sig, name = "atypicalmem_sig")
FeaturePlot(object = Bcell_integrated, 
            features = "atypicalmem_sig1", 
            raster = "true", 
            split.by = "cancer_type",
            order = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(object = Bcell_integrated, 
            features = "atypicalmem_sig1", 
            raster = "true", 
            split.by = "study",
            order = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


#Transitional/Proliferative 2
sig = list(c('SUGCT', 'BCL6', 'SSBP2', 'STAG3', 'RGS13', 'LCK', 'MYO1E', 'IRAG2', 'LAT2', 'RRAS2', 'MARCKSL1', 'SERPINA9', 'LMO2', 'CD22'))
Bcell_integrated <- AddModuleScore(object = Bcell_integrated, features = sig, name = "transprolif2_sig")
FeaturePlot(object = Bcell_integrated, 
            features = "transprolif2_sig1", 
            raster = "true", 
            split.by = "cancer_type",
            order = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(object = Bcell_integrated, 
            features = "transprolif2_sig1", 
            raster = "true", 
            split.by = "study",
            order = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


#Transitional/Proliferative 1
sig = list(c('H4C3', 'VPREB3', 'STMN1', 'HMGB2', 'EIF6', 'PSMA7', 'CRIP1', 'PTMA', 'TYMS', 'TUBB', 'AURKB', 'UBE2C', 'MKI67', 'TUBA1B'))
Bcell_integrated <- AddModuleScore(object = Bcell_integrated, features = sig, name = "transprolif1_sig")
FeaturePlot(object = Bcell_integrated, 
            features = "transprolif1_sig1", 
            raster = "true", 
            split.by = "cancer_type",
            order = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(object = Bcell_integrated, 
            features = "transprolif1_sig1", 
            raster = "true", 
            split.by = "study",
            order = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))



#MT1X-high Plasma/Plasmablasts
sig = list(c('MT1X', 'MT2A', 'MT1E', 'MT1F'))
Bcell_integrated <- AddModuleScore(object = Bcell_integrated, features = sig, name = "MT1Xplasma_sig")
FeaturePlot(object = Bcell_integrated, 
            features = "MT1Xplasma_sig1", 
            raster = "true", 
            split.by = "cancer_type",
            order = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(object = Bcell_integrated, 
            features = "MT1Xplasma_sig1", 
            raster = "true", 
            split.by = "study",
            order = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


#IGKC-high Plasma/Plasmablasts
sig = list(c('IGKC'))
Bcell_integrated <- AddModuleScore(object = Bcell_integrated, features = sig, name = "IGKCplasma_sig")
FeaturePlot(object = Bcell_integrated, 
            features = "IGKC", 
            raster = "true", 
            split.by = "cancer_type",
            order = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(object = Bcell_integrated, 
            features = "IGKC", 
            raster = "true", 
            split.by = "study",
            order = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


#Conventional Plasma cells
sig = list(c('DERL3', 'FKBP11', 'PRDX4', 'TMED2', 'DPP7', 'FIS1', 'ISG20', 'EZR', 'PSME2', 'PRDM1', 'MCL1', 'MZB1', 'XBP1', 'VIM', 'SSR4', 'JCHAIN', 'JUN', 'CD38', 'SDC1', 'CD27'))
Bcell_integrated <- AddModuleScore(object = Bcell_integrated, features = sig, name = "Convplasma_sig")
FeaturePlot(object = Bcell_integrated, 
            features = "Convplasma_sig1", 
            raster = "true", 
            split.by = "cancer_type",
            order = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(object = Bcell_integrated, 
            features = "Convplasma_sig1", 
            raster = "true", 
            split.by = "study",
            order = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


#Stressed Plasma cells
sig = list(c('FOS', 'DNAJB1', 'DUSP1', 'ZFAND2A', 'HSPA1B', 'HSPA1A', 'HSPB1', 'HSPA6'))
Bcell_integrated <- AddModuleScore(object = Bcell_integrated, features = sig, name = "Stressplasma_sig")
FeaturePlot(object = Bcell_integrated, 
            features = "Stressplasma_sig1", 
            raster = "true", 
            split.by = "cancer_type",
            order = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

FeaturePlot(object = Bcell_integrated, 
            features = "Stressplasma_sig1", 
            raster = "true", 
            split.by = "study",
            order = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
