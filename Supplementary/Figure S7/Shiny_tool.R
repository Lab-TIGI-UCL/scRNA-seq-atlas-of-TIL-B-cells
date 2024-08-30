#Creation of Shiny Tool


reqPkg = c("data.table", "Matrix", "hdf5r", "reticulate", "ggplot2", 
           "gridExtra", "glue", "readr", "RColorBrewer", "R.utils", "Seurat")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}


reqPkg = c("shiny", "shinyhelper", "data.table", "Matrix", "DT", "hdf5r", 
           "reticulate", "ggplot2", "gridExtra", "magrittr", "ggdendro")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}

devtools::install_github("SGDDNB/ShinyCell")

library(ShinyCell)
library(Seurat)

Bcell_integrated = readRDS("~/Tumourlowres_subsetted_integrated.rds")

bcell_meta = Bcell_integrated@meta.data
head(rownames(bcell_meta))
bcell_meta$cellids_rownames = rownames(bcell_meta)
merged = merge(x = bcell_meta, by.x = "cellids_rownames",
               y = subsetted_meta[, c("cells", "cancer_type", "patient", "treatment", "treatment_point", "platform", "tissue")], by.y = "cells")
Bcell_integrated@meta.data = merged

orderings = match(colnames(Bcell_integrated), subsetted_meta$cells)
new_order = Bcell_integrated@meta.data[orderings, ]
all(colnames(Bcell_integrated) == new_order$cells)#should be true
rownames(new_order) = new_order$cellids_rownames
Bcell_integrated@meta.data = new_order

Bcell_integrated@meta.data$res.0.2 <- NULL
Bcell_integrated@meta.data$res.0.3 <- NULL
Bcell_integrated@meta.data$res.0.4 <- NULL
Bcell_integrated@meta.data$res.0.5 <- NULL
Bcell_integrated@meta.data$res.0.6 <- NULL
Bcell_integrated@meta.data$res.0.7 <- NULL
Bcell_integrated@meta.data$res.0.8 <- NULL
Bcell_integrated@meta.data$res.0.9 <- NULL
Bcell_integrated@meta.data$cell <- NULL
Bcell_integrated@meta.data$region <- NULL
Bcell_integrated@meta.data$type <- NULL
Bcell_integrated@meta.data$Sample <- NULL
Bcell_integrated@meta.data$Sample2 <- NULL
Bcell_integrated@meta.data$b1 <- NULL
Bcell_integrated@meta.data$b2 <- NULL
Bcell_integrated@meta.data$b3 <- NULL
Bcell_integrated@meta.data$b4 <- NULL
Bcell_integrated@meta.data$b5 <- NULL
Bcell_integrated@meta.data$b6 <- NULL
Bcell_integrated@meta.data$cluster <- NULL
Bcell_integrated@meta.data$well <- NULL
Bcell_integrated@meta.data$plate <- NULL
Bcell_integrated@meta.data$cell_id <- NULL
Bcell_integrated@meta.data$sample_name <- NULL
Bcell_integrated@meta.data$patient_id <- NULL
Bcell_integrated@meta.data$smokingHx <- NULL
Bcell_integrated@meta.data$histolgy <- NULL
Bcell_integrated@meta.data$Notes <- NULL
Bcell_integrated@meta.data$driver_gene <- NULL
Bcell_integrated@meta.data$driver_mutation <- NULL
Bcell_integrated@meta.data$stage.at.dx <- NULL
Bcell_integrated@meta.data$pathlogy_review <- NULL
Bcell_integrated@meta.data$biopsy_type <- NULL
Bcell_integrated@meta.data$biopsy_site <- NULL
Bcell_integrated@meta.data$primary_or_metastaic <- NULL
Bcell_integrated@meta.data$early_treatment_status <- NULL
Bcell_integrated@meta.data$best_rxn_status <- NULL
Bcell_integrated@meta.data$biopsy_time_status <- NULL
Bcell_integrated@meta.data$biopsy_timing <- NULL
Bcell_integrated@meta.data$treatment_history <- NULL
Bcell_integrated@meta.data$analysis <- NULL
Bcell_integrated@meta.data$treatment_history_detail <- NULL
Bcell_integrated@meta.data$line_of_therapy <- NULL
Bcell_integrated@meta.data$treatment_type <- NULL
Bcell_integrated@meta.data$treatment.x <- NULL
Bcell_integrated@meta.data$percent_PFS_ref_values <- NULL
Bcell_integrated@meta.data$percent_PFS_reference_values <- NULL
Bcell_integrated@meta.data$early_bx_day <- NULL
Bcell_integrated@meta.data$infections <- NULL
Bcell_integrated@meta.data$gender <- NULL
Bcell_integrated@meta.data$race <- NULL
Bcell_integrated@meta.data$secondary_mutation <- NULL
Bcell_integrated@meta.data$percent.PFS.reference.values <- NULL
Bcell_integrated@meta.data$pfs_over_under <- NULL
Bcell_integrated@meta.data$pfs_day <- NULL
Bcell_integrated@meta.data$pfs_month <- NULL
Bcell_integrated@meta.data$ca_dx_OS <- NULL
Bcell_integrated@meta.data$percent.ercc <- NULL
Bcell_integrated@meta.data$percent.ribo <- NULL
Bcell_integrated@meta.data$sample_type <- NULL
Bcell_integrated@meta.data$sort_plate_number <- NULL
Bcell_integrated@meta.data$Sequence_Run1 <- NULL
Bcell_integrated@meta.data$stage <- NULL
Bcell_integrated@meta.data$treatment_status <- NULL
Bcell_integrated@meta.data$RNA_snn_res.0.1 <- NULL
Bcell_integrated@meta.data$RNA_snn_res.0.3 <- NULL
Bcell_integrated@meta.data$RNA_snn_res.0.5 <- NULL
Bcell_integrated@meta.data$RNA_snn_res.0.7 <- NULL
Bcell_integrated@meta.data$RNA_snn_res.0.9 <- NULL
Bcell_integrated@meta.data$main_seurat_cluster <- NULL
Bcell_integrated@meta.data$immune_annotation <- NULL
Bcell_integrated@meta.data$pANN_0.25_0.09_218 <- NULL
Bcell_integrated@meta.data$DF.classifications_0.25_0.09_218 <- NULL
Bcell_integrated@meta.data$nonimmune_general_annotation <- NULL
Bcell_integrated@meta.data$general_annotation <- NULL
Bcell_integrated@meta.data$general_annotation1 <- NULL
Bcell_integrated@meta.data$Final_annotation <- NULL
Bcell_integrated@meta.data$integrated_snn_res.0.1 <- NULL
Bcell_integrated@meta.data$integrated_snn_res.0.4 <- NULL
Bcell_integrated@meta.data$integrated_snn_res.0.5 <- NULL
Bcell_integrated@meta.data$integrated_snn_res.0.6 <- NULL
Bcell_integrated@meta.data$integrated_snn_res.0.8 <- NULL
Bcell_integrated@meta.data$integrated_snn_res.0.9 <- NULL
Bcell_integrated@meta.data$integrated_snn_res.1.1 <- NULL
Bcell_integrated@meta.data$integrated_snn_res.1.3 <- NULL
Bcell_integrated@meta.data$integrated_snn_res.1.4 <- NULL
Bcell_integrated@meta.data$integrated_snn_res.1.6 <- NULL
Bcell_integrated@meta.data$integrated_snn_res.1.8 <- NULL
Bcell_integrated@meta.data$integrated_snn_res.1.9 <- NULL
Bcell_integrated@meta.data$integrated_snn_res.2.5 <- NULL
Bcell_integrated@meta.data$seurat_clusters <- NULL

Idents(object = Bcell_integrated) <- "integrated_snn_res.0.2"
Bcell_integrated <- RenameIdents(object = Bcell_integrated, '3' = "Naive B cells", '1' = "Activated B cells", '0' = "Resting Memory B cells", '6' = "IGKC-high Plasma/Plasmablasts", '2' = "Conventional Plasma cells", '4' = "Stressed Plasma cells", '5' = "Atypical Memory B cells", '7' = "Proliferative B cells", '8' = "GC B cells", '9' = "MT1X-high Plasma/Plasmablasts")

saveRDS(Bcell_integrated, "Tidied_Bcell.rds")

scConf = createConfig(seu)

scConf <- scConf[c(1,2,3,4,5,6,82, 94,95,96,97,98),]

makeShinyApp(seu, scConf, gene.mapping = TRUE, gex.assay = "SCT",
             shiny.title = "ShinyCell Quick Start") 
