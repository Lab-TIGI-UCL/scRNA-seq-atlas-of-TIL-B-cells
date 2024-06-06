#Immune cell types in each study from full metadata

Full_metadata<-readRDS("~/combined_metadata2.rds")

#Subset to studies used
table(Full_metadata$study)
library(dplyr)    
subset_meta <- Full_metadata %>% filter(study %in% c("Azizi", "Bi", "BRAUN", "Chan", "KIM", "KRISHNA", "Li", "Maynard", "Pelka", "QIAN", "VISHWAKARMA", "Wu", "YOST", "ZhangYY", "Zilionis"))
table(subset_meta$study)

#Plot immune cell types in each study

library(pals)

study <- subset_meta$study
cells <- subset_meta$cluster_group
value <- abs(rnorm(length(study) , 0 , 15))
data <- data.frame(study, cells, value)


ggplot(data, aes(fill=cells, y=value, x=study)) + 
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values = c("B_cells" = "#B71C1C","DCs" = "#EF6C00", "Endothelial" = "#F9A825",
                               "Enteric_glia" = "#7CB342","Epithelial"="#689F38",
                               "Erythroblast" = "#4CAF50", "Fibroblasts" = "#81C784", "Granulocyte" = "#4DB6AC", 
                               "Hepatocyte" = "#26A69A", "ILC" = "#29B6F6", "Macrophages" = "#039BE5", 
                               "Mast" = "#3F51B5", "Megakaryocyte" = "#7E57C2", "Melanocytes" = "#AB47BC",
                               "Mesothelial" = "#8E24AA", "Monocytes" = "#D81B60", "Myeloid" = "#F06292", "ND" = "#F48FB1",
                               "NK_cells" = "#FFD54F", "NKT" = "#FFB74D", "osteoclast" = "#FF8A65", "Patient_specific" = "#A1887F",
                               "Plasma_cells" = "#FFD700", "Platelets" = "#B0BEC5", "PTPRC- population" = "#90A4AE", "RBC" = "#78909C",
                               "T_cell" = "#607D8B", "Tumor" = "#455A64", "Undetermined" = "#37474F", "Doublets" = "#263238"))
