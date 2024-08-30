#Pathway radar plot

library(data.table)
library(fgsea)
library(ggplot2)
library(msigdbr)
library(dplyr)
library(tidyverse)
library(fmsb)

m_df <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:KEGG")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

markers <- read.csv("~/tumoursubset__0.2_RNA_rPCA.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()


markers$cluster[markers$cluster=="0"]<-"Resting Memory B cells"
markers$cluster[markers$cluster=="1"]<-"Activated Memory B cells"
markers$cluster[markers$cluster=="2"]<-"Conventional Plasma cells"
markers$cluster[markers$cluster=="3"]<-"Naive B cells"
markers$cluster[markers$cluster=="4"]<-"Stressed Plasma cells"
markers$cluster[markers$cluster=="5"]<-"Atypical Memory B cells"
markers$cluster[markers$cluster=="6"]<-"IGKC-high Plasma/Plasmablasts"
markers$cluster[markers$cluster=="7"]<-"Transitional/Proliferative 1"
markers$cluster[markers$cluster=="8"]<-"Transitional/Proliferative 2"
markers$cluster[markers$cluster=="9"]<-"MT1X-high Plasma/Plasmablasts"

colnames(markers) [1] <- "gene1" 

fgseaRes<-data.frame()
for(i in c("Naive B cells", "Activated Memory B cells", "Resting Memory B cells", "IGKC-high Plasma/Plasmablasts", "Conventional Plasma cells", "Stressed Plasma cells", "Atypical Memory B cells", "Transitional/Proliferative 1", "Transitional/Proliferative 2", "MT1X-high Plasma/Plasmablasts")){
  cluster.genes<- markers %>% dplyr::filter(cluster == i) %>%arrange(desc(avg_log2FC)) %>% dplyr::select(gene,avg_log2FC) 
  ranks<- deframe(cluster.genes)
  fgseaRes1<- fgsea(fgsea_sets, stats = ranks, nperm = 1000) 
  fgseaRes1$cluster<-i
  fgseaRes<-bind_rows(fgseaRes, fgseaRes1)
}



top5<- fgseaRes %>%
  filter(pval < 0.05)%>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = NES)



pathway<-as.data.frame(unique(top5$pathway))

colnames(pathway)<-"pathway"

for(i in c("Naive B cells", "Activated Memory B cells", "Resting Memory B cells", "IGKC-high Plasma/Plasmablasts", "Conventional Plasma cells", "Stressed Plasma cells", "Atypical Memory B cells", "Transitional/Proliferative 1", "Transitional/Proliferative 2", "MT1X-high Plasma/Plasmablasts")){
  cluster.pathway<- fgseaRes %>% dplyr::filter(cluster == i)%>%dplyr::select(pathway,NES)
  pathway<-left_join(pathway, cluster.pathway,by ="pathway")
}


colnames(pathway)[-1]<-c("Naive B cells", "Activated Memory B cells", "Resting Memory B cells", "IGKC-high Plasma/Plasmablasts", "Conventional Plasma cells", "Stressed Plasma cells", "Atypical Memory B cells", "Transitional/Proliferative 1", "Transitional/Proliferative 2", "MT1X-high Plasma/Plasmablasts")


rownames(pathway)<-pathway$pathway

pathway[is.na(pathway)] <- 0

pathway <- data.frame(t(pathway[-1]))

#Add rows with highest and lowest values
pathway <- rbind(rep(2.2374599,25) , rep(-2.9921065,25) , pathway)

#Rename pathways
pathway <- pathway %>% 
  rename(
    TCR_Signalling = KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY,
    NK_mediated_cytotoxicity = KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY,
    VEGF_Signalling= KEGG_VEGF_SIGNALING_PATHWAY,
    IgA_production= KEGG_INTESTINAL_IMMUNE_NETWORK_FOR_IGA_PRODUCTION,
    GVHD = KEGG_GRAFT_VERSUS_HOST_DISEASE,
    Cell_adhesion= KEGG_CELL_ADHESION_MOLECULES_CAMS,
    T1D = KEGG_TYPE_I_DIABETES_MELLITUS,
    Asthma = KEGG_ASTHMA,
    Protein_export = KEGG_PROTEIN_EXPORT,
    Vibrio_Cholerae_Infection = KEGG_VIBRIO_CHOLERAE_INFECTION,
    Antigen_presentation = KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION,
    N_Glycan_Biosynthesis = KEGG_N_GLYCAN_BIOSYNTHESIS,
    Cardiac_muscle_contraction = KEGG_CARDIAC_MUSCLE_CONTRACTION,
    Insulin_signalling =KEGG_INSULIN_SIGNALING_PATHWAY,
    MTOR_signalling =KEGG_MTOR_SIGNALING_PATHWAY,
    Ribosome = KEGG_RIBOSOME,
    MAPK_signalling = KEGG_MAPK_SIGNALING_PATHWAY,
    Endocytosis = KEGG_ENDOCYTOSIS,
    Prion_Diseases = KEGG_PRION_DISEASES,
    Spliceosome = KEGG_SPLICEOSOME,
    Lupus =KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS,
    Oxidative_Phosphorylation = KEGG_OXIDATIVE_PHOSPHORYLATION,
    Parkinsons = KEGG_PARKINSONS_DISEASE,
    Viral_Mycarditis = KEGG_VIRAL_MYOCARDITIS,
    Leishmania = KEGG_LEISHMANIA_INFECTION,
    Allograft_rejection =KEGG_ALLOGRAFT_REJECTION,
    Autoimmune_thyroid = KEGG_AUTOIMMUNE_THYROID_DISEASE
  )

# Select relevant columns
pathway <- pathway[,c(1,2,3,4,7,10,12,13,17,18,19,22)] 

#Plot Radar Plot
radarchart(pathway, plwd = 5, pcol = colors_border) + 
  legend(x=1.5, y=1, legend = rownames(pathway[-c(1,2),]), bty = "n", pch=20, text.col = "black", col=colors_border, cex=1, pt.cex=3)

colors_border= c("#9E9E9E","#D64550",
                 "#FF977E","#A43B76",
                 "#9A64A0","#750985",
                 "#5ECBC8", "#28788D", "#37A794", 
                 "#4C5D8A", "#039BE5", "#BF360C")

