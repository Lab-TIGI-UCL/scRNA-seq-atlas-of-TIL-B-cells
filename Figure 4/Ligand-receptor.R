########### R script for running Ligand-receptor analysis ###########
########### Created by Danwen Qian on 09-05-2024. ###########
########### Last modified by Danwen Qian on 09-05-2024. ###########

library(tidyverse)
library(magrittr)
library(liana)
library(CellChat)
library(Seurat)
library(RColorBrewer)

BCELL<-readRDS("~/BCELL.rds")
TCELL<-readRDS("~/combined_TCELL_integrated_16dataset.rds")

BCELL_proportion<-as.data.frame.array(prop.table(table(BCELL$patient_ID,BCELL$clusters),margin=1))
TCELL_proportion<-as.data.frame.array(prop.table(table(TCELL$patient_ID,TCELL$clusters),margin=1))
BCELL_proportion$patient_ID<-rownames(BCELL_proportion)
TCELL_proportion$patient_ID<-rownames(TCELL_proportion)
proportion<-merge(BCELL_proportion,TCELL_proportion,by="patient_ID")
write.csv(proportion,"~/proportion.csv")

p<-ggplot(proportion, aes_string("CD4_Treg","Naive.B.cells")) +geom_point(aes())+theme_classic()+geom_smooth(method = 'lm', se = T)+ggpubr::stat_cor()+ylab("Naive B cells")+coord_cartesian(ylim =c(0,1))
ggsave("Naive.B.cells.pdf",p,height=4, width=4,dpi=600)

TCELL_patient<-TCELL@meta.data[,c(6,14)] %>% distinct(patient_ID, .keep_all = TRUE)
BCELL_patient<-BCELL@meta.data[,23:22] %>% distinct(patient_ID, .keep_all = TRUE)
patient<-merge(TCELL_patient,BCELL_patient,by="patient_ID")
TCELL<-subset(TCELL, sample_ID %in% patient$patient_ID)
TCELL<-subset(TCELL,clusters %in% c("CD8_Tem","CD4/CD8_Tnaive","CD4_Treg","CD8_Teff","CD4/CD8_terminal_Tex","CD8_GZMK_Tex","CD4 Tcm","CD4_Th","CD4/CD8_Tisg","NK_like_T" ))
BCELL<-subset(BCELL, patient_ID %in% patient$patient_ID)
DefaultAssay(BCELL) <- "RNA"

data<-merge(TCELL,BCELL)
saveRDS(data,"~/TCELL_BCELL_interation.rds")



liana_results <- liana_wrap(data,method = c("natmi", "connectome", "logfc", "sca", "cellphonedb","call_cellchat"),resource = c("Consensus"))

saveRDS(liana_results,"~/TCELL_BCELL_interation_results.rds")

liana_results <- liana_results %>%
  liana_aggregate()

write.csv(liana_results,"~/TCELL_BCELL_liana_results.csv")

pdf("~/chord_Bcells.pdf", width = 6, height = 6)
liana_results %>%
     filter(aggregate_rank <= 0.05) %>%
     chord_freq(source_groups = c("CD8_Tem","CD4/CD8_Tnaive","CD4_Treg","CD8_Teff","CD4/CD8_terminal_Tex","CD8_GZMK_Tex","CD4 Tcm","CD4_Th","CD4/CD8_Tisg","NK_like_T"),
                target_groups = c("Naive B cells", "Activated B cells", "Resting Memory B cells", "IGKC-high Plasma/Plasmablasts", "Conventional Plasma cells", "Stressed Plasma cells", "Atypical Memory B cells", "Proliferative B cells", "GC B cells", "MT1X-high Plasma/Plasmablasts"))
dev.off()

liana_dotplot1<-function (liana_res, source_groups = NULL, target_groups = NULL,
    ntop = NULL, specificity = "natmi.edge_specificity", magnitude = "sca.LRscore",
    y.label = "Interactions (Ligand -> Receptor)", size.label = "Interaction\nSpecificity",
    colour.label = "Expression\nMagnitude", show_complex = TRUE,
    size_range = c(2, 10))
{
    if (show_complex) {
        entities <- c("ligand.complex", "receptor.complex")
    }
    else {
        entities <- c("ligand", "receptor")
    }
    if (!is.null(source_groups)) {
        liana_mod <-liana_res %>% dplyr::filter(source %in% source_groups)
    } else {
        liana_mod <-liana_res
    }
    if (!is.null(target_groups)) {
        liana_mod <- liana_mod %>% dplyr::filter(target %in% target_groups)
    } else {
        liana_mod <- liana_mod
    }
    
    if (!is.null(ntop)) {
        top_int <- liana_mod %>% distinct_at(entities) %>% head(ntop)
        liana_mod %<>% inner_join(top_int, by = entities)
    }
    liana_mod %<>% rename(magnitude = !!magnitude) %>% rename(specificity = !!specificity) %>%
        unite(entities, col = "interaction", sep = " -> ") %>%
        unite(c("source", "target"), col = "source_target", remove = FALSE)
    cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
        "#0072B2", "#D55E00", "#CC79A7", "#DF69A7")
    suppressWarnings(ggplot(liana_mod, aes(x = source, y = interaction,
        colour = magnitude, size = specificity, group = source)) +
        geom_point() + scale_color_gradientn(colours = viridis::viridis(20)) +
        scale_size_continuous(range = size_range) + facet_grid(. ~
        target, space = "free", scales = "free", switch = "y") +
        labs(y = y.label, colour = colour.label, size = size.label,
            x = "Source", title = "Target") + theme_bw(base_size = 20) +
        theme(legend.text = element_text(size = 16), axis.text.x = element_text(colour = cbPalette[1:length(unique(liana_mod$target))],
            face = "bold", size = 23), axis.title.x = element_text(colour = "gray6"),
            axis.text.y = element_text(size = 18, vjust = 0.5),
            legend.title = element_text(size = 18), panel.spacing = unit(0.1,
                "lines"), strip.background = element_rect(fill = NA),
            plot.title = element_text(vjust = 0, hjust = 0.5,
                colour = "gray6"), strip.text = element_text(size = 24,
                colour = "gray6")))
}

for(i in c("Naive B cells", "B cells", "Resting Memory B cells", "IGKC-high Plasma/Plasmablasts", "Conventional Plasma cells", "Stressed Plasma cells", "Atypical Memory B cells", "Proliferative B cells", "GC B cells", "MT1X-high Plasma/Plasmablasts")){
    p<-liana_results %>%
        filter(aggregate_rank <= 0.05) %>%
        liana_dotplot1(source_groups =  c("CD8_Tem","CD4/CD8_Tnaive","CD4_Treg","CD8_Teff","CD4/CD8_terminal_Tex","CD8_GZMK_Tex","CD4 Tcm","CD4_Th","CD4/CD8_Tisg","NK_like_T"),
                      target_groups = c(i),
                      show_complex = F)+ theme(axis.text.x = element_text(size=16,hjust = 1,vjust = 0.5,angle = 90,color = "black"))+scale_colour_gradientn(colours = colorRampPalette(brewer.pal(9, "OrRd"))(100))
    ggsave(paste0("~//",i,"_Tcells.pdf"), p, height=10, width=11, dpi=600)
}



