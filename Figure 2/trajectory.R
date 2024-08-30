########### R script for trajectory analysis. ###########
########### Created by Danwen Qian on 29-04-2024. ###########
########### Last modified by Danwen Qian on 29-04-2024. ###########
###monocle3

library(monocle3)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
library(ggsci)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)

options(Seurat.object.assay.version = "v3")



combined_integrated<-readRDS("BCELL.rds")
Idents(object = combined_integrated) <- "integrated_snn_res.0.2"
combined_integrated <- RenameIdents(object = combined_integrated, '3' = "Naive B cells", '1' = "Activated B cells", '0' = "Resting Memory B cells", '6' = "IGKC-high Plasma/Plasmablasts", '2' = "Conventional Plasma cells", '4' = "Stressed Plasma cells", '5' = "Atypical Memory B cells", '7' = "Proliferative B cells", '8' = "GC B cells", '9' = "MT1X-high Plasma/Plasmablasts")
combined_integrated$clusters<-combined_integrated@active.ident

DefaultAssay(combined_integrated) <- "RNA"
combined_integrated <- NormalizeData(combined_integrated, normalization.method = "LogNormalize", scale.factor = 10000)

cds <- as.cell_data_set(combined_integrated)
cds <- cluster_cells(cds)
rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name

cds <- learn_graph(cds)

custom_color = c("#E66C37","Resting Memory B cells" = "#D64550", "Activated B cells" = "#FF977E","Atypical Memory B cells" = "#A43B76","Proliferative B cells" = "#9A64A0","GC B cells"="#750985", "Conventional Plasma cells" = "#5ECBC8", "Stressed Plasma cells" = "#28788D", "IGKC-high Plasma/Plasmablasts" = "#37A794",  "MT1X-high Plasma/Plasmablasts" = "#4C5D8A")
p<-plot_cells(cds,
color_cells_by = "clusters",
label_groups_by_cluster=FALSE,
label_cell_groups = TRUE,
label_leaves=FALSE,
label_branch_points=FALSE,
group_label_size = 4)+scale_color_manual(values = custom_color)
ggsave("trajectory.pdf", p, height=5, width=5, dpi=600)

cds <- order_cells(cds)
colData(cds)$pseudotime <- pseudotime(cds)
combined_integrated$pseudotime <- pseudotime(cds)

p<-plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=3)
ggsave("pseudotime.pdf", p, height=5, width=6, dpi=600)

cds <- estimate_size_factors(cds)

mature_genes<-c("CD19", "MS4A1", "JCHAIN","IGHG1","IGHA1","IGHM")

mature_cds <- cds[rowData(cds)$gene_short_name %in% mature_genes,]

mature_cds <- order_cells(mature_cds)
p<-plot_genes_in_pseudotime(mature_cds,
                         color_cells_by="clusters",
                         min_expr=0.5,ncol=3)+scale_color_manual(values = custom_color)
ggsave("genes_in_pseudotime.pdf", p, height=5, width=10, dpi=600)

#branches analysis
cds_sub_Aty <- choose_graph_segments(cds, clear_cds = F, return_list = T)
colData(cds)$chosen<-"Unchosen"
colData(cds)$chosen[colnames(cds) %in% cds_sub_Aty$cells]<-"Chosen"
p<-plot_cells(cds,
              color_cells_by = "chosen",
              label_groups_by_cluster=FALSE,
              label_cell_groups = F,
              label_leaves=FALSE,
              label_branch_points=FALSE,
              group_label_size = 5,trajectory_graph_segment_size= 0.5)+scale_color_manual(values=c("purple","grey"))
ggsave("branch_plasma.pdf", p, height=5, width=6, dpi=600)

cds_sub_Aty <- cds[,cds_sub_Aty$cells]
subset_pr_test_res <- graph_test(cds_sub_Aty,neighbor_graph="principal_graph", cores=4)
subset_pr_test_res<- subset(subset_pr_test_res, q_value < 0.05)
write.csv(subset_pr_test_res,"Atp_plasma.csv")

high_affinity=c("BATF","GARS1","GART","MIF","MYC","SPP1","UCK2",
                "CD320","TNFRSF8","AURKA","BUB1","CCNA2","CCNB1","CCNB2","CCND2","CDC20","CDC25C","KIF22","PLK1","NFIL3","PML")
Low_affinity=c("PPP1R15A","RGS1","CCR6","CD22","CD38","CD72","FCER2","ICOSLG","TNFRSF18","BACH2","EGR3","ELK4","FOXP1","JUN","NR4A1","REL")
exhaustion_genes = c('PDCD1','CD160','FASLG','CD244','LAG3','TNFRSF1B','CCR5','CCL3',
                     'CCL4','CXCL10','CTLA4','LGALS1','LGALS3','PTPN13','RGS16','ISG20',
                     'MX1','IRF4','EOMES','PBX3','NFATC1','NR4A2','CKS2','GAS2',
                     'ENTPD1','CA2',"CD52", "APOE",  "PTGDS", "PIM2", "DERL3")


Bactivated_genes=c("CD69","CD83","IER2","DUSP2","IL6","NR4A2","JUN","CCR7","GPR183")
BCSR_genes=c("APEX1","APEX2","XRCC5","XRCC6","POLD2","AICDA")
CSR_m=c("APEX1","XRCC5","XRCC6","POLD2","POLE3", #CSR machinery
        "NCL","NME2","DDX21",#IgH locus
        "NPM1","SERBP1",#CSR interactors
        "MIR155HG","HSP90AB1",#AICDA/AICDA stability
        "BATF","HIVEP3","BHLHE40","IRF4")#TF

feature <- list(high_affinity=high_affinity,
                Low_affinity=Low_affinity,
                exhaustion_genes=exhaustion_genes,
                Bactivated_genes=Bactivated_genes,
                BCSR_genes=BCSR_genes,
                CSR_m=CSR_m)

cds_subset_Aty <- subset(combined_integrated,cells=colnames(cds_sub_Aty))

cds_subset_Aty <- AddModuleScore(cds_subset_Aty,features = feature,name=c("high_affinity","Low_affinity","exhaustion_genes","Bactivated_genes","BCSR_genes","CSR_m"))


cds_subset_Aty$path<-"Aty"

cds_subset<-merge(cds_subset_GC,y=list(cds_subset_Aty,cds_subset_plamsa))

saveRDS(cds_subset,"BCELL_branches.rds")

cds_subset$pseudotime <- (cds_subset$pseudotime-min(cds_subset$pseudotime))/(max(cds_subset$pseudotime)-min(cds_subset$pseudotime))


for(i in 33:38){
  p <- ggplot(cds_subset@meta.data) +
  stat_smooth(aes_string(x = "pseudotime", y = colnames(cds_subset@meta.data)[i],color="path"),method = "lm", formula = y ~ poly(x, 3),se=TRUE)+
  theme_classic()
  ggsave(paste0(colnames(cds_subset@meta.data)[i],".pdf"), p, height=4, width=10, dpi=600)
}

Idents(cds_subset) <- cds_subset$path
markers <- FindAllMarkers(cds_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top20<-markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)

write.csv(top20,file="top20_branches.csv")

p<-plot_cells(cds, genes=c("MARCKSL1","MS4A1","CD22","MX1","ISG15","CAPG","IGHA1","IGHG1","JCHAIN"),
show_trajectory_graph=FALSE,
label_cell_groups=FALSE,
label_leaves=FALSE)
ggsave("branches_gene.pdf", p, height=10, width=12, dpi=600)



other_genes<-c("S100A9","C1QA","CD1A","VCAN","ISG15","CXCL9","AXL","SIGLEC6")
other_cds <- cds_sub[rowData(cds_sub)$gene_short_name %in% other_genes,]
other_cds <- order_cells(other_cds)
p<-plot_genes_in_pseudotime(other_cds,
                         color_cells_by="clusters_0.4",
                         min_expr=0.5,ncol=4)+scale_color_d3("category20c")
ggsave("monocle3/other_genes_in_pseudotime.pdf", p, height=4, width=10, dpi=600)

##Analyzing branches in single-cell trajectories,run for direct and indirect
#cds_subset <- choose_graph_segments(cds_sub,clear_cds=FALSE)
cds_subset <- choose_graph_segments(cds_sub, clear_cds = F, return_list = T)
colData(cds_sub)$chosen<-"Unchosen"
colData(cds_sub)$chosen[colnames(cds_sub) %in% cds_subset$cells]<-"Chosen"
p<-plot_cells(cds_sub,
              color_cells_by = "chosen",
              label_groups_by_cluster=FALSE,
              label_cell_groups = T,
              label_leaves=FALSE,
              label_branch_points=FALSE,
              group_label_size = 5,trajectory_graph_segment_size= 0.5)+scale_color_manual(values=c("purple","grey"))
ggsave("monocle3/branch_directly.pdf", p, height=5, width=5, dpi=600)


cds_subset <- cds_sub[,cds_subset$cells]

subset_pr_test_res <- graph_test(cds_subset,neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))
cds_subset <- preprocess_cds(cds_subset)
gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=0.01,random_seed=1)
write.csv(gene_module_df,"monocle3/gene_module_directly.csv")

cell_group_df <- tibble::tibble(cell=row.names(colData(cds_subset)),
                                cell_group=colData(cds_subset)$clusters_0.4)
agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
write.csv(agg_mat,"monocle3/modules_directly.csv")

agg_mat <-agg_mat[,c(1:2,5:6)]
p<-pheatmap::pheatmap(agg_mat,
                      scale="column", clustering_method="ward.D2")
ggsave("monocle3/module_directly.pdf", p, height=10, width=3, dpi=600)

#annotation of module
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
m_t2g <- msigdbr(species = "Homo sapiens") %>% filter(gs_cat %in% c("C5","H"))
m_t2g<-m_t2g %>% dplyr::select(gs_name, gene_symbol)
results<-list()

for (i in 1:63){
    genelist<-gene_module_df[gene_module_df$module==i,]$id
    em<-enricher(genelist, TERM2GENE=m_t2g)
    pathway<-em@result %>% filter(p.adjust<0.05)
    results[[paste0("Module_",i)]]<-pathway$Description
}

n.obs <- sapply(results, length)
seq.max <- seq_len(max(n.obs))
results <- sapply(results, "[", i = seq.max)
write.csv(results,"monocle3/module_directly_annotation.csv",row.names = F)

CCR7_dir<-(gene_module_df %>% subset(module %in% c("21","52","17","54","20","2")))$id
em<-enricher(CCR7_dir, TERM2GENE=m_t2g)
p<-barplot(em,showCategory = 20)
ggsave("monocle3/CCR7_dir_pathway.pdf", p, height=8, width=7, dpi=600)
pathway<-em@result %>% filter(p.adjust<0.05)
write.csv(pathway,"monocle3/CCR7_dir_pathway.csv",row.names = F)
ego2<-pairwise_termsim(em)
p<-emapplot(ego2, repel=TRUE,showCategory=20,cex.params = list(category_label =0.6,line = 0.5))
ggsave("monocle3/CCR7_dir_pathway_emapplot.pdf", p, height=8, width=8, dpi=600)

dys_dir<-(gene_module_df %>% subset(module %in% c("46","49","34","18","6","30")))$id
em<-enricher(dys_dir, TERM2GENE=m_t2g)
p<-barplot(em,showCategory = 20)
ggsave("monocle3/dys_dir_pathway.pdf", p, height=8, width=7, dpi=600)
pathway<-em@result %>% filter(p.adjust<0.05)
write.csv(pathway,"monocle3/CCR7_dir_pathway.csv",row.names = F)
ego2<-pairwise_termsim(em)
p<-emapplot(ego2, repel=TRUE,showCategory=20,cex.params = list(category_label =0.6,line = 0.5))
ggsave("monocle3/CCR7_dir_pathway_emapplot.pdf", p, height=8, width=8, dpi=600)

###indirectly way
cds_subset <- choose_graph_segments(cds_sub, clear_cds = F, return_list = T)
colData(cds_sub)$chosen<-"Unchosen"
colData(cds_sub)$chosen[colnames(cds_sub) %in% cds_subset$cells]<-"Chosen"
p<-plot_cells(cds_sub,
              color_cells_by = "chosen",
              label_groups_by_cluster=FALSE,
              label_cell_groups = T,
              label_leaves=FALSE,
              label_branch_points=FALSE,
              group_label_size = 5,trajectory_graph_segment_size= 0.5)+scale_color_manual(values=c("purple","grey"))
ggsave("monocle3/branch_indirectly.pdf", p, height=5, width=5, dpi=600)

cds_subset <- cds_sub[,cds_subset$cells]

subset_pr_test_res <- graph_test(cds_subset,neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))
cds_subset <- preprocess_cds(cds_subset)
gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=0.01,random_seed=1)
write.csv(gene_module_df,"monocle3/gene_module_indirectly.csv")

cell_group_df <- tibble::tibble(cell=row.names(colData(cds_subset)),
                                cell_group=colData(cds_subset)$clusters_0.4)
agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
write.csv(agg_mat,"monocle3/modules_indirectly.csv")

agg_mat <-agg_mat[,c(1:2,4:5,7)]
p<-pheatmap::pheatmap(agg_mat,
                      scale="column", clustering_method="ward.D2")
ggsave("monocle3/module_indirectly.pdf", p, height=10, width=4, dpi=600)

#annotation of module
results<-list()

for (i in 1:50){
    genelist<-gene_module_df[gene_module_df$module==i,]$id
    em<-enricher(genelist, TERM2GENE=m_t2g)
    pathway<-em@result %>% filter(p.adjust<0.05)
    results[[paste0("Module_",i)]]<-pathway$Description
}

n.obs <- sapply(results, length)
seq.max <- seq_len(max(n.obs))
results <- sapply(results, "[", i = seq.max)
write.csv(results,"monocle3/module_indirectly_annotation.csv",row.names = F)

CCR7_indir<-(gene_module_df %>% subset(module %in% c("43","3","36","32","5","23","29")))$id
em<-enricher(CCR7_indir, TERM2GENE=m_t2g)
p<-barplot(em,showCategory = 20)
ggsave("monocle3/CCR7_indir_pathway.pdf", p, height=8, width=7, dpi=600)
pathway<-em@result %>% filter(p.adjust<0.05)
write.csv(pathway,"monocle3/CCR7_indir_pathway.csv",row.names = F)
ego2<-pairwise_termsim(em)
p<-emapplot(ego2, repel=TRUE,showCategory=20,cex.params = list(category_label =0.6,line = 0.5))
ggsave("monocle3/CCR7_dir_pathway_emapplot.pdf", p, height=8, width=8, dpi=600)

dys_indir<-(gene_module_df %>% subset(module %in% c("22","24","4","20")))$id
em<-enricher(dys_indir, TERM2GENE=m_t2g)
p<-barplot(em,showCategory = 20)
ggsave("monocle3/dys_indir_pathway.pdf", p, height=8, width=7, dpi=600)
pathway<-em@result %>% filter(p.adjust<0.05)
write.csv(pathway,"monocle3/CCR7_dir_pathway.csv",row.names = F)
ego2<-pairwise_termsim(em)
p<-emapplot(ego2, repel=TRUE,showCategory=20,cex.params = list(category_label =0.6,line = 0.5))
ggsave("monocle3/CCR7_dir_pathway_emapplot.pdf", p, height=8, width=8, dpi=600)

##DF
Branch1<-choose_graph_segments(cds_sub, clear_cds = F, return_list = T)
Branch2<-choose_graph_segments(cds_sub, clear_cds = F, return_list = T)

overlap<-Branch1$cells[Branch1$cells %in% Branch2$cells]
Branch1<-Branch1$cells[!(Branch1$cells %in% overlap)]
Branch2<-Branch2$cells[!(Branch2$cells %in% overlap)]

DC_branch<-subset(combined_DC_integrated, cells = c(Branch1,Branch2))
DC_branch$Branch<-"Direct_branch"
DC_branch$Branch[colnames(DC_branch)  %in% Branch2]<-"Indirect_branch"

DC_branch_dys<-subset(DC_branch,clusters_0.4=="cDC2_CD1C+_B")
Idents(object=DC_branch_dys)<-"Branch"
de.markers <- FindMarkers(DC_branch_dys, ident.1 = "Direct_branch", ident.2 = "Indirect_branch",logfc.threshold=0)
write.csv(de.markers,file="monocle3/Branch_dys_DE.csv")
de.markers$Group<-"not-significant"
de.markers$Group[which((de.markers$p_val_adj < 0.05) & (abs(de.markers$avg_log2FC) > 1))]="significant"
de.markers$logP<--log10(de.markers$p_val_adj)
de.markers<-de.markers[de.markers$logP!="Inf",]
de.markers$Lable<-""
up.gene<-rownames(de.markers[((de.markers$p_val_adj < 0.05) & de.markers$avg_log2FC > 1),] %>% slice_max(n = 5, order_by = avg_log2FC))
down.gene<-rownames(de.markers[((de.markers$p_val_adj < 0.05) & de.markers$avg_log2FC < -1),] %>% slice_min(n = 5, order_by = avg_log2FC))
top10<-c(as.character(up.gene),as.character(down.gene))
de.markers$Lable[match(top10,rownames(de.markers))]<-top10
p<-ggscatter(de.markers,x="avg_log2FC", y="logP", color="Group",palette=c("#BBBBBB","#CC0000"),size=2, label=de.markers$Lable, font.label=8,repel=T,
xlab="log2FoldChange",
ylab="-log10 (P_value_adjust)")+
geom_hline(yintercept=1.30, linetype="dashed")+
geom_vline(xintercept=c(-1,1), linetype="dashed")+
theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12))
ggsave("monocle3/Branch_dys_volcano.pdf", p, height=5, width=5, dpi=600)

DC_branch_CCR7<-subset(DC_branch,clusters_0.4=="CCR7+_DC")
Idents(object=DC_branch_CCR7)<-"Branch"
de.markers <- FindMarkers(DC_branch_CCR7, ident.1 = "Direct_branch", ident.2 = "Indirect_branch",logfc.threshold=0)
write.csv(de.markers,file="monocle3/Branch_CCR7_DE.csv")
de.markers$Group<-"not-significant"
de.markers$Group[which((de.markers$p_val_adj < 0.05) & (abs(de.markers$avg_log2FC) > 1))]="significant"
de.markers$logP<--log10(de.markers$p_val_adj)
de.markers<-de.markers[de.markers$logP!="Inf",]
de.markers$Lable<-""
up.gene<-rownames(de.markers[((de.markers$p_val_adj < 0.05) & de.markers$avg_log2FC > 1),] %>% slice_max(n = 5, order_by = avg_log2FC))
down.gene<-rownames(de.markers[((de.markers$p_val_adj < 0.05) & de.markers$avg_log2FC < -1),] %>% slice_min(n = 5, order_by = avg_log2FC))
top10<-c(as.character(up.gene),as.character(down.gene))
de.markers$Lable[match(top10,rownames(de.markers))]<-top10
p<-ggscatter(de.markers,x="avg_log2FC", y="logP", color="Group",palette=c("#BBBBBB","#CC0000"),size=2, label=de.markers$Lable, font.label=8,repel=T,
xlab="log2FoldChange",
ylab="-log10 (P_value_adjust)")+
geom_hline(yintercept=1.30, linetype="dashed")+
geom_vline(xintercept=c(-1,1), linetype="dashed")+
theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12))
ggsave("monocle3/Branch_CCR7_volcano.pdf", p, height=6, width=6, dpi=600)

###non-skin
non_skin<-subset(combined_DC_integrated,cancer_type %in% c("BC","CRC","LC","Ova","RCC")&clusters_0.4 %in% c("cDC2_CD1C+_A","cDC2_CD1C+_B","cDC2_C1Q+","cDC2_HSP+","CCR7+_DC","cDC2_CD207+","cDC2_ISG15+","cDC2_CXCL9+","AS_DC","cDC2_S100B+"))
cds_sub <- as.cell_data_set(non_skin)
cds_sub <- cluster_cells(cds_sub,resolution =2*10^-4,random_seed = 3)
cds_sub <- learn_graph(cds_sub)
p<-plot_cells(cds_sub,
              color_cells_by = "clusters_0.4",
              label_groups_by_cluster=FALSE,
              label_cell_groups = TRUE,
              label_leaves=FALSE,
              label_branch_points=FALSE,
              group_label_size = 5)
ggsave("monocle3/non_skin_tractrpy.png", p, height=12, width=12, dpi=600)
cds_sub <- order_cells(cds_sub)
p<-plot_cells(cds_sub,
              color_cells_by = "pseudotime",
              label_cell_groups=FALSE,
              label_leaves=FALSE,
              label_branch_points=FALSE,
              graph_label_size=3,
              label_principal_points = FALSE)
ggsave("monocle3/non_skin_pseudotime.png", p, height=9, width=12, dpi=600)
cds_subset <- choose_graph_segments(cds_sub,clear_cds=FALSE)









