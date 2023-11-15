#CPI Data Processing


# Loading libraries 
library(tidyverse) 
library(reshape2)
library(ggplot2)
library(ggpubr)



#Load DEGs for all clusters
bcell_signatures <- read.csv("~/tumoursubset__0.2_RNA_rPCA.csv", header = TRUE, check.names = FALSE) %>% as.data.frame()
bcell_signatures$cluster = paste0("cluster_",bcell_signatures$cluster) 
bcell_signatures <- bcell_signatures %>% as.data.frame()

#Load CPI expression matrix
library(data.table)
exp_matrix3 <- fread("~excl_outliers_lm_study_tissueType_resid.csv", header = TRUE) %>% as.data.frame()

#Rework data
rownames(exp_matrix3) <- exp_matrix3$Patient 
exp_matrix3 <- exp_matrix3 %>% select(-Patient)
exp_matrix3 <- exp_matrix3 %>% t()%>%as.data.frame()
exp_matrix3$Gene <- rownames(exp_matrix3)

#Load survival data
survival_dat <- read.csv("~/master_sample_sheet6.csv")
survival_dat <- survival_dat %>% select(Patient, Response_responder, Tum_Type)

##########
# Step 1: Large dataframe of patients with expressed means per cluster 
######### 

bcell_genes <- unique(bcell_signatures$gene) %>% as.list()  
exp_matrix3 <- exp_matrix3 %>% filter(Gene %in% bcell_genes) %>% as.data.frame()  
patients <- exp_matrix3 %>% select(-Gene)
patients <- colnames(patients) %>% unique() 
clusters <- unique(bcell_signatures$cluster) 


# formation of dataframe 
exp_means_df <- matrix(nrow = 938, ncol = 10) %>% as.data.frame()
rownames(exp_means_df) <- patients
colnames(exp_means_df) <- clusters
exp_means_df <- exp_means_df %>% as.data.frame()


# For loop to fill dataframe
for (i in patients) {
  
  tmp <- exp_means_df[i, ] %>% as.data.frame()
  full_mat <- exp_matrix3 %>% select(Gene, i) %>% as.data.frame()
  
  print(paste0(i, "_started"))
  
  for (j in clusters) {
    
    j <- j %>% as.character()
    
    new_tmp <- bcell_signatures %>% as.data.frame()
    new_tmp <- new_tmp %>% select(cluster, gene) %>% as.data.frame()
    
    gene_list <- new_tmp %>% filter(cluster == j)
    gene_list <- gene_list$gene
    
    tmp_mat <- full_mat %>% dplyr::filter(Gene %in% gene_list)
    tmp_mat <- tmp_mat %>% select(-Gene) %>% as.data.frame()
    colnames(tmp_mat) <- "temporary"
    exp_means <- mean(tmp_mat$temporary)
    
    exp_means_df[i, j] <- exp_means
    
    print(paste0(j, "_done"))
  }
  
  print(paste0(i, "_ended"))
}

exp_means_df <- exp_means_df %>% as.data.frame()

exp_means_df$survival <- NA
exp_means_df$cancer_type <- NA

for (i in patients) {
  if (i %in% survival_dat$Patient) {
    tmp <- survival_dat %>% filter(Patient == i)
    tt <- tmp$Tum_Type %>% as.character()
    s <- tmp$Response_responder %>% as.character()
    
    exp_means_df[i, "survival"] <- s
    exp_means_df[i, "cancer_type"] <- tt
  }
}

exp_means_df <- exp_means_df %>% as.data.frame()
saved_means_df <- exp_means_df

#To remove blanks
exp_means_df <- exp_means_df[!(exp_means_df$survival == ""), ]

