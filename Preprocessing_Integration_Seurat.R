#refinedBcellmeta subsetted to exclude cells in clusters identified as artifacts, technical, T/B doublets and MNP/B doublets


### Load datasets, select B cells, preprocess and integrate them 

library(Seurat)
library(SeuratDisk)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)


# Initialization

meta <- read.csv("/SAN/colcc/tigilab-general/evie/Bcellatlas_files_subsetted_meta.csv")
setwd ("/SAN/colcc/NMD_inputs/scRNA_datasets")
source ("/SAN/colcc/NMD_analysis/github_scripts/NextGenTargets_WorkStream2/scRNA_tools/tools_v1.2.R")

studies <- c(
  "Azizi",
  "BRAUN",
  "Chan",
  "KRISHNA",
  "KIM",
  "YOST",
  "Li",
  "Maynard",
  "Bi",
  "Pelka",
  "QIAN",
  "VISHWAKARMA",
  "Wu",
  "Zilionis",
  "ZhangYY")

files <- c(
  "/SAN/colcc/NMD_inputs/scRNA_datasets/AZIZI_CELL_2018/download/Azizi_human.h5Seurat", 
  "/SAN/colcc/NMD_inputs/scRNA_datasets/BRAUN_CANCERCELL_2021/download/BRAUN.h5Seurat", 
  "/SAN/colcc/NMD_inputs/scRNA_datasets/Chan_CancerCell_2021/Chan_combined.h5seurat", 
  "/SAN/colcc/NMD_inputs/scRNA_datasets/KRISHNA_CANCERCELL_2021/KRISHNA_RCC.h5seurat",
  "/SAN/colcc/NMD_inputs/scRNA_datasets/KIM_NCOMMS_2020/KIM_counts.h5Seurat",
  "/SAN/colcc/NMD_inputs/scRNA_datasets/YOST_NMED_2019/BASAL_CELL_CARCINOMA_raw_counts.h5Seurat", # failed to get raw count - we need it
  "/SAN/colcc/NMD_inputs/scRNA_datasets/LI_CELL_2019/Li.h5Seurat",
  "/SAN/colcc/NMD_inputs/scRNA_datasets/Maynard_CELL_2020/Maynard.h5seurat", 
  "/SAN/colcc/NMD_inputs/scRNA_datasets/BI_CANCERCELL_2021/download/Bi.h5seurat",
  "/SAN/colcc/NMD_inputs/scRNA_datasets/Pelka_CELL_2021/Pelka_CRC.h5Seurat",
  "/SAN/colcc/NMD_inputs/scRNA_datasets/QIAN_NATURECELL_2020/QIAN_NATURECELL_2020_combined4_raw_counts.h5Seurat",
  "/SAN/colcc/NMD_inputs/scRNA_datasets/VISHWAKARMA_COMMBIOL_2021/VISHWAKARMA_COMMBIOL_2021_combined6_raw_counts.h5Seurat",
  "/SAN/colcc/NMD_inputs/scRNA_datasets/Wu_NatureComm_2021/Wu.h5Seurat",
  "/SAN/colcc/NMD_inputs/scRNA_datasets/ZILIONIS_CELL_2019/Zilionis_rawcounts.h5Seurat", 
  "/SAN/colcc/NMD_inputs/scRNA_datasets/ZhangYY_CancerCell_2021/ZhangYY.h5seurat")

exp <- list()

# Step 1 : loading data 
#  Step 1' : build cell id lists from each experiment

for (i in 1:15) {  
  if (length(grep("h5.eurat", files[i]))) {
    data <- LoadH5Seurat(files[i])
  } else {
    data <- readRDS(files[i])
  }
  cells <- meta %>% filter(study == studies[i])
  print(paste("filtered cells in iteration ", i, sep = ""))
  exp[i] <- subset(data, cells = cells$cells) 
  cat(paste0("### Study : ", studies[i], " (", files[i], ")\n"))
  print (data)
  print (exp[[i]])
  
}

#saveRDS(data_cells, "cellids_in_rawdata.rds")
gs2new <- readRDS("/SAN/colcc/NMD_inputs/DC_results/tools/GeneSymbol2Current.rds")

#Set working directory
setwd ("/SAN/colcc/tigilab-general/evie")
studies <- studies[1:15]
nstudy <- length(studies)

getStats <- function (meta, exp) {
  sizes <- c()
  real.sizes <- c()
  gene.sizes <- c()
  for (i in 1:15) {
    cells <- meta %>% filter(study == studies[i])
    study  <- studies[i]
    sizes[i]      <- dim(cells)[1]  
    real.sizes[i] <- dim(exp[[i]])[2]
    gene.sizes[i] <- dim(exp[[i]])[1]
  }
  stats <- data.frame(id = studies, genes = gene.sizes, cell_meta = sizes, cell_real = real.sizes)
  return(stats)
}
stats <- getStats(meta, exp)

names(exp) <- studies
saveRDS(exp, "Tumour_subsetted_dataset.rds")

# Step 2 : adding Gene Symbols   - exp 5

exp$KRISHNA <- convertEnsemblToGeneSymbol_KeepAll2(exp$KRISHNA)

# Step 3 : Gene Symbol matching  

for (i in 1:nstudy) {
  exp[[i]] <- doConvertGeneName(exp[[i]], gs2new, shrink=T)
}

saveRDS(exp, "Tumour_subsetted_genename_updated.rds")



# Step 4 : remove low quality genes and do SCTransform

for (i in 1:nstudy) {
  study <- studies[i]
  cat ("processing ", study, "...\n")
  exp[[i]]$study <- study
  exp[[i]]$percent.mt <- PercentageFeatureSet(exp[[i]], pattern = "^MT-")
  exp[[i]] <- subset(exp[[i]], subset= nFeature_RNA > 200&nCount_RNA>300)
  exp[[i]] <- SCTransform(exp[[i]], method = "glmGamPoi", vars.to.regress = "percent.mt")
}

saveRDS(exp, "Tumour_subsetted_sctransformed.rds")

# *** refer to /home/murai/DC/scRNA_genelists

DONTRUN = FALSE

if (!DONTRUN) {
  
  
  ## Step 5: Preparing for Integration
  
  features <- SelectIntegrationFeatures(object.list = exp, nfeatures = 2000)     
  exp      <- PrepSCTIntegration(object.list = exp, anchor.features = features) 
  anchors  <- FindIntegrationAnchors(object.list = exp, normalization.method = "SCT", anchor.features = features) 
  saveRDS(anchors, "Tumour_subsetted_integrate_anchors_2800.rds")
}
  #Step 6
  if (TRUE) {
    source("/SAN/colcc/NMD_inputs/DC_results/tools/integrate_functions_mod.R")    
    exp <- IntegrateDataW(anchors, normalization.method = "SCT")
    saveRDS(exp, "Tumour_subsetted_integrated_beforeUmap.rds")
    
    exp <- RunPCA(exp, verbose = FALSE)
    exp <- RunUMAP(exp, reduction = "pca", dims = 1:50)
    exp <- FindNeighbors(exp, reduction = "pca", dims = 1:50)
    exp <- FindClusters(exp, resolution = 0.5)
    
    ## Step 7: Build an Elbow Plot for residual SSE of SNN Clustering
    
    sse_results <- data.frame()
    ncell <- dim(exp@reductions$pca@cell.embeddings)[1]
    DefaultAssay(exp) <- "integrated"
    # running
    for (r in c(0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 0.9, 1.1, 1.3, 1.4, 1.6, 1.8, 1.9, 2.5)) {
      exp <- FindClusters(exp, resolution = r)
      ncl <- length(unique(Idents(exp)))
      s   <- getSSEfromExp(exp)
      v   <- s/(ncell-1)
      sse_results <- rbind(sse_results, data.frame(res=r, nclust=ncl, sse=s, var=v))
      cat (paste(r, ncl, s, v), "\n")
    }
    
    # Prepare for gene expression analysis
    DefaultAssay(exp) <- "SCT"
    saveRDS(exp, "Tumourlowres_subsetted_integrated.rds")
  }

