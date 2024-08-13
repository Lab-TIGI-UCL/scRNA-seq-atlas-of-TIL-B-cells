library(Giotto)
checkGiottoEnvironment()
library(data.table)
library(rcartocolor)
library(reticulate)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(ggpubr)
library(UCell)


#Fig5A - proportion of tumour cells in 100 nearest neighbours

gem = load("~/SMI_Giotto_Object.RData") #publicly available Giotto object with metadata included
cosmx_meta = gem@cell_metadata$rna
subset_df <- subset(cosmx_meta, cell_type %in% c("B-cell", "plasmablast"),
                    select = c("cell_type", "prop_tumor_in_100_neighbors"))

ggplot(data=subset_df, mapping = aes(x = cell_type, y = prop_tumor_in_100_neighbors)) +
  geom_violin(aes(fill= cell_type), size=0.8) +
  geom_boxplot(color = "black", alpha = 0.9, outlier.shape = NA, size=1.7, fill = NA) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 28, hjust = 0.5),
    axis.text.x = element_text(size =24, angle = 45, hjust =1, color="black"),
    axis.text.y = element_text(size =24, color="black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(50, 250, 50, 250),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black"),
    strip.text.x = element_text(size = 20),
    legend.position = "none"
  ) +
  scale_y_continuous(limits=c(0,1.35),breaks = seq(0,1,0.25), expand = expansion(mult = c(0,0.01))) +
  labs(
    x = "Cell Type",
    y="Proportion of tumor cells in 100 neighbors"
  ) +
  scale_fill_manual(values = c("forestgreen", "mediumorchid1")) +
  stat_compare_means(label = "p.format", paired = F, size=6, label.y.npc = 0.98)


#Fig5B - distribution of B cells and plasma cells in spatial niches

load("~/SMI_Giotto_Object.RData")
metadata <- gem@cell_metadata$rna


### percentage of B and plasma in each niche out of all cells in corresponding niche
proportions <- metadata %>%
  filter(cell_type %in% c("B-cell", "plasmablast")) %>%  # Filter for B-cell and plasmablast
  group_by(Run_Tissue_name, niche, cell_type) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%  # Count cells for each cell_type
  left_join(metadata %>%                                 # Join with total cell count in each niche
              group_by(Run_Tissue_name, niche) %>%
              summarise(total_cells = n(), .groups = 'drop'), 
            by = c("Run_Tissue_name", "niche")) %>%
  mutate(proportion = round((cell_count / total_cells) * 100, 2)) %>%  # Calculate proportion
  ungroup()

proportions = proportions[proportions$total_cells > 50,]

compute_p_values <- function(data) {
  data %>%
    group_by(niche) %>%
    do(tidy(t.test(proportion ~ cell_type, data = .)))
}

# Get p-values
p_values <- compute_p_values(proportions)
format_pvalue <- function(pvalue) {
  if (pvalue < 0.001) {
    return(sprintf("p = %.1e", pvalue))
  } else {
    return(sprintf("p = %.3f", pvalue))
  }
}
p_values <- p_values %>%
  mutate(label = sapply(p.value, format_pvalue))

ggplot(proportions, aes(x = niche, y = proportion, fill = cell_type)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.5), 
              size = 0.5, alpha = 1) +
  geom_boxplot(width = 0.5, alpha = 0.8, outlier.shape = NA) + 
  geom_vline(xintercept = seq(1.5, length(unique(proportions$niche)) - 0.5, by = 1), 
             linetype = "dotted", color = "gray") +
  scale_fill_manual(values = c("#D7282C", "#4C8BBD")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = 'white'),
        axis.line = element_line(colour = "black"),
        #legend.position = "none"
  ) +
  geom_text(data = p_values, aes(x = niche, y = max(proportions$proportion) + 0.1, label = label),
            vjust = 0, size = 3.5, inherit.aes = FALSE)

#Fig5C - Cell proximity enrichment score between all cell types in the CosMx dataset for the Lung 5 Rep 1 sample.

#Giotto instructions
results_folder = '/CosMx/Lung5_Rep1/'
my_python_path = NULL
instrs = createGiottoInstructions(save_dir = results_folder, save_plot = TRUE,
                                  show_plot = FALSE, return_plot = FALSE,
                                  python_path = my_python_path)
data_path = '/CosMx/All+SMI+Flat+data/Lung5_Rep1/Lung5_Rep1-Flat_files_and_images'

# CosMx loading function
fov_join = createGiottoCosMxObject(cosmx_dir = data_path,
                                   data_to_use = 'subcellular', # only subcellular
                                   FOVs = 1:30,
                                   instructions = instrs)

# Set up vector of image names
id_set = sprintf("%02d", 1:30)
new_names = paste0("fov0", id_set)
image_names = paste0(new_names, '-image')

# Aggregate subcellular features
fov_join = calculateOverlapRaster(fov_join, feat_info = 'rna')
fov_join = overlapToMatrix(fov_join, feat_info = 'rna')

# Filtering and normalisation
fov_join <- filterGiotto(gobject = fov_join,
                         feat_type = 'rna',
                         expression_threshold = 1,
                         feat_det_in_min_cells = 5,
                         min_det_feats_per_cell = 20)
fov_join <- normalizeGiotto(gobject = fov_join,
                            feat_type = 'rna',
                            norm_methods = 'standard',
                            verbose = TRUE)
fov_join <- normalizeGiotto(gobject = fov_join,
                            feat_type = 'rna',
                            scalefactor = 5000,
                            verbose = TRUE,
                            norm_methods = 'pearson_resid',
                            update_slot = 'pearson')
fov_join = addStatistics(gobject = fov_join,
                         expression_values = 'normalized',
                         feat_type = 'rna')

# Dimension reduction
fov_join = calculateHVF(fov_join,
                        method = 'var_p_resid',
                        expression_values = 'pearson')
fov_join = runPCA(fov_join,
                  scale_unit = FALSE,
                  center = FALSE,
                  expression_values = 'pearson')
fov_join <- runUMAP(fov_join,
                    dimensions_to_use = 1:10,
                    n_threads = 8)

# Clustering 
fov_join <- createNearestNetwork(gobject = fov_join,
                                 dimensions_to_use = 1:10,
                                 k = 10)
fov_join <- doLeidenCluster(gobject = fov_join, resolution = 0.8,
                            name = "leiden_clus_res0.8", n_iterations = 1000)

# Spatial expression patterns
fov_join <- createSpatialNetwork(fov_join, method = 'kNN', k = 10, name = 'k10_network')
cell_proximities = cellProximityEnrichment(fov_join,
                                           cluster_column = 'cell_type',
                                           spatial_network_name = 'k10_network',
                                           adjust_method = 'fdr',
                                           number_of_simulations = 2000)
cellProximityHeatmap(fov_join,
                     CPscore = cell_proximities,
                     order_cell_types = T, scale = T,
                     color_breaks = c(-1.5, 0, 1.5),
                     color_names = c('blue', 'white', 'red'),
                     save_param = list(base_height = 2000, base_width = 2000,
                                       save_name = 'Fig3C'))




#Fig5D - Visualise all cell types spatially
spatPlot2D(gobject = fov_join,
           cell_color = "cell_type",
           cell_color_code = cell_type_colors,
           show_image = F,
           point_size = 0.3,
           background_color = "black",
           show_legend = T,
           save_param = list(base_height = 8,
                             save_name = 'Fig5D1'))

#subset region of interest
gem_meta_lung5 = metadata[metadata$Run_Tissue_name == "Lung5_Rep1",]
gem_meta_lung5 = gem_meta_lung5[gem_meta_lung5$fov %in% c(1:2, 6:7, 11:12),]
fov_join = subsetGiotto(fov_join, cell_ids = gem_meta_lung5$cell_ID)

cell_type_names = sort(unique(gem_meta_lung5$cell_type))
cell_type_colors = c("#ed2224", "#f6989a", "pink", "#9b4e43", "#f78e45", "#bfbebe", "#1f5429", "#fed307", 
                     "#6abd45", "#b1cc71", "#724fa0", "#0f1131", "#1f9799", "#4a97d2", 
                     "#a3519f", "#fce5d8")
cell_type_colors = setNames(cell_type_colors, cell_type_names)

spatPlot2D(gobject = fov_join,
           cell_color = "niche",
           #cell_color_code = cell_type_colors,
           #select_cell_groups = c("B_cell", "plasmablast", "tumour",
           #"CD4_T_cell", "CD8_T_cell"),
           #show_other_cells = F,
           show_image = F,
           image_name = image_names,
           point_size = 1,
           background_color = "black",
           show_legend = T,
           save_plot = T,
           return_plot = F,
           save_param = list(base_height = 8,
                             save_name = 'Fig5D2'))


#Fig5E - atlas-derived signatures

gem_meta_lung5 = gem_meta_lung5[gem_meta_lung5$fov == 6,]
fov_join = subsetGiotto(fov_join, cell_ids = gem_meta_lung5$cell_ID)
cells_to_keep = gem_meta_lung5[gem_meta_lung5$cell_type %in% c("B_cell")]$cell_ID
mat = fov_join@expression[["cell"]][["rna"]][["raw"]]@exprMat

b.gene.sets = list(naiveBcell_sig = c("CD19", "MS4A1", "SELL", "TCL1A"),
                   memoryBcell_sig = c("CD19", "MS4A1", "CD27", "HLA-DPB1"),
                   mt1x_plasma_sig = c("JCHAIN", "MT1X", "MT2A", "IGHG1"))

scores <- ScoreSignatures_UCell(mat[,gem_meta_lung5[gem_meta_lung5$cell_type %in% c("B_cell")]$cell_ID], features = b.gene.sets)
fov_join_b = subsetGiotto(fov_join, cell_ids = cells_to_keep)
fov_join_b@cell_metadata$cell$rna@metaDT$sig_naiveB = scores[,1]
fov_join_b@cell_metadata$cell$rna@metaDT$sig_memB = scores[,2]

spatInSituPlotPoints(fov_join_b,
                     show_image = F,
                     image_name = image_names,
                     #feats = list('rna' = c("MS4A1", "SELL", "TNFSF13B", "CD40")),
                     spat_unit = 'cell',
                     point_size = 0.3,
                     show_polygon = TRUE,
                     use_overlap = F,
                     polygon_feat_type = 'cell',
                     polygon_color = 'black',
                     polygon_line_size = 0.03,
                     polygon_fill = 'sig_memB',
                     polygon_fill_as_factor = F,
                     polygon_fill_gradient = viridis(20, option = "A"),
                     #polygon_fill_gradient_midpoint = 0.5,
                     #polygon_fill_code = cell_type_colors,
                     #polygon_alpha = 0.6,
                     save_plot = F,
                     return_plot = T,
                     save_param = list(base_height = 8, base_width = 8,
                                       save_name = 'Fig5E1'))

scores <- ScoreSignatures_UCell(mat[,gem_meta_lung5[gem_meta_lung5$cell_type %in% c("plasmablast")]$cell_ID], features = b.gene.sets)
fov_join_plasma = subsetGiotto(fov_join, cell_ids = gem_meta_lung5[gem_meta_lung5$cell_type %in% c("plasmablast")]$cell_ID)
fov_join_plasma@cell_metadata$cell$rna@metaDT$sig_mt1x = scores[,3]


#Fig 5F - Ligand receptor spatial mapping

mat2 = mat[c("SELL", "TNFSF13B", "CD40"), cells_to_keep]
b_id = gem_meta_lung5[gem_meta_lung5$cell_type %in% c("B_cell")]$cell_ID
mat_b = mat2[c("SELL", "CD40"),cells_to_keep]
non_zero_both_rows <- which(mat_b[1, ] != 0 & mat_b[2, ] != 0)
cd40_sell_b <- colnames(mat_b[, non_zero_both_rows])

mat_mem = mat[c("CD27", "TNFRSF13B"),b_id]
non_zero_both_rows <- which(mat_mem[1, ] != 0 & mat_mem[2, ] != 0)
cd27_tnfrsf13b_b <- colnames(mat_b[, non_zero_both_rows])

mat_mt1x = mat[c("MT1X", "CD63"), gem_meta_lung5[gem_meta_lung5$cell_type %in% c("plasmablast")]$cell_ID]
non_zero_both_rows <- which(mat_mt1x[1, ] != 0 & mat_mt1x[2, ] != 0)
cd63_mt1x_b <- colnames(mat_mt1x[, non_zero_both_rows])

cd4_id = gem_meta_lung5[gem_meta_lung5$cell_type %in% c("CD4_T_cell")]$cell_ID
mat_cd4 = mat["TIMP1",cd4_id]
tnfsf13b_cd4 = names(mat_cd4[which(mat_cd4 != 0)])
timp1_cd4 = names(mat_cd4[which(mat_cd4 != 0)])

fov_join_naive = subsetGiotto(fov_join, cell_ids = c(cd40_sell_b, tnfsf13b_cd4))
fov_join_mem = subsetGiotto(fov_join, cell_ids = c(cd27_tnfrsf13b_b, tnfsf13b_cd4))
fov_join_mt1x = subsetGiotto(fov_join, cell_ids = c(cd63_mt1x_b, timp1_cd4))

spatPlot2D(gobject = fov_join_naive,
           cell_color = "cell_type",
           show_image = F,
           image_name = image_names,
           point_size = 1,
           background_color = "black",
           cell_color_code = cell_type_colors,
           select_cell_groups = c("B_cell", "CD4_T_cell"),
           show_other_cells = F,
           show_legend = T,
           show_plot = T,
           show_network = T,
           spatial_network_name = "k10_network",
           network_color = "white",
           network_alpha = 0.6,
           save_plot = F,
           return_plot = T,
           save_param = list(base_height = 8,
                             base_width = 20,
                             save_name = 'Fig5F1'))
