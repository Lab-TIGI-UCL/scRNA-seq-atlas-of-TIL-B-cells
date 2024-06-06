library(Giotto)
checkGiottoEnvironment()
library(data.table)
library(rcartocolor)
library(reticulate)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(ggpubr)


#Fig3A - proportion of tumour cells in 100 nearest neighbours

gem = load("C:/CosMx/SMI_Giotto_Object.RData") #publicly available Giotto object with metadata included
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


#Fig3B - distribution of B cells and plasma cells in spatial niches

load("~/SMI_Giotto_Object.RData")
metadata <- gem@cell_metadata$rna


### percentage of B and plasma in each niche out of all cells in corresponding niche
head(metadata)

cell_niche = metadata %>% group_by(cell_type, niche) %>% 
  summarise(count = n()) %>% as.data.frame()
head(cell_niche)


niche_tot = metadata %>% group_by(niche) %>% summarise(total_count = n()) %>% as.data.frame()
head(niche_tot)

cell_niche = merge(x = cell_niche, y = niche_tot, by = "niche")


### percentage of cells in each niche
cell_niche$count_pct = 100*(cell_niche$count / cell_niche$total_count)
head(cell_niche)


### Plotting
cell_niche$cell_clean = "NA"
cell_niche[cell_niche$cell_type == "B-cell", "cell_clean"] = "B-cell"
cell_niche[cell_niche$cell_type == "plasmablast", "cell_clean"] = "plasmablast"

cell_niche$niche <- factor(cell_niche$niche, levels = c("immune","macrophages", "neutrophils", "lymphoid structure", "myeloid-enriched stroma", "plasmablast-enriched stroma", "stroma", "tumor-stroma boundary", "tumor interior"))

ggplot(data=cell_niche,
       mapping=aes(x=niche, y=count_pct, alluvium = cell_type, stratum = cell_type)) +
  geom_alluvium(aes(fill=cell_clean), color="black", alpha = 1, position = position_dodge(width=0)) +
  geom_stratum(width = 1/12) + ggtitle("Distribution of B & Plasma cells out of all cells in each niche")


#Fig3C - Cell proximity enrichment score between all cell types in the CosMx dataset for the Lung 5 Rep 1 sample.

#Giotto instructions
results_folder = 'C:/CosMx/Lung5_Rep1/'
my_python_path = NULL
instrs = createGiottoInstructions(save_dir = results_folder, save_plot = TRUE,
                                  show_plot = FALSE, return_plot = FALSE,
                                  python_path = my_python_path)
data_path = 'C:/CosMx/All+SMI+Flat+data/Lung5_Rep1/Lung5_Rep1-Flat_files_and_images'

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




#Fig3D - Visualise all cell types spatially
spatPlot2D(gobject = fov_join,
           cell_color = "cell_type",
           cell_color_code = cell_type_colors,
           show_image = F,
           point_size = 0.3,
           background_color = "black",
           show_legend = T,
           save_param = list(base_height = 8,
                             save_name = 'Fig3D'))


#Fig3E - Visualise cell types of interest spatially
  
spatPlot2D(gobject = fov_join,
           cell_color = "cell_type",
           select_cell_groups = c("B-cell", "plasmablast"),
           show_image = F,
           point_size = 0.3,
           background_color = "black",
           show_legend = T,
           save_param = list(base_height = 8,
                             save_name = 'Fig3D'))

#Fig3F - visualisation of transcripts of interest in FOV2

spatInSituPlotPoints(fov_join,
                     show_image = TRUE,
                     image_name = image_names,
                     feats = list('rna' = c('CD19', 'EPCAM', 'MS4A1', 'TCL1A', 'SELL', 'HLA-DPB1', 'CD27', 'CD69', 'CD38', 'JCHAIN', 'HSPA1A', "TNFRSF13B", 'IGKC', 'MT1X', 'MKI67')),
                     feats_color_code = pal15,
                     spat_unit = 'cell',
                     point_size = 0.8,
                     show_polygon = TRUE,
                     use_overlap = FALSE,
                     polygon_feat_type = 'cell',
                     polygon_color = 'white',
                     polygon_line_size = 0.03,
                     save_param = list(base_height = 5,
                                       save_name = 'All_markers'))




