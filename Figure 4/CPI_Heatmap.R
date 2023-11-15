#CPI Heatmap

library(pheatmap)

master_sample_sheet6 <- read.csv("~/master_sample_sheet6.csv")

#Convert row names into column
d <- cbind(rownames(exp_means_df), data.frame(exp_means_df, row.names=NULL))
colnames(d)[1] = "Patient"

#Add in Studies
df = merge(x = d,
           y = master_sample_sheet6[, c("Patient", "Study")], by = "Patient")


#Rename studies to include cancer type
df$Study[df$Study == 'SNYDER_PLOSMED_2017'] <- 'BLADDER_SNYDER_PLOSMED_2017'
df$Study[df$Study == 'VANALLEN_SCIENCE_2015'] <- 'MELANOMA_VANALLEN_SCIENCE_2015'
df$Study[df$Study == 'SNYDER_NEJM_2014'] <- 'MELANOMA_SNYDER_NEJM_2014'
df$Study[df$Study == 'Liu_NatureMedicine_2019'] <- 'MELANOMA_LIU_NATUREMEDICINE_2019'
df$Study[df$Study == 'HUGO_CELL_2016'] <- 'MELANOMA_HUGO_CELL_2016'
df$Study[df$Study == 'KIM_NMED_2018'] <- 'GASTRIC_KIM_NMED_2018'
df$Study[df$Study == 'MARIATHASAN_NATURE_2018'] <- 'BLADDER_MARIATHASAN_NATURE_2018'
df$Study[df$Study == 'MCDERMOT_NMED_2018'] <- 'RENAL_MARIATHASAN_NATURE_2018'
df$Study[df$Study == 'RIAZ_CELL_2017'] <- 'MELANOMA_RIAZ_CELL_2017'
df$Study[df$Study == 'SHIM_AOO_2020'] <- 'LUNG_SHIM_AOO_2020'


#Rename clusters to identities
colnames(df)[2] = "Resting Memory B cells"
colnames(df)[3] = "Activated Memory B cells"
colnames(df)[4] = "Conventional Plasma cells"
colnames(df)[5] = "Naive B cells"
colnames(df)[6] = "Stressed Plasma cells"
colnames(df)[7] = "Atypical Memory B cells"
colnames(df)[8] = "IGKC-high Plasma/Plasmablasts"
colnames(df)[9] = "Transitional/Proliferative 1"
colnames(df)[10] = "Transitional/Proliferative 2"
colnames(df)[11] = "MT1X-high Plasma/Plasmablasts"


shap_subclass = function(x) {s
  return (shapiro.test(x)$p.value) }


head(df)
a = df %>% group_by(Study, survival) %>% 
  select(-Patient, -cancer_type) %>%
  summarise_all(.funs = shap_subclass) %>%
  as.data.frame()

# define vector of celltypes
celltypes = c("Naive B cells", "Resting Memory B cells", 
              "Activated Memory B cells", "Atypical Memory B cells", 
              "Transitional/Proliferative 1", "Transitional/Proliferative 2", 
              "MT1X-high Plasma/Plasmablasts", "IGKC-high Plasma/Plasmablasts", 
              "Conventional Plasma cells", "Stressed Plasma cells")

# define vector of studies
studies = c("BLADDER_SNYDER_PLOSMED_2017", "MELANOMA_VANALLEN_SCIENCE_2015", 
            "MELANOMA_SNYDER_NEJM_2014", "MELANOMA_LIU_NATUREMEDICINE_2019", 
            "MELANOMA_HUGO_CELL_2016", "GASTRIC_KIM_NMED_2018", 
            "BLADDER_MARIATHASAN_NATURE_2018", "RENAL_MARIATHASAN_NATURE_2018",
            "MELANOMA_RIAZ_CELL_2017", "LUNG_SHIM_AOO_2020"
)



out = as.data.frame(matrix(NA, nrow = length(celltypes), ncol = length(studies)))
rownames(out) = celltypes
colnames(out) = studies

for (i in 1:length(celltypes)) {
  curr_cell = celltypes[i]
  
  for (j in 1:length(studies)) {
    curr_study = studies[j]
    
    htest =  wilcox.test(x = df[(df$Study == studies[j]) & (df$survival == "response"), curr_cell],
                         y = df[(df$Study == studies[j]) & (df$survival == "no_response"), curr_cell])
    
    out[curr_cell, curr_study] = htest$p.value
  }
}


#Meta-analysis across all cohorts (switch out cluster name)
tmp <- df %>% select("survival", "Stressed Plasma cells")
colnames(tmp) <- c("survival", "cell")
resp <- tmp %>% filter(survival == "response")
resp <- resp$cell

nresp <- tmp %>% filter(survival == "no_response")
nresp <- nresp$cell

wilcox.test(x = resp, y = nresp)

p <- c(0.07402, 0.9041, 0.05551, 0.1875, 0.3567, 0.3898, 0.9352, 0.04483, 0.5092, 0.2907)
p <- p.adjust(p)



#Pheatmap

df2 <- df %>% 
  select(-Patient, -cancer_type) %>% 
  group_by(Study, survival) %>% summarise_all(.funs = median) %>% ungroup() %>% 
  pivot_longer(!Study:survival, names_to='cluster', values_to='mean_exprs') %>%   
  pivot_wider(names_from='survival', values_from = 'mean_exprs') %>%  
  mutate(ratio = response / no_response)


#Cap all ratios
df2$ratio[df2$ratio >= 20] <- 20
df2$ratio[df2$ratio <= -20] <- -20

paletteLength <- 400
myColor <- colorRampPalette(c("darkred", "white", "darkblue"))(paletteLength)
myBreaks <- c(seq(min(df2$ratio), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(df2$ratio)/paletteLength, max(df2$ratio),
                  length.out=floor(paletteLength/2)))


mat = reshape2::dcast(data = df2, formula = Study ~ cluster) %>% as.data.frame()
mat <- mat[, c("Study", "Naive B cells", "Resting Memory B cells", "Activated Memory B cells", 
               "Atypical Memory B cells", "Transitional/Proliferative 1", "Transitional/Proliferative 2", 
               "MT1X-high Plasma/Plasmablasts", "IGKC-high Plasma/Plasmablasts", "Conventional Plasma cells", 
               "Stressed Plasma cells")] %>% as.data.frame()
rownames(mat) = mat$Study
mat$Study = NULL
mat = mat %>% t() %>% as.data.frame()
mat = as.matrix(mat)


pheatmap(mat, color=myColor, breaks=myBreaks, border_colour = "white", cluster_cols = F, cluster_rows = F, fontsize_row=10, fontsize_col=8, 
         display_numbers = matrix(ifelse(out < 0.05, "*", ""), nrow(out)))

