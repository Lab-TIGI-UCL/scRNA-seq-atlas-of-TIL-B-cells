#Pie Charts - Doughnuts

#CANCER TYPE by number of cells
df <- data.frame(
  group = c("Breast Cancer", "Basal Cell Carcinoma", "Colorectal Cancer", "Lung Cancer", "Melanoma", "Ovarian Cancer", "Renal Cell Carcinoma"),
  cells = c(37820, 6508, 32016, 38664, 6743, 735, 3615)
)
head(df)
df$fraction = df$cells / sum(df$cells)
df$ymax = cumsum(df$fraction)
df$ymin = c(0, head(df$ymax, n=-1))
df$labelPosition <- (df$ymax + df$ymin) / 2
df$label <- paste0(df$group, "\n cells: ", df$cells)



ggplot(df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=group)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=3) +
  scale_fill_manual(values = c("Breast Cancer" = "#E66C37","Basal Cell Carcinoma" = "#28788D",
                               "Colorectal Cancer" = "#D64550","Lung Cancer" = "#37A794",
                               "Melanoma" = "#4C5D8A","Renal Cell Carcinoma"="#A43B76",
                               "Ovarian Cancer" = "#5ECBC8")) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

ggsave(
  "pie_cancertype_paper.pdf",
  plot = last_plot(),
  width = 250,
  height = 250,
  units = "mm",
  dpi = 300) 


#CANCER TYPE by number of samples
df <- data.frame(
  group = c("Breast Cancer", "Basal Cell Carcinoma", "Colorectal Cancer", "Lung Cancer", "Melanoma", "Ovarian Cancer", "Renal Cell Carcinoma"),
  cells = c(51, 18, 71, 159, 24, 5, 30)
)
head(df)
df$fraction = df$cells / sum(df$cells)
df$ymax = cumsum(df$fraction)
df$ymin = c(0, head(df$ymax, n=-1))
df$labelPosition <- (df$ymax + df$ymin) / 2
df$label <- paste0(df$group, "\n cells: ", df$cells)



ggplot(df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=group)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=3) +
  scale_fill_manual(values = c("Breast Cancer" = "#E66C37","Basal Cell Carcinoma" = "#28788D",
                               "Colorectal Cancer" = "#D64550","Lung Cancer" = "#37A794",
                               "Melanoma" = "#4C5D8A","Renal Cell Carcinoma"="#A43B76",
                               "Ovarian Cancer" = "#5ECBC8")) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

ggsave(
  "pie_cancertype_samples_paper.pdf",
  plot = last_plot(),
  width = 250,
  height = 250,
  units = "mm",
  dpi = 300) 


#TREATMENT
df <- data.frame(
  group = c("Naive", "Chemo", "CPI", "Target", "Chemo+CPI", "Chemo+Target", "CPI+Target", "Chemo+CPI+Target"),
  cells = c(77389, 11088, 8761, 1130, 25972, 453, 1113, 195))

head(df)


df$fraction = df$cells / sum(df$cells)
df$ymax = cumsum(df$fraction)
df$ymin = c(0, head(df$ymax, n=-1))
df$labelPosition <- (df$ymax + df$ymin) / 2
df$label <- paste0(df$group, "\n cells: ", df$cells)



ggplot(df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=group)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=3) +
  scale_fill_manual(values = c("Naive" = "#37A794","Chemo" = "#A43B76",
                               "CPI" = "#E66C37","Target" = "#D64550",
                               "Chemo+CPI" = "#4C5D8A","Chemo+Target"="#28788D",
                               "CPI+Target" = "#5ECBC8", "Chemo+CPI+Target" = "#FF977E")) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

ggsave(
  "pie_treatment_paper.pdf",
  plot = last_plot(),
  width = 250,
  height = 250,
  units = "mm",
  dpi = 300) 



#TISSUE
df <- data.frame(
  group = c("Tumor", "Metastasis"),
  cells = c(108757, 17344))

head(df)


df$fraction = df$cells / sum(df$cells)
df$ymax = cumsum(df$fraction)
df$ymin = c(0, head(df$ymax, n=-1))
df$labelPosition <- (df$ymax + df$ymin) / 2
df$label <- paste0(df$group, "\n cells: ", df$cells)



ggplot(df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=group)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=3) +
  scale_fill_manual(values = c("Tumor" = "#28788D","Metastasis" = "#5ECBC8")) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

ggsave(
  "pie_tissue_paper.pdf",
  plot = last_plot(),
  width = 250,
  height = 250,
  units = "mm",
  dpi = 300) 


#Study
df <- data.frame(
  group = c("Azizi", "Bi", "Braun", "Chan", "Kim", "Krishna", "Li", "Maynard", "Pelka", "Qian", "Vishwakarma", "Wu", "Yost", "Zhang", "Zilionis" ),
  cells = c(883, 1335, 1459, 5825, 15186, 659, 6743, 1706, 28334, 13258, 162, 3180, 6508, 34394, 6469))

head(df)


df$fraction = df$cells / sum(df$cells)
df$ymax = cumsum(df$fraction)
df$ymin = c(0, head(df$ymax, n=-1))
df$labelPosition <- (df$ymax + df$ymin) / 2
df$label <- paste0(df$group, "\n cells: ", df$cells)



ggplot(df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=group)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=3) +
  scale_fill_manual(values = c("Azizi" = "#9E9E9E","Bi" = "#D64550",
                               "Braun" = "#FF977E","Chan" = "#A43B76",
                               "Kim" = "#9A64A0","Krishna"="#750985",
                               "Li" = "#5ECBC8", "Maynard" = "#28788D", "Pelka" = "#37A794", 
                               "Qian" = "#4C5D8A", "Vishwakarma" = "#039BE5", "Wu" = "#BF360C", 
                               "Yost" = "#8E24AA", "Zhang" = "#E66C37", "Zilionis" = "#BA68C8")) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

ggsave(
  "pie_study_paper.pdf",
  plot = last_plot(),
  width = 250,
  height = 250,
  units = "mm",
  dpi = 300) 


#Platform
df <- data.frame(
  group = c("10X", "GEXSCOPE", "inDrop", "MARS-seq", "Smart-seq2"),
  cells = c(107120, 3180, 7352, 6743, 1706))

head(df)


df$fraction = df$cells / sum(df$cells)
df$ymax = cumsum(df$fraction)
df$ymin = c(0, head(df$ymax, n=-1))
df$labelPosition <- (df$ymax + df$ymin) / 2
df$label <- paste0(df$group, "\n cells: ", df$cells)



ggplot(df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=group)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=3) +
  scale_fill_manual(values = c("10X" = "#37A794","inDrop" = "#D64550",
                               "MARS-seq" = "#FF977E","Smart-seq2" = "#A43B76",
                               "GEXSCOPE" = "#9A64A0")) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

ggsave(
  "pie_platform_paper.pdf",
  plot = last_plot(),
  width = 250,
  height = 250,
  units = "mm",
  dpi = 300) 