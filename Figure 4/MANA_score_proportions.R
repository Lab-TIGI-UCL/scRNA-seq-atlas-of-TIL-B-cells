#MANA score proportion analysis


df <- data.frame(
  cluster = c("CD4 Tregs", "Naive B cells", "Resting Memory B cells", "Activated Memory B cells", "Atypical Memory B cells", "Transitional/Proliferative 1", "Transitional/Proliferative 2", "MT1X-high Plasma/Plasmablasts", "IGKC-high Plasma/Plasmablasts", "Conventional Plasma cells", "Stressed Plasma cells", "CD4 Tregs", "Naive B cells", "Resting Memory B cells", "Activated Memory B cells", "Atypical Memory B cells", "Transitional/Proliferative 1", "Transitional/Proliferative 2", "MT1X-high Plasma/Plasmablasts", "IGKC-high Plasma/Plasmablasts", "Conventional Plasma cells", "Stressed Plasma cells"),
  proportion = c(0.0979624731, 0.11735226, 0.25012658, 0.18206614, 0.06734406, 0.03965345, 0.01600902, 0.02346403, 0.06032059, 0.16208643, 0.08157743, 0.021691897, 0.10347634, 0.19901096, 0.20670528, 0.07140136, 0.02915075, 0.01028003, 0.02251695, 0.06869792, 0.17981057, 0.10894984),
  MANA = c("Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "Low", "High", "High", "High", "High", "High", "High", "High", "High", "High", "High", "High"))

df$cluster <- factor(df$cluster,levels = c("CD4 Tregs", "Naive B cells", "Resting Memory B cells", "Activated Memory B cells", "Atypical Memory B cells", "Transitional/Proliferative 1", "Transitional/Proliferative 2", "MT1X-high Plasma/Plasmablasts", "IGKC-high Plasma/Plasmablasts", "Conventional Plasma cells", "Stressed Plasma cells"))

ggplot(data=df, aes(x=cluster, y=proportion, fill=MANA)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values = c("High" = "#28788D","Low" = "#D64550")) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1)
  )


#Adding proportion ratio (Low/High, calculated manually)

df <- data.frame(
  cluster = c("CD4 Tregs", "Naive B cells", "Resting Memory B cells", "Activated Memory B cells", "Atypical Memory B cells", "Transitional/Proliferative 1", "Transitional/Proliferative 2", "MT1X-high Plasma/Plasmablasts", "IGKC-high Plasma/Plasmablasts", "Conventional Plasma cells", "Stressed Plasma cells"),
  ratio = c(4.5160860343, 1.1340975145, 1.256848266, 0.8808006259, 0.9431761524, 0.1360289186, 1.5572931207, 1.0420607587, 0.878055551, 0.9014288203, 0.7487613566)
)

df$cluster <- factor(df$cluster,levels = c("CD4 Tregs", "Naive B cells", "Resting Memory B cells", "Activated Memory B cells", "Atypical Memory B cells", "Transitional/Proliferative 1", "Transitional/Proliferative 2", "MT1X-high Plasma/Plasmablasts", "IGKC-high Plasma/Plasmablasts", "Conventional Plasma cells", "Stressed Plasma cells"))

ggplot(data=df, aes(x=cluster, y=ratio)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1)
  )

ggplot(data=df, aes(x=cluster, y=ratio)) + geom_point()