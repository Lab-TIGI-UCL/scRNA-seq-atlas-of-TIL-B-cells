#Correlation plot

metrics_output<-read.table("~/meta_analysis_input_data.txt",sep="\t",header=T,stringsAsFactors=F)

coverage<-read.table("~/coverage_purity_stats.csv",header=T,sep=",",stringsAsFactors=F)
metrics_output2 <- metrics_output[!(is.na(metrics_output$TMB)) & !(is.na(metrics_output$CD8A)),]
metrics_output2 <- merge(metrics_output2,coverage,by.x="case",by.y="Patient")
exp_means_df$case <- rownames(exp_means_df)
metrics_output2 <- merge(metrics_output2, exp_means_df, by.x = "case")
metrics_output2$coverage<-NULL

only_data <- full_join(exp_means_df, metrics_output2, by = "case") %>% select(-c(case))
only_data<-select_if(only_data, is.numeric)

only_data <- only_data[,c(1:30)] 
only_data <- only_data[ -c(19, 21:23, 25) ]

#Reorder columns 
library(tidyverse)
library(Hmisc)
library(rstatix)
library(ggsignif)
library(ggpubr)
library(ggrepel)
library(reshape2)

only_data <- only_data[, c(4, 1, 2, 6, 8, 9, 10, 7, 3, 5, 11:25)]

cors <- rcorr(as.matrix(only_data), type = "spearman") 

rho <- as.data.frame(cors$r)
rho$names <- colnames(only_data)

rho_melt <- reshape2::melt(rho) %>% dplyr::rename(rho_x = names,
                                                  rho_y = variable,
                                                  rho = value)

p <- as.data.frame(cors$P)
p$names <- colnames(only_data)

p_melt <- reshape2::melt(p) %>% dplyr::rename(p_x = names,
                                              p_y = variable,
                                              p = value)
p_melt$rounded_p <- round(p_melt$p, digits = 5)

p_melt1 <- p_melt[order(p_melt$rounded_p),]
p_melt1$q <- p.adjust(p_melt1$rounded_p, method = "fdr")  

p_sig <- p_melt1 %>% add_significance(p.col = "q")
p_sig$q.signif[p_sig$q.signif == "ns"] <- ""

clusters <- colnames(only_data)

ggplot(data = rho_melt, aes(x=rho_x, y=rho_y, fill = rho)) +
  geom_tile(stat = "identity", color = "black") +
  coord_fixed() +
  theme_minimal() +
  scale_fill_gradient2(mid="white",low="darkblue",high="darkred", limits=c(-1,1), space = "Lab",
                       midpoint = 0, name = "Rho") +
  geom_text(data=p_sig, inherit.aes = F, mapping = aes(x=p_x, y=p_y, label = q.signif), size = 4) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 10, color = "black", angle = 60, hjust = 1),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 15, color = "black"),
    legend.text = element_text(size = 12, color = "black")) +
  scale_x_discrete(limits = clusters) +
  scale_y_discrete(limits = clusters)   
