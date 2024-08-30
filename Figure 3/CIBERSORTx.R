#CIBERSORTx

########### R script for CPI1000. ###########
########### Created by Danwen Qian on 17-05-2024. ###########
########### Last modified by Danwen Qian on 17-05-2024. ###########


BCELL<-readRDS("~/BCELL.rds")
sampled.cells <- sample(x =colnames(BCELL) , size =10000, replace = F)
BCELL.downsample<-subset(BCELL,cells=sampled.cells)
counts<-as.matrix(BCELL.downsample@assays$RNA@counts)
colnames(counts)<-as.vector(BCELL.downsample$clusters)
keep_feature <- rowSums(counts > 0) > 300
counts <- counts[keep_feature, ]
write.table(counts,sep="\t","~/BCELL_reference.txt",col.names=NA)


library(ggplot2)
library(dplyr)
library(ggsci)
library(forestplot)
library(ggbeeswarm)
library(ggpubr)
library(scales)

BCELL<-read.csv("CIBERSORTx_BCELL.csv",row.names = 1)
BCELL_col<-read.csv("CIBERSORTx_BCELL.csv",row.names = 1,header  = F)
colnames(BCELL)<-BCELL_col[1,]
BCELL$response[BCELL$response=="non_responder"]<-0
BCELL$response[BCELL$response=="responder"]<-1
BCELL$response<-as.numeric(BCELL$response)

###Forest plot
results<-data.frame(matrix(ncol=6,nrow=10))
results[,1]<-colnames(BCELL)[1:10]
colnames(results)<-c("Var","OR_mean","OR_1","OR_2","Pvalue","OR")
for(i in 1:10){
  this_dat<-BCELL[,c(i,14,15)]
  names(this_dat)[1]<-"BCELL"
  fit_results<-glm(response ~ BCELL,
                   data = this_dat,
                   family = binomial(link = "logit"))
  fit.result<-summary(fit_results)
  results[i,2]<-fit.result$coefficients[2,1]
  results[i,3:4]<-confint(fit_results)[2,]
  results[i,5]<-fit.result$coefficients[2,4]
}
results$OR<-paste0(round(results$OR_mean,2),
                   " (",
                   round(results$OR_1,2),
                   "-",
                   round(results$OR_2,2),
                   ")")
results$Pvalue<-round(results$Pvalue,3)
results[2:11,]<-results[1:10,]
results[1,]<-c("Variable","","","","P-value","OR (95% CI)")

write.csv(results,"Cibersortx_BCELL_OR.csv",row.names = F)

pdf("forestplot_BCELL.pdf",  height=6, width=6)
forestplot(labeltext=as.matrix(results[,c(1,6,5)]),
           mean=as.numeric(results$OR_mean) ,
           lower=as.numeric(results$OR_1),
           upper=as.numeric(results$OR_2),
           zero=0,
           boxsize=0.2,
           graph.pos=2,clip=c(-10,20),
           lineheight = unit(7,'mm'),
           colgap=unit(2,'mm'),
           lwd.zero=1.5,
           lwd.ci=2,
           col=fpColors(box='#458B00',
                        summary='#8B008B',
                        lines = 'black',
                        zero = '#7AC5CD'),
           xlab="OR",
           lwd.xaxis =1,
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.85),
                            xlab  = gpar(cex = 0.8),
                            cex = 0.9),
           lty.ci = "solid",
           title = "Forestplot",
           line.margin = 0.08)
dev.off()


BCELL<-read.csv("CIBERSORTx_BCELL.csv",row.names = 1)

cbPalette2 <- c("#377EB8","#E41A1C")
for(i in c(1:10)){
  this_data<-BCELL[,c(i,14,15)]
  y<-colnames(this_data)[1]
  colnames(this_data)[1]<-"BCELL"
  
  p1<-ggplot(data=this_data,aes(x=response, y=BCELL,color=response))+theme_classic()+geom_boxplot(lwd=0.2,outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(strip.background = element_rect(size = 0.4),axis.line = element_line(size = 0.2),legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means(aes(label = sprintf("P = %.3f", as.numeric(..p.format..))),hjust = 0)+scale_fill_manual(values=cbPalette2)+scale_color_manual(values = cbPalette2)+theme(plot.margin = unit(c(1,1,1,1), "lines"))+ ggtitle("")+xlab("")+scale_y_continuous(trans=pseudo_log_trans(base = 10))+ylab(BCELL_col[1,i])
  ggsave(paste0("Cibersortx_BCELL/",y,".pdf"),plot = p1, height=4, width=4,dpi=600)
  
  p1<-ggplot(data=this_data,aes(x=response, y=BCELL,color=response))+theme_classic()+geom_boxplot(lwd=0.2,outlier.shape = NA) +geom_quasirandom(size = 2, shape=18,alpha=0.5)+ theme(strip.background = element_rect(size = 0.4),axis.line = element_line(size = 0.2),legend.position = "none",text = element_text(size=10),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+stat_compare_means(aes(label = sprintf("P = %.3f", as.numeric(..p.format..))),hjust = 0)+scale_fill_manual(values=cbPalette2)+scale_color_manual(values = cbPalette2)+theme(plot.margin = unit(c(1,1,1,1), "lines"))+ ggtitle("")+xlab("")+facet_wrap(~cancer_type,nrow = 2,scales = "free_y")+scale_y_continuous(trans=pseudo_log_trans(base = 10))+ylab(BCELL_col[1,i])
  ggsave(paste0("Cibersortx_BCELL/",y,"_by_cancer.pdf"),plot = p1,height=6, width=8,dpi=600)
  
}