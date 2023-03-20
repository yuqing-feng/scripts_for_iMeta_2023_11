# Metabolome

# Step 1: Load the data 
# Step 2: filt no hit data
# Step 3: replace 1 with 0
# Step 4: remove less than 30%
# Step 5: format the data with METLIN and HMDB:
#       METLIN for PCoA
#       HMDB for Sig and Pathway after glm

library(openxlsx)
m1 <- read.xlsx('~/metabolome/raw_2nd.xlsx',sheet = 1)
m2 <- read.xlsx('~/metabolome/raw_2nd.xlsx',sheet = 2)
dim(m1)
colnames(m1)
colnames(m2)
colnames(m1)[12:14] <- c('KEGG','METLIN','HMP')
m1_ids <- m1[,c(1,12,14,13)]
m2_ids <- m2[,c(1,72:74)]
m_ids <- rbind(m1_ids,m2_ids)
head(m_ids)
m_ids_uniq <- m_ids[!duplicated(m_ids$METLIN),]
dim(m_ids)
dim(m_ids_uniq)
# 3790 metlin ID
m1_2nd <- m1[,c(2:11,13)] 
m2_2nd <- m2[,c(2:71,74)]
head(m2_2nd)
dim(m1_2nd)
dim(m2_2nd)
# Define the function to convert less or equal 1 to 0.

replace_1_to_0 <- function(data){
  nC <- ncol(data)
  data1 <- data[,1:(nC-1)]
  data1[data1<=1] <- 0
  data1 <- cbind(data1,data[,nC])
  colnames(data1)[nC] <- 'METLIN'
  return(data1)
}
m1_3rd <- replace_1_to_0(m1_2nd)
head(m1_3rd)
m2_3rd <- replace_1_to_0(m2_2nd)
head(m2_3rd)
m1_4th <- aggregate(m1_3rd[,1:10], list(m1_3rd$METLIN),sum)
m2_4th <- aggregate(m2_3rd[,1:70], list(m2_3rd$METLIN),sum)
colnames(m1_4th)[1] <- 'Metlin'
colnames(m2_4th)[1] <- 'Metlin'
m_all <- merge(m1_4th,m2_4th,by = 'Metlin',all = T)
dim(m_all)
dim(m1_4th)
dim(m2_4th)
m_all[is.na(m_all)] <- 0
# View(m_all)
rownames(m_all) <- m_all[,1]
m_all <- m_all[,-1]
dim(m_all)
# filter metabolites present in less than 30% of them
m_all_metlin <- m_all[rowSums(m_all!=0)>=24,]
dim(m_all)
dim(m_all_metlin)
# 1109 metlin left
length(sort(table(m_ids_uniq[m_ids_uniq$METLIN%in%rownames(m_all_metlin),3]),decreasing = F))
# 692 HMDB
dim(m_all_metlin)

m_1109 <- c()
return_half_min <- function(data){
  nR <- nrow(data)
  nC <- ncol(data)
  temp <- c(rep(0,nR))
  for(k in 1:nR){
    temp[k] <- max(data[k,])
  }
  for(i in 1:nR){
    for(j in 1:nC){
      if(data[i,j]!=0&data[i,j]<temp[i]){
        temp[i] <- data[i,j]
      }
    }
  }
  return(temp/2)
}
m_1109 <- return_half_min(m_all_metlin)
m_1109

for(i in 1:nrow(m_all_metlin)){
  for(j in 1:ncol(m_all_metlin)){
    if(m_all_metlin[i,j]==0){
      m_all_metlin[i,j] <- m_1109[i]
    }
  }
}
min(m_all_metlin)
m_all_metlin_log <- log2(m_all_metlin)
dim(m_all_metlin)
m_all_metlin_log_scale <- scale(t(m_all_metlin_log))

# PCoA
# library(vegan)
# vegdist_met <- vegdist(as.matrix(t(m_all_metlin_log+1)))
# pcoa_met <- cmdscale(vegdist_met, eig = TRUE)
# df_pcoa_met <- as.data.frame(pcoa_met$points)
# colnames(df_pcoa_met) <- c('PCoA1', 'PCoA2')
# PCo1_met <- round(eigenvals(pcoa_met)[1]/sum(eigenvals(pcoa_met)) * 100, 1)
# PCo2_met <- round(eigenvals(pcoa_met)[2]/sum(eigenvals(pcoa_met)) * 100, 1)
# print(PCo1_met)#54.7
# print(PCo2_met)#15.8
# head(df_pcoa_met)
# df_pcoa_met$DPH <- rep(c(1,4,7,14,21,28,35,42),each = 10)
# df_pcoa_met$DPH <- as.factor(df_pcoa_met$DPH)
# library(ggplot2)
# f_pcoa <- ggplot(df_pcoa_met, aes(x = PCoA1, y = PCoA2, color = DPH)) +
#   geom_point() +
#   theme_bw() +
#   xlab('PCoA 1 (54.7%)') +
#   scale_color_manual(values = c("#DEEBF7","#C6DBEF","#9ECAE1","#6BAED6","#4292C6","#2171B5","#08519C","#08306B")) +
#   ylab('PCoA 2 (15.8%)')
# f_pcoa
# f_pcoa1 <- ggplot(df_pcoa_met, aes(x = DPH, y = PCoA1, fill = DPH)) +
#   geom_boxplot() + 
#   theme_bw() +
#   scale_fill_manual(values = c("#DEEBF7","#C6DBEF","#9ECAE1","#6BAED6","#4292C6","#2171B5","#08519C","#08306B")) +
#   theme(panel.grid.minor = element_line(color = NA), legend.position = 'none',
#         axis.text.x = element_text(angle = 90,vjust = 1,hjust = 1))
# f_pcoa1
# f_pcoa2 <- ggplot(df_pcoa_met, aes(x = DPH, y = PCoA2, fill = DPH)) +
#   geom_boxplot() + 
#   theme_bw() +
#   scale_fill_manual(values = c("#DEEBF7","#C6DBEF","#9ECAE1","#6BAED6","#4292C6","#2171B5","#08519C","#08306B")) +
#   theme(panel.grid.minor = element_line(color = NA), legend.position = 'none',
#         axis.text.x = element_text(angle = 90,vjust = 1,hjust = 1))
# f_pcoa2
# library(ggpubr)
# pdf("~/Desktop/01maturity/figures/P6/0pcoa.pdf",width = 7,height = 2)
# ggarrange(f_pcoa,f_pcoa1,f_pcoa2,align = 'h',ncol = 3,widths = c(2,1,1))
# dev.off()
# DPH <- rep(c(1,4,7,14,21,28,35,42),each = 10)
# set.seed(1234)
# adonis(t(m_all_metlin_log)~ DPH)
# 12.1%
# pairwise adonis
# library(pairwiseAdonis)
# data(iris)
# data_for_adonis <- data.frame(t(m_all_metlin_log),DPH)
# dim(data_for_adonis)
# set.seed(1234)
# pairwise.adonis(data_for_adonis[,1:1108],data_for_adonis$DPH)

# head(data_for_adonis)

# if (F) {
#   if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#   
#   BiocManager::install("ropls")
# }
library(ropls)
library(ggplot2)
library(ggsci)
library(Cairo)
library(tidyverse)
m_all_metlin_log_scale
DPH = rep(c(1,4,7,14,21,28,35,42),each = 10)
summary(glm(m_all_metlin_log_scale[,1]~DPH))$coefficients[2,4]
dim(m_all_metlin_log_scale)
sig_glm <- c()
for(i in 1:1109){
  sig_glm[i] <- summary(glm(m_all_metlin_log_scale[,i]~DPH))$coefficients[2,4]
}
SIG_GLM <- data.frame(sig_glm)
rownames(SIG_GLM) <- colnames(m_all_metlin_log_scale)
SIG_GLM$Q <- p.adjust(SIG_GLM$sig_glm,method = 'BH')
sum(SIG_GLM<0.05)

SIG_GLM
m_ids2 <- m_ids[!duplicated(m_ids$METLIN),]
SIG_GLM2 <- merge(SIG_GLM, m_ids2,by.x = 'row.names',by.y = 'METLIN',all.x = T)
head(SIG_GLM2)
SIG_GLM3 <- SIG_GLM2[SIG_GLM2$Q<0.05,]
dim(SIG_GLM3)
# 641 Q < 0.05
# All 1109
SIG_GLM3
dim(SIG_GLM3)

#pls-da
DPH = rep(c(1,4,7,14,21,28,35,42),each = 10)
plsda = opls(t(m_all_metlin_log),DPH)
# variation: plsda@modelDF
# importance: sort(plsda@vipVn)
# points: plsda@scoreMN

# sample scores plot
sample.score = plsda@scoreMN %>% 
  as.data.frame() 
sample.score
sample.score$DPH <- DPH
 
  
col_pal <-  colorRampPalette(c("white", "#3C5488"))(9)[2:9]
f6A <- ggplot(sample.score, aes(p1, p2, color = as.factor(DPH))) +
  geom_point() +
  stat_ellipse(level = 0.95, linetype = 1, 
               size = 0.3, show.legend = FALSE)+
  scale_color_manual('DPH',values = col_pal) +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        legend.position = c(0.2,0.4)) +
  xlab('P1 30.3%') +
  ylab('P2 21.0%')
f6A

f6B <- ggplot(sample.score, aes(as.factor(DPH),p1)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(color = as.factor(DPH))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none') +
  xlab('DPH') +
  ylab('Value of P1') +
  coord_flip() +
  scale_color_manual(values = col_pal)
f6B
library(reshape2)
pva_plsda <- data.frame(melt(pairwise.wilcox.test(sample.score$p1,sample.score$DPH,'BH')$p.value))%>%
  na.omit()
colnames(pva_plsda) <- c("T1",'T2','FDR')
pva_plsda

# VIP scores plot
vip.score = as.data.frame(plsda@vipVn)
colnames(vip.score) = 'vip'
head(vip.score)
vip.score2 <- merge(vip.score,SIG_GLM2,by.x = 'row.names',by.y = 'Row.names',all.x = T) %>%
  select(c('Compound','vip')) %>%
  filter(vip > 1.5)
head(vip.score2)
dim(vip.score2)
vip.score3 <- vip.score2[order(vip.score2$vip)>31,]
dim(vip.score3)
vip.score3$Compound[7] <- 'gamma-Glutamyl-beta-cyanoalanine'
vip.score3$Compound <- factor(vip.score3$Compound,levels = vip.score3$Compound[order(vip.score3$vip,decreasing = F)])
f6C <-  ggplot(vip.score3, aes(Compound, vip)) +
  geom_segment(aes(x = Compound, xend = Compound,
                   y = 0, yend = vip)) +
  geom_point(shape = 21, size = 3, color = '#3C5488' ,fill = '#3C5488') +
  geom_point(aes(1,2.5), color = 'white') +
  #  geom_hline(yintercept = 1, linetype = 'dashed') +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = '', y = 'VIP value') +
  theme_bw() +
  theme(legend.position = 'none',
        legend.text = element_text(color = 'black', family = 'Arial', face = 'plain'),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black',size = 8, family = 'Arial', face = 'plain'),
        axis.text.x = element_text(angle = 90),
        axis.title = element_text(color = 'black',size = 8, family = 'Arial', face = 'plain'),
        axis.ticks = element_line(color = 'black'),
        axis.ticks.x = element_blank()) +
  coord_flip()
f6C

m_all_metlin_log_scale
summary(glm(m_all_metlin_log_scale[,1]~DPH))$coefficients[2,4]
dim(m_all_metlin_log_scale)
sig_glm <- c()
for(i in 1:1109){
  sig_glm[i] <- summary(glm(m_all_metlin_log_scale[,i]~DPH))$coefficients[2,4]
}
SIG_GLM <- data.frame(sig_glm)
rownames(SIG_GLM) <- colnames(m_all_metlin_log_scale)
SIG_GLM$Q <- p.adjust(SIG_GLM$sig_glm,method = 'BH')
sum(SIG_GLM<0.05)

SIG_GLM
m_ids2 <- m_ids[!duplicated(m_ids$METLIN),]
SIG_GLM2 <- merge(SIG_GLM, m_ids2,by.x = 'row.names',by.y = 'METLIN',all.x = T)
head(SIG_GLM2)
SIG_GLM3 <- SIG_GLM2[SIG_GLM2$Q<0.05,]
dim(SIG_GLM3)
# 641 Q < 0.05
# All 1109
SIG_GLM3
dim(SIG_GLM3)
# write.xlsx(SIG_GLM3,'~/metabolome/sig_Q_0.05_metabolites_all.xlsx')

m_all_metlin_log2 <- m_all_metlin_log_scale[,colnames(m_all_metlin_log_scale)%in%SIG_GLM3$Row.names]
dim(m_all_metlin_log2)
Group = rep(c('T1','T2','T3','T4','T5','T6','T7','T8'),each = 10)
m_all_metlin_log3 <- data.frame(m_all_metlin_log2,Group)
head(m_all_metlin_log3)
dim(m_all_metlin_log3)
m_all_metlin_log4 <- aggregate(m_all_metlin_log3[,1:641],list(Group),mean)
head(m_all_metlin_log4)
dim(m_all_metlin_log4)

m_all_metlin_plot <- t(m_all_metlin_log4[,2:642])
colnames(m_all_metlin_plot) <- c(1,4,7,14,21,28,35,42)
library(RColorBrewer)
library(ggplotify)

f6D <- as.ggplot(pheatmap::pheatmap(m_all_metlin_plot,cluster_cols = F,color = colorRampPalette(c('#3C5488','white','#E64B35'))(18),
                             cutree_rows = 6,show_rownames = F))
out <- pheatmap::pheatmap(m_all_metlin_plot,cluster_cols = F,color = colorRampPalette(c('#3C5488','white','#E64B35'))(18),
                cutree_rows = 6,show_rownames = F)
# Six clusters
out1 <- data.frame(sort(cutree(out$tree_row, k = 6)))
dim(out1)
colnames(out1) <- 'Cluster'
out1$Metlin <- gsub('X','',rownames(out1))
head(out1)
out2 <- merge(out1,m_ids_uniq,by.x = 'Metlin',by.y = 'METLIN',all.x = T)
head(out1)
head(out2)
# write.xlsx(out2,'~/metabolome/metabolites_cluster.xlsx')
table(out2$Cluster)
out2$Cluster_new[out2$Cluster==6] <- 'C1'  
out2$Cluster_new[out2$Cluster==4] <- 'C2'
out2$Cluster_new[out2$Cluster==1] <- 'C3'
out2$Cluster_new[out2$Cluster==2] <- 'C4'
out2$Cluster_new[out2$Cluster==3] <- 'C5'
out2$Cluster_new[out2$Cluster==5] <- 'C6'
table(out2$Cluster)
table(out2$Cluster_new)

m_all_metlin_p2 <- m_all_metlin_plot
rownames(m_all_metlin_p2) <- gsub('X','',rownames(m_all_metlin_plot))
head(m_all_metlin_p2)
m_all_metlin_p3 <- merge(m_all_metlin_p2, out2[,c(1,6)], by.x = 'row.names', by.y = 'Metlin',all.x = T )
head(m_all_metlin_p3)
library(dplyr)
library(reshape2)

Met_C1 <- m_all_metlin_p3 %>% filter(Cluster_new == 'C1') %>% select(-c('Cluster_new')) %>% melt(.)
Met_C2 <- m_all_metlin_p3 %>% filter(Cluster_new == 'C2') %>% select(-c('Cluster_new')) %>% melt(.)
Met_C3 <- m_all_metlin_p3 %>% filter(Cluster_new == 'C3') %>% select(-c('Cluster_new')) %>% melt(.)
Met_C4 <- m_all_metlin_p3 %>% filter(Cluster_new == 'C4') %>% select(-c('Cluster_new')) %>% melt(.)
Met_C5 <- m_all_metlin_p3 %>% filter(Cluster_new == 'C5') %>% select(-c('Cluster_new')) %>% melt(.)
Met_C6 <- m_all_metlin_p3 %>% filter(Cluster_new == 'C6') %>% select(-c('Cluster_new')) %>% melt(.)
head(Met_C1)
colnames(Met_C1) <- colnames(Met_C2) <- colnames(Met_C3) <- colnames(Met_C4) <- colnames(Met_C5) <- colnames(Met_C6) <- c('Metlin','DPH','Value')
library(ggplot2)
f6E1 <- ggplot(Met_C1, aes(x = DPH, y = Value, fill = DPH)) +
  geom_violin() +
  theme_bw() +
  scale_fill_manual(values = col_pal) +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),legend.position = 'none') +
  stat_summary(fun=mean, geom="point", shape=16, size=2,col = '#E64B35')+
  ggtitle('Cluster 1')
f6E1
f6E2 <- ggplot(Met_C2, aes(x = DPH, y = Value, fill = DPH)) +
  geom_violin() +
  theme_bw() +
  scale_fill_manual(values = col_pal) +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),legend.position = 'none') +
  stat_summary(fun=mean, geom="point", shape=16, size=2,col = '#E64B35')+
  ggtitle('Cluster 2')
f6E2
f6E3 <- ggplot(Met_C3, aes(x = DPH, y = Value, fill = DPH)) +
  geom_violin() +
  theme_bw() +
  scale_fill_manual(values = col_pal) +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),legend.position = 'none') +
  stat_summary(fun=mean, geom="point", shape=16, size=2,col = '#E64B35')+
  ggtitle('Cluster 3')
f6E3
f6E4 <- ggplot(Met_C4, aes(x = DPH, y = Value, fill = DPH)) +
  geom_violin() +
  theme_bw() +
  scale_fill_manual(values = col_pal) +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),legend.position = 'none') +
  stat_summary(fun=mean, geom="point", shape=16, size=2,col = '#E64B35')+
  ggtitle('Cluster 4')
f6E4
f6E5 <- ggplot(Met_C5, aes(x = DPH, y = Value, fill = DPH)) +
  geom_violin() +
  theme_bw() +
  scale_fill_manual(values = col_pal) +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),legend.position = 'none') +
  stat_summary(fun=mean, geom="point", shape=16, size=2,col = '#E64B35')+
  ggtitle('Cluster 5')
f6E5
f6E6 <- ggplot(Met_C6, aes(x = DPH, y = Value, fill = DPH)) +
  geom_violin() +
  theme_bw() +
  scale_fill_manual(values = col_pal) +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),legend.position = 'none') +
  stat_summary(fun=mean, geom="point", shape=16, size=2,col = '#E64B35') +
  ggtitle('Cluster 6')
f6E6
library(ggpubr)
library(openxlsx)
# write.xlsx(out2,'~/metabolome/clusters_metabolites.xlsx')


# replace the original picture of part F
require(dplyr)
cat_cluster0 <- read.xlsx('clusters_for_all_metabolites_from_six_clusters.xlsx',sheet = 1) 
ids_keep <- names(sort(table(cat_cluster0$Super_class),decreasing = T))[c(1:6)]
head(cat_cluster0)
for(i in 1:nrow(cat_cluster0)){
  if(cat_cluster0$Super_class[i]%in% ids_keep){
    cat_cluster0$SC[i] <- cat_cluster0$Super_class[i]
  }
  else{cat_cluster0$SC[i] <- 'Others'
  }
}
table(cat_cluster0$SC)
head(cat_cluster0)
cat_cluster1 <- table(cat_cluster0[,c(9:10)]) %>%
  melt() %>%
  data.frame()
colnames(cat_cluster1) <- c('Cluster','Super_class','Proportion')
cat_cluster1$Super_class <- factor(cat_cluster1$Super_class, levels = c(names(sort(table(cat_cluster0$Super_class),decreasing = T))[c(1:8)],'Others'))
fF7 <- ggplot(cat_cluster1, aes(x = Cluster, y = Proportion, fill = Super_class)) +
  geom_bar(stat = 'identity',position = 'fill') +
  scale_fill_manual(values = c('#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#DC0000',
                               'gray')) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

SIG_GLM4 <- SIG_GLM3[order(SIG_GLM3$Q,decreasing = F),][1:30,c(1,2,4)]
SIG_GLM4
SIG_GLM5 <- merge(SIG_GLM4, out2, by.x = 'Row.names',by.y = 'Metlin',all.x = T)
dim(SIG_GLM5)
SIG_GLM5$Value <- -log10(SIG_GLM5$sig_glm)
head(SIG_GLM5)
SIG_GLM5 <- SIG_GLM5[,c(1,3,8,9)]
head(SIG_GLM5)
colnames(SIG_GLM5) <- c('Metlin','Compound','Cluster','Value')
SIG_GLM5$Cluster <- gsub('C','Cluster',SIG_GLM5$Cluster)
SIG_GLM5$Compound <- factor(SIG_GLM5$Compound,levels = SIG_GLM5$Compound[order(SIG_GLM5$Value,decreasing = F)])
fG1 <- ggplot(SIG_GLM5,aes(x = Compound, y = Value, fill = Cluster)) +
  geom_bar(stat = 'identity') +
  xlab('Metabolite') +
  ylab('-log10(P)') +
  coord_flip() +
  theme_bw() +
  theme(legend.position = 'left',panel.grid.minor = element_line(color = NA),legend.key.size = unit(0.1,'cm'),
        axis.text.x = element_text(angle = 90,vjust = 1,hjust = 1),axis.text.y = element_text(size = 6))+
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F')) 
fG1
metlin_id <- SIG_GLM5$Metlin[order(SIG_GLM5$Value,decreasing = T)]
m_all_metlin_plot2 <- m_all_metlin_plot[gsub('X','',rownames(m_all_metlin_plot))%in%metlin_id,]
m_all_metlin_plot2
rownames(m_all_metlin_plot2) <- gsub('X','',rownames(m_all_metlin_plot2))
library(reshape2)
m_plot <- m_all_metlin_plot2
for(i in 1:30){
  m_plot[i,] <- rank(m_all_metlin_plot2[i,])
}
m_plot
m_plot2 <- melt(m_plot)
head(m_plot2)
colnames(m_plot2) <- c('Metlin','DPH','Rank')
m_plot2$Metlin <- factor(m_plot2$Metlin,levels = rev(metlin_id))
m_plot2$Rank <- as.factor(m_plot2$Rank)
m_plot2$DPH <- as.factor(m_plot2$DPH)
fG2 <- ggplot(m_plot2,aes(x = DPH,y = Metlin,fill = Rank)) +
  geom_tile() +
  theme_bw() + 
  ylab('') +
  scale_fill_manual(values = colorRampPalette(c("white", "#3C5488"))(11)[2:11]) +
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),legend.key.size = unit(0.1,'cm'),
        axis.text.x = element_text(angle = 90,vjust = 1,hjust = 1)) 
fG2

wilcox.test(m1$value[1:320],m2$value[1:320])
wilcox.test(m1$value[321:640],m2$value[321:640])

