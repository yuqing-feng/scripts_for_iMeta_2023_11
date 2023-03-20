# Multi-omics correlation analyses
# Jun 15th, 2022

# Step 1: Prepare the data of metabolome and amplication sequencing (n = 2)
# Step 2: PA to determine the key communities
# Step 3: mmvec to detect the correlations
# Step 4: mixOmics
# Step 5: Pearson correlation directly

# Step 1: Load the data
library(dplyr)
metab_log <- read.csv('log_transformed_metab.csv',header = T,row.names = 1)
AA_16S <- read.csv('L6_AA_16S.csv',header = T,row.names = 1)
AA_ITS <- read.csv('L6_AA_ITS.csv',header = T,row.names = 1)
Groups <- read.csv('Groups.csv',header = T)
AA_16S_all <- cbind(AA_16S,Groups$Segment)
AA_ITS_all <- cbind(AA_ITS,Groups$Segment)

colnames(AA_16S_all)[533] <- 'Segment'
colnames(AA_ITS_all)[604] <- 'Segment'

AA_16S_C <- AA_16S_all %>%
  filter(Segment=='Cecum') %>%
  select(-c(Segment))
AA_16S_D <- AA_16S_all %>%
  filter(Segment=='Duodenum') %>%
  select(-c(Segment))
AA_16S_I <- AA_16S_all %>%
  filter(Segment=='Ileum') %>%
  select(-c(Segment))
AA_16S_J <- AA_16S_all %>%
  filter(Segment=='Jejunum')%>%
  select(-c(Segment))
AA_ITS_C <- AA_ITS_all %>%
  filter(Segment=='Cecum') %>%
  select(-c(Segment))
AA_ITS_D <- AA_ITS_all %>%
  filter(Segment=='Duodenum') %>%
  select(-c(Segment))
AA_ITS_I <- AA_ITS_all %>%
  filter(Segment=='Ileum') %>%
  select(-c(Segment))
AA_ITS_J <- AA_ITS_all %>%
  filter(Segment=='Jejunum')%>%
  select(-c(Segment))

library(vegan)
dist_16S_C <- vegdist(log2(AA_16S_C+1))
dist_16S_D <- vegdist(log2(AA_16S_D+1))
dist_16S_I <- vegdist(log2(AA_16S_I+1))
dist_16S_J <- vegdist(log2(AA_16S_J+1))
dist_ITS_C <- vegdist(log2(AA_ITS_C+1))
dist_ITS_D <- vegdist(log2(AA_ITS_D+1))
dist_ITS_I <- vegdist(log2(AA_ITS_I+1))
dist_ITS_J <- vegdist(log2(AA_ITS_J+1))

metab0 <- read.csv('log_transformed_metab.csv',header = T,row.names = 1)
metab <- metab0[,c(1,8,9,10,2:7,11:80)]
metab_dist <- vegdist(t(metab+1))

# protest(x,y)$signif -- P values
# protest(x,y)$ss -- Procrustes Sum of Squares
resid_16S_C <- data.frame(residuals(procrustes(dist_16S_C,metab_dist,scale = F)))
median(resid_16S_C[,1]) # 1.13
set.seed(1234)
cor_16S_C <- protest(dist_16S_C,metab_dist)
cor_16S_C$signif # 0.001
cor_16S_C$ss # 0.3913997

resid_16S_D <- data.frame(residuals(procrustes(dist_16S_D,metab_dist,scale = F)))
median(resid_16S_D[,1]) # 0.62
set.seed(1234)
cor_16S_D <- protest(dist_16S_D,metab_dist)
cor_16S_D$signif # 0.001
cor_16S_D$ss # 0.399745

resid_16S_I <- data.frame(residuals(procrustes(dist_16S_I,metab_dist,scale = F)))
median(resid_16S_I[,1]) # 0.96
set.seed(1234)
cor_16S_I <- protest(dist_16S_I,metab_dist)
cor_16S_I$signif # 0.001
cor_16S_I$ss # 0.4913162
resid_16S_J <- data.frame(residuals(procrustes(dist_16S_J,metab_dist,scale = F)))
median(resid_16S_J[,1]) # 0.91
set.seed(1234)
cor_16S_J <- protest(dist_16S_J,metab_dist)
cor_16S_J$signif # 0.001
cor_16S_J$ss # 0.3957452

resid_ITS_C <- data.frame(residuals(procrustes(dist_ITS_C,metab_dist,scale = F)))
median(resid_ITS_C[,1]) # 1.08
set.seed(1234)
cor_ITS_C <- protest(dist_ITS_C,metab_dist)
cor_ITS_C$signif # 0.001
cor_ITS_C$ss # 0.5672115
resid_ITS_D <- data.frame(residuals(procrustes(dist_ITS_D,metab_dist,scale = F)))
median(resid_ITS_D[,1]) # 0.77
set.seed(1234)
cor_ITS_D <- protest(dist_ITS_D,metab_dist)
cor_ITS_D$signif # 0.001
cor_ITS_D$ss # 0.5349073
resid_ITS_I <- data.frame(residuals(procrustes(dist_ITS_I,metab_dist,scale = F)))
median(resid_ITS_I[,1]) # 0.67
set.seed(1234)
cor_ITS_I <- protest(dist_ITS_I,metab_dist)
cor_ITS_I$signif # 0.001
cor_ITS_I$ss # 0.5121084
set.seed(1234)
cor_ITS_J <- protest(dist_ITS_J,metab_dist)
cor_ITS_J$signif # 0.001
cor_ITS_J$ss # 0.5197852
resid_ITS_J <- data.frame(residuals(procrustes(dist_ITS_J,metab_dist,scale = F)))
median(resid_ITS_J[,1]) # 0.86

all_residuals <- data.frame(resid_16S_C,resid_16S_D,resid_16S_I,resid_16S_J,resid_ITS_C,resid_ITS_D,resid_ITS_I,resid_ITS_J)
colnames(all_residuals) <- c('16S_Cecum','16S_Duodenum','16S_Ileum','16S_Jejunum','ITS_Cecum','ITS_Duodenum','ITS_Ileum','ITS_Jejunum')
library(reshape2)
all_res2 <- melt(all_residuals)
colnames(all_res2) <- c('Source','Residual')
all_res2$Sequencing <- rep(c('Bacteria','Fungi'),each = 320)
all_res2$Segment <- c(rep(c('Cecum','Duodenum','Ileum','Jejunum'),each = 80),rep(c('Cecum','Duodenum','Ileum','Jejunum'),each = 80))

library(ggplot2)
all_res2$Segment <- factor(all_res2$Segment,levels = c('Duodenum','Jejunum','Ileum','Cecum'))

f7A <- ggplot(all_res2, aes(x = Segment, y = Residual, color = Segment)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(alpha = 0.5,width = 0.2) +
  facet_grid(.~Sequencing) +
  theme_bw() +
  theme(legend.position = 'none',axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1)) +
  scale_color_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488'))+
  ylim(c(0,6))
f7A

library(ggpubr)

all_res2_16S <- all_res2[all_res2$Sequencing=='Bacteria',]
all_res2_ITS <- all_res2[all_res2$Sequencing=='Fungi',]
pairwise.t.test(all_res2_16S$Residual,all_res2_16S$Segment,p.adjust.method = 'BH')
pva_16S <- data.frame(melt(pairwise.wilcox.test(all_res2_16S$Residual,all_res2_16S$Segment,p.adjust.method ='BH')$p.value))%>%
  na.omit()
pva_16S
pva_ITS <- data.frame(melt(pairwise.wilcox.test(all_res2_ITS$Residual,all_res2_ITS$Segment,p.adjust.method ='BH')$p.value))%>%
  na.omit()
pva_ITS

# mmvec
# step biplot
library(openxlsx)
library(dplyr)
d_16S_C <- read.xlsx('~/C_16S/beta.xlsx',sheet = 1) %>%
  select(-c('PC3'))
d_16S_D <- read.xlsx('~/D_16S/beta.xlsx',sheet = 1) %>%
  select(-c('PC3'))
d_16S_I <- read.xlsx('~/I_16S/beta.xlsx',sheet = 1) %>%
  select(-c('PC3'))
d_16S_J <- read.xlsx('~/J_16S/beta.xlsx',sheet = 1) %>%
  select(-c('PC3'))
d_ITS_C <- read.xlsx('~/C_ITS/beta.xlsx',sheet = 1) %>%
  select(-c('PC3'))
d_ITS_D <- read.xlsx('~/D_ITS/beta.xlsx',sheet = 1) %>%
  select(-c('PC3'))
d_ITS_I <- read.xlsx('~/I_ITS/beta.xlsx',sheet = 1) %>%
  select(-c('PC3'))
d_ITS_J <- read.xlsx('~/J_ITS/beta.xlsx',sheet = 1) %>%
  select(-c('PC3'))
d_16S_C
bac_id <- read.xlsx('taxa1.xlsx',sheet = 1)%>%
  select(c('ID','Genus','Phylum'))
fun_id <- read.xlsx('taxa2.xlsx', sheet = 2) %>%
  select(c('ID','Genus','Phylum'))
head(bac_id)
head(fun_id)
d_16S_C2 <- merge(d_16S_C, bac_id, by.x = 'Taxa', by.y = 'ID',all.x = T)
d_16S_D2 <- merge(d_16S_D, bac_id, by.x = 'Taxa', by.y = 'ID',all.x = T)
d_16S_I2 <- merge(d_16S_I, bac_id, by.x = 'Taxa', by.y = 'ID',all.x = T)
d_16S_J2 <- merge(d_16S_J, bac_id, by.x = 'Taxa', by.y = 'ID',all.x = T)
# write.xlsx(d_16S_C2,'~/out_reorder/C_16S/taxa.xlsx', sheet = 1)
# write.xlsx(d_16S_D2,'~/out_reorder/D_16S/taxa.xlsx', sheet = 1)
# write.xlsx(d_16S_I2,'~/out_reorder/I_16S/taxa.xlsx', sheet = 1)
# write.xlsx(d_16S_J2,'~/out_reorder/J_16S/taxa.xlsx', sheet = 1)
d_ITS_C2 <- merge(d_ITS_C, fun_id, by.x = 'Taxa', by.y = 'ID',all.x = T)
d_ITS_D2 <- merge(d_ITS_D, fun_id, by.x = 'Taxa', by.y = 'ID',all.x = T)
d_ITS_I2 <- merge(d_ITS_I, fun_id, by.x = 'Taxa', by.y = 'ID',all.x = T)
d_ITS_J2 <- merge(d_ITS_J, fun_id, by.x = 'Taxa', by.y = 'ID',all.x = T)
# write.xlsx(d_ITS_C2, '~/out_reorder/C_ITS/taxa.xlsx', sheet = 1)
# write.xlsx(d_ITS_D2, '~/out_reorder/D_ITS/taxa.xlsx', sheet = 1)
# write.xlsx(d_ITS_I2, '~/out_reorder/I_ITS/taxa.xlsx', sheet = 1)
# write.xlsx(d_ITS_J2, '~/out_reorder/J_ITS/taxa.xlsx', sheet = 1)

metlin_cluster <- read.xlsx('clusters_metabolites.xlsx',sheet = 1) %>%
  select(c('Metlin','Cluster_new'))
head(metlin_cluster)

m_16S_C <- read.xlsx('~/out_reorder/C_16S/beta.xlsx',sheet = 2) %>%
  select(-c('PC3'))
m_16S_D <- read.xlsx('~/out_reorder/D_16S/beta.xlsx',sheet = 2) %>%
  select(-c('PC3'))
m_16S_I <- read.xlsx('~/out_reorder/I_16S/beta.xlsx', sheet = 2) %>%
  select(-c('PC3'))
m_16S_J <- read.xlsx('~/out_reorder/J_16S/beta.xlsx', sheet = 2) %>%
  select(-c('PC3'))
m_ITS_C <- read.xlsx('~/out_reorder/C_ITS/beta.xlsx',sheet = 2) %>%
  select(-c('PC3'))
m_ITS_D <- read.xlsx('~/out_reorder/D_ITS/beta.xlsx',sheet = 2) %>%
  select(-c('PC3'))
m_ITS_I <- read.xlsx('~/out_reorder/I_ITS/beta.xlsx', sheet = 2) %>%
  select(-c('PC3'))
m_ITS_J <- read.xlsx('~/out_reorder/J_ITS/beta.xlsx', sheet = 2) %>%
  select(-c('PC3'))

m_16S_C2 <- merge(m_16S_C, metlin_cluster, by.x = 'Metabolite', by.y = 'Metlin',all.x = T)
m_16S_D2 <- merge(m_16S_D, metlin_cluster, by.x = 'Metabolite', by.y = 'Metlin',all.x = T)
m_16S_I2 <- merge(m_16S_I, metlin_cluster, by.x = 'Metabolite', by.y = 'Metlin',all.x = T)
m_16S_J2 <- merge(m_16S_J, metlin_cluster, by.x = 'Metabolite', by.y = 'Metlin',all.x = T)
m_ITS_C2 <- merge(m_ITS_C, metlin_cluster, by.x = 'Metabolite', by.y = 'Metlin',all.x = T)
m_ITS_D2 <- merge(m_ITS_D, metlin_cluster, by.x = 'Metabolite', by.y = 'Metlin',all.x = T)
m_ITS_I2 <- merge(m_ITS_I, metlin_cluster, by.x = 'Metabolite', by.y = 'Metlin',all.x = T)
m_ITS_J2 <- merge(m_ITS_J, metlin_cluster, by.x = 'Metabolite', by.y = 'Metlin',all.x = T)

# write.xlsx(m_16S_C2,'~/out_reorder/C_16S/metab.xlsx',sheet = 1)
# write.xlsx(m_16S_D2,'~/out_reorder/D_16S/metab.xlsx',sheet = 1)
# write.xlsx(m_16S_I2,'~/out_reorder/I_16S/metab.xlsx',sheet = 1)
# write.xlsx(m_16S_J2,'~/out_reorder/J_16S/metab.xlsx',sheet = 1)
# write.xlsx(m_ITS_C2,'~/out_reorder/C_ITS/metab.xlsx',sheet = 1)
# write.xlsx(m_ITS_D2,'~/out_reorder/D_ITS/metab.xlsx',sheet = 1)
# write.xlsx(m_ITS_I2,'~/out_reorder/I_ITS/metab.xlsx',sheet = 1)
# write.xlsx(m_ITS_J2,'~/out_reorder/J_ITS/metab.xlsx',sheet = 1)


format_taxa <- function(d_16S_C2){
  d_16S_C2$Length <- sqrt((d_16S_C2$PC1)^2+ (d_16S_C2$PC2)^2)
  d_16S_C2$Rank <- rank(d_16S_C2$Length)
  temp1 <- nrow(d_16S_C2)
  temp1
  temp2 <- temp1 - 5
  for(i in 1:nrow(d_16S_C2)){
    if(d_16S_C2$Rank[i]>temp2){d_16S_C2$Genus2[i] <- d_16S_C2$Genus[i]}
    else(d_16S_C2$Genus2[i] <- '')
  }
  return(d_16S_C2)
}

d_16S_C2 <- format_taxa(d_16S_C2)
d_16S_D2 <- format_taxa(d_16S_D2)
d_16S_I2 <- format_taxa(d_16S_I2)
d_16S_J2 <- format_taxa(d_16S_J2)
d_ITS_C2 <- format_taxa(d_ITS_C2)
d_ITS_D2 <- format_taxa(d_ITS_D2)
d_ITS_I2 <- format_taxa(d_ITS_I2)
d_ITS_J2 <- format_taxa(d_ITS_J2)


p_16S_D <- ggplot(m_16S_D2,aes(x = PC1,y = PC2)) + geom_point(aes(color = Cluster_new),size = 0.5,data = m_16S_D2) + #1
  geom_segment(aes(x = 0, y = 0, xend = PC1/200, yend = PC2/200,col = Genus2), 
               arrow = arrow(length = unit(0.2,'cm')),data = d_16S_D2)+ # 1
  geom_text(data = d_16S_D2,aes(x = PC1/200,y = PC2/200,label = Genus2)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks.length = unit(0,'cm'),legend.position = 'none') +
  scale_color_manual(values = c('gray','#4DBBD5','#4DBBD5','#4DBBD5','#4DBBD5','#4DBBD5','#4DBBD5',
                                '#E64B35','#E64B35','#E64B35','#E64B35','#E64B35')) +
  ggtitle('Duodenum (Bacteria)')
p_16S_D

p_16S_J <- ggplot(m_16S_J2,aes(x = PC1,y = PC2)) + geom_point(aes(color = Cluster_new),size = 0.5,data = m_16S_J2) + #3
  geom_segment(aes(x = 0, y = 0, xend = PC1/200, yend = PC2/200,col = Genus2), 
               arrow = arrow(length = unit(0.2,'cm')),data = d_16S_J2)+ # 1
  geom_text(data = d_16S_J2,aes(x = PC1/200,y = PC2/200,label = Genus2)) + # 1
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks.length = unit(0,'cm'),legend.position = 'none') +
  scale_color_manual(values = c('gray','#E64B35','#E64B35','#4DBBD5','#4DBBD5','#4DBBD5','#4DBBD5','#4DBBD5','#4DBBD5',
                                '#E64B35','#E64B35','#E64B35')) +
  ggtitle('Jejunum (Bacteria)') # 1
p_16S_J

p_16S_I <- ggplot(m_16S_I2,aes(x = PC1,y = PC2)) + geom_point(aes(color = Cluster_new),size = 0.5,data = m_16S_I2) + #3
  geom_segment(aes(x = 0, y = 0, xend = PC1/200, yend = PC2/200,col = Genus2), 
               arrow = arrow(length = unit(0.2,'cm')),data = d_16S_I2)+ # 1
  geom_text(data = d_16S_I2,aes(x = PC1/200,y = PC2/200,label = Genus2)) + # 1
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks.length = unit(0,'cm'),legend.position = 'none') +
  scale_color_manual(values = c('gray','#4DBBD5','#4DBBD5','#4DBBD5','#4DBBD5','#4DBBD5','#4DBBD5',
                                '#E64B35','#E64B35','#E64B35','#E64B35','#E64B35')) +
  ggtitle('Ileum (Bacteria)') # 1
p_16S_I

p_16S_C <- ggplot(m_16S_C2,aes(x = PC1,y = PC2)) + geom_point(aes(color = Cluster_new),size = 0.5,data = m_16S_C2) + #3
  geom_segment(aes(x = 0, y = 0, xend = PC1/200, yend = PC2/200,col = Genus2), 
               arrow = arrow(length = unit(0.2,'cm')),data = d_16S_C2)+ # 1
  geom_text(data = d_16S_C2,aes(x = PC1/200,y = PC2/200,label = Genus2)) + # 1
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks.length = unit(0,'cm'),legend.position = 'none') +
  scale_color_manual(values = c('gray','#E64B35','#E64B35','#E64B35','#4DBBD5','#4DBBD5','#4DBBD5','#4DBBD5','#4DBBD5','#4DBBD5',
                                '#E64B35','#E64B35')) +
  ggtitle('Cecum (Bacteria)') # 1
p_16S_C


p_ITS_D <- ggplot(m_ITS_D2,aes(x = PC1,y = PC2)) + geom_point(aes(color = Cluster_new),size = 0.5,data = m_ITS_D2) + #1
  geom_segment(aes(x = 0, y = 0, xend = PC1/200, yend = PC2/200,col = Genus2), 
               arrow = arrow(length = unit(0.2,'cm')),data = d_ITS_D2)+ # 1
  geom_text(data = d_ITS_D2,aes(x = PC1/200,y = PC2/200,label = Genus2)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks.length = unit(0,'cm'),legend.position = 'none') +
  scale_color_manual(values = c('gray','#4DBBD5','#4DBBD5','#4DBBD5','#4DBBD5','#4DBBD5','#4DBBD5',
                                '#E64B35','#E64B35','#E64B35','#E64B35','#E64B35')) +
  ggtitle('Duodenum (Fungi)')
p_ITS_D

p_ITS_J <- ggplot(m_ITS_J2,aes(x = PC1,y = PC2)) + geom_point(aes(color = Cluster_new),size = 0.5,data = m_ITS_J2) + #3
  geom_segment(aes(x = 0, y = 0, xend = PC1/200, yend = PC2/200,col = Genus2), 
               arrow = arrow(length = unit(0.2,'cm')),data = d_ITS_J2)+ # 1
  geom_text(data = d_ITS_J2,aes(x = PC1/200,y = PC2/200,label = Genus2)) + # 1
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks.length = unit(0,'cm'),legend.position = 'none') +
  scale_color_manual(values = c('gray','#4DBBD5','#4DBBD5','#4DBBD5','#4DBBD5','#4DBBD5','#4DBBD5',
                                '#E64B35','#E64B35','#E64B35','#E64B35','#E64B35')) +
  ggtitle('Jejunum (Fungi)') # 1
p_ITS_J

p_ITS_I <- ggplot(m_ITS_I2,aes(x = PC1,y = PC2)) + geom_point(aes(color = Cluster_new),size = 0.5,data = m_ITS_I2) + #3
  geom_segment(aes(x = 0, y = 0, xend = PC1/200, yend = PC2/200,col = Genus2), 
               arrow = arrow(length = unit(0.2,'cm')),data = d_ITS_I2)+ # 1
  geom_text(data = d_ITS_I2,aes(x = PC1/200,y = PC2/200,label = Genus2)) + # 1
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks.length = unit(0,'cm'),legend.position = 'none') +
  scale_color_manual(values = c('gray','#E64B35','#4DBBD5','#4DBBD5','#4DBBD5','#4DBBD5','#4DBBD5','#4DBBD5',
                                '#E64B35','#E64B35','#E64B35','#E64B35')) +
  ggtitle('Ileum (Fungi)') # 1
p_ITS_I

p_ITS_C <- ggplot(m_ITS_C2,aes(x = PC1,y = PC2)) + geom_point(aes(color = Cluster_new),size = 0.5,data = m_ITS_C2) + #3
  geom_segment(aes(x = 0, y = 0, xend = PC1/200, yend = PC2/200,col = Genus2), 
               arrow = arrow(length = unit(0.2,'cm')),data = d_ITS_C2)+ # 1
  geom_text(data = d_ITS_C2,aes(x = PC1/200,y = PC2/200,label = Genus2)) + # 1
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks.length = unit(0,'cm'),legend.position = 'none') +
  scale_color_manual(values = c('gray','#E64B35','#4DBBD5','#4DBBD5','#4DBBD5','#4DBBD5','#4DBBD5','#4DBBD5',
                                '#E64B35','#E64B35','#E64B35','#E64B35')) +
  ggtitle('Cecum (Fungi)')  
p_ITS_C
library(patchwork)
library(ComplexHeatmap)
library(reshape2)
library(ggplot2)
library(tibble)
library(dplyr)
library(circlize)
h_16S_D <- read.table('~/out_reorder/D_16S/ranks.txt',header = T,row.names = 1) %>%
  rownames_to_column() %>%
  melt() %>%
  # slice_max(value, n = 100) %>%
  filter(value > quantile(value, 0.99)) %>%
  mutate(Annotation = 'Duodenum_16S') %>%
  rename('Metlin' = 'rowname') %>%
  rename('Taxa' = 'variable') 
h_16S_J <- read.table('~/out_reorder/J_16S/ranks.txt',header = T,row.names = 1) %>%
  rownames_to_column() %>%
  melt() %>%
  # slice_max(value, n = 100) %>%
  filter(value > quantile(value, 0.99)) %>%
  mutate(Annotation = 'Jejunum_16S')%>%
  rename('Metlin' = 'rowname') %>%
  rename('Taxa' = 'variable') 
h_16S_I <- read.table('~/out_reorder/I_16S/ranks.txt',header = T, row.names = 1) %>%
  rownames_to_column() %>%
  melt() %>%
  # slice_max(value, n = 100) %>%
  filter(value > quantile(value, 0.99)) %>%
  mutate(Annotation = 'Ileum_16S')%>%
  rename('Metlin' = 'rowname') %>%
  rename('Taxa' = 'variable') 
h_16S_C <- read.table('~/out_reorder/C_16S/ranks.txt',header = T, row.names = 1) %>%
  rownames_to_column() %>%
  melt() %>%
  # slice_max(value, n = 100) %>%
  filter(value > quantile(value, 0.99)) %>%
  mutate(Annotation = 'Cecum_16S') %>%
  rename('Metlin' = 'rowname') %>%
  rename('Taxa' = 'variable')
h_ITS_D <- read.table('~/out_reorder/D_ITS/ranks.txt',header = T,row.names = 1) %>%
  rownames_to_column() %>%
  melt() %>%
  # slice_max(value, n = 100) %>%
  filter(value > quantile(value, 0.99)) %>%
  mutate(Annotation = 'Duodenum_ITS') %>%
  rename('Metlin' = 'rowname') %>%
  rename('Taxa' = 'variable') 
h_ITS_J <- read.table('~/out_reorder/J_ITS/ranks.txt',header = T,row.names = 1) %>%
  rownames_to_column() %>%
  melt() %>%
  # slice_max(value, n = 100) %>%
  filter(value > quantile(value, 0.99)) %>%
  mutate(Annotation = 'Jejunum_ITS')%>%
  rename('Metlin' = 'rowname') %>%
  rename('Taxa' = 'variable') 
h_ITS_I <- read.table('~/out_reorder/I_ITS/ranks.txt',header = T, row.names = 1) %>%
  rownames_to_column() %>%
  melt() %>%
  # slice_max(value, n = 100) %>%
  filter(value > quantile(value, 0.99)) %>%
  mutate(Annotation = 'Ileum_ITS')%>%
  rename('Metlin' = 'rowname') %>%
  rename('Taxa' = 'variable') 
h_ITS_C <- read.table('~/out_reorder/C_ITS/ranks.txt',header = T, row.names = 1) %>%
  rownames_to_column() %>%
  melt() %>%
  # slice_max(value, n = 100) %>%
  filter(value > quantile(value, 0.99)) %>%
  mutate(Annotation = 'Cecum_ITS') %>%
  rename('Metlin' = 'rowname') %>%
  rename('Taxa' = 'variable')
all_metab_taxa <- rbind(h_16S_D,h_16S_J,h_16S_I,h_16S_C,
                        h_ITS_D,h_ITS_J,h_ITS_I,h_ITS_C)
le_metlin <- unique(all_metab_taxa$Metlin)
length(le_metlin)
all_metab_taxa2 <- all_metab_taxa
all_metab_taxa2$value <- 1
all_metab_taxa2$Taxa2 <- paste(all_metab_taxa2$Annotation,all_metab_taxa2$Taxa,sep = '-')
head(all_metab_taxa2)
all_metab_taxa2

all_metab_taxa3 <- dcast(all_metab_taxa2[,c(1,3,5)],Metlin~Taxa2) 
all_metab_taxa3[is.na(all_metab_taxa3)] <- 0
rownames(all_metab_taxa3) <- all_metab_taxa3$Metlin
all_metab_taxa3 <- all_metab_taxa3[,-1]
head(all_metab_taxa3)

taxa_seg <- strsplit(colnames(all_metab_taxa3),'-')
taxa_seg
Taxa_Seg <- c()
for(i in 1:length(taxa_seg)){
  Taxa_Seg[i] <- taxa_seg[[i]][1]
}
Taxa_Seg

data_for_order <- data.frame(colnames(all_metab_taxa3),colSums(all_metab_taxa3))
data_for_order
colnames(data_for_order) <- c('Genus','Number')
data_for_order$original <- 1:nrow(data_for_order)
for(i in 1:nrow(data_for_order)){
  data_for_order$cat1[i] <- strsplit(data_for_order$Genus,'-')[[i]][1]
}
data_for_order$cat1[data_for_order$cat1=='Duodenum_16S'] <- 'A'
data_for_order$cat1[data_for_order$cat1=='Duodenum_ITS'] <- 'B'
data_for_order$cat1[data_for_order$cat1=='Jejunum_16S'] <- 'C'
data_for_order$cat1[data_for_order$cat1=='Jejunum_ITS'] <- 'D'
data_for_order$cat1[data_for_order$cat1=='Ileum_16S'] <- 'E'
data_for_order$cat1[data_for_order$cat1=='Ileum_ITS'] <- 'F'
data_for_order$cat1[data_for_order$cat1=='Cecum_16S'] <- 'G'
data_for_order$cat1[data_for_order$cat1=='Cecum_ITS'] <- 'H'
head(data_for_order)
data_for_order$cat1 <- factor(data_for_order$cat1,levels = c('A','B','C','D','E','F','G','H'),ordered = TRUE)
data_for_order2 <- data_for_order[order(data_for_order[,4], -data_for_order[,2]), ]
all_metab_taxa4 <- all_metab_taxa3[,data_for_order2$Genus]

data_for_order3 <- data.frame(rownames(all_metab_taxa3),rowSums(all_metab_taxa3))
colnames(data_for_order3) <- c('Metabolite','Number')
head(data_for_order3)
metlin_to_cluster <- openxlsx::read.xlsx('clusters_metabolites.xlsx',sheet = 1) %>%
  select(c('Metlin','Cluster_new'))
head(metlin_to_cluster)
data_for_order3_2 <- merge(data_for_order3,metlin_to_cluster, by.x = 'Metabolite',by.y = 'Metlin', all.x = T)
head(data_for_order3_2)
data_for_order4 <- data_for_order3_2[order(data_for_order3_2$Cluster_new ,-data_for_order3_2[,2]),]
all_metab_taxa5 <- all_metab_taxa4[data_for_order4$Metabolite,]
Source1 <- c()
for(i in 1:nrow(data_for_order2)){
  Source1[i] <- strsplit(data_for_order2$Genus,'-')[[i]][1]
}
Source1
Number_Taxa <- colSums(all_metab_taxa5)
Number_Taxa
col_fun1 <- colorRamp2(c(1,72),c('white','#7E6148FF'))
col_fun2 <- colorRamp2(c(1,25),c('white','#8491B4FF'))
column_ha <- HeatmapAnnotation(Source = Source1,Number1 = Number_Taxa,
                               col = list(Number1 = col_fun1,
                                          Source = c('Duodenum_16S'='#E64B35','Jejunum_16S'='#E64B35','Ileum_16S'='#E64B35','Cecum_16S'='#E64B35',
                                                     'Duodenum_ITS'='#4DBBD5','Jejunum_ITS'='#4DBBD5','Ileum_ITS'='#4DBBD5','Cecum_ITS'='#4DBBD5')))
column_ha
Number_metab <- rowSums(all_metab_taxa5)
rowname_ha <- HeatmapAnnotation( Number2 = Number_metab, Cluster = data_for_order4$Cluster_new,which = 'row',
                                 col = list(Cluster = c('C1'='#00A087','C2'='#3C5488','C3'='#F39B7F','C4'='#8491B4','C5'='#91D1C2','C6'='#B09C85'),Number = col_fun2))
pd1 <- Heatmap(all_metab_taxa5,name = 'Co-occurrence', top_annotation = column_ha,
               show_row_names = FALSE,show_column_names = FALSE,
               col = c('white','#E64B35'),
               cluster_columns = FALSE,
               cluster_rows = FALSE,
               border = TRUE,
               column_title = 'Genus',
               column_title_side = 'bottom',
               row_title = 'Metabolite',
               show_row_dend = FALSE, column_dend_reorder = unique(all_metab_taxa2$Taxa2),
               show_heatmap_legend = FALSE,
               right_annotation = rowname_ha) 
pd1

# Fisher exast test

all_metab_taxa
all_metab_taxa_clu <- merge(all_metab_taxa, metlin_to_cluster, by = 'Metlin', all.x = T)
data_for_exact <- data.frame(table(all_metab_taxa_clu[,4:5]))
data_for_exact
merge1 <- data.frame(table(metlin_to_cluster$Cluster_new)) %>%
  rename('Cluster' = 'Var1') %>%
  rename('Freq_me' = 'Freq')
merge2 <- data.frame(Source = c('Duodenum_16S','Jejunum_16S','Ileum_16S','Cecum_16S',
                                'Duodenum_ITS','Jejunum_ITS','Ileum_ITS','Cecum_ITS'),
                     Number = c(65,45,20,93,46,65,56,91)) %>%
  rename('Taxa' = 'Number')
head(data_for_exact)
head(merge1)
head(merge2)
data_for_exact2 <- merge(data_for_exact,merge1,by.x = 'Cluster_new',by.y = 'Cluster',all.x = T)
data_for_exact3 <- merge(data_for_exact2,merge2,by.x = 'Annotation',by.y = 'Source',all.x = T)
data_for_exact3$Ratio <- data_for_exact3$Freq/(data_for_exact3$Freq_me*data_for_exact3$Taxa) 
head(data_for_exact3)
data_for_exact3$Annotation <- factor(data_for_exact3$Annotation,levels = rev(c('Cecum_ITS','Cecum_16S','Ileum_ITS','Ileum_16S',
                                                                               'Jejunum_ITS','Jejunum_16S','Duodenum_ITS','Duodenum_16S')))
head(data_for_exact3)

f7C <- ggplot(data_for_exact3,aes(y = Cluster_new,x = Annotation,fill = Ratio * 100)) +
  geom_tile() +
  scale_fill_gradient2(low = '#3C5488',mid =  'white',high = '#E64B35',midpoint = 1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  xlab('Source') +
  ylab('Cluster') +
  geom_text(data = data_for_exact3, aes(y = Cluster_new,x = Annotation, label = round(Ratio*100,2))) +
  scale_y_discrete(labels = c(1,2,3,4,5,6)) +
  scale_x_discrete(labels = rev(c("Cecum (F)","Cecum (B)",
                                  "Ileum (F)","Ileum (B)",
                                  "Jejunum (F)",'Jejunum (B)',
                                  'Duodenum (F)','Duodenum (B)')))

f7C
library(networkD3)
h_D <- rbind(h_16S_D[,c(1,2)],h_ITS_D[,c(1,2)])
h_J <- rbind(h_16S_J[,c(1,2)],h_ITS_J[,c(1,2)])
h_I <- rbind(h_16S_I[,c(1,2)],h_ITS_I[,c(1,2)])
h_C <- rbind(h_16S_C[,c(1,2)],h_ITS_C[,c(1,2)])
metlin_cluster
taxa_id <- rbind(bac_id,fun_id)
convert_id_to_cat <- function(corr,metlin_cluster,taxa_id){
  corr1 <- merge(corr, metlin_cluster, by.x = 'Metlin', by.y = 'Metlin',all.x = T)
  corr2 <- merge(corr1, taxa_id[,c(1,3)], by.x = 'Taxa', by.y = 'ID',all.x = T)
  corr3 <- corr2[,c(3,4)]
  colnames(corr3) <- c('Target','Source')
  corr3$Count <- paste(corr3$Target,corr3$Source,sep = '-')
  require(dplyr)
  corr4 <- count(corr3, Count)
  require(stringr)
  corr5 <- data.frame(str_split_fixed(corr4$Count, '-', 2))
  corr5$number <- corr4$n
  colnames(corr5) <- c('Target','Source','Number')
  return(corr5)
}
hD2 <- convert_id_to_cat(h_D, metlin_cluster, taxa_id)
hJ2 <- convert_id_to_cat(h_J, metlin_cluster, taxa_id)
hI2 <- convert_id_to_cat(h_I, metlin_cluster, taxa_id)
hC2 <- convert_id_to_cat(h_C, metlin_cluster, taxa_id)
node_D <- data.frame(name = c(as.character(hD2$Source),
                              as.character(hD2$Target)) %>% unique())
node_J <- data.frame(name = c(as.character(hJ2$Source),
                              as.character(hJ2$Target)) %>% unique())
node_I <- data.frame(name = c(as.character(hI2$Source),
                              as.character(hI2$Target)) %>% unique())
node_C <- data.frame(name = c(as.character(hC2$Source),
                              as.character(hC2$Target)) %>% unique())
hD2$IDsource <- match(hD2$Source, node_D$name)-1 
hD2$IDtarget <- match(hD2$Target, node_D$name)-1
hJ2$IDsource <- match(hJ2$Source, node_J$name)-1 
hJ2$IDtarget <- match(hJ2$Target, node_J$name)-1
hI2$IDsource <- match(hI2$Source, node_I$name)-1 
hI2$IDtarget <- match(hI2$Target, node_I$name)-1
hC2$IDsource <- match(hC2$Source, node_C$name)-1 
hC2$IDtarget <- match(hC2$Target, node_C$name)-1
color_D <- 'd3.scaleOrdinal().domain(["Firmicutes", "Bacteroidetes","Proteobacteria",
"Ascomycota", "Basidiomycota", "Actinobacteria", "Verrucomicrobia", "Mortierellomycota",
"Epsilonbacteraeota","C1","C2","C3","C4","C5","C6"]).
range(["#8491B4", "#00A087" , "#91D1C2", "#4DBBD5", "#F39B7F", "#B09C85", "#B09C85", "#B09C85",
"#B09C85","#3C5488","#3C5488","#3C5488","#3C5488","#3C5488","#3C5488"])'
color_J <- 'd3.scaleOrdinal().domain(["Firmicutes", "Bacteroidetes","Proteobacteria",
"Ascomycota", "Basidiomycota", "Actinobacteria", "Verrucomicrobia", "Mortierellomycota",
"Epsilonbacteraeota","Mucoromycota","Glomeromycota","C1","C2","C3","C4","C5","C6"]).
range(["#8491B4", "#00A087" , "#91D1C2", "#4DBBD5", "#F39B7F", "#B09C85", "#B09C85", "#B09C85",
"#B09C85","#B09C85","#B09C85","#3C5488","#3C5488","#3C5488","#3C5488","#3C5488","#3C5488"])'
color_I <- 'd3.scaleOrdinal().domain(["Firmicutes", "Bacteroidetes","Proteobacteria",
"Ascomycota", "Basidiomycota", "Actinobacteria", "Verrucomicrobia", "Mortierellomycota",
"Epsilonbacteraeota","Mucoromycota","Glomeromycota","C1","C2","C3","C4","C5","C6"]).
range(["#8491B4", "#00A087" , "#91D1C2", "#4DBBD5", "#F39B7F", "#B09C85", "#B09C85", "#B09C85",
"#B09C85","#B09C85","#B09C85","#3C5488","#3C5488","#3C5488","#3C5488","#3C5488","#3C5488"])'
color_C <- 'd3.scaleOrdinal().domain(["Firmicutes", "Bacteroidetes","Proteobacteria",
"Ascomycota", "Basidiomycota", "Actinobacteria", "Verrucomicrobia", "Mortierellomycota",
"Epsilonbacteraeota","Mucoromycota","Glomeromycota","C1","C2","C3","C4","C5","C6"]).
range(["#8491B4", "#00A087" , "#91D1C2", "#4DBBD5", "#F39B7F", "#B09C85", "#B09C85", "#B09C85",
"#B09C85","#B09C85","#B09C85","#3C5488","#3C5488","#3C5488","#3C5488","#3C5488","#3C5488"])'
P_C1 <- sankeyNetwork(Links = hD2, Nodes = node_D,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "Number", NodeID = "name", colourScale = color_D,
                      sinksRight=FALSE)
P_C1
P_C2 <- sankeyNetwork(Links = hJ2, Nodes = node_J,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "Number", NodeID = "name",  colourScale = color_J,
                      sinksRight=FALSE)
P_C2
P_C3 <- sankeyNetwork(Links = hI2, Nodes = node_I,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "Number", NodeID = "name", colourScale = color_I,
                      sinksRight=FALSE)
P_C3
P_C4 <- sankeyNetwork(Links = hC2, Nodes = node_C,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "Number", NodeID = "name", colourScale = color_C,
                      sinksRight=FALSE)

P_C1
P_C2
P_C3 
P_C4
# saveNetwork(P_C1,'~/Desktop/sankey_D.html')
# saveNetwork(P_C2,'~/Desktop/sankey_J.html')
# saveNetwork(P_C3,'~/Desktop/sankey_I.html')
# saveNetwork(P_C4,'~/Desktop/sankey_C.html')

library(data.table)
h_4S <- rbind(h_D, h_J, h_I, h_C)
h_4S_1 <- setDT(h_4S)[,list(Count = .N), names(h_4S)]
dim(h_4S)
dim(h_4S_1)
head(h_4S_1)
head(taxa_id)
head(metlin_2_hmdb)
h_4S_2 <- merge(h_4S_1, taxa_id, by.x = 'Taxa', by.y = 'ID', all.x = T)
h_4S_3 <- merge(h_4S_2, metlin_2_hmdb[,c(1,2,5)], by.x = 'Metlin', by.y = 'Metlin', all.x = T)
ID_network1 <- h_4S_3[,c(1,6,7)]
dim(ID_network1)
ID_network1 <- ID_network1[!duplicated(ID_network1),]
dim(ID_network1)
colnames(ID_network1) <- c('ID','Lable','Group')
head(ID_network1)
ID_network2 <- h_4S_3[,c(2,4,5)]
dim(ID_network2)
ID_network2 <- ID_network2[!duplicated(ID_network2),]
colnames(ID_network2) <- c('ID','Lable','Group')
dim(ID_network2)
head(ID_network2)
ID_network <- rbind(ID_network1, ID_network2)
dim(ID_network)
require(openxlsx)
write.xlsx(ID_network,'ID_4_network.xlsx')
write.xlsx(h_4S_1,'edges.xlsx')

hm_metori_p1 <- read.xlsx('hmdb_from_metori_and_source.xlsx', sheet = 1)%>%
  filter(Source %in% c('Co-Metabolism','Host','Microbiota'))
hm_metori_p1_met <- merge(hm_metori_p1[,c(1,3)], metlin_2_hmdb[,c(1,2,4)], by.x = 'HMDB.ID',by.y = 'HMP',all.x = T)
h_p1_sub <- h_4S_1[h_4S_1$Metlin%in%hm_metori_p1_met$Metlin,]
h_p1_sub_1 <- merge(h_p1_sub, taxa_id, by.x = 'Taxa', by.y = 'ID', all.x = T)
h_p1_sub_2 <- merge(h_p1_sub_1, hm_metori_met[,c(3,4,2)], by.x = 'Metlin', by.y = 'Metlin', all.x = T)


ID_network3 <- h_p1_sub_2[,c(1,6,7)]
ID_network3 <- ID_network3[!duplicated(ID_network3),]
dim(ID_network3)
colnames(ID_network3) <- c('ID','Lable','Group')
head(ID_network3)
ID_network4 <- h_p1_sub_2[,c(2,4,5)]
dim(ID_network4)
ID_network4 <- ID_network4[!duplicated(ID_network4),]
colnames(ID_network4) <- c('ID','Lable','Group')
ID_network5 <- rbind(ID_network3, ID_network4)
write.xlsx(ID_network5,'ID_4_network_p1_sub.xlsx')
write.xlsx(h_p1_sub,'edges_p1_sub.xlsx')
 

hm_metori_p2 <- read.xlsx('hmdb_from_metori_and_source.xlsx', sheet = 1)%>%
  filter(Source %in% c('Unknown'))
table(hm_metori_p2$Source)
hm_metori_p2_met <- merge(hm_metori_p2[,c(1,3)], metlin_2_hmdb[,c(1,2,4)], by.x = 'HMDB.ID',by.y = 'HMP',all.x = T)
h_p2_sub <- h_4S_1[h_4S_1$Metlin%in%hm_metori_p2_met$Metlin,]
h_p2_sub_1 <- merge(h_p2_sub, taxa_id, by.x = 'Taxa', by.y = 'ID', all.x = T)
h_p2_sub_2 <- merge(h_p2_sub_1, hm_metori_met[,c(3,4,2)], by.x = 'Metlin', by.y = 'Metlin', all.x = T)


ID_network6 <- h_p2_sub_2[,c(1,6,7)]
ID_network6 <- ID_network6[!duplicated(ID_network6),]
dim(ID_network6)
colnames(ID_network6) <- c('ID','Lable','Group')
head(ID_network6)
ID_network7 <- h_p2_sub_2[,c(2,4,5)]
dim(ID_network7)
ID_network7 <- ID_network7[!duplicated(ID_network7),]
colnames(ID_network7) <- c('ID','Lable','Group')
ID_network8 <- rbind(ID_network6, ID_network7)
write.xlsx(ID_network8,'ID_4_network_p2_sub.xlsx')
write.xlsx(h_p2_sub,'edges_p2_sub.xlsx')


enrich0 <- read.xlsx('enrichment.xlsx', sheet = 1)
enrich0$log <- log(enrich0$P.value, base = 0.05)
enrich0
enrich0$Pathway <- factor(enrich0$Pathway,levels = enrich0$Pathway[order(enrich0$log)])
ggplot(enrich0, aes(x = Pathway, y = log, fill = Class)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,vjust = 1, hjust = 1),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  ylab('log0.05(P value)') +
  scale_fill_manual(values = c('#00A087','#F39B7F'))
hm_metori_p2_met2 <- merge(hm_metori_p2_met[,c(1,4)], hmdb_db[,c(1,2,3,5,6,7,8,9,10)], by.x = 'HMDB.ID', by.y = 'ID', all.x = T)
head(hm_metori_p2_met2)
# write.csv(hm_metori_p2_met2,'~/Desktop/temp.csv')

# output degree of all network
h_4S_1_taxa <- merge(h_4S_1, taxa_id, by.x = 'Taxa', by.y = 'ID', all.x = T)
h_4S_1_taxa_met <- merge(h_4S_1_taxa, metlin_2_hmdb, by.x = 'Metlin', by.y = 'Metlin', all.x = T)
sort(table(h_4S_1_taxa$Genus))
h_4S_1_met <- merge(h_4S_1, metlin_2_hmdb, by.x = 'Metlin', by.y = 'Metlin', all.x = T)
h_4S_1_met$New <- paste(h_4S_1_met$Compound,h_4S_1_met$HMP, h_4S_1_met$Metlin,sep = '#')
data.frame(sort(table(h_4S_1_met$New),decreasing = T))
# write.csv(data.frame(sort(table(h_4S_1_met$New),decreasing = T)), '~/Desktop/temp1.csv')
# output degree of sub network with certain categories
h_p1_sub_taxa <- merge(h_p1_sub, taxa_id, by.x = 'Taxa', by.y = 'ID',all.x = T)
h_p1_sub_taxa_met <- merge(h_p1_sub_taxa, metlin_2_hmdb, by.x = 'Metlin', by.y = 'Metlin', all.x = T)  
data.frame(sort(table(h_p1_sub_taxa_met$Compound), decreasing = T))
unique(h_p1_sub_taxa[,c(4,5)])


# sankey plot for the last one
require(openxlsx)
require(dplyr)
sankey_all <- read.xlsx('sankey.xlsx',sheet = 2)
head(sankey_all)
sankey1 <- sankey_all[sankey_all$Cat=='Lipids and lipid-like molecules',]
sankey2 <- sankey_all[sankey_all$Cat=='Organoheterocyclic compounds',]
sankey3 <- sankey_all[sankey_all$Cat=='Organic acids and derivatives',]
node_san1 <- data.frame(name = c(as.character(sankey1$Source),
                              as.character(sankey1$Target)) %>% unique())
node_san2 <- data.frame(name = c(as.character(sankey2$Source),
                                 as.character(sankey2$Target)) %>% unique())
node_san3 <- data.frame(name = c(as.character(sankey3$Source),
                                 as.character(sankey3$Target)) %>% unique())
sankey1$IDsource <- match(sankey1$Source, node_san1$name)-1 
sankey1$IDtarget <- match(sankey1$Target, node_san1$name)-1
sankey2$IDsource <- match(sankey2$Source, node_san2$name)-1 
sankey2$IDtarget <- match(sankey2$Target, node_san2$name)-1
sankey3$IDsource <- match(sankey3$Source, node_san3$name)-1 
sankey3$IDtarget <- match(sankey3$Target, node_san3$name)-1
require(networkD3)
color_san1 <- 'd3.scaleOrdinal().domain(["Endocannabinoids","Fatty Acyls","Glycerolipids",
"Glycerophospholipids","Prenol lipids","Sphingolipids","Steroids and steroid derivatives",
"Unknown_1","Eicosanoids","Fatty acid esters","Fatty acids and conjugates","Fatty acyl glycosides",
"Fatty acyl thioesters","Fatty alcohol esters","Fatty alcohols","Fatty amides","Lineolic acids and derivatives",
"Glycerol vinyl ethers","Monoradylglycerols","Glycerophosphocholines","Glycerophosphoethanolamines",
"Glycerophosphoglycerols","Glycerophosphoinositols","Glycerophosphoserines","Diterpenoids","Monoterpenoids",
"Polyprenols","Retinoids","Sesquiterpenoids","Sesterterpenoids","Terpene glycosides","Terpene lactones",
"Triterpenoids", "Phosphosphingolipids","Glycosphingolipids","Bile acids, alcohols and derivatives","Cholestane steroids",
"Hydroxysteroids","Pyranoisoflavonoids","Steroid lactones","Steroidal glycosides","Stigmastanes and derivatives",
"Vitamin D and derivatives","Unknown_2"]).
range([rep("#00A087",7),rep("gray",37)])'

sankeyNetwork(Links = sankey1, Nodes = node_san1,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "Number", NodeID = "name", 
                      sinksRight=FALSE)
sankeyNetwork(Links = sankey2, Nodes = node_san2,
              Source = "IDsource", Target = "IDtarget",
              Value = "Number", NodeID = "name",  
              sinksRight=FALSE)
sankeyNetwork(Links = sankey3, Nodes = node_san3,
              Source = "IDsource", Target = "IDtarget",
              Value = "Number", NodeID = "name",  
              sinksRight=FALSE)
