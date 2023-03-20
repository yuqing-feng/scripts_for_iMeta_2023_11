# DMM analyses of the data on 16S amplicon sequencing and ITS amplicon sequencing
rm(list = ls())
library(DirichletMultinomial) #1.34.0
library(dplyr) 
library(reshape2)
library(tibble)
# 16S-based analyses
# Load data generated from qiime directly
raw_16S_L6 <- read.csv('level-6_no_archaea.csv')
format_16S_L6 <- raw_16S_L6 %>%
  remove_rownames %>% 
  column_to_rownames(var="index") %>%
  select(-c('Group1','Group2'))%>%
  as.matrix()
format_16S_L6

raw_ITS_L6 <- read.csv('level-6_for_bar.csv')
format_ITS_L6 <- raw_ITS_L6 %>%
  remove_rownames %>% 
  column_to_rownames(var="X") %>%
  select(-c('Group1','Group2')) %>%
  as.matrix()
format_ITS_L6

# step 1
set.seed(1234)
fit_16S <- lapply(1:6, dmn, count = format_16S_L6, verbose=TRUE)
fit_16S
# step 2
lplc_16S <- sapply(fit_16S, laplace) 
best_16S <- fit_16S[[which.min(lplc_16S)]]
best_16S
cluster_16S <- data.frame(apply(mixture(best_16S), 1, which.max))
head(cluster_16S)
colnames(cluster_16S) <- 'Cluster' 
importance_16S <- melt(fitted(best_16S))
head(importance_16S)
head(cluster_16S)
cluster_16S$Segment <- rep(c('Cecum','Duodenum','Ileum','Jejunum'),80)
head(cluster_16S)
table(cluster_16S)
# write.csv(cluster_16S,'cluster_16S.csv')
head(importance_16S)
colnames(importance_16S) <- c('Genus','Cluster','Importance')
# write.csv(importance_16S,'importance_16S.csv')

# step 1
set.seed(1234)
fit_ITS <- lapply(1:6, dmn, count = format_ITS_L6, verbose=TRUE)

# step 2
lplc_ITS <- sapply(fit_ITS, laplace)
best_ITS <- fit_ITS[[which.min(lplc_ITS)]]
best_ITS
d_ITS <- melt(fitted(best_ITS))
cluster_ITS <- data.frame(apply(mixture(best_ITS), 1, which.max)) 
colnames(cluster_ITS) <- 'Cluster'
importance_ITS <- melt(fitted(best_ITS))


head(importance_ITS)
colnames(importance_ITS) <- c("Genus",'Cluster','Importance')
head(cluster_ITS)
cluster_ITS$Segment <- rep(c('Cecum','Duodenum','Ileum','Jejunum'),80)
head(cluster_ITS)
table(cluster_ITS)

# write.csv(cluster_ITS,'cluster_ITS.csv')
head(importance_16S)
# write.csv(importance_ITS,'importance_ITS.csv')



library(ggplot2)
library(ggpubr)
library(reshape2)
cluster_16S <- read.csv('cluster_16S.csv',row.names = 1)
cluster_ITS <- read.csv('cluster_ITS.csv',row.names = 1)
head(cluster_16S)
cluster_16S$DPH <- rep(c(1,4,7,14,21,28,35,42),each=40)
cluster_16S
cluster_ITS$DPH <- rep(c(1,4,7,14,21,28,35,42),each=40)
cluster_ITS
cluster_16S2 <- melt(table(cluster_16S[,1:3]))
cluster_16S2$Cluster <- as.factor(cluster_16S2$Cluster)
cluster_16S2$Segment <- factor(cluster_16S2$Segment,levels = c('Duodenum','Jejunum','Ileum','Cecum'))
f14_1 <- ggplot(cluster_16S2,aes(x = '', y = value, fill = Cluster)) + 
  geom_bar(stat = 'identity',width = 1) +
  coord_polar('y',start = 0) + 
  facet_grid(Segment~DPH) + 
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F')) +
  theme_bw() + 
  xlab('Segment') + ylab("DPH") +
  theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid  = element_blank(),
        panel.border = element_blank())
f14_1
cluster_ITS2 <- melt(table(cluster_ITS[,1:3]))
cluster_ITS2$Cluster <- as.factor(cluster_ITS2$Cluster)
cluster_ITS2$Segment <- factor(cluster_ITS2$Segment,levels = c('Duodenum','Jejunum','Ileum','Cecum'))
f14_2 <- ggplot(cluster_ITS2,aes(x = '', y = value, fill = Cluster)) + 
  geom_bar(stat = 'identity',width = 1) +
  coord_polar('y',start = 0) + 
  facet_grid(Segment~DPH) + 
  scale_fill_manual(values = c('#8491B4','#91D1C2','#DC0000','#7E6148')) +
  theme_bw() + 
  xlab('Segment') + ylab("DPH") +
  theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid  = element_blank(),
        panel.border = element_blank())
f14_1 / f14_2


# prepare the contribute of the enterotypes
bac_contri <- read.table('importance_16S.csv',header = T,row.names = 1,sep = ',')
fun_contri <- read.table('importance_ITS.csv',header = T,row.names = 1, sep = ',')

bac1_contri <- bac_contri %>%
  filter(Cluster =='1') %>%
  slice_max(Importance, n = 5)
bac2_contri <- bac_contri %>%
  filter(Cluster =='2') %>%
  slice_max(Importance, n = 5)
bac3_contri <- bac_contri %>%
  filter(Cluster =='3') %>%
  slice_max(Importance, n = 5)
bac4_contri <- bac_contri %>%
  filter(Cluster =='4') %>%
  slice_max(Importance, n = 5)
bac5_contri <- bac_contri %>%
  filter(Cluster =='5') %>%
  slice_max(Importance, n = 5)
bac_all <- rbind(bac1_contri,bac2_contri,bac3_contri,bac4_contri,bac5_contri)
bac_all
bac_all$Genus[11] <- bac_all$Genus[8] <- 'f__Brocadiaceae'
f14_3 <- ggplot(bac_all,aes(x = Cluster,y = Genus,fill = log2(Importance+1))) +
  geom_tile() +
  scale_fill_continuous(low = '#E64B351F',high = '#E64B35') + 
  theme_bw() + 
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor = element_line(colour=NA)) 
f14_3
fun1_contri <- fun_contri %>%
  filter(Cluster =='1') %>%
  slice_max(Importance, n = 5)
fun2_contri <- fun_contri %>%
  filter(Cluster =='2') %>%
  slice_max(Importance, n = 5)
fun3_contri <- fun_contri %>%
  filter(Cluster =='3') %>%
  slice_max(Importance, n = 5)
fun4_contri <- fun_contri %>%
  filter(Cluster =='4') %>%
  slice_max(Importance, n = 5)

fun_all <- rbind(fun1_contri,fun2_contri,fun3_contri,fun4_contri)
fun_all$Genus[15] <- 'f__Nectriaceae'
f14_4 <- ggplot(fun_all,aes(x = Cluster,y = Genus,fill = log2(Importance+1))) +
  geom_tile() +
  scale_fill_continuous(low = '#4DBBD51F',high = '#4DBBD5') + 
  theme_bw() + 
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor = element_line(colour=NA)) 
f14_4
f14_3 + f14_4
# FEAST
# calculate the Transfer efficiency of the data from 
detach("package:DirichletMultinomial", unload = TRUE)
library(FEAST)

meta1 <- Load_metadata(metadata_path = '~/material/meta_jejunum.txt')
meta2 <- Load_metadata(metadata_path = '~/material/meta_ileum.txt')
meta3 <- Load_metadata(metadata_path = '~/material/meta_cecum.txt')
L6_bac <- Load_CountMatrix(CountMatrix_path = '~/material/bacteria.txt')
L6_fun <- Load_CountMatrix(CountMatrix_path = '~/material/fungi.txt')

out_bac_DJ <- FEAST(C = L6_bac, 
                    metadata = meta1,
                    different_sources_flag = 1,
                    dir_path = '~/material/04feast2/',
                    outfile = '~/material/DJ_bac_source_contributions_matrix.txt')
out_bac_JI <- FEAST(C = L6_bac, 
                    metadata = meta2,
                    different_sources_flag = 1,
                    dir_path = '~/material/04feast2/',
                    outfile = '~/material/JI_bac_source_contributions_matrix.txt')
out_bac_IC <- FEAST(C = L6_bac, 
                    metadata = meta3,
                    different_sources_flag = 1,
                    dir_path = '~/material/04feast2/',
                    outfile = '~/material/IC_bac_source_contributions_matrix.txt')
out_fun_DJ <- FEAST(C = L6_fun, 
                    metadata = meta1,
                    different_sources_flag = 1,
                    dir_path = '~/material/04feast2/',
                    outfile = '~/material/DJ_fun_source_contributions_matrix.txt')
out_fun_JI <- FEAST(C = L6_fun, 
                    metadata = meta2,
                    different_sources_flag = 1,
                    dir_path = '~/material/04feast2/',
                    outfile = '~/material/JI_fun_source_contributions_matrix.txt')
out_fun_IC <- FEAST(C = L6_fun, 
                    metadata = meta3,
                    different_sources_flag = 1,
                    dir_path = '~/material/04feast2/',
                    outfile = '~/material/IC_fun_source_contributions_matrix.txt')


for(i in 1:80){
  eff_b1[i] <- out_bac_DJ[[i]][i]
}
for(i in 1:80){
  eff_b2[i] <- out_bac_JI[[i*2]][i]
}
for(i in 1:80){
  eff_b3[i] <- out_bac_IC[[160+i]][i]
}
for(i in 1:80){
  eff_f1[i] <- out_fun_DJ[[i]][i]
}
for(i in 1:80){
  eff_f2[i] <- out_fun_JI[[i*2]][i]
}
for(i in 1:80){
  eff_f3[i] <- out_fun_IC[[160+i]][i]
}
out_bac_DJ
eff_all <- data.frame('Efficiency' = c(eff_b1,eff_b2,eff_b3,eff_f1,eff_f2,eff_f3),
                      'Segment' = c(rep('DJ',80),rep('JI',80),rep('IC',80),rep('DJ',80),rep('JI',80),rep('IC',80)),
                      'Group' = c(rep('Bacteria',240),rep('Fungi',240)),
                      'DPH' = c(rep(c(1,4,7,14,21,28,35,42),each=10),rep(c(1,4,7,14,21,28,35,42),each=10),
                                rep(c(1,4,7,14,21,28,35,42),each=10),rep(c(1,4,7,14,21,28,35,42),each=10),
                                rep(c(1,4,7,14,21,28,35,42),each=10),rep(c(1,4,7,14,21,28,35,42),each=10)))
eff_all$Segment <- factor(eff_all$Segment, levels = c('DJ','JI','IC'))
f14_11 <- ggplot(eff_all,aes(x = as.factor(DPH),y = Efficiency,fill = Segment)) + 
  geom_boxplot(outlier.shape = NA) +
  theme_bw() + 
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087')) + 
  facet_grid(Group~Segment) + 
  xlab('DPH') + 
  theme(legend.position = 'none')

eff_all %>%
  filter(Group == 'Bacteria') %>%
  select(c(Segment, Efficiency)) %>%
  rownames_to_column() %>%
  rename(group = Segment) %>%
  KwWlx(i = 3)

eff_all %>%
  filter(Group == 'Fungi') %>%
  select(c(Segment, Efficiency)) %>%
  rownames_to_column() %>%
  rename(group = Segment) %>%
  KwWlx(i = 3)

cluster_16S2
cluster_ITS2
library(openxlsx)
ab_temp <- read.xlsx('01data_for_quantitive.xlsx',sheet = 11)
rownames(cluster_16S) == ab_temp$ID
head(cluster_16S)
head(ab_temp)
clus_16S_ab <- cbind(cluster_16S,ab_temp[,4:5])
clus_ITS_ab <- cbind(cluster_ITS,ab_temp[,4:5])

library(EasyAovWlxPlot)
head(clus_16S_ab)
test1 <- data.frame(rownames(clus_16S_ab),clus_16S_ab$Cluster,clus_16S_ab$`log10(16S)`)
test1           
colnames(test1)[2] <- 'group'
KwWlx(test1,i=3) 

clus_16S_ab$Cluster <- as.factor(clus_16S_ab$Cluster)
f14_5 <- ggplot(clus_16S_ab,aes(x = Cluster, y = `log10(16S)`, fill = Cluster)) + 
  geom_boxplot(outlier.shape = NA) +
  ylab('Number of copy (log10)') +
  xlab('Cluster') +
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor = element_line(colour=NA),
        legend.position = 'none') +
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F')) +
  annotate('text',x = 1,y = 12,label = 'b') +
  annotate('text',x = 2,y = 12,label = 'c') +
  annotate('text',x = 3,y = 12,label = 'a') +
  annotate('text',x = 4,y = 12,label = 'd') + 
  annotate('text',x = 5,y = 12,label = 'cd') +
  ggtitle('Bacteria')



test2 <- data.frame(rownames(clus_ITS_ab),clus_ITS_ab$Cluster,clus_16S_ab$`log10(ITS)`)
test2           
colnames(test2)[2] <- 'group'
KwWlx(test2,i=3) 

clus_ITS_ab$Cluster <- as.factor(clus_ITS_ab$Cluster)
f14_6 <- ggplot(clus_ITS_ab,aes(x = Cluster, y = `log10(ITS)`, fill = Cluster)) + 
  geom_boxplot(outlier.shape = NA) +
  ylab('Number of copy (log10)') +
  xlab('Cluster') +
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor = element_line(colour=NA),
        legend.position = 'none') +
  scale_fill_manual(values = c('#8491B4','#91D1C2','#DC0000','#7E6148')) +
  annotate('text',x = 1,y = 12,label = 'a') +
  annotate('text',x = 2,y = 12,label = 'b') +
  annotate('text',x = 3,y = 12,label = 'a') +
  annotate('text',x = 4,y = 12,label = 'b') +
  ggtitle('Fungi')
f14_5 + f14_6

raw_16S_L6 <- read.csv('level-6_no_archaea.csv')
raw_ITS_L6 <- read.csv('level-6_for_bar.csv')
library(vegan)
alpha_bacteria <- diversity(raw_16S_L6[,2:530])
alpha_fungi <- diversity(raw_ITS_L6[,2:604])
a_bac <- data.frame('ID'= rownames(cluster_16S),
               'group' = as.factor(cluster_16S$Cluster),
               'Shannon' = alpha_bacteria)
a_bac
a_fun <- data.frame('ID' = rownames(cluster_ITS),
                    'group' = as.factor(cluster_ITS$Cluster),
                    'Shannon' = alpha_fungi)

KwWlx(a_bac,i=3) 

KwWlx(a_fun,i=3) 
f14_7 <- ggplot(a_bac, aes(x = group,y = Shannon,fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  ylab('Shannon index') +
  xlab('Cluster') +
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor = element_line(colour=NA),
        legend.position = 'none') +
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F')) +
  annotate('text',x = 1,y = 4,label = 'c') +
  annotate('text',x = 2,y = 4,label = 'c') +
  annotate('text',x = 3,y = 4,label = 'a') +
  annotate('text',x = 4,y = 4,label = 'b') + 
  annotate('text',x = 5,y = 4,label = 'b') +
  ggtitle('Bacteria')
f14_8 <- ggplot(a_fun, aes(x = group, y = Shannon, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  ylab('Shannon index') +
  xlab('Cluster') +
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor = element_line(colour=NA),
        legend.position = 'none') +
  scale_fill_manual(values = c('#8491B4','#91D1C2','#DC0000','#7E6148')) +
  annotate('text',x = 1,y = 4.5,label = 'b') +
  annotate('text',x = 2,y = 4.5,label = 'c') +
  annotate('text',x = 3,y = 4.5,label = 'a') +
  annotate('text',x = 4,y = 4.5,label = 'd') +
  ggtitle('Fungi')
# f14_9 f14_10 beta
library(tibble)
format_16S_L6 <- raw_16S_L6 %>%
  remove_rownames %>% 
  column_to_rownames(var="index") %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  mutate(across(Granulicella:uncultured.bacterium.40, ~ . / sum)) %>%
  select(-c('sum')) %>%
  mutate(Group1 = factor(Group1, levels = c('Duodenum','Jejunum','Ileum','Cecum'))) %>%
  select(-c('Group2'))
format_16S_L6
format_ITS_L6 <- raw_ITS_L6 %>%
  remove_rownames %>% 
  column_to_rownames(var="X") %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  mutate(across(X__:unidentified.53, ~ . / sum)) %>%
  select(-c('sum')) %>%
  mutate(Group1 = factor(Group1, levels = c('Duodenum','Jejunum','Ileum','Cecum'))) %>%
  select(-c('Group2'))
format_ITS_L6

format_16S_L6_AA <- format_16S_L6
format_ITS_L6_AA <- format_ITS_L6

for(i in 1:nrow(ab_temp)){
  format_16S_L6_AA[i,1:529] <- format_16S_L6[i,1:529] * 10^(ab_temp$`log10(16S)`[i])
  format_ITS_L6_AA[i,1:603] <- format_ITS_L6[i,1:603] * 10^(ab_temp$`log10(ITS)`[i])
}

vegdist_16S_AA <- vegdist(as.matrix(format_16S_L6_AA[,1:529]))
pcoa_16S_AA <- cmdscale(vegdist_16S_AA, eig = TRUE)
df_pcoa_16S_AA <- as.data.frame(pcoa_16S_AA$points)
colnames(df_pcoa_16S_AA) <- c('PCoA1', 'PCoA2')
PCo1_16S_AA <- round(eigenvals(pcoa_16S_AA)[1]/sum(eigenvals(pcoa_16S_AA)) * 100, 1)
PCo2_16S_AA <- round(eigenvals(pcoa_16S_AA)[2]/sum(eigenvals(pcoa_16S_AA)) * 100, 1)
print(PCo1_16S_AA)#23.5
print(PCo2_16S_AA)#16.7
df_pcoa_16S_AA$Cluster <- as.factor(clus_16S_ab$Cluster)
f14_9 <- ggplot(df_pcoa_16S_AA, aes(x = PCoA1, y = PCoA2, color = Cluster)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F')) +
  xlab('PCoA 1 (23.5%)') +
  ylab('PCoA 2 (16.7%)') +
  ggtitle('Bacteria') +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor = element_line(colour=NA))

vegdist_ITS_AA <- vegdist(as.matrix(format_ITS_L6_AA[,1:603]))
pcoa_ITS_AA <- cmdscale(vegdist_ITS_AA, eig = TRUE)
df_pcoa_ITS_AA <- as.data.frame(pcoa_ITS_AA$points)
colnames(df_pcoa_ITS_AA) <- c('PCoA1', 'PCoA2')
PCo1_ITS_AA <- round(eigenvals(pcoa_ITS_AA)[1]/sum(eigenvals(pcoa_ITS_AA)) * 100, 1)
PCo2_ITS_AA <- round(eigenvals(pcoa_ITS_AA)[2]/sum(eigenvals(pcoa_ITS_AA)) * 100, 1)
print(PCo1_ITS_AA)#19.6
print(PCo2_ITS_AA)#14.2
df_pcoa_ITS_AA$Cluster <- as.factor(cluster_ITS$Cluster)
f14_10 <- ggplot(df_pcoa_ITS_AA, aes(x = PCoA1, y = PCoA2, color = Cluster)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) +
  xlab('PCoA 1 (19.6%)') +
  ylab('PCoA 2 (14.2%)') +
  ggtitle('Fungi') +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor = element_line(colour=NA))
f14_10
library(patchwork)
svg('~/Desktop/figure materials/2nd/figureS6A.svg',width = 10.5,height = 3)
f14_3 + f14_4
dev.off()
svg('~/Desktop/figure materials/2nd/figureS6B.svg',width = 10.5, height = 9)
f14_5 + f14_6 + f14_7 + f14_8 + f14_9 + f14_10 + plot_layout(nrow =3)
dev.off()
library(pairwiseAdonis)
pairwise.adonis(format_16S_L6_AA[,1:529],clus_16S_ab$Cluster)
pairwise.adonis(format_ITS_L6_AA[,1:603],clus_ITS_ab$Cluster)



# BetaNTI and RCbray
# Calculate the betaNTI and RCbray of 16S and ITS
# betaNTI: ape and iCAMP
# RCbray: vegan and script online
# Calculated on server (64 threads and 384 GB memory in total)

# betaNTI
library(ape)
library(iCAMP)
# comp_16S_0 <- read.table('04table/table.tsv',header=T,row.names = 1,sep = '\t')
# cS_16S_1 <- colSums(comp_16S_0)
# comp_16S_1 <- comp_16S_0
# for(i in 1:320){
#   comp_16S_1[,i] <- comp_16S_0[,i]/cS_16S_1[i]
# }
# comp_16S_2 <- t(comp_16S_1)
# comp_16S_2_c <- comp_16S_2[seq(from = 1, to = 320,by = 4),]
# root_tree_16S <- read.tree('03rooted_tree/tree.nwk')
# nti_16S_C <- ses.comdistnt(comp_16S_2_c, cophenetic(root_tree_16S), abundance.weighted = T,cores = 20)
# write.table(nti_16S_C$comdistnt.obs.z,'NTI_16S_C.txt',quote = F)
# 
# comp_16S_2_d <- comp_16S_2[seq(from = 2, to = 320,by = 4),]
# root_tree_16S <- read.tree('03rooted_tree/tree.nwk')
# nti_16S_D <- ses.comdistnt(comp_16S_2_d, cophenetic(root_tree_16S), abundance.weighted = T,cores = 20)
# write.table(nti_16S_D$comdistnt.obs.z,'NTI_16S_D.txt',quote = F)
# 
# comp_16S_2_i <- comp_16S_2[seq(from = 3, to = 320,by = 4),]
# root_tree_16S <- read.tree('03rooted_tree/tree.nwk')
# nti_16S_I <- ses.comdistnt(comp_16S_2_i, cophenetic(root_tree_16S), abundance.weighted = T,cores = 20)
# write.table(nti_16S_I$comdistnt.obs.z,'NTI_16S_I.txt',quote = F)
# 
# comp_16S_2_j <- comp_16S_2[seq(from = 4, to = 320,by = 4),]
# root_tree_16S <- read.tree('03rooted_tree/tree.nwk')
# nti_16S_J <- ses.comdistnt(comp_16S_2_j, cophenetic(root_tree_16S), abundance.weighted = T,cores = 20)
# write.table(nti_16S_J$comdistnt.obs.z,'NTI_16S_J.txt',quote = F)


# betaNTI ITS
# comp_ITS_0 <- read.table('04table/table.tsv',header=T,row.names = 1,sep = '\t')
# cS_ITS_1 <- colSums(comp_ITS_0)
# comp_ITS_1 <- comp_ITS_0
# for(i in 1:320){
#   comp_ITS_1[,i] <- comp_ITS_0[,i]/cS_ITS_1[i]
# }
# comp_ITS_2 <- t(comp_ITS_1)
# comp_ITS_2_c <- comp_ITS_2[seq(from = 1, to = 320,by = 4),]
# root_tree_ITS <- read.tree('03rooted_tree/tree.nwk')
# nti_ITS_C <- ses.comdistnt(comp_ITS_2_c, cophenetic(root_tree_ITS), abundance.weighted = T,cores = 20)
# write.table(nti_ITS_C$comdistnt.obs.z,'NTI_ITS_C.txt',quote = F)
# 
# comp_ITS_2_d <- comp_ITS_2[seq(from = 2, to = 320,by = 4),]
# root_tree_ITS <- read.tree('03rooted_tree/tree.nwk')
# nti_ITS_D <- ses.comdistnt(comp_ITS_2_d, cophenetic(root_tree_ITS), abundance.weighted = T,cores = 20)
# write.table(nti_ITS_D$comdistnt.obs.z,'NTI_ITS_D.txt',quote = F)
# 
# comp_ITS_2_i <- comp_ITS_2[seq(from = 3, to = 320,by = 4),]
# root_tree_ITS <- read.tree('03rooted_tree/tree.nwk')
# nti_ITS_I <- ses.comdistnt(comp_ITS_2_i, cophenetic(root_tree_ITS), abundance.weighted = T,cores = 20)
# write.table(nti_ITS_I$comdistnt.obs.z,'NTI_ITS_I.txt',quote = F)
# 
# comp_ITS_2_j <- comp_ITS_2[seq(from = 4, to = 320,by = 4),]
# root_tree_ITS <- read.tree('03rooted_tree/tree.nwk')
# nti_ITS_J <- ses.comdistnt(comp_ITS_2_j, cophenetic(root_tree_ITS), abundance.weighted = T,cores = 20)
# write.table(nti_ITS_J$comdistnt.obs.z,'NTI_ITS_J.txt',quote = F)

# RCbray
# raup_crick=function(spXsite, plot_names_in_col1=FALSE, classic_metric=FALSE, split_ties=TRUE, reps=9999, set_all_species_equal=FALSE, as.distance.matrix=TRUE, report_similarity=FALSE){
#   ##expects a species by site matrix for spXsite, with row names for plots, or optionally plots named in column 1. 
#   # By default calculates a modification of the Raup-Crick metric (standardizing the metric to range from -1 to 1 instead of 0 to 1). 
#   # Specifying classic_metric=TRUE instead calculates the original Raup-Crick metric that ranges from 0 to 1. 
#   # The option split_ties (defaults to TRUE) adds half of the number of null observations that are equal to the observed number of shared species to the calculation- this is highly recommended.
#   # The argument report_similarity defaults to FALSE so the function reports a dissimilarity (which is appropriate as a measure of beta diversity).
#   # Setting report_similarity=TRUE returns a measure of similarity, as Raup and Crick originally specified.  
#   # If ties are split (as we recommend) the dissimilarity (default) and similarity (set report_similarity=TRUE) calculations can be flipped by multiplying by -1 (for our modification, which ranges from -1 to 1) or by subtracting the metric from 1 (for the classic metric which ranges from 0 to 1). 
#   # If ties are not split (and there are ties between the observed and expected shared number of species) this conversion will not work.
#   # The argument reps specifies the number of randomizations (a minimum of 999 is recommended- default is 9999).  
#   # set_all_species_equal weights all species equally in the null model instead of weighting species by frequency of occupancy.  
#   ##Note that the choice of how many plots (rows) to include has a real impact on the metric, as species and their occurrence frequencies across the set of plots is used to determine gamma and the frequency with which each species is drawn from the null model	
#   ##this section moves plot names in column 1 (if specified as being present) into the row names of the matrix and drops the column of names
#   if(plot_names_in_col1){
#     row.names(spXsite)<-spXsite[,1]
#     spXsite<-spXsite[,-1]
#   }
#   ## count number of sites and total species richness across all plots (gamma)
#   n_sites<-nrow(spXsite)
#   gamma<-ncol(spXsite)
#   ##make the spXsite matrix into a pres/abs. (overwrites initial spXsite matrix):
#   ceiling(spXsite/max(spXsite))->spXsite
#   ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
#   occur<-apply(spXsite, MARGIN=2, FUN=sum)
#   ##NOT recommended- this is a non-trivial change to the metric:
#   ##sets all species to occur with equal frequency in the null model
#   ##e.g.- discards any occupancy frequency information
#   if(set_all_species_equal){
#     occur<-rep(1,gamma)
#   }
#   ## determine how many unique species richness values are in the dataset
#   ##this is used to limit the number of null communities that have to be calculated
#   alpha_levels<-sort(unique(apply(spXsite, MARGIN=1, FUN=sum)))
#   ##make_null:
#   ##alpha_table is used as a lookup to help identify which null distribution to use for the tests later.  It contains one row for each combination of alpha richness levels. 
#   alpha_table<-data.frame(c(NA), c(NA))
#   names(alpha_table)<-c("smaller_alpha", "bigger_alpha")
#   col_count<-1
#   ##null_array will hold the actual null distribution values.  Each element of the array corresponds to a null distribution for each combination of alpha values.  The alpha_table is used to point to the correct null distribution- the row numbers of alpha_table correspond to the [[x]] indices of the null_array.  Later the function will find the row of alpha_table with the right combination of alpha values.  That row number is used to identify the element of null_array that contains the correct null distribution for that combination of alpha levels. 
#   null_array<-list()
#   ##looping over each combination of alpha levels:
#   for(a1 in 1:length(alpha_levels)){
#     for(a2 in a1:length(alpha_levels)){
#       ##build a null distribution of the number of shared species for a pair of alpha values:
#       null_shared_spp<-NULL
#       for(i in 1:reps){
#         ##two empty null communities of size gamma:
#         com1<-rep(0,gamma)
#         com2<-rep(0,gamma)
#         ##add alpha1 number of species to com1, weighting by species occurrence frequencies:
#         com1[sample(1:gamma, alpha_levels[a1], replace=FALSE, prob=occur)]<-1
#         ##same for com2:
#         com2[sample(1:gamma, alpha_levels[a2], replace=FALSE, prob=occur)]<-1
#         ##how many species are shared in common?
#         null_shared_spp[i]<-sum((com1+com2)>1)
#       }
#       ##store null distribution, record values for alpha 1 and 2 in the alpha_table to help find the correct null distribution later:
#       null_array[[col_count]]<-null_shared_spp
#       alpha_table[col_count, which(names(alpha_table)=="smaller_alpha")]<-alpha_levels[a1]
#       alpha_table[col_count, which(names(alpha_table)=="bigger_alpha")]<-alpha_levels[a2]
#       #increment the counter for the columns of the alpha table/ elements of the null array
#       col_count<-col_count+1
#     }
#   }
#   ##create a new column with both alpha levels to match on:
#   alpha_table$matching<-paste(alpha_table[,1], alpha_table[,2], sep="_")
#   #####################
#   ##do the test:
#   ##build a site by site matrix for the results, with the names of the sites in the row and col names:
#   results<-matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(spXsite), row.names(spXsite)))
#   ##for each pair of sites (duplicates effort now to make a full matrix instead of a half one- but this part should be minimal time as compared to the null model building)
#   for(i in 1:n_sites){
#     for(j in 1:n_sites){
#       ##how many species are shared between the two sites:
#       n_shared_obs<-sum((spXsite[i,]+spXsite[j,])>1)
#       ## what was the observed richness of each site?
#       obs_a1<-sum(spXsite[i,])
#       obs_a2<-sum(spXsite[j,])
#       ##place these alphas into an object to match against alpha_table (sort so smaller alpha is first)
#       obs_a_pair<-sort(c(obs_a1, obs_a2))
#       ##match against the alpha table- row index identifies which element of the null array contains the correct null distribution for the observed combination of alpha values:
#       null_index<-which(alpha_table$matching==paste(obs_a_pair[1], obs_a_pair[2], sep="_"))
#       ##how many null observations is the observed value tied with?
#       num_exact_matching_in_null<-sum(null_array[[null_index]]==n_shared_obs)
#       ##how many null values are bigger than the observed value?
#       num_greater_in_null<-sum(null_array[[null_index]]>n_shared_obs)
#       rc<-(num_greater_in_null)/reps
#       if(split_ties){
#         rc<-((num_greater_in_null+(num_exact_matching_in_null)/2)/reps)
#       }
#       if(!classic_metric){
#         ##our modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
#         rc<-(rc-.5)*2
#       }
#       ## at this point rc represents an index of dissimilarity- multiply by -1 to convert to a similarity as specified in the original 1979 Raup Crick paper
#       if(report_similarity & !classic_metric){
#         rc<- rc*-1
#       }
#       ## the switch to similarity is done differently if the original 0 to 1 range of the metric is used:
#       if(report_similarity & classic_metric){
#         rc<- 1-rc
#       }
#       ##store the metric in the results matrix:
#       results[i,j]<-round(rc, digits=2)
#     }
#   }
#   if(as.distance.matrix){
#     results<-as.dist(results)
#   }	
#   return(results)
# }

# function 
# raup_crick()

# RCbray 16S
library(phyloseq)
library(dplyr)
# comp_16S_0 <- read.table('table.tsv',header=T,row.names = 1, sep = '\t')
# dim(comp_16S_0)
# comp_16S_c <- comp_16S_0[,seq(from = 1, to = 320, by = 4)]
# comp_16S_d <- comp_16S_0[,seq(from = 2, to = 320, by = 4)]
# comp_16S_i <- comp_16S_0[,seq(from = 3, to = 320, by = 4)]
# comp_16S_j <- comp_16S_0[,seq(from = 4, to = 320, by = 4)]
# 
# set.seed(1234)
# comp_16S_c2 <- phyloseq(otu_table(comp_16S_c, taxa_are_rows = TRUE)) %>%
#   rarefy_even_depth(., min(colSums(comp_16S_c)))
# set.seed(1234)
# comp_16S_d2 <- phyloseq(otu_table(comp_16S_d, taxa_are_rows = TRUE)) %>%
#   rarefy_even_depth(., min(colSums(comp_16S_d)))
# set.seed(1234)
# comp_16S_i2 <- phyloseq(otu_table(comp_16S_i, taxa_are_rows = TRUE)) %>%
#   rarefy_even_depth(., min(colSums(comp_16S_i)))
# set.seed(1234)
# comp_16S_j2 <- phyloseq(otu_table(comp_16S_j, taxa_are_rows = TRUE)) %>%
#   rarefy_even_depth(., min(colSums(comp_16S_j)))
# 
# RC_16S_C <- raup_crick(t(comp_16S_c2@.Data))
# RC_16S_D <- raup_crick(t(comp_16S_d2@.Data))
# RC_16S_I <- raup_crick(t(comp_16S_i2@.Data))
# RC_16S_J <- raup_crick(t(comp_16S_j2@.Data))
# RC_16S_C2 <- as.matrix(RC_16S_C)
# RC_16S_D2 <- as.matrix(RC_16S_D)
# RC_16S_I2 <- as.matrix(RC_16S_I)
# RC_16S_J2 <- as.matrix(RC_16S_J)
# write.table(RC_16S_C2,'RCbray_16S_C.txt',quote=F)
# write.table(RC_16S_D2,'RCbray_16S_D.txt',quote=F)
# write.table(RC_16S_I2,'RCbray_16S_I.txt',quote=F)
# write.table(RC_16S_J2,'RCbray_16S_J.txt',quote=F)
# 
# RCbray ITS
# library(phyloseq)
# library(dplyr)
# comp_ITS_0 <- read.table('table.tsv',header=T,row.names = 1, sep = '\t')
# dim(comp_ITS_0)
# comp_ITS_c <- comp_ITS_0[,seq(from = 1, to = 320, by = 4)]
# comp_ITS_d <- comp_ITS_0[,seq(from = 2, to = 320, by = 4)]
# comp_ITS_i <- comp_ITS_0[,seq(from = 3, to = 320, by = 4)]
# comp_ITS_j <- comp_ITS_0[,seq(from = 4, to = 320, by = 4)]
# 
# set.seed(1234)
# comp_ITS_c2 <- phyloseq(otu_table(comp_ITS_c, taxa_are_rows = TRUE)) %>%
#   rarefy_even_depth(., min(colSums(comp_ITS_c)))
# set.seed(1234)
# comp_ITS_d2 <- phyloseq(otu_table(comp_ITS_d, taxa_are_rows = TRUE)) %>%
#   rarefy_even_depth(., min(colSums(comp_ITS_d)))
# set.seed(1234)
# comp_ITS_i2 <- phyloseq(otu_table(comp_ITS_i, taxa_are_rows = TRUE)) %>%
#   rarefy_even_depth(., min(colSums(comp_ITS_i)))
# set.seed(1234)
# comp_ITS_j2 <- phyloseq(otu_table(comp_ITS_j, taxa_are_rows = TRUE)) %>%
#   rarefy_even_depth(., min(colSums(comp_ITS_j)))
# 
# RC_ITS_C <- raup_crick(t(comp_ITS_c2@.Data))
# RC_ITS_D <- raup_crick(t(comp_ITS_d2@.Data))
# RC_ITS_I <- raup_crick(t(comp_ITS_i2@.Data))
# RC_ITS_J <- raup_crick(t(comp_ITS_j2@.Data))
# RC_ITS_C2 <- as.matrix(RC_ITS_C)
# RC_ITS_D2 <- as.matrix(RC_ITS_D)
# RC_ITS_I2 <- as.matrix(RC_ITS_I)
# RC_ITS_J2 <- as.matrix(RC_ITS_J)
# write.table(RC_ITS_C2,'RCbray_ITS_C.txt',quote=F)
# write.table(RC_ITS_D2,'RCbray_ITS_D.txt',quote=F)
# write.table(RC_ITS_I2,'RCbray_ITS_I.txt',quote=F)
# write.table(RC_ITS_J2,'RCbray_ITS_J.txt',quote=F)


# Calculate the stochasticity of the gut microbiota
# two parts:
# First: the stochasticity of each segment
# Second: the stochasticity at each time point of each segment
library(dplyr)
library(reshape2)
library(tidyr)
betaNTI_16S_C0 <- read.table('NTI_16S_C.txt', header = T, row.names = 1) 
betaNTI_16S_C0[upper.tri(betaNTI_16S_C0)] <- NA
betaNTI_16S_C <- betaNTI_16S_C0 %>%
  rownames_to_column() %>%
  gather(ID2, Value, -rowname) %>%
  mutate(values = -Value) %>%
  select(-c(Value)) %>%
  rename(NTI = values) %>%
  rename(R1 = rowname) %>%
  rename(ID1 = ID2) %>%
  drop_na()

betaNTI_16S_D0 <- read.table('NTI_16S_D.txt', header = T, row.names = 1)
betaNTI_16S_D0[upper.tri(betaNTI_16S_D0)] <- NA
betaNTI_16S_D <- betaNTI_16S_D0 %>%
  rownames_to_column() %>%
  gather(ID2, Value, -rowname) %>%
  mutate(values = -Value) %>%
  select(-c(Value)) %>%
  rename(NTI = values) %>%
  rename(R1 = rowname) %>%
  rename(ID1 = ID2) %>%
  drop_na()

betaNTI_16S_I0 <- read.table('NTI_16S_I.txt', header = T,row.names = 1)
betaNTI_16S_I0[upper.tri(betaNTI_16S_I0)] <- NA
betaNTI_16S_I <- betaNTI_16S_I0 %>%
  rownames_to_column() %>%
  gather(ID2, Value, -rowname) %>%
  mutate(values = -Value) %>%
  select(-c(Value)) %>%
  rename(NTI = values) %>%
  rename(R1 = rowname) %>%
  rename(ID1 = ID2) %>%
  drop_na()

betaNTI_16S_J0 <- read.table('NTI_16S_J.txt', header = T,row.names = 1)
betaNTI_16S_J0[upper.tri(betaNTI_16S_J0)] <- NA
betaNTI_16S_J <- betaNTI_16S_J0 %>%
  rownames_to_column() %>%
  gather(ID2, Value, -rowname) %>%
  mutate(values = -Value) %>%
  select(-c(Value)) %>%
  rename(NTI = values) %>%
  rename(R1 = rowname) %>%
  rename(ID1 = ID2) %>%
  drop_na()

betaNTI_ITS_C0 <- read.table('NTI_ITS_C.txt', header = T, row.names = 1) 
betaNTI_ITS_C0[upper.tri(betaNTI_ITS_C0)] <- NA
betaNTI_ITS_C <- betaNTI_ITS_C0 %>%
  rownames_to_column() %>%
  gather(ID2, Value, -rowname) %>%
  mutate(values = -Value) %>%
  select(-c(Value)) %>%
  rename(NTI = values) %>%
  rename(R1 = rowname) %>%
  rename(ID1 = ID2) %>%
  drop_na()

betaNTI_ITS_D0 <- read.table('NTI_ITS_D.txt', header = T, row.names = 1)
betaNTI_ITS_D0[upper.tri(betaNTI_ITS_D0)] <- NA
betaNTI_ITS_D <- betaNTI_ITS_D0 %>%
  rownames_to_column() %>%
  gather(ID2, Value, -rowname) %>%
  mutate(values = -Value) %>%
  select(-c(Value)) %>%
  rename(NTI = values) %>%
  rename(R1 = rowname) %>%
  rename(ID1 = ID2) %>%
  drop_na()

betaNTI_ITS_I0 <- read.table('NTI_ITS_I.txt', header = T,row.names = 1)
betaNTI_ITS_I0[upper.tri(betaNTI_ITS_I0)] <- NA
betaNTI_ITS_I <- betaNTI_ITS_I0 %>%
  rownames_to_column() %>%
  gather(ID2, Value, -rowname) %>%
  mutate(values = -Value) %>%
  select(-c(Value)) %>%
  rename(NTI = values) %>%
  rename(R1 = rowname) %>%
  rename(ID1 = ID2) %>%
  drop_na()

betaNTI_ITS_J0 <- read.table('NTI_ITS_J.txt', header = T,row.names = 1)
betaNTI_ITS_J0[upper.tri(betaNTI_ITS_J0)] <- NA
betaNTI_ITS_J <- betaNTI_ITS_J0 %>%
  rownames_to_column() %>%
  gather(ID2, Value, -rowname) %>%
  mutate(values = -Value) %>%
  select(-c(Value)) %>%
  rename(NTI = values) %>%
  rename(R1 = rowname) %>%
  rename(ID1 = ID2) %>%
  drop_na()


RCbray_16S_C0 <- read.table('RCbray_16S_C.txt', header = T, row.names = 1) 
RCbray_16S_C0[upper.tri(RCbray_16S_C0)] <- NA
for(i in 1:80){
  RCbray_16S_C0[i,i] <- NA
}
RCbray_16S_C <- RCbray_16S_C0 %>%
  rownames_to_column() %>%
  gather(ID2, Value, -rowname) %>%
  rename(RC = Value) %>%
  rename(R2 = rowname) %>%
  drop_na()

RCbray_16S_D0 <- read.table('RCbray_16S_D.txt', header = T, row.names = 1)
RCbray_16S_D0[upper.tri(RCbray_16S_D0)] <- NA
for(i in 1:80){
  RCbray_16S_D0[i,i] <- NA
}
RCbray_16S_D <- RCbray_16S_D0 %>%
  rownames_to_column() %>%
  gather(ID2, Value, -rowname) %>%
  rename(RC = Value) %>%
  rename(R2 = rowname) %>%
  drop_na()

RCbray_16S_I0 <- read.table('RCbray_16S_I.txt', header = T,row.names = 1)
RCbray_16S_I0[upper.tri(RCbray_16S_I0)] <- NA
for(i in 1:80){
  RCbray_16S_I0[i,i] <- NA
}
RCbray_16S_I <- RCbray_16S_I0 %>%
  rownames_to_column() %>%
  gather(ID2, Value, -rowname) %>%
  rename(RC = Value) %>%
  rename(R2 = rowname) %>%
  drop_na()

RCbray_16S_J0 <- read.table('RCbray_16S_J.txt', header = T,row.names = 1)
RCbray_16S_J0[upper.tri(RCbray_16S_J0)] <- NA
for(i in 1:80){
  RCbray_16S_J0[i,i] <- NA
}
RCbray_16S_J <- RCbray_16S_J0 %>%
  rownames_to_column() %>%
  gather(ID2, Value, -rowname) %>%
  rename(RC = Value) %>%
  rename(R2 = rowname) %>%
  drop_na()

RCbray_ITS_C0 <- read.table('RCbray_ITS_C.txt', header = T, row.names = 1) 
RCbray_ITS_C0[upper.tri(RCbray_ITS_C0)] <- NA
for(i in 1:80){
  RCbray_ITS_C0[i,i] <- NA
}
RCbray_ITS_C <- RCbray_ITS_C0 %>%
  rownames_to_column() %>%
  gather(ID2, Value, -rowname) %>%
  rename(RC = Value) %>%
  rename(R2 = rowname) %>%
  drop_na()

RCbray_ITS_D0 <- read.table('RCbray_ITS_D.txt', header = T, row.names = 1)
RCbray_ITS_D0[upper.tri(RCbray_ITS_D0)] <- NA
for(i in 1:80){
  RCbray_ITS_D0[i,i] <- NA
}
RCbray_ITS_D <- RCbray_ITS_D0 %>%
  rownames_to_column() %>%
  gather(ID2, Value, -rowname) %>%
  rename(RC = Value) %>%
  rename(R2 = rowname) %>%
  drop_na()

RCbray_ITS_I0 <- read.table('RCbray_ITS_I.txt', header = T,row.names = 1)
RCbray_ITS_I0[upper.tri(RCbray_ITS_I0)] <- NA
for(i in 1:80){
  RCbray_ITS_I0[i,i] <- NA
}
RCbray_ITS_I <- RCbray_ITS_I0 %>%
  rownames_to_column() %>%
  gather(ID2, Value, -rowname) %>%
  rename(RC = Value) %>%
  rename(R2 = rowname) %>%
  drop_na()

RCbray_ITS_J0 <- read.table('RCbray_ITS_J.txt', header = T,row.names = 1)
RCbray_ITS_J0[upper.tri(RCbray_ITS_J0)] <- NA
for(i in 1:80){
  RCbray_ITS_J0[i,i] <- NA
}
RCbray_ITS_J <- RCbray_ITS_J0 %>%
  rownames_to_column() %>%
  gather(ID2, Value, -rowname) %>%
  rename(RC = Value) %>%
  rename(R2 = rowname) %>%
  drop_na()


ls()
head(betaNTI_16S_C)
head(betaNTI_ITS_C)
head(RCbray_16S_C)
head(RCbray_ITS_C)
NTI_RC_16S_C <- cbind(betaNTI_16S_C,RCbray_16S_C) %>%
  select(-c('R2','ID2'))
head(NTI_RC_16S_C)
NTI_RC_16S_D <- cbind(betaNTI_16S_D,RCbray_16S_D) %>%
  select(-c('R2','ID2'))
NTI_RC_16S_I <- cbind(betaNTI_16S_I,RCbray_16S_I) %>%
  select(-c('R2','ID2'))
NTI_RC_16S_J <- cbind(betaNTI_16S_J,RCbray_16S_J) %>%
  select(-c('R2','ID2'))
NTI_RC_ITS_C <- cbind(betaNTI_ITS_C,RCbray_ITS_C) %>%
  select(-c('R2','ID2'))
NTI_RC_ITS_D <- cbind(betaNTI_ITS_D,RCbray_ITS_D) %>%
  select(-c('R2','ID2'))
NTI_RC_ITS_I <- cbind(betaNTI_ITS_I,RCbray_ITS_I) %>%
  select(-c('R2','ID2'))
NTI_RC_ITS_J <- cbind(betaNTI_ITS_J,RCbray_ITS_J) %>%
  select(-c('R2','ID2'))
ls()

ls()
# 16S C
for(i in 1:nrow(NTI_RC_16S_C)){
  if(NTI_RC_16S_C$NTI[i]>=2){
    NTI_RC_16S_C$Cat[i] <- 'Homogeneous_selection'
  }
  else if(NTI_RC_16S_C$NTI[i]<= -2){
    NTI_RC_16S_C$Cat[i] <- 'Heterogenous_selection'
  }
  else if(NTI_RC_16S_C$RC[i] >= 0.95){
    NTI_RC_16S_C$Cat[i] <- 'Disperal_limitation'
  }
  else if(NTI_RC_16S_C$RC[i] <= -0.95){
    NTI_RC_16S_C$Cat[i] <- 'Homogenizing_disperal'
  }
  else{NTI_RC_16S_C$Cat[i] <- 'Undominated'}
}
# 16S D
for(i in 1:nrow(NTI_RC_16S_D)){
  if(NTI_RC_16S_D$NTI[i]>=2){
    NTI_RC_16S_D$Cat[i] <- 'Homogeneous_selection'
  }
  else if(NTI_RC_16S_D$NTI[i]<= -2){
    NTI_RC_16S_D$Cat[i] <- 'Heterogenous_selection'
  }
  else if(NTI_RC_16S_D$RC[i] >= 0.95){
    NTI_RC_16S_D$Cat[i] <- 'Disperal_limitation'
  }
  else if(NTI_RC_16S_D$RC[i] <= -0.95){
    NTI_RC_16S_D$Cat[i] <- 'Homogenizing_disperal'
  }
  else{NTI_RC_16S_D$Cat[i] <- 'Undominated'}
}
# 16S I
for(i in 1:nrow(NTI_RC_16S_I)){
  if(NTI_RC_16S_I$NTI[i]>=2){
    NTI_RC_16S_I$Cat[i] <- 'Homogeneous_selection'
  }
  else if(NTI_RC_16S_I$NTI[i]<= -2){
    NTI_RC_16S_I$Cat[i] <- 'Heterogenous_selection'
  }
  else if(NTI_RC_16S_I$RC[i] >= 0.95){
    NTI_RC_16S_I$Cat[i] <- 'Disperal_limitation'
  }
  else if(NTI_RC_16S_I$RC[i] <= -0.95){
    NTI_RC_16S_I$Cat[i] <- 'Homogenizing_disperal'
  }
  else{NTI_RC_16S_I$Cat[i] <- 'Undominated'}
}
# 16S J
for(i in 1:nrow(NTI_RC_16S_J)){
  if(NTI_RC_16S_J$NTI[i]>=2){
    NTI_RC_16S_J$Cat[i] <- 'Homogeneous_selection'
  }
  else if(NTI_RC_16S_J$NTI[i]<= -2){
    NTI_RC_16S_J$Cat[i] <- 'Heterogenous_selection'
  }
  else if(NTI_RC_16S_J$RC[i] >= 0.95){
    NTI_RC_16S_J$Cat[i] <- 'Disperal_limitation'
  }
  else if(NTI_RC_16S_J$RC[i] <= -0.95){
    NTI_RC_16S_J$Cat[i] <- 'Homogenizing_disperal'
  }
  else{NTI_RC_16S_J$Cat[i] <- 'Undominated'}
}
# ITS C
for(i in 1:nrow(NTI_RC_ITS_C)){
  if(NTI_RC_ITS_C$NTI[i]>=2){
    NTI_RC_ITS_C$Cat[i] <- 'Homogeneous_selection'
  }
  else if(NTI_RC_ITS_C$NTI[i]<= -2){
    NTI_RC_ITS_C$Cat[i] <- 'Heterogenous_selection'
  }
  else if(NTI_RC_ITS_C$RC[i] >= 0.95){
    NTI_RC_ITS_C$Cat[i] <- 'Disperal_limitation'
  }
  else if(NTI_RC_ITS_C$RC[i] <= -0.95){
    NTI_RC_ITS_C$Cat[i] <- 'Homogenizing_disperal'
  }
  else{NTI_RC_ITS_C$Cat[i] <- 'Undominated'}
}
# ITS D
for(i in 1:nrow(NTI_RC_ITS_D)){
  if(NTI_RC_ITS_D$NTI[i]>=2){
    NTI_RC_ITS_D$Cat[i] <- 'Homogeneous_selection'
  }
  else if(NTI_RC_ITS_D$NTI[i]<= -2){
    NTI_RC_ITS_D$Cat[i] <- 'Heterogenous_selection'
  }
  else if(NTI_RC_ITS_D$RC[i] >= 0.95){
    NTI_RC_ITS_D$Cat[i] <- 'Disperal_limitation'
  }
  else if(NTI_RC_ITS_D$RC[i] <= -0.95){
    NTI_RC_ITS_D$Cat[i] <- 'Homogenizing_disperal'
  }
  else{NTI_RC_ITS_D$Cat[i] <- 'Undominated'}
}
# ITS I
for(i in 1:nrow(NTI_RC_ITS_I)){
  if(NTI_RC_ITS_I$NTI[i]>=2){
    NTI_RC_ITS_I$Cat[i] <- 'Homogeneous_selection'
  }
  else if(NTI_RC_ITS_I$NTI[i]<= -2){
    NTI_RC_ITS_I$Cat[i] <- 'Heterogenous_selection'
  }
  else if(NTI_RC_ITS_I$RC[i] >= 0.95){
    NTI_RC_ITS_I$Cat[i] <- 'Disperal_limitation'
  }
  else if(NTI_RC_ITS_I$RC[i] <= -0.95){
    NTI_RC_ITS_I$Cat[i] <- 'Homogenizing_disperal'
  }
  else{NTI_RC_ITS_I$Cat[i] <- 'Undominated'}
}
# ITS J
for(i in 1:nrow(NTI_RC_ITS_J)){
  if(NTI_RC_ITS_J$NTI[i]>=2){
    NTI_RC_ITS_J$Cat[i] <- 'Homogeneous_selection'
  }
  else if(NTI_RC_ITS_J$NTI[i]<= -2){
    NTI_RC_ITS_J$Cat[i] <- 'Heterogenous_selection'
  }
  else if(NTI_RC_ITS_J$RC[i] >= 0.95){
    NTI_RC_ITS_J$Cat[i] <- 'Disperal_limitation'
  }
  else if(NTI_RC_ITS_J$RC[i] <= -0.95){
    NTI_RC_ITS_J$Cat[i] <- 'Homogenizing_disperal'
  }
  else{NTI_RC_ITS_J$Cat[i] <- 'Undominated'}
}

NTI_RC_16S_C1 <- NTI_RC_16S_C %>%
  mutate(ID = paste0(substr(R1,1,3),substr(ID1,1,3))) %>%
  filter(ID %in% c('D01D01','D04D04','D07D07','D14D14','D21D21','D28D28','D35D35','D42D42')) %>%
  mutate(Category = gsub("_", " ", Cat)) %>%
  mutate(ID_mask = paste0('ID',1:360)) %>%
  mutate(Date = substr(ID,1,3)) %>%
  select(c(Date,ID_mask,Category,NTI,RC))
NTI_RC_16S_D1 <- NTI_RC_16S_D %>%
  mutate(ID = paste0(substr(R1,1,3),substr(ID1,1,3))) %>%
  filter(ID %in% c('D01D01','D04D04','D07D07','D14D14','D21D21','D28D28','D35D35','D42D42')) %>%
  mutate(Category = gsub("_", " ", Cat)) %>%
  mutate(ID_mask = paste0('ID',1:360)) %>%
  mutate(Date = substr(ID,1,3)) %>%
  select(c(Date,ID_mask,Category,NTI,RC))
NTI_RC_16S_I1 <- NTI_RC_16S_I %>%
  mutate(ID = paste0(substr(R1,1,3),substr(ID1,1,3))) %>%
  filter(ID %in% c('D01D01','D04D04','D07D07','D14D14','D21D21','D28D28','D35D35','D42D42')) %>%
  mutate(Category = gsub("_", " ", Cat)) %>%
  mutate(ID_mask = paste0('ID',1:360)) %>%
  mutate(Date = substr(ID,1,3)) %>%
  select(c(Date,ID_mask,Category,NTI,RC))
NTI_RC_16S_J1 <- NTI_RC_16S_J %>%
  mutate(ID = paste0(substr(R1,1,3),substr(ID1,1,3))) %>%
  filter(ID %in% c('D01D01','D04D04','D07D07','D14D14','D21D21','D28D28','D35D35','D42D42')) %>%
  mutate(Category = gsub("_", " ", Cat)) %>%
  mutate(ID_mask = paste0('ID',1:360)) %>%
  mutate(Date = substr(ID,1,3)) %>%
  select(c(Date,ID_mask,Category,NTI,RC))

NTI_RC_16S_C2 <- table(NTI_RC_16S_C1[,c('Date','Category')]) %>%
  melt() %>%
  mutate(Segment = 'Cecum') %>%
  mutate(Ratio = value/45 * 100)

NTI_RC_16S_D2 <- table(NTI_RC_16S_D1[,c('Date','Category')]) %>%
  melt() %>%
  mutate(Segment = 'Duodenum') %>%
  mutate(Ratio = value/45 * 100)

NTI_RC_16S_I2 <- table(NTI_RC_16S_I1[,c('Date','Category')]) %>%
  melt() %>%
  mutate(Segment = 'Ileum') %>%
  mutate(Ratio = value/45 * 100)
NTI_RC_16S_J2 <- table(NTI_RC_16S_J1[,c('Date','Category')]) %>%
  melt()%>%
  mutate(Segment = 'Jejunum') %>%
  mutate(Ratio = value/45 * 100)

NTI_RC_16S <- rbind(NTI_RC_16S_D2, NTI_RC_16S_J2, NTI_RC_16S_I2, NTI_RC_16S_C2)
head(NTI_RC_16S)
NTI_RC_16S$Segment <- factor(NTI_RC_16S$Segment, levels = c('Duodenum','Jejunum','Ileum','Cecum'))
NTI_RC_16S$Category <- factor(NTI_RC_16S$Category, levels = rev(c('Undominated','Homogeneous selection','Heterogenous selection',
                                                                  'Disperal limitation','Homogenizing disperal')))
NTI_RC_bar1 <- ggplot(NTI_RC_16S, aes(x = Date, y = Ratio, fill = Category)) +
  geom_bar(stat = 'identity',position = 'stack') +
  facet_grid(.~Segment) +
  xlab('DPH') +
  scale_fill_manual(values = rev(c('#E64B35','#4DBBD5','#3C5488','#F39B7F'))) +
  theme_bw() +
  ylab("Ratio (%)") +
  scale_x_discrete(label = c(1,4,7,14,21,28,35,42)) +
  theme(axis.text.x = element_text(angle = 90,hjust = 1)) +
  ggtitle('Bacteria')
NTI_RC_ITS_C1 <- NTI_RC_ITS_C %>%
  mutate(ID = paste0(substr(R1,1,3),substr(ID1,1,3))) %>%
  filter(ID %in% c('D01D01','D04D04','D07D07','D14D14','D21D21','D28D28','D35D35','D42D42')) %>%
  mutate(Category = gsub("_", " ", Cat)) %>%
  mutate(ID_mask = paste0('ID',1:360)) %>%
  mutate(Date = substr(ID,1,3)) %>%
  select(c(Date,ID_mask,Category,NTI,RC))
NTI_RC_ITS_D1 <- NTI_RC_ITS_D %>%
  mutate(ID = paste0(substr(R1,1,3),substr(ID1,1,3))) %>%
  filter(ID %in% c('D01D01','D04D04','D07D07','D14D14','D21D21','D28D28','D35D35','D42D42')) %>%
  mutate(Category = gsub("_", " ", Cat)) %>%
  mutate(ID_mask = paste0('ID',1:360)) %>%
  mutate(Date = substr(ID,1,3)) %>%
  select(c(Date,ID_mask,Category,NTI,RC))
NTI_RC_ITS_I1 <- NTI_RC_ITS_I %>%
  mutate(ID = paste0(substr(R1,1,3),substr(ID1,1,3))) %>%
  filter(ID %in% c('D01D01','D04D04','D07D07','D14D14','D21D21','D28D28','D35D35','D42D42')) %>%
  mutate(Category = gsub("_", " ", Cat)) %>%
  mutate(ID_mask = paste0('ID',1:360)) %>%
  mutate(Date = substr(ID,1,3)) %>%
  select(c(Date,ID_mask,Category,NTI,RC))
NTI_RC_ITS_J1 <- NTI_RC_ITS_J %>%
  mutate(ID = paste0(substr(R1,1,3),substr(ID1,1,3))) %>%
  filter(ID %in% c('D01D01','D04D04','D07D07','D14D14','D21D21','D28D28','D35D35','D42D42')) %>%
  mutate(Category = gsub("_", " ", Cat)) %>%
  mutate(ID_mask = paste0('ID',1:360)) %>%
  mutate(Date = substr(ID,1,3)) %>%
  select(c(Date,ID_mask,Category,NTI,RC))

NTI_RC_ITS_C2 <- table(NTI_RC_ITS_C1[,c('Date','Category')]) %>%
  melt() %>%
  mutate(Segment = 'Cecum') %>%
  mutate(Ratio = value/45 * 100)

NTI_RC_ITS_D2 <- table(NTI_RC_ITS_D1[,c('Date','Category')]) %>%
  melt() %>%
  mutate(Segment = 'Duodenum') %>%
  mutate(Ratio = value/45 * 100)

NTI_RC_ITS_I2 <- table(NTI_RC_ITS_I1[,c('Date','Category')]) %>%
  melt() %>%
  mutate(Segment = 'Ileum') %>%
  mutate(Ratio = value/45 * 100)
NTI_RC_ITS_J2 <- table(NTI_RC_ITS_J1[,c('Date','Category')]) %>%
  melt()%>%
  mutate(Segment = 'Jejunum') %>%
  mutate(Ratio = value/45 * 100)

NTI_RC_ITS <- rbind(NTI_RC_ITS_D2, NTI_RC_ITS_J2, NTI_RC_ITS_I2, NTI_RC_ITS_C2)
head(NTI_RC_ITS)
NTI_RC_ITS$Segment <- factor(NTI_RC_ITS$Segment, levels = c('Duodenum','Jejunum','Ileum','Cecum'))
NTI_RC_ITS$Category <- factor(NTI_RC_ITS$Category, levels = rev(c('Undominated','Homogeneous selection','Heterogenous selection',
                                                                  'Disperal limitation','Homogenizing disperal')))
NTI_RC_bar2 <- ggplot(NTI_RC_ITS, aes(x = Date, y = Ratio, fill = Category)) +
  geom_bar(stat = 'identity',position = 'stack') +
  facet_grid(.~Segment) +
  xlab('DPH') +
  scale_fill_manual(values = rev(c('#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F'))) +
  theme_bw() +
  ylab("Ratio (%)") +
  scale_x_discrete(label = c(1,4,7,14,21,28,35,42)) +
  theme(axis.text.x = element_text(angle = 90,hjust = 1)) +
  ggtitle('Fungi')
library(ggpubr)
NTI_RC_16S_RA <- NTI_RC_16S %>%
  filter(Category %in% 'Undominated') %>%
  group_by(Segment)%>%
  summarise_at(vars(Ratio), funs(mean(., na.rm=TRUE))) %>%
  data.frame() %>%
  mutate(Sequence = '16S')
NTI_RC_ITS_RA <- NTI_RC_ITS %>%
  filter(Category %in% 'Undominated') %>%
  group_by(Segment)%>%
  summarise_at(vars(Ratio), funs(mean(., na.rm=TRUE))) %>%
  data.frame() %>%
  mutate(Sequence = 'ITS')
rbind(NTI_RC_16S_RA, NTI_RC_ITS_RA) %>%
  ggplot(aes(x = Segment, y = Ratio, fill = Sequence)) +
  geom_bar(stat = 'identity',position = 'dodge')+
  scale_fill_manual(values = c('#426E9D','#DD4D4F')) +
  theme_bw() +
  theme(legend.position = c(0.2,0.2)) +
  ylab('Ratio (%)')


NTI_RC_bar1 <- ggplot(NTI_RC_16S, aes(x = Date, y = Ratio, fill = Category)) +
  geom_bar(stat = 'identity',position = 'stack') +
  facet_grid(.~Segment) +
  xlab('DPH') +
  scale_fill_manual(values = rev(c('#E64B35','#4DBBD5','#3C5488','#F39B7F'))) +
  theme_bw() +
  ylab("Ratio (%)") +
  scale_x_discrete(label = c(1,4,7,14,21,28,35,42)) +
  theme(axis.text.x = element_text(angle = 90,hjust = 1)) +
  ggtitle('Bacteria')
NTI_RC_bar2 <- ggplot(NTI_RC_ITS, aes(x = Date, y = Ratio, fill = Category)) +
  geom_bar(stat = 'identity',position = 'stack') +
  facet_grid(.~Segment) +
  xlab('DPH') +
  scale_fill_manual(values = rev(c('#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F'))) +
  theme_bw() +
  ylab("Ratio (%)") +
  scale_x_discrete(label = c(1,4,7,14,21,28,35,42)) +
  theme(axis.text.x = element_text(angle = 90,hjust = 1)) +
  ggtitle('Fungi')

f14_1 / f14_2 | f14_11 
NTI_RC_bar1 + NTI_RC_bar2 + plot_layout(guides = 'collect')
