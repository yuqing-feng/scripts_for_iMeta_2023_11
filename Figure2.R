# Maturaty of the gut microbiota 
# Part2: Comparision of the results from relative abundance and absosulte abundance
# 2.1 dual y axis to present the tendency of them.
library(openxlsx)
library(dplyr)
library(tibble)
library(reshape2)
library(ggplot2)
library(patchwork)
ab_temp <- read.xlsx('01data_for_quantitive.xlsx',sheet = 11)
raw_16S <- read.csv('level-2.csv')
format_16S <- raw_16S %>%
  remove_rownames %>% 
  column_to_rownames(var="index") %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  select(-c('D_0__Archaea.D_1__Euryarchaeota')) %>%
  mutate(across(D_0__Bacteria.D_1__Acidobacteria:D_0__Bacteria.D_1__Zixibacteria, ~ . / sum)) %>%
  select(-c('sum')) %>%
  setNames(gsub("D_0__Bacteria.D_1__", "", names(.))) %>%
  mutate(DPH = rep(c(1,4,7,14,21,28,35,42),each = 40)) %>%
  mutate(Group1 = factor(Group1, levels = c('Duodenum','Jejunum','Ileum','Cecum'))) %>%
  select(c('Firmicutes','DPH','Group1')) 
format_16S

raw_ITS <- read.csv('level-2.csv')
format_ITS <- raw_ITS %>%
  remove_rownames %>% 
  column_to_rownames(var="index") %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  select(-c('k__Fungi.p__unidentified')) %>%
  mutate(across(k__Fungi.p__Ascomycota:k__Fungi.p__Rozellomycota, ~ . / sum)) %>%
  select(-c('sum')) %>%
  setNames(gsub("k__Fungi.p__", "", names(.))) %>%
  mutate(DPH = rep(c(1,4,7,14,21,28,35,42),each = 40)) %>%
  mutate(Group1 = factor(Group1, levels = c('Duodenum','Jejunum','Ileum','Cecum'))) %>%
  select(c('Ascomycota','DPH','Group1')) 
format_ITS

for(i in 1:nrow(ab_temp)){
  format_16S$Firmicutes2[i] <- format_16S$Firmicutes[i] * 10^(ab_temp$`log10(16S)`[i])
  format_ITS$Ascomycota2[i] <- format_ITS$Ascomycota[i] * 10^(ab_temp$`log10(ITS)`[i])
}
head(format_16S)


f12_1 <- ggplot(format_16S, aes(x=DPH)) +
  geom_smooth(aes(y=Firmicutes),color = '#4DBBD5') +
  geom_smooth(aes(y = log10(Firmicutes2+1)/10), color = '#E64B35') + 
  facet_wrap(.~Group1,nrow = 1) +
  scale_y_continuous(name = "Relative abundance", sec.axis = sec_axis(~.*10, name="Number of copies (log10)"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
        panel.grid.minor = element_line(color = NA)) +
  ggtitle('Firmicutes') +
  scale_x_continuous(breaks = c(1,7,14,21,28,35,42)) 
f12_1

f12_2 <- ggplot(format_ITS, aes(x=DPH)) +
  geom_smooth(aes(y=Ascomycota),color = '#4DBBD5') +
  geom_smooth(aes(y = log10(Ascomycota2+1)/10), color = '#E64B35') + 
  facet_wrap(.~Group1,nrow = 1) +
  scale_y_continuous(
    name = "Relative abundance", 
    sec.axis = sec_axis(~.*10, name="Number of copies (log10)"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
        panel.grid.minor = element_line(color = NA),
        legend.position = 'none') +
  ggtitle('Ascomycota') +
  scale_x_continuous(breaks = c(1,7,14,21,28,35,42)) 
f12_2

f12_1 + f12_2 + plot_layout(nrow = 2)



raw_16S_L6 <- read.csv('level-6_no_archaea.csv')
format_16S_L6 <- raw_16S_L6 %>%
  remove_rownames %>% 
  column_to_rownames(var="index") %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  mutate(across(Granulicella:uncultured.bacterium.40, ~ . / sum)) %>%
  select(-c('sum')) %>%
  mutate(DPH = rep(c(1,4,7,14,21,28,35,42), each = 40)) %>%
  mutate(Group1 = factor(Group1, levels = c('Duodenum','Jejunum','Ileum','Cecum'))) %>%
  select(c('Lactobacillus','DPH','Group1'))
format_16S_L6

raw_ITS_L6 <- read.csv('level-6_for_bar.csv')
format_ITS_L6 <- raw_ITS_L6 %>%
  remove_rownames %>% 
  column_to_rownames(var="X") %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  mutate(across(X__:unidentified.53, ~ . / sum)) %>%
  select(-c('sum')) %>%
  mutate(DPH = rep(c(1,4,7,14,21,28,35,42), each = 40)) %>%
  mutate(Group1 = factor(Group1, levels = c('Duodenum','Jejunum','Ileum','Cecum'))) %>%
  select(c('Candida','DPH','Group1'))
format_ITS_L6

for(i in 1:nrow(ab_temp)){
  format_16S_L6$Lactobacillus2[i] <- format_16S_L6$Lactobacillus[i] * 10^(ab_temp$`log10(16S)`[i])
  format_ITS_L6$Candida2[i] <- format_ITS_L6$Candida[i] * 10^(ab_temp$`log10(ITS)`[i])
}

f12_3 <- ggplot(format_16S_L6, aes(x=DPH)) +
  geom_smooth(aes(y=Lactobacillus),color = '#4DBBD5') +
  geom_smooth(aes(y = log10(Lactobacillus2+1)/10), color = '#E64B35') + 
  facet_wrap(.~Group1,nrow = 1) +
  scale_y_continuous(name = "Relative abundance", sec.axis = sec_axis(~.*10, name="Number of copies (log10)"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
        panel.grid.minor = element_line(color = NA)) +
  ggtitle('Lactobacillus') +
  scale_x_continuous(breaks = c(1,7,14,21,28,35,42)) 
f12_3
f12_4 <- ggplot(format_ITS_L6, aes(x=DPH)) +
  geom_smooth(aes(y=Candida),color = '#4DBBD5') +
  geom_smooth(aes(y = log10(Candida2+1)/10), color = '#E64B35') + 
  facet_wrap(.~Group1,nrow = 1) +
  scale_y_continuous(name = "Relative abundance", sec.axis = sec_axis(~.*10, name="Number of copies (log10)"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
        panel.grid.minor = element_line(color = NA)) +
  ggtitle('Candida') +
  scale_x_continuous(breaks = c(1,7,14,21,28,35,42)) 
f12_4
svg('~/Desktop/figure materials/2nd/figure 2a.svg',width = 10.5,height = 5)
f12_1 + f12_2 + f12_3 + f12_4 + plot_layout(nrow = 2)
dev.off()
# part2 PA analysis and the others

ab_temp
raw_16S_L6 <- read.csv('level-6_no_archaea.csv')
raw_ITS_L6 <- read.csv('level-6_for_bar.csv')
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
library(vegan)
# first PA analysis based on 16S RA and AA
dist_16S_RA <- vegdist(format_16S_L6[,1:529])
dist_ITS_RA <- vegdist(format_ITS_L6[,1:603])
dist_16S_AA <- vegdist(format_16S_L6_AA[,1:529])
dist_ITS_AA <- vegdist(format_ITS_L6_AA[,1:603])
PA_16S <- procrustes(dist_16S_AA,dist_16S_RA)
PA_ITS <- procrustes(dist_ITS_AA,dist_ITS_RA)
residual_16S <- data.frame(residuals(PA_16S))
residual_ITS <- data.frame(residuals(PA_ITS))

PA_16S_XY <- data.frame(AA_x = PA_16S$X[,1], AA_y = PA_16S$X[,2],
                        RA_x = PA_16S$Yrot[,1], RA_y = PA_16S$Yrot[,2])
PA_ITS_XY <- data.frame(AA_x = PA_ITS$X[,1], AA_y = PA_ITS$X[,2],
                        RA_x = PA_ITS$Yrot[,1], RA_y = PA_ITS$Yrot[,2])

protest(dist_16S_AA,dist_16S_RA)
protest(dist_ITS_AA,dist_ITS_RA)

f13_1 <- ggplot(PA_16S_XY) +
  geom_segment(aes(x = RA_x, y = RA_y, xend = AA_x, yend = AA_y, color = 'gray'), alpha = 0.3) +
  geom_point(aes(x = RA_x, y = RA_y, color = '#4DBBD5')) +
  geom_point(aes(x = AA_x, y = AA_y, color = '#E64B35')) +
  theme_bw() +
  xlab('PCoA 1') + 
  ylab('PCoA 2') +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor = element_line(colour=NA),
        legend.position = c(0.2,0.3)) +
  annotate('text', x = -0.7, y = -0.25, label = 'P = 0.001\nCorrelation = 0.622') +
  scale_color_manual(name = 'Method',values = c('#4DBBD5','#E64B35','gray')) + 
  ggtitle('Bacteria') 
f13_1
f13_2 <- ggplot(PA_ITS_XY) +
  geom_segment(aes(x = RA_x, y = RA_y, xend = AA_x, yend = AA_y, color = 'gray'), alpha = 0.3) +
  geom_point(aes(x = RA_x, y = RA_y, color = '#4DBBD5')) +
  geom_point(aes(x = AA_x, y = AA_y, color = '#E64B35')) +
  theme_bw() +
  xlab('PCoA 1') + 
  ylab('PCoA 2') +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor = element_line(colour=NA),
        legend.position = c(0.2,0.3))  +
  scale_color_manual(name = 'Method',values = c('#4DBBD5','#E64B35','gray')) +
  annotate('text', x = -0.6, y = 0, label = 'P = 0.001\nCorrelation = 0.459') +
  ggtitle('Fungi')
f13_2

dph <- rep(c(1,4,7,14,21,28,35,42), each = 40)
set.seed(1234)
adonis(format_16S_L6[,1:529]~ format_16S_L6$Group1 + dph )
set.seed(1234)
adonis(format_16S_L6_AA[,1:529]~ format_16S_L6_AA$Group1 + dph)
set.seed(1234)
adonis(format_ITS_L6[,1:603] ~ format_ITS_L6$Group1 + dph)
set.seed(1234)
adonis(format_ITS_L6_AA[,1:603] ~ format_ITS_L6_AA$Group1 + dph)

explains = data.frame(Method = c('RA','AA','RA','AA'),
                      Individual = c(0.31110,0.21089,0.29979,0.24770) * 100,
                      Segment = c(0.32671,0.27479,0.18011,0.11674) * 100,
                      DPH = c(0.08565,0.03220,0.07234,0.02779) * 100,
                      Kingdom = c('Bacteria','Bacteria','Fungi','Fungi')) 
explains
library(reshape2)
exp_melt <- melt(explains)
colnames(exp_melt)[3:4] <- c('Factor','Explanation') 
exp_melt$Factor <- factor(exp_melt$Factor, levels = c('Individual', 'Segment', 'DPH'))
f13_3 <- ggplot(exp_melt, aes(x = Factor, y = Explanation, fill = Method)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = c('#E64B35','#4DBBD5')) +
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor = element_line(colour=NA),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ylab('Explanation (%)') +
  facet_wrap(.~Kingdom)
f13_3
svg('~/Desktop/figure materials/2nd/figure 2b.svg',width = 10.5,height = 4)
f13_1 + f13_2 + f13_3
dev.off()

#########
# Figure 2 D
Genus_16S_RA_ID <- read.table('~/Desktop/01maturity/material/02micro/output/L6_RA_16S.csv',header = T,sep = ',',row.names = 1)%>%
  colSums()
Genus_16S_RA0 <- read.table('~/Desktop/01maturity/material/02micro/output/L6_RA_16S.csv',header = T,sep = ',',row.names = 1) 
head(Genus_16S_RA0)
Genus_16S_RA <- Genus_16S_RA0[,Genus_16S_RA_ID>0.001*280]
Genus_16S_AA0 <- read.table('~/Desktop/01maturity/material/02micro/output/L6_AA_16S.csv',header = T,sep = ',',row.names = 1)
Genus_16S_AA <- Genus_16S_AA0[,Genus_16S_RA_ID>0.001*280]
rownames(Genus_16S_AA) == rownames(Genus_16S_RA)
colnames(Genus_16S_AA) == colnames(Genus_16S_RA)
rownames(Genus_16S_AA)
Genus_16S_AA$Seg <- c(rep(c('Cecum','Duodenum','Ileum','Jejunum'),80))
Genus_16S_RA$Seg <- c(rep(c('Cecum','Duodenum','Ileum','Jejunum'),80))
rm(Genus_16S_RA_ID)
rm(Genus_16S_AA0)
rm(Genus_16S_RA0)
dim(Genus_16S_AA)
taxa_16S <- colnames(Genus_16S_RA)[1:58]
colnames(Genus_16S_AA)[1:58] <- colnames(Genus_16S_RA)[1:58] <- paste0('Genus',1:58)

Genus_16S_RA_C <- Genus_16S_RA[Genus_16S_RA$Seg=='Cecum',] %>%
  select(-c('Seg'))
Genus_16S_RA_D <- Genus_16S_RA[Genus_16S_RA$Seg=='Duodenum',] %>%
  select(-c('Seg'))
Genus_16S_RA_I <- Genus_16S_RA[Genus_16S_RA$Seg=='Ileum',] %>%
  select(-c('Seg'))
Genus_16S_RA_J <- Genus_16S_RA[Genus_16S_RA$Seg=='Jejunum',] %>%
  select(-c('Seg'))
Genus_16S_AA_C <- Genus_16S_AA[Genus_16S_AA$Seg=='Cecum',] %>%
  select(-c('Seg'))
Genus_16S_AA_D <- Genus_16S_AA[Genus_16S_AA$Seg=='Duodenum',] %>%
  select(-c('Seg'))
Genus_16S_AA_I <- Genus_16S_AA[Genus_16S_AA$Seg=='Ileum',] %>%
  select(-c('Seg'))
Genus_16S_AA_J <- Genus_16S_AA[Genus_16S_AA$Seg=='Jejunum',] %>%
  select(-c('Seg'))

ab_2_cor <- function(x_RA, y_AA){
  p_RA <- matrix(, nrow = ncol(x_RA), ncol = ncol(x_RA))
  cor_RA <- matrix(, nrow = ncol(x_RA), ncol = ncol(x_RA))
  for(i in 1:ncol(x_RA)){
    for(j in 1:ncol(x_RA)){
      p_RA[i,j] <- cor.test(x_RA[,i], x_RA[,j], method = 'spearman')$p.value
      cor_RA[i,j] <- cor.test(x_RA[,i], x_RA[,j], method = 'spearman')$estimate
    }
  }
  cor_RA[is.na(cor_RA)] <- 0
  colnames(p_RA) <- rownames(p_RA) <- colnames(cor_RA) <- rownames(cor_RA) <- paste0('G',1:ncol(x_RA))
  
  p_AA <- matrix(, nrow = ncol(y_AA), ncol = ncol(y_AA))
  cor_AA <- matrix(, nrow = ncol(y_AA), ncol = ncol(y_AA))
  for(i in 1:ncol(y_AA)){
    for(j in 1:ncol(y_AA)){
      p_AA[i,j] <- cor.test(y_AA[,i], y_AA[,j], method = 'spearman')$p.value
      cor_AA[i,j] <- cor.test(y_AA[,i], y_AA[,j], method = 'spearman')$estimate
    }
  }
  colnames(p_AA) <- rownames(p_AA) <- colnames(cor_AA) <- rownames(cor_AA) <- paste0('G',1:ncol(y_AA))
  
  # calculate the 
  ind_RA <- which(lower.tri(cor_RA, diag = F), arr.ind = TRUE)
  Col_RA = dimnames(cor_RA)[[2]][ind_RA[,2]]
  Row_RA = dimnames(cor_RA)[[1]][ind_RA[,1]]
  Cor_RA = cor_RA[ind_RA]
  Pva_RA = p_RA[ind_RA]
  Q_RA = p.adjust(Pva_RA, method = 'BH')
  Cor_RA[Q_RA >= 0.05] <- 0
  out_RA <- data.frame(Col_RA,Row_RA,Cor_RA)
  
  ind_AA <- which(upper.tri(cor_AA, diag = F), arr.ind = TRUE)
  Col_AA = dimnames(cor_AA)[[2]][ind_AA[,2]]
  Row_AA = dimnames(cor_AA)[[1]][ind_AA[,1]]
  Cor_AA = cor_AA[ind_AA]
  Pva_AA = p_AA[ind_AA]
  Q_AA = p.adjust(Pva_AA, method = 'BH')
  Cor_AA[Q_AA >= 0.05] <- 0
  out_AA <- data.frame(Col_AA,Row_AA,Cor_AA)
  colnames(out_RA) <- colnames(out_AA) <- c('Col','Row','Cor')
  out_RA$Cat <- 'RA'
  out_AA$Cat <- 'AA'
  return(rbind(out_RA, out_AA))
}
library(ggplot2)
Cecum <- ab_2_cor(Genus_16S_RA_C, Genus_16S_AA_C)
Cecum$Row <- factor(Cecum$Row, levels = paste0('G',1:58))
Cecum$Col <- factor(Cecum$Col, levels = paste0('G',1:58))
f14_C <- ggplot(Cecum, aes(x = Col, y = Row, color = Cor)) +
  geom_point(aes(size = abs(Cor))) +
  theme_bw() +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
        panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_color_gradient2(low="#3C5488", mid="white", high="#E64B35", midpoint=0)  +
  scale_size(range=c(0,2),
             breaks=c(0,0.25,0.5,0.75,1),
             labels=c("0","0.25","0.50","0.75","1.00"),
             name = "",
             guide="legend") +
  guides(size = 'none') +
  xlab('Absolute abundance') +
  ylab('Relative abundance') +
  ggtitle('Cecum')
f14_C
Ileum <- ab_2_cor(Genus_16S_RA_I, Genus_16S_AA_I)
Ileum$Row <- factor(Ileum$Row, levels = paste0('G',1:58))
Ileum$Col <- factor(Ileum$Col, levels = paste0('G',1:58))
f14_I <- ggplot(Ileum, aes(x = Col, y = Row, color = Cor)) +
  geom_point(aes(size = abs(Cor))) +
  theme_bw() +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
        panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_color_gradient2(low="#3C5488", mid="white", high="#E64B35", midpoint=0)  +
  scale_size(range=c(0,2),
             breaks=c(0,0.25,0.5,0.75,1),
             labels=c("0","0.25","0.50","0.75","1.00"),
             name = "",
             guide="legend") +
  guides(size = 'none') +
  xlab('Absolute abundance') +
  ylab('Relative abundance') +
  ggtitle('Ileum')
f14_I
Duodenum <- ab_2_cor(Genus_16S_RA_D, Genus_16S_AA_D)
Duodenum$Row <- factor(Duodenum$Row, levels = paste0('G',1:58))
Duodenum$Col <- factor(Duodenum$Col, levels = paste0('G',1:58))
f14_D <- ggplot(Duodenum, aes(x = Col, y = Row, color = Cor)) +
  geom_point(aes(size = abs(Cor))) +
  theme_bw() +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
        panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_color_gradient2(low="#3C5488", mid="white", high="#E64B35", midpoint=0)  +
  scale_size(range=c(0,2),
             breaks=c(0,0.25,0.5,0.75,1),
             labels=c("0","0.25","0.50","0.75","1.00"),
             name = "",
             guide="legend") +
  guides(size = 'none') +
  xlab('Absolute abundance') +
  ylab('Relative abundance') +
  ggtitle('Duodenum')

Jejunum <- ab_2_cor(Genus_16S_RA_J, Genus_16S_AA_J)
Jejunum$Row <- factor(Jejunum$Row, levels = paste0('G',1:58))
Jejunum$Col <- factor(Jejunum$Col, levels = paste0('G',1:58))
f14_J <- ggplot(Jejunum, aes(x = Col, y = Row, color = Cor)) +
  geom_point(aes(size = abs(Cor))) +
  theme_bw() +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
        panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_color_gradient2(low="#3C5488", mid="white", high="#E64B35", midpoint=0)  +
  scale_size(range=c(0,2),
             breaks=c(0,0.25,0.5,0.75,1),
             labels=c("0","0.25","0.50","0.75","1.00"),
             name = "",
             guide="legend") +
  guides(size = 'none') +
  xlab('Absolute abundance') +
  ylab('Relative abundance') +
  ggtitle('Jejunum')
f14_J
f14_D


cor2confusion <- function(correlation){
  CR1 <- correlation[correlation$Cat=='RA',]
  CR2 <- correlation[correlation$Cat=='AA',]
  CR1$Cat2 <- paste(CR1$Col,CR1$Row,sep = '-')
  CR1_1 <- CR1[,c(3,5)]
  CR2$Cat2 <- paste(CR2$Row,CR2$Col,sep = '-')
  CR2_1 <- CR2[,c(3,5)]
  colnames(CR1_1) <- c('CRA','Cat_R')
  colnames(CR2_1) <- c('CAA','Cat_A')
  CR_all <- merge(CR1_1,CR2_1,by.x = 'Cat_R',by.y = 'Cat_A')
  CR_all$CRA[CR_all$CRA>0] <- 'Pos'
  CR_all$CAA[CR_all$CAA>0] <- 'Pos'
  CR_all$CRA[CR_all$CRA<0] <- 'Neg'
  CR_all$CAA[CR_all$CAA<0] <- 'Neg'
  CR_all$CRA[CR_all$CRA==0] <- 'NS'
  CR_all$CAA[CR_all$CAA==0] <- 'NS'
  print(table(CR_all$CRA))
  print(table(CR_all$CAA))
  CR_final <- melt(table(CR_all$CRA,CR_all$CAA))
  colnames(CR_final) <- c("RA","AA",'Number')
  CR_final$Ratio <- round(CR_final$Number/sum(CR_final$Number) * 100, 1)
  CR_final$Cons[CR_final$RA=='Pos'&CR_final$AA=='Pos'] <- 'Consist'
  CR_final$Cons[CR_final$RA=='Neg'&CR_final$AA=='Neg'] <- 'Consist'
  CR_final$Cons[CR_final$RA=='NS'&CR_final$AA=='NS'] <- 'Consist'
  CR_final$Cons[is.na(CR_final$Cons)] <- 'Diff'
  return(CR_final[,c(1,2,4,5)])
}

Genus_ITS_RA_ID <- read.table('L6_RA_ITS.csv',header = T,sep = ',',row.names = 1)%>%
  colSums()
Genus_ITS_RA0 <- read.table('L6_RA_ITS.csv',header = T,sep = ',',row.names = 1) 
head(Genus_ITS_RA0)
Genus_ITS_RA <- Genus_ITS_RA0[,Genus_ITS_RA_ID>0.001*280]
Genus_ITS_AA0 <- read.table('L6_AA_ITS.csv',header = T,sep = ',',row.names = 1)
Genus_ITS_AA <- Genus_ITS_AA0[,Genus_ITS_RA_ID>0.001*280]
rownames(Genus_ITS_AA) == rownames(Genus_ITS_RA)
colnames(Genus_ITS_AA) == colnames(Genus_ITS_RA)
rownames(Genus_ITS_AA)
Genus_ITS_AA$Seg <- c(rep(c('Cecum','Duodenum','Ileum','Jejunum'),80))
Genus_ITS_RA$Seg <- c(rep(c('Cecum','Duodenum','Ileum','Jejunum'),80))
rm(Genus_ITS_RA_ID)
rm(Genus_ITS_AA0)
rm(Genus_ITS_RA0)
dim(Genus_ITS_AA)
taxa_ITS <- colnames(Genus_ITS_RA)[1:82]
colnames(Genus_ITS_AA)[1:82] <- colnames(Genus_ITS_RA)[1:82] <- paste0('Genus',1:82)

Genus_ITS_RA_C <- Genus_ITS_RA[Genus_ITS_RA$Seg=='Cecum',] %>%
  select(-c('Seg'))
Genus_ITS_RA_D <- Genus_ITS_RA[Genus_ITS_RA$Seg=='Duodenum',] %>%
  select(-c('Seg'))
Genus_ITS_RA_I <- Genus_ITS_RA[Genus_ITS_RA$Seg=='Ileum',] %>%
  select(-c('Seg'))
Genus_ITS_RA_J <- Genus_ITS_RA[Genus_ITS_RA$Seg=='Jejunum',] %>%
  select(-c('Seg'))
Genus_ITS_AA_C <- Genus_ITS_AA[Genus_ITS_AA$Seg=='Cecum',] %>%
  select(-c('Seg'))
Genus_ITS_AA_D <- Genus_ITS_AA[Genus_ITS_AA$Seg=='Duodenum',] %>%
  select(-c('Seg'))
Genus_ITS_AA_I <- Genus_ITS_AA[Genus_ITS_AA$Seg=='Ileum',] %>%
  select(-c('Seg'))
Genus_ITS_AA_J <- Genus_ITS_AA[Genus_ITS_AA$Seg=='Jejunum',] %>%
  select(-c('Seg'))


Cecum_ITS <- ab_2_cor(Genus_ITS_RA_C, Genus_ITS_AA_C)
Cecum_ITS[is.na(Cecum_ITS)] <- 0
Cecum_ITS$Row <- factor(Cecum_ITS$Row, levels = paste0('G',1:82))
Cecum_ITS$Col <- factor(Cecum_ITS$Col, levels = paste0('G',1:82))
p15_C <- ggplot(Cecum_ITS, aes(x = Col, y = Row, color = Cor)) +
  geom_point(aes(size = abs(Cor))) +
  theme_bw() +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
        panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_color_gradient2(low="#3C5488", mid="white", high="#E64B35", midpoint=0)  +
  scale_size(range=c(0,2),
             breaks=c(0,0.25,0.5,0.75,1),
             labels=c("0","0.25","0.50","0.75","1.00"),
             name = "",
             guide="legend") +
  guides(size = 'none') +
  xlab('Absolute abundance') +
  ylab('Relative abundance') +
  ggtitle('Cecum')

Ileum_ITS <- ab_2_cor(Genus_ITS_RA_I, Genus_ITS_AA_I)
Ileum_ITS[is.na(Ileum_ITS)] <- 0

Ileum_ITS$Row <- factor(Ileum_ITS$Row, levels = paste0('G',1:82))
Ileum_ITS$Col <- factor(Ileum_ITS$Col, levels = paste0('G',1:82))
p15_I <- ggplot(Ileum_ITS, aes(x = Col, y = Row, color = Cor)) +
  geom_point(aes(size = abs(Cor))) +
  theme_bw() +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
        panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_color_gradient2(low="#3C5488", mid="white", high="#E64B35", midpoint=0)  +
  scale_size(range=c(0,2),
             breaks=c(0,0.25,0.5,0.75,1),
             labels=c("0","0.25","0.50","0.75","1.00"),
             name = "",
             guide="legend") +
  guides(size = 'none') +
  xlab('Absolute abundance') +
  ylab('Relative abundance') +
  ggtitle('Ileum')

p15_I
Duodenum_ITS <- ab_2_cor(Genus_ITS_RA_D, Genus_ITS_AA_D)
Duodenum_ITS[is.na(Duodenum_ITS)] <- 0

Duodenum_ITS$Row <- factor(Duodenum_ITS$Row, levels = paste0('G',1:82))
Duodenum_ITS$Col <- factor(Duodenum_ITS$Col, levels = paste0('G',1:82))
p15_D <- ggplot(Duodenum_ITS, aes(x = Col, y = Row, color = Cor)) +
  geom_point(aes(size = abs(Cor))) +
  theme_bw() +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
        panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_color_gradient2(low="#3C5488", mid="white", high="#E64B35", midpoint=0)  +
  scale_size(range=c(0,2),
             breaks=c(0,0.25,0.5,0.75,1),
             labels=c("0","0.25","0.50","0.75","1.00"),
             name = "",
             guide="legend") +
  guides(size = 'none') +
  xlab('Absolute abundance') +
  ylab('Relative abundance') +
  ggtitle('Duodenum')
p15_D
Jejunum_ITS <- ab_2_cor(Genus_ITS_RA_J, Genus_ITS_AA_J)
Jejunum_ITS[is.na(Jejunum_ITS)] <- 0
Jejunum_ITS$Row <- factor(Jejunum_ITS$Row, levels = paste0('G',1:82))
Jejunum_ITS$Col <- factor(Jejunum_ITS$Col, levels = paste0('G',1:82))
p15_J <- ggplot(Jejunum_ITS, aes(x = Col, y = Row, color = Cor)) +
  geom_point(aes(size = abs(Cor))) +
  theme_bw() +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
        panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_color_gradient2(low="#3C5488", mid="white", high="#E64B35", midpoint=0)  +
  scale_size(range=c(0,2),
             breaks=c(0,0.25,0.5,0.75,1),
             labels=c("0","0.25","0.50","0.75","1.00"),
             name = "",
             guide="legend") +
  guides(size = 'none') +
  xlab('Absolute abundance') +
  ylab('Relative abundance') +
  ggtitle('Jejunum')

table(Duodenum$Cor > 0 & Duodenum$Cat == 'AA')[2]
table(Duodenum$Cor > 0 & Duodenum$Cat == 'RA')[2]
table(Jejunum$Cor > 0 & Jejunum$Cat == 'AA')[2]
table(Jejunum$Cor > 0 & Jejunum$Cat == 'RA')[2]
table(Ileum$Cor > 0 & Ileum$Cat == 'AA')[2]
table(Ileum$Cor > 0 & Ileum$Cat == 'RA')[2]
table(Cecum$Cor > 0 & Cecum$Cat == 'AA')[2]
table(Cecum$Cor > 0 & Cecum$Cat == 'RA')[2]

table(Duodenum_ITS$Cor > 0 & Duodenum_ITS$Cat == 'AA')[2]
table(Duodenum_ITS$Cor > 0 & Duodenum_ITS$Cat == 'RA')[2]
table(Jejunum_ITS$Cor > 0 & Jejunum_ITS$Cat == 'AA')[2]
table(Jejunum_ITS$Cor > 0 & Jejunum_ITS$Cat == 'RA')[2]
table(Ileum_ITS$Cor > 0 & Ileum_ITS$Cat == 'AA')[2]
table(Ileum_ITS$Cor > 0 & Ileum_ITS$Cat == 'RA')[2]
table(Cecum_ITS$Cor > 0 & Cecum_ITS$Cat == 'AA')[2]
table(Cecum_ITS$Cor > 0 & Cecum_ITS$Cat == 'RA')[2]

confu_C <- cor2confusion(Cecum)
confu_D <- cor2confusion(Duodenum) # AA Neg
confu_I <- cor2confusion(Ileum) # AA Neg
confu_J <- cor2confusion(Jejunum) 
confu_C_ITS <- cor2confusion(Cecum_ITS)
confu_D_ITS <- cor2confusion(Duodenum_ITS)
confu_I_ITS <- cor2confusion(Ileum_ITS) # AA Neg
confu_J_ITS <- cor2confusion(Jejunum_ITS) # AA Neg
add_confu <- data.frame(RA = c('Pos','Neg','NS'),
                        AA = c('Neg','Neg','Neg'),
                        Ratio = c(0.0,0.0,0.0),
                        Cons = c('Diff','Consist','Diff'))
confu_D <- rbind(confu_D, add_confu)
confu_I <- rbind(confu_I, add_confu)
confu_I_ITS <- rbind(confu_I_ITS, add_confu)
confu_J_ITS <- rbind(confu_J_ITS, add_confu)

relevel_matrix <- function(confusion_matrix){
  confusion_matrix$RA <- factor(confusion_matrix$RA, levels = c('Pos','Neg','NS'))
  confusion_matrix$AA <- factor(confusion_matrix$AA, levels = c('NS','Neg','Pos'))
  return(confusion_matrix)
}
confu_D <- relevel_matrix(confu_D)
confu_J <- relevel_matrix(confu_J)
confu_I <- relevel_matrix(confu_I)
confu_C <- relevel_matrix(confu_C)
confu_D_ITS <- relevel_matrix(confu_D_ITS)
confu_J_ITS <- relevel_matrix(confu_J_ITS)
confu_I_ITS <- relevel_matrix(confu_I_ITS)
confu_C_ITS <- relevel_matrix(confu_C_ITS)

p16_1 <- ggplot(confu_D,aes(x = RA, y = AA)) +
  geom_tile(aes(fill = Cons)) +
  geom_text(aes(label = Ratio,size = 4)) +
  theme_bw() + 
  theme(legend.position = 'none',panel.grid.major = element_line(color = NA), 
        panel.grid.minor = element_line(color = NA),panel.border = element_blank(),
        axis.ticks = element_blank()) +
  scale_fill_manual(values = c('#E64B35','#4DBBD5')) +
  xlab('Relative abundance') + ylab('Absolute abundance') +
  ggtitle('Duodenum') + scale_x_discrete(position = "top")
p16_1
p16_2 <- ggplot(confu_J,aes(x = RA, y = AA)) +
  geom_tile(aes(fill = Cons)) +
  geom_text(aes(label = Ratio,size = 4)) +
  theme_bw() + 
  theme(legend.position = 'none',panel.grid.major = element_line(color = NA), 
        panel.grid.minor = element_line(color = NA),panel.border = element_blank(),
        axis.ticks = element_blank()) +
  scale_fill_manual(values = c('#E64B35','#4DBBD5')) +
  xlab('Relative abundance') + ylab('Absolute abundance') +
  ggtitle('Jejunum') + scale_x_discrete(position = "top")
p16_2
p16_3 <- ggplot(confu_I,aes(x = RA, y = AA)) +
  geom_tile(aes(fill = Cons)) +
  geom_text(aes(label = Ratio,size = 4)) +
  theme_bw() + 
  theme(legend.position = 'none',panel.grid.major = element_line(color = NA), 
        panel.grid.minor = element_line(color = NA),panel.border = element_blank(),
        axis.ticks = element_blank()) +
  scale_fill_manual(values = c('#E64B35','#4DBBD5')) +
  xlab('Relative abundance') + ylab('Absolute abundance') +
  ggtitle('Ileum') + scale_x_discrete(position = "top")
p16_3
p16_4 <- ggplot(confu_C,aes(x = RA, y = AA)) +
  geom_tile(aes(fill = Cons)) +
  geom_text(aes(label = Ratio,size = 4)) +
  theme_bw() + 
  theme(legend.position = 'none',panel.grid.major = element_line(color = NA), 
        panel.grid.minor = element_line(color = NA),panel.border = element_blank(),
        axis.ticks = element_blank()) +
  scale_fill_manual(values = c('#E64B35','#4DBBD5')) +
  xlab('Relative abundance') + ylab('Absolute abundance') +
  ggtitle('Cecum') + scale_x_discrete(position = "top")
p16_4
p17_1 <- ggplot(confu_D_ITS,aes(x = RA, y = AA)) +
  geom_tile(aes(fill = Cons)) +
  geom_text(aes(label = Ratio,size = 4)) +
  theme_bw() + 
  theme(legend.position = 'none',panel.grid.major = element_line(color = NA), 
        panel.grid.minor = element_line(color = NA),panel.border = element_blank(),
        axis.ticks = element_blank()) +
  scale_fill_manual(values = c('#E64B35','#4DBBD5')) +
  xlab('Relative abundance') + ylab('Absolute abundance') +
  ggtitle('Duodenum') + scale_x_discrete(position = "top")
p17_1
p17_2 <- ggplot(confu_J_ITS,aes(x = RA, y = AA)) +
  geom_tile(aes(fill = Cons)) +
  geom_text(aes(label = Ratio,size = 4)) +
  theme_bw() + 
  theme(legend.position = 'none',panel.grid.major = element_line(color = NA), 
        panel.grid.minor = element_line(color = NA),panel.border = element_blank(),
        axis.ticks = element_blank()) +
  scale_fill_manual(values = c('#E64B35','#4DBBD5')) +
  xlab('Relative abundance') + ylab('Absolute abundance') +
  ggtitle('Jejunum') + scale_x_discrete(position = "top")
p17_2
p17_3 <- ggplot(confu_I_ITS,aes(x = RA, y = AA)) +
  geom_tile(aes(fill = Cons)) +
  geom_text(aes(label = Ratio,size = 4)) +
  theme_bw() + 
  theme(legend.position = 'none',panel.grid.major = element_line(color = NA), 
        panel.grid.minor = element_line(color = NA),panel.border = element_blank(),
        axis.ticks = element_blank()) +
  scale_fill_manual(values = c('#E64B35','#4DBBD5')) +
  xlab('Relative abundance') + ylab('Absolute abundance') +
  ggtitle('Ileum') + scale_x_discrete(position = "top")
p17_3
p17_4 <- ggplot(confu_C_ITS,aes(x = RA, y = AA)) +
  geom_tile(aes(fill = Cons)) +
  geom_text(aes(label = Ratio,size = 4)) +
  theme_bw() + 
  theme(legend.position = 'none',panel.grid.major = element_line(color = NA), 
        panel.grid.minor = element_line(color = NA),panel.border = element_blank(),
        axis.ticks = element_blank()) +
  scale_fill_manual(values = c('#E64B35','#4DBBD5')) +
  xlab('Relative abundance') + ylab('Absolute abundance') +
  ggtitle('Cecum') + 
  scale_x_discrete(position = "top") 
p17_4
svg('~/Desktop/figure materials/2nd/figure 2C.svg',width = 10.5,height = 5)
p16_1 + p16_2 + p16_3 + p16_4 + p17_1 + p17_2 + p17_3 + p17_4 + plot_layout(nrow = 2)
dev.off()
