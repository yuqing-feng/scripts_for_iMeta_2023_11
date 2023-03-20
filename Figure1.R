# Maturaty of the gut microbiota 

# Step1: Calculate the standard curve of qPCR
# part one Figure 001
# Standard curve
x1 <- c(9.942367777,13.16022936,16.5225298,19.86947911,22.81842211)
y1 <- c(61000000,6100000,610000,61000,6100)
#ITS Pichia pastoris
x2 <- c(14.0476737,17.35078621,21.22638702,24.95912552,28.33756542)
y2 <- c(270000000,27000000,2700000,270000,27000)
ec <- lm(log10(y1)~x1)
ec
pp <- lm(log10(y2)~x2)
pp

# step 2 calculate the relative abundance of Archeae, Bacteria and Fungi
rm(list = ls())
library(openxlsx)
library(dplyr)
library(reshape2)
# laod the absolute abundance of the three kingdomss
ab0 <- read.xlsx('~/01qpcr/01data_for_quantitive.xlsx',sheet = 10)
ab0
colnames(ab0)[4:6] <- c('Bacteria','Fungi','Archaea')
ab1 <- ab0 %>% select(-c('Date','Number')) %>%
  melt() 
colnames(ab1)[2:3] <- c('Kingdom','Value')

library(ggplot2)
ab1_bf <- ab1[ab1$Kingdom!='Archaea',]
ab1_ar <- ab1[ab1$Kingdom=='Archaea',]
ab1_bf$Segment <- factor(ab1_bf$Segment,levels = c('Duodenum','Jejunum','Ileum','Cecum'))
ab1_bf$Kingdom <- factor(ab1_bf$Kingdom,levels = c('Bacteria','Fungi'))
p1_1 <- ggplot(ab1_bf, aes(x = Kingdom, y = Value, fill = Segment)) +
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) +
  theme_bw() + 
  ylab('Log10(copies per gram)')
p1_1
library(ggpubr)
ab1_ar$Archaea[ab1_ar$Value==0] <- 'Absent'
ab1_ar$Archaea[ab1_ar$Value!=0] <- 'Present'
head(ab1_ar)

p1_2 <- table(ab1_ar[,c(1,4)]) %>% 
  melt() %>%
  mutate(Segment = factor(Segment, levels = c('Duodenum','Jejunum','Ileum','Cecum')))%>%
  ggplot(aes(x = "", y = value, fill = Archaea)) + 
  geom_bar(stat = 'identity',width = 1, position = position_fill()) +
  coord_polar(theta = "y") +
  annotate('text', x = 0.01, y = 0.01, label = c('00.0%'))+
  facet_wrap( ~ Segment)+
  xlab('') + ylab('') +
  theme_bw() +
  theme(legend.position = "bottom",panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        axis.ticks = element_blank(),axis.text = element_blank()) +
  scale_fill_manual(values = c('#3C5488','#00A087')) 
p1_2

p1_3 <- ab0[ab0$Archaea!=0,] %>%
  select(c('Segment','Bacteria','Fungi','Archaea')) %>%
  melt() %>%
  mutate(Segment = factor(Segment, levels = c('Duodenum','Jejunum','Cecum')))%>%
  ggplot(aes(x = variable,y = value,fill = variable)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(Segment~.) +
  theme_bw() +
  theme(legend.position = 'none',panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)) +
  scale_y_continuous(breaks = c(2,3,4,5,6,7,8,9,10,11,12)) +
  xlab('Kingdom') + ylab('Log10(copies per gram)') +
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087'))
p1_3

ab_bf <- ab0[,c(1,2,4,5)] %>%
  melt()
head(ab_bf)
svg('')

ab_bf$Segment <- factor(ab_bf$Segment, levels = c('Duodenum','Jejunum','Ileum','Cecum'))
ab_bf$DPH <- rep(c(rep(c(1,4,7,14,21,28,35,42),each = 10)),8)

f1_1 <- ggplot(ab_bf, aes(x = DPH, y = value,  color = Segment, group = Segment)) +
  geom_smooth(size = 1) +
  geom_jitter(width = 0.1,alpha = .5,size = 0.5) +
  scale_color_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) +
  theme_bw() +
  ylab('Log10(copies per gram)') +
  xlab("DPH") +
  scale_x_continuous(breaks = c(1,7,14,21,28,35,42)) +
  theme(legend.position = "right",panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  facet_wrap(.~variable)
f1_1

head(ab0)
ab0$Ratio <- ab0$Bacteria-ab0$Fungi
range(ab0$Ratio)
ab_r <- aggregate(Ratio~Segment + Date,ab0,median)
ab_r
ab_r$Segment <- factor(ab_r$Segment, levels = c('Duodenum','Jejunum','Ileum','Cecum'))
ab_r$Date <- factor(ab_r$Date, levels = c('D42','D35','D28','D21','D14','D07','D04','D01'))
f1_2 <- ggplot(ab_r,aes(x = Segment,y = Date,color = Ratio)) +
  geom_point(size = 4) +
  scale_color_gradient(low = 'white',high = '#3C5488') +
  theme_bw() +
  scale_y_discrete(labels= rev(c(1,4,7,14,21,28,35,42))) +
  xlab('Segment') +
  theme(legend.position = "right",panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_color_continuous(low = '#3C54881f', high = '#3C5488',
                         breaks = c(1, 2, 3, 4), labels = c(expression(10^1), expression(10^2), expression(10^3),expression(10^4)))
f1_2

library(patchwork)
layout <- "AAAB"
f1_1 + f1_2 + plot_layout(design = layout) 
svg(f1_1,'~/Desktop/01maturity/figure materials/3rd_svg/Figure1a.svg',width = 13,height = 5)
f1_1
dev.off()

library(tibble)
raw_16S <- read.csv('~/Desktop/01maturity/material/02micro/16S/level-2.csv')
dim(raw_16S)
format_16S <- raw_16S %>%
  remove_rownames %>% 
  column_to_rownames(var="index") %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  select(-c('D_0__Archaea.D_1__Euryarchaeota')) %>%
  mutate(across(D_0__Bacteria.D_1__Acidobacteria:D_0__Bacteria.D_1__Zixibacteria, ~ . / sum)) %>%
  select(-c('sum','Group2')) %>%
  setNames(gsub("D_0__Bacteria.D_1__", "", names(.))) %>%
  select(c('Firmicutes','Proteobacteria','Epsilonbacteraeota','Bacteroidetes','Actinobacteria',
           'Verrucomicrobia','Group1')) %>%
  mutate(Others = 1- rowSums(across(where(is.numeric)))) %>%
  relocate(Group1, .after = Others)
format_16S

raw_ITS <- read.csv('~/Desktop/01maturity/material/02micro/ITS/level-2.csv')
dim(raw_ITS)
format_ITS <- raw_ITS %>%
  remove_rownames %>% 
  column_to_rownames(var="index") %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  select(-c('k__Fungi.p__unidentified','Group2')) %>%
  mutate(across(k__Fungi.p__Ascomycota:k__Fungi.p__Rozellomycota, ~ . / sum)) %>%
  select(-c('sum')) %>%
  setNames(gsub("k__Fungi.p__", "", names(.))) %>%
  select(c('Ascomycota','Basidiomycota','Mortierellomycota','Mucoromycota','Glomeromycota','Chytridiomycota','Group1')) %>%
  mutate(Others = 1- rowSums(across(where(is.numeric)))) %>%
  relocate(Group1, .after = Others)
format_ITS
ab_temp <- read.xlsx('~/Desktop/01maturity/material/01qpcr/01data_for_quantitive.xlsx',sheet = 11)
head(ab_temp)
rownames(format_ITS) == ab_temp$ID
rownames(format_16S) == ab_temp$ID
format_ab_ITS <- format_ITS
format_ab_16S <- format_16S
for(i in 1:nrow(ab_temp)){
  format_ab_16S[i,1:7] <- format_16S[i,1:7] * 10^(ab_temp$`log10(16S)`[i])
  format_ab_ITS[i,1:7] <- format_ITS[i,1:7] * 10^(ab_temp$`log10(ITS)`[i])
}
format_ab_16S$DPH <- substring(rownames(format_ab_16S),1,3)
format_ab_ITS$DPH <- substring(rownames(format_ab_ITS),1,3)

raw_16S_L6 <- read.csv('~/Desktop/01maturity/material/02micro/16S/level-6_no_archaea.csv')
dim(raw_16S_L6)
format_16S_L6 <- raw_16S_L6 %>%
  remove_rownames %>% 
  column_to_rownames(var="index") %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  mutate(across(Granulicella:uncultured.bacterium.40, ~ . / sum)) %>%
  select(-c('sum','Group2')) %>%
  select(c('Lactobacillus','Streptococcus','Escherichia.Shigella','Helicobacter','Enterococcus',
           'Faecalibacterium','Candidatus.Arthromitus','X.Ruminococcus..torques.group','Group1')) %>%
  mutate(Others = 1- rowSums(across(where(is.numeric)))) %>%
  relocate(Group1, .after = Others)
format_16S_L6


raw_ITS_L6 <- read.csv('~/Desktop/01maturity/material/02micro/ITS/level-6_for_bar.csv')
dim(raw_ITS_L6)
format_ITS_L6 <- raw_ITS_L6 %>%
  remove_rownames %>% 
  column_to_rownames(var="X") %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  select(-c('Group2')) %>%
  mutate(across(X__:unidentified.53, ~ . / sum)) %>%
  select(-c('sum')) %>%
  select(c('Candida','Fusarium','Trichosporon','Plectosphaerella','Aspergillus','Cladosporium','Tetracladium','Talaromyces','Group1')) %>%
  mutate(Others = 1- rowSums(across(where(is.numeric)))) %>%
  relocate(Group1, .after = Others)
format_ITS
ab_temp <- read.xlsx('~/Desktop/01maturity/material/01qpcr/01data_for_quantitive.xlsx',sheet = 11)
head(ab_temp)
rownames(format_ITS) == ab_temp$ID
rownames(format_16S) == ab_temp$ID
format_ab_ITS_L6 <- format_ITS_L6
format_ab_16S_L6 <- format_16S_L6
for(i in 1:nrow(ab_temp)){
  format_ab_16S_L6[i,1:7] <- format_16S_L6[i,1:7] * 10^(ab_temp$`log10(16S)`[i])
  format_ab_ITS_L6[i,1:7] <- format_ITS_L6[i,1:7] * 10^(ab_temp$`log10(ITS)`[i])
}
format_ab_16S_L6$DPH <- substring(rownames(format_ab_16S_L6),1,3)
format_ab_ITS_L6$DPH <- substring(rownames(format_ab_ITS_L6),1,3)


# phylum level relative abundance
format_16S_RA2 <- data.frame(format_16S,DPH = rep(c('D01','D04','D07','D14','D21','D28','D35','D42'), each = 40))
format_ITS_RA2 <- data.frame(format_ITS,DPH = rep(c('D01','D04','D07','D14','D21','D28','D35','D42'), each = 40))
f6_1 <- format_16S_RA2 %>%
  melt() %>%
  group_by(DPH, Group1, variable)%>% 
  summarise(mean = mean(value)) %>%
  filter(Group1=='Duodenum') %>%
  ggplot(aes(x = DPH, y = mean,fill = variable)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_x_discrete(labels = c(1,4,7,14,21,28,35,42)) +
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#DC0000','#7E6148')) +
  ylab('Relative abundnce') +
  guides(fill=guide_legend(title="Phylum")) +
  ggtitle('Duodenum')
f6_1
f6_2 <- format_16S_RA2 %>%
  melt() %>%
  group_by(DPH, Group1, variable)%>% 
  summarise(mean = mean(value)) %>%
  filter(Group1=='Jejunum') %>%
  ggplot(aes(x = DPH, y = mean,fill = variable)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_x_discrete(labels = c(1,4,7,14,21,28,35,42)) +
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#DC0000','#7E6148')) +
  ylab('Relative abundnce') +
  guides(fill=guide_legend(title="Phylum")) +
  ggtitle('Jejunum')
f6_2
f6_3 <- format_16S_RA2 %>%
  melt() %>%
  group_by(DPH, Group1, variable)%>% 
  summarise(mean = mean(value)) %>%
  filter(Group1=='Ileum') %>%
  ggplot(aes(x = DPH, y = mean,fill = variable)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_x_discrete(labels = c(1,4,7,14,21,28,35,42)) +
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#DC0000','#7E6148')) +
  ylab('Relative abundnce') +
  guides(fill=guide_legend(title="Phylum")) +
  ggtitle('Ileum')
f6_3
f6_4 <- format_16S_RA2 %>%
  melt() %>%
  group_by(DPH, Group1, variable)%>% 
  summarise(mean = mean(value)) %>%
  filter(Group1=='Cecum') %>%
  ggplot(aes(x = DPH, y = mean,fill = variable)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_x_discrete(labels = c(1,4,7,14,21,28,35,42)) +
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#DC0000','#7E6148')) +
  ylab('Relative abundnce') +
  guides(fill=guide_legend(title="Phylum")) +
  ggtitle('Cecum')
f6_4

f7_1 <- format_ITS_RA2 %>%
  melt() %>%
  group_by(DPH, Group1, variable)%>% 
  summarise(mean = mean(value)) %>%
  filter(Group1=='Duodenum') %>%
  ggplot(aes(x = DPH, y = mean,fill = variable)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_x_discrete(labels = c(1,4,7,14,21,28,35,42)) +
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#DC0000','#7E6148')) +
  ylab('Relative abundnce') +
  guides(fill=guide_legend(title="Phylum")) +
  ggtitle('Duodenum')
f7_1
f7_2 <- format_ITS_RA2 %>%
  melt() %>%
  group_by(DPH, Group1, variable)%>% 
  summarise(mean = mean(value)) %>%
  filter(Group1=='Jejunum') %>%
  ggplot(aes(x = DPH, y = mean,fill = variable)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_x_discrete(labels = c(1,4,7,14,21,28,35,42)) +
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#DC0000','#7E6148')) +
  ylab('Relative abundnce') +
  guides(fill=guide_legend(title="Phylum")) +
  ggtitle('Jejunum')
f7_2
f7_3 <- format_ITS_RA2 %>%
  melt() %>%
  group_by(DPH, Group1, variable)%>% 
  summarise(mean = mean(value)) %>%
  filter(Group1=='Ileum') %>%
  ggplot(aes(x = DPH, y = mean,fill = variable)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_x_discrete(labels = c(1,4,7,14,21,28,35,42)) +
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#DC0000','#7E6148')) +
  ylab('Relative abundnce') +
  guides(fill=guide_legend(title="Phylum")) +
  ggtitle('Ileum')
f7_3
f7_4 <- format_ITS_RA2 %>%
  melt() %>%
  group_by(DPH, Group1, variable)%>% 
  summarise(mean = mean(value)) %>%
  filter(Group1=='Cecum') %>%
  ggplot(aes(x = DPH, y = mean,fill = variable)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_x_discrete(labels = c(1,4,7,14,21,28,35,42)) +
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#DC0000','#7E6148')) +
  ylab('Relative abundnce') +
  guides(fill=guide_legend(title="Phylum")) +
  ggtitle('Cecum')
f7_4
svg('~/Desktop/figure materials/2nd/figure S3.svg',width = 10.5,height = 7)
{f6_1 + f6_2 + f6_3 + f6_4 + plot_layout(guides = 'collect',nrow = 1)}/{f7_1 + f7_2 + f7_3 + f7_4 + plot_layout(guides = 'collect',nrow = 1)}
dev.off()
# genus level relative abundance
format_16S_L6_RA2 <- data.frame(format_16S_L6,DPH = rep(c('D01','D04','D07','D14','D21','D28','D35','D42'), each = 40))
format_ITS_L6_RA2 <- data.frame(format_ITS_L6,DPH = rep(c('D01','D04','D07','D14','D21','D28','D35','D42'), each = 40))
f8_1 <- format_16S_L6_RA2 %>%
  melt() %>%
  group_by(DPH, Group1, variable)%>% 
  summarise(mean = mean(value)) %>%
  filter(Group1=='Duodenum') %>%
  ggplot(aes(x = DPH, y = mean,fill = variable)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_x_discrete(labels = c(1,4,7,14,21,28,35,42)) +
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#DC0000','#7E6148')) +
  ylab('Relative abundnce') +
  guides(fill=guide_legend(title="Genus")) +
  ggtitle('Duodenum')
f8_1
f8_2 <- format_16S_L6_RA2 %>%
  melt() %>%
  group_by(DPH, Group1, variable)%>% 
  summarise(mean = mean(value)) %>%
  filter(Group1=='Jejunum') %>%
  ggplot(aes(x = DPH, y = mean,fill = variable)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_x_discrete(labels = c(1,4,7,14,21,28,35,42)) +
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#DC0000','#7E6148')) +
  ylab('Relative abundnce') +
  guides(fill=guide_legend(title="Genus")) +
  ggtitle('Jejunum')
f8_2
f8_3 <- format_16S_L6_RA2 %>%
  melt() %>%
  group_by(DPH, Group1, variable)%>% 
  summarise(mean = mean(value)) %>%
  filter(Group1=='Ileum') %>%
  ggplot(aes(x = DPH, y = mean,fill = variable)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_x_discrete(labels = c(1,4,7,14,21,28,35,42)) +
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#DC0000','#7E6148')) +
  ylab('Relative abundnce') +
  guides(fill=guide_legend(title="Genus")) +
  ggtitle('Ileum')
f8_3
f8_4 <- format_16S_L6_RA2 %>%
  melt() %>%
  group_by(DPH, Group1, variable)%>% 
  summarise(mean = mean(value)) %>%
  filter(Group1=='Cecum') %>%
  ggplot(aes(x = DPH, y = mean,fill = variable)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_x_discrete(labels = c(1,4,7,14,21,28,35,42)) +
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#DC0000','#7E6148')) +
  ylab('Relative abundnce') +
  guides(fill=guide_legend(title="Genus")) +
  ggtitle('Cecum')
f8_4

f9_1 <- format_ITS_L6_RA2 %>%
  melt() %>%
  group_by(DPH, Group1, variable)%>% 
  summarise(mean = mean(value)) %>%
  filter(Group1=='Duodenum') %>%
  ggplot(aes(x = DPH, y = mean,fill = variable)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_x_discrete(labels = c(1,4,7,14,21,28,35,42)) +
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#DC0000','#7E6148')) +
  ylab('Relative abundnce') +
  guides(fill=guide_legend(title="Genus")) +
  ggtitle('Duodenum')
f9_1
f9_2 <- format_ITS_L6_RA2 %>%
  melt() %>%
  group_by(DPH, Group1, variable)%>% 
  summarise(mean = mean(value)) %>%
  filter(Group1=='Jejunum') %>%
  ggplot(aes(x = DPH, y = mean,fill = variable)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_x_discrete(labels = c(1,4,7,14,21,28,35,42)) +
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#DC0000','#7E6148')) +
  ylab('Relative abundnce') +
  guides(fill=guide_legend(title="Genus")) +
  ggtitle('Jejunum')
f9_2
f9_3 <- format_ITS_L6_RA2 %>%
  melt() %>%
  group_by(DPH, Group1, variable)%>% 
  summarise(mean = mean(value)) %>%
  filter(Group1=='Ileum') %>%
  ggplot(aes(x = DPH, y = mean,fill = variable)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_x_discrete(labels = c(1,4,7,14,21,28,35,42)) +
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#DC0000','#7E6148')) +
  ylab('Relative abundnce') +
  guides(fill=guide_legend(title="Genus")) +
  ggtitle('Ileum')
f9_3
f9_4 <- format_ITS_L6_RA2 %>%
  melt() %>%
  group_by(DPH, Group1, variable)%>% 
  summarise(mean = mean(value)) %>%
  filter(Group1=='Cecum') %>%
  ggplot(aes(x = DPH, y = mean,fill = variable)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_x_discrete(labels = c(1,4,7,14,21,28,35,42)) +
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#DC0000','#7E6148')) +
  ylab('Relative abundnce') +
  guides(fill=guide_legend(title="Genus")) +
  ggtitle('Cecum')
f9_4

f6_1 + f6_2 + f6_3 + f6_4 + f7_1 + f7_2 + f7_3 + f7_4 + plot_layout(design = layout3,guides = 'collect') 

f8_1 + f8_2 + f8_3 + f8_4 + f9_1 + f9_2 + f9_3 + f9_4 + plot_layout(design = layout3,guides = 'collect') 

# alpha diversity
library(vegan)
alpha_bacteria <- diversity(raw_16S_L6[,2:530])
alpha_fungi <- diversity(raw_ITS_L6[,2:604])
alpha_bacteria
alpha_all <- data.frame(alpha_bacteria,alpha_fungi,ab_temp[,c(1,2)])
head(alpha_all)
alpha_all$Segment <- factor(alpha_all$Segment, levels = c('Duodenum','Jejunum','Ileum','Cecum'))
alpha_all$DPH <- rep(c(1,4,7,14,21,28,35,42),each = 40)
f10_1 <- ggplot(alpha_all, aes(x = DPH, y = alpha_bacteria,  color = Segment, group = Segment)) +
  geom_smooth(size = 1) +
  geom_jitter(width = 0.1,alpha = .5,size = 0.5) +
  scale_color_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) +
  theme_bw() +
  ylab('Shannon index') +
  xlab("DPH") +
  scale_x_continuous(breaks = c(1,7,14,21,28,35,42)) +
  theme(legend.position = "right",panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_color_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) +
  ggtitle('Bacteria')
f10_1
head(alpha_all)
f10_2 <- ggplot(alpha_all, aes(x = DPH, y = alpha_fungi,  color = Segment, group = Segment)) +
  geom_smooth(size = 1) +
  geom_jitter(width = 0.1,alpha = .5,size = 0.5) +
  scale_color_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) +
  theme_bw() +
  ylab('Shannon index') +
  xlab("DPH") +
  scale_x_continuous(breaks = c(1,7,14,21,28,35,42)) +
  theme(legend.position = "right",panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_color_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) +
  ggtitle('Fungi') 
f10_1
f10_2
f10_1 + f10_2 + plot_layout(guides = 'collect')


# bray-curtis dissimarity
library(tidyverse)
raw_16S_ra <-  raw_16S_L6 %>%
  remove_rownames %>% 
  column_to_rownames(var="index") %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  mutate(across(Granulicella:uncultured.bacterium.40, ~ . / sum)) %>%
  select(-c('sum','Group1')) 

raw_ITS_ra <- raw_ITS_L6 %>%
  remove_rownames %>% 
  column_to_rownames(var="X") %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  mutate(across(X__:unidentified.53, ~ . / sum)) %>%
  select(-c('sum','Group1'))


raw_16S_ab <- raw_16S_ra
raw_ITS_ab <- raw_ITS_ra
for(i in 1:nrow(ab_temp)){
  raw_16S_ab[i,1:529] <- raw_16S_ab[i,1:529] * 10^(ab_temp$`log10(16S)`[i])
  raw_ITS_ab[i,1:603] <- raw_ITS_ab[i,1:603] * 10^(ab_temp$`log10(ITS)`[i])
}

library(vegan)
dist_16S <- as.matrix(vegdist(raw_16S_ab[,1:c(ncol(raw_16S_ab)-1)]))
dist_ITS <- as.matrix(vegdist(raw_ITS_ab[,1:c(ncol(raw_ITS_ab)-1)]))
# dist_16S[raw_16S_ab$ID=='Cecum_D01C',raw_16S_ab$ID=='Cecum_D04C']

dist_16S_C1 <- dist_16S[raw_16S_ab$Group2=='D01C',raw_16S_ab$Group2=='D04C'] %>% melt() %>% mutate(Group = 'Cecum_1vs4') %>% select(-c('Var1','Var2'))
dist_16S_C2 <- dist_16S[raw_16S_ab$Group2=='D04C',raw_16S_ab$Group2=='D07C'] %>% melt() %>% mutate(Group = 'Cecum_4vs7') %>% select(-c('Var1','Var2'))
dist_16S_C3 <- dist_16S[raw_16S_ab$Group2=='D07C',raw_16S_ab$Group2=='D14C'] %>% melt() %>% mutate(Group = 'Cecum_7vs14') %>% select(-c('Var1','Var2'))
dist_16S_C4 <- dist_16S[raw_16S_ab$Group2=='D14C',raw_16S_ab$Group2=='D21C'] %>% melt() %>% mutate(Group = 'Cecum_14vs21') %>% select(-c('Var1','Var2'))
dist_16S_C5 <- dist_16S[raw_16S_ab$Group2=='D21C',raw_16S_ab$Group2=='D28C'] %>% melt() %>% mutate(Group = 'Cecum_21vs28') %>% select(-c('Var1','Var2'))
dist_16S_C6 <- dist_16S[raw_16S_ab$Group2=='D28C',raw_16S_ab$Group2=='D35C'] %>% melt() %>% mutate(Group = 'Cecum_28vs35') %>% select(-c('Var1','Var2'))
dist_16S_C7 <- dist_16S[raw_16S_ab$Group2=='D35C',raw_16S_ab$Group2=='D42C'] %>% melt() %>% mutate(Group = 'Cecum_35vs42') %>% select(-c('Var1','Var2'))

dist_16S_D1 <- dist_16S[raw_16S_ab$Group2=='D01D',raw_16S_ab$Group2=='D04D'] %>% melt() %>% mutate(Group = 'Duodenum_1vs4') %>% select(-c('Var1','Var2'))
dist_16S_D2 <- dist_16S[raw_16S_ab$Group2=='D04D',raw_16S_ab$Group2=='D07D'] %>% melt() %>% mutate(Group = 'Duodenum_4vs7') %>% select(-c('Var1','Var2'))
dist_16S_D3 <- dist_16S[raw_16S_ab$Group2=='D07D',raw_16S_ab$Group2=='D14D'] %>% melt() %>% mutate(Group = 'Duodenum_7vs14') %>% select(-c('Var1','Var2'))
dist_16S_D4 <- dist_16S[raw_16S_ab$Group2=='D14D',raw_16S_ab$Group2=='D21D'] %>% melt() %>% mutate(Group = 'Duodenum_14vs21') %>% select(-c('Var1','Var2'))
dist_16S_D5 <- dist_16S[raw_16S_ab$Group2=='D21D',raw_16S_ab$Group2=='D28D'] %>% melt() %>% mutate(Group = 'Duodenum_21vs28') %>% select(-c('Var1','Var2'))
dist_16S_D6 <- dist_16S[raw_16S_ab$Group2=='D28D',raw_16S_ab$Group2=='D35D'] %>% melt() %>% mutate(Group = 'Duodenum_28vs35') %>% select(-c('Var1','Var2'))
dist_16S_D7 <- dist_16S[raw_16S_ab$Group2=='D35D',raw_16S_ab$Group2=='D42D'] %>% melt() %>% mutate(Group = 'Duodenum_35vs42') %>% select(-c('Var1','Var2'))

dist_16S_I1 <- dist_16S[raw_16S_ab$Group2=='D01I',raw_16S_ab$Group2=='D04I'] %>% melt() %>% mutate(Group = 'Ileum_1vs4') %>% select(-c('Var1','Var2'))
dist_16S_I2 <- dist_16S[raw_16S_ab$Group2=='D04I',raw_16S_ab$Group2=='D07I'] %>% melt() %>% mutate(Group = 'Ileum_4vs7') %>% select(-c('Var1','Var2'))
dist_16S_I3 <- dist_16S[raw_16S_ab$Group2=='D07I',raw_16S_ab$Group2=='D14I'] %>% melt() %>% mutate(Group = 'Ileum_7vs14') %>% select(-c('Var1','Var2'))
dist_16S_I4 <- dist_16S[raw_16S_ab$Group2=='D14I',raw_16S_ab$Group2=='D21I'] %>% melt() %>% mutate(Group = 'Ileum_14vs21') %>% select(-c('Var1','Var2'))
dist_16S_I5 <- dist_16S[raw_16S_ab$Group2=='D21I',raw_16S_ab$Group2=='D28I'] %>% melt() %>% mutate(Group = 'Ileum_21vs28') %>% select(-c('Var1','Var2'))
dist_16S_I6 <- dist_16S[raw_16S_ab$Group2=='D28I',raw_16S_ab$Group2=='D35I'] %>% melt() %>% mutate(Group = 'Ileum_28vs35') %>% select(-c('Var1','Var2'))
dist_16S_I7 <- dist_16S[raw_16S_ab$Group2=='D35I',raw_16S_ab$Group2=='D42I'] %>% melt() %>% mutate(Group = 'Ileum_35vs42') %>% select(-c('Var1','Var2'))

dist_16S_J1 <- dist_16S[raw_16S_ab$Group2=='D01J',raw_16S_ab$Group2=='D04J'] %>% melt() %>% mutate(Group = 'Jejunum_1vs4') %>% select(-c('Var1','Var2'))
dist_16S_J2 <- dist_16S[raw_16S_ab$Group2=='D04J',raw_16S_ab$Group2=='D07J'] %>% melt() %>% mutate(Group = 'Jejunum_4vs7') %>% select(-c('Var1','Var2'))
dist_16S_J3 <- dist_16S[raw_16S_ab$Group2=='D07J',raw_16S_ab$Group2=='D14J'] %>% melt() %>% mutate(Group = 'Jejunum_7vs14') %>% select(-c('Var1','Var2'))
dist_16S_J4 <- dist_16S[raw_16S_ab$Group2=='D14J',raw_16S_ab$Group2=='D21J'] %>% melt() %>% mutate(Group = 'Jejunum_14vs21') %>% select(-c('Var1','Var2'))
dist_16S_J5 <- dist_16S[raw_16S_ab$Group2=='D21J',raw_16S_ab$Group2=='D28J'] %>% melt() %>% mutate(Group = 'Jejunum_21vs28') %>% select(-c('Var1','Var2'))
dist_16S_J6 <- dist_16S[raw_16S_ab$Group2=='D28J',raw_16S_ab$Group2=='D35J'] %>% melt() %>% mutate(Group = 'Jejunum_28vs35') %>% select(-c('Var1','Var2'))
dist_16S_J7 <- dist_16S[raw_16S_ab$Group2=='D35J',raw_16S_ab$Group2=='D42J'] %>% melt() %>% mutate(Group = 'Jejunum_35vs42') %>% select(-c('Var1','Var2'))

dist_16S_combined <- rbind(dist_16S_C1,dist_16S_C2,dist_16S_C3,dist_16S_C4,dist_16S_C5,dist_16S_C6,dist_16S_C7,
                           dist_16S_D1,dist_16S_D2,dist_16S_D3,dist_16S_D4,dist_16S_D5,dist_16S_D6,dist_16S_D7,
                           dist_16S_I1,dist_16S_I2,dist_16S_I3,dist_16S_I4,dist_16S_I5,dist_16S_I6,dist_16S_I7,
                           dist_16S_J1,dist_16S_J2,dist_16S_J3,dist_16S_J4,dist_16S_J5,dist_16S_J6,dist_16S_J7) %>% 
  separate(Group, c("Segment","DPH"), sep = "([_])") %>%
  mutate(DPH = factor(DPH, levels = c('1vs4','4vs7','7vs14','14vs21','21vs28','28vs35','35vs42'))) %>%
  mutate(Segment = factor(Segment, levels = c('Duodenum','Jejunum','Ileum','Cecum')))



dist_ITS_C1 <- dist_ITS[raw_ITS_ab$Group2=='D01C',raw_ITS_ab$Group2=='D04C'] %>% melt() %>% mutate(Group = 'Cecum_1vs4') %>% select(-c('Var1','Var2'))
dist_ITS_C2 <- dist_ITS[raw_ITS_ab$Group2=='D04C',raw_ITS_ab$Group2=='D07C'] %>% melt() %>% mutate(Group = 'Cecum_4vs7') %>% select(-c('Var1','Var2'))
dist_ITS_C3 <- dist_ITS[raw_ITS_ab$Group2=='D07C',raw_ITS_ab$Group2=='D14C'] %>% melt() %>% mutate(Group = 'Cecum_7vs14') %>% select(-c('Var1','Var2'))
dist_ITS_C4 <- dist_ITS[raw_ITS_ab$Group2=='D14C',raw_ITS_ab$Group2=='D21C'] %>% melt() %>% mutate(Group = 'Cecum_14vs21') %>% select(-c('Var1','Var2'))
dist_ITS_C5 <- dist_ITS[raw_ITS_ab$Group2=='D21C',raw_ITS_ab$Group2=='D28C'] %>% melt() %>% mutate(Group = 'Cecum_21vs28') %>% select(-c('Var1','Var2'))
dist_ITS_C6 <- dist_ITS[raw_ITS_ab$Group2=='D28C',raw_ITS_ab$Group2=='D35C'] %>% melt() %>% mutate(Group = 'Cecum_28vs35') %>% select(-c('Var1','Var2'))
dist_ITS_C7 <- dist_ITS[raw_ITS_ab$Group2=='D35C',raw_ITS_ab$Group2=='D42C'] %>% melt() %>% mutate(Group = 'Cecum_35vs42') %>% select(-c('Var1','Var2'))

dist_ITS_D1 <- dist_ITS[raw_ITS_ab$Group2=='D01D',raw_ITS_ab$Group2=='D04D'] %>% melt() %>% mutate(Group = 'Duodenum_1vs4') %>% select(-c('Var1','Var2'))
dist_ITS_D2 <- dist_ITS[raw_ITS_ab$Group2=='D04D',raw_ITS_ab$Group2=='D07D'] %>% melt() %>% mutate(Group = 'Duodenum_4vs7') %>% select(-c('Var1','Var2'))
dist_ITS_D3 <- dist_ITS[raw_ITS_ab$Group2=='D07D',raw_ITS_ab$Group2=='D14D'] %>% melt() %>% mutate(Group = 'Duodenum_7vs14') %>% select(-c('Var1','Var2'))
dist_ITS_D4 <- dist_ITS[raw_ITS_ab$Group2=='D14D',raw_ITS_ab$Group2=='D21D'] %>% melt() %>% mutate(Group = 'Duodenum_14vs21') %>% select(-c('Var1','Var2'))
dist_ITS_D5 <- dist_ITS[raw_ITS_ab$Group2=='D21D',raw_ITS_ab$Group2=='D28D'] %>% melt() %>% mutate(Group = 'Duodenum_21vs28') %>% select(-c('Var1','Var2'))
dist_ITS_D6 <- dist_ITS[raw_ITS_ab$Group2=='D28D',raw_ITS_ab$Group2=='D35D'] %>% melt() %>% mutate(Group = 'Duodenum_28vs35') %>% select(-c('Var1','Var2'))
dist_ITS_D7 <- dist_ITS[raw_ITS_ab$Group2=='D35D',raw_ITS_ab$Group2=='D42D'] %>% melt() %>% mutate(Group = 'Duodenum_35vs42') %>% select(-c('Var1','Var2'))

dist_ITS_I1 <- dist_ITS[raw_ITS_ab$Group2=='D01I',raw_ITS_ab$Group2=='D04I'] %>% melt() %>% mutate(Group = 'Ileum_1vs4') %>% select(-c('Var1','Var2'))
dist_ITS_I2 <- dist_ITS[raw_ITS_ab$Group2=='D04I',raw_ITS_ab$Group2=='D07I'] %>% melt() %>% mutate(Group = 'Ileum_4vs7') %>% select(-c('Var1','Var2'))
dist_ITS_I3 <- dist_ITS[raw_ITS_ab$Group2=='D07I',raw_ITS_ab$Group2=='D14I'] %>% melt() %>% mutate(Group = 'Ileum_7vs14') %>% select(-c('Var1','Var2'))
dist_ITS_I4 <- dist_ITS[raw_ITS_ab$Group2=='D14I',raw_ITS_ab$Group2=='D21I'] %>% melt() %>% mutate(Group = 'Ileum_14vs21') %>% select(-c('Var1','Var2'))
dist_ITS_I5 <- dist_ITS[raw_ITS_ab$Group2=='D21I',raw_ITS_ab$Group2=='D28I'] %>% melt() %>% mutate(Group = 'Ileum_21vs28') %>% select(-c('Var1','Var2'))
dist_ITS_I6 <- dist_ITS[raw_ITS_ab$Group2=='D28I',raw_ITS_ab$Group2=='D35I'] %>% melt() %>% mutate(Group = 'Ileum_28vs35') %>% select(-c('Var1','Var2'))
dist_ITS_I7 <- dist_ITS[raw_ITS_ab$Group2=='D35I',raw_ITS_ab$Group2=='D42I'] %>% melt() %>% mutate(Group = 'Ileum_35vs42') %>% select(-c('Var1','Var2'))

dist_ITS_J1 <- dist_ITS[raw_ITS_ab$Group2=='D01J',raw_ITS_ab$Group2=='D04J'] %>% melt() %>% mutate(Group = 'Jejunum_1vs4') %>% select(-c('Var1','Var2'))
dist_ITS_J2 <- dist_ITS[raw_ITS_ab$Group2=='D04J',raw_ITS_ab$Group2=='D07J'] %>% melt() %>% mutate(Group = 'Jejunum_4vs7') %>% select(-c('Var1','Var2'))
dist_ITS_J3 <- dist_ITS[raw_ITS_ab$Group2=='D07J',raw_ITS_ab$Group2=='D14J'] %>% melt() %>% mutate(Group = 'Jejunum_7vs14') %>% select(-c('Var1','Var2'))
dist_ITS_J4 <- dist_ITS[raw_ITS_ab$Group2=='D14J',raw_ITS_ab$Group2=='D21J'] %>% melt() %>% mutate(Group = 'Jejunum_14vs21') %>% select(-c('Var1','Var2'))
dist_ITS_J5 <- dist_ITS[raw_ITS_ab$Group2=='D21J',raw_ITS_ab$Group2=='D28J'] %>% melt() %>% mutate(Group = 'Jejunum_21vs28') %>% select(-c('Var1','Var2'))
dist_ITS_J6 <- dist_ITS[raw_ITS_ab$Group2=='D28J',raw_ITS_ab$Group2=='D35J'] %>% melt() %>% mutate(Group = 'Jejunum_28vs35') %>% select(-c('Var1','Var2'))
dist_ITS_J7 <- dist_ITS[raw_ITS_ab$Group2=='D35J',raw_ITS_ab$Group2=='D42J'] %>% melt() %>% mutate(Group = 'Jejunum_35vs42') %>% select(-c('Var1','Var2'))

dist_ITS_combined <- rbind(dist_ITS_C1,dist_ITS_C2,dist_ITS_C3,dist_ITS_C4,dist_ITS_C5,dist_ITS_C6,dist_ITS_C7,
                           dist_ITS_D1,dist_ITS_D2,dist_ITS_D3,dist_ITS_D4,dist_ITS_D5,dist_ITS_D6,dist_ITS_D7,
                           dist_ITS_I1,dist_ITS_I2,dist_ITS_I3,dist_ITS_I4,dist_ITS_I5,dist_ITS_I6,dist_ITS_I7,
                           dist_ITS_J1,dist_ITS_J2,dist_ITS_J3,dist_ITS_J4,dist_ITS_J5,dist_ITS_J6,dist_ITS_J7) %>% 
  separate(Group, c("Segment","DPH"), sep = "([_])") %>%
  mutate(DPH = factor(DPH, levels = c('1vs4','4vs7','7vs14','14vs21','21vs28','28vs35','35vs42'))) %>%
  mutate(Segment = factor(Segment, levels = c('Duodenum','Jejunum','Ileum','Cecum')))
  

dist_16S_combined 
head(dist_16S_combined)
dist_16S_combined$Compare <- rep(c(rep(c(1,2,3,4,5,6,7), each = 100)),4)
dist_ITS_combined$Compare <- rep(c(rep(c(1,2,3,4,5,6,7), each = 100)),4)
f11_1 <- ggplot(dist_16S_combined,aes(x = Compare,y = value, color = Segment, group = Segment)) +
  geom_smooth() +
  theme_bw() +
  ylab('Bray-Curtis dissimilarity') +
  theme(legend.position = "right",panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)) +
  scale_color_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) +
  scale_x_discrete(labels = c('1vs4','4vs7','7vs14','14vs21','21vs28','28vs35','35vs42')) + 
  ggtitle('Bacteria') +
  xlab('Comparision')
f11_1
f11_2 <- ggplot(dist_ITS_combined,aes(x = Compare,y = value, color = Segment,group = Segment)) +
  geom_smooth() +
  theme_bw() +
  ylab('Bray-Curtis dissimilarity') +
  theme(legend.position = "right",panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)) +
  scale_color_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) +
  scale_x_discrete(labels = c('1vs4','4vs7','7vs14','14vs21','21vs28','28vs35','35vs42')) + 
  ggtitle('Fungi') +
  xlab('Comparision')
f11_2
f11_3 <- ggplot(dist_16S_combined,aes(x = as.factor(Compare),y = value, color = Segment, group = Segment)) +
  geom_smooth() +
  theme_bw() +
  ylab('Bray-Curtis dissimilarity') +
  theme(legend.position = "right",panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)) +
  scale_color_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) +
  scale_x_discrete(labels = c('1vs4','4vs7','7vs14','14vs21','21vs28','28vs35','35vs42')) + 
  ggtitle('Bacteria') +
  xlab('Comparision')
f11_3
f11_4 <- ggplot(dist_ITS_combined,aes(x = as.factor(Compare),y = value, color = Segment,group = Segment)) +
  geom_smooth() +
  theme_bw() +
  ylab('Bray-Curtis dissimilarity') +
  theme(legend.position = "right",panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)) +
  scale_color_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) +
  scale_x_discrete(labels = c('1vs4','4vs7','7vs14','14vs21','21vs28','28vs35','35vs42')) + 
  ggtitle('Fungi') +
  xlab('Comparision')
f11_4
# PCoA
raw_16S_L6 <- read.csv('~/Desktop/01maturity/material/02micro/16S/level-6_no_archaea.csv')
format_16S_L6 <- raw_16S_L6 %>%
  remove_rownames %>% 
  column_to_rownames(var="index") %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  mutate(across(Granulicella:uncultured.bacterium.40, ~ . / sum)) %>%
  select(-c('sum')) %>%
  mutate(DPH = rep(c(1,4,7,14,21,28,35,42), each = 40)) %>%
  mutate(Group1 = factor(Group1, levels = c('Duodenum','Jejunum','Ileum','Cecum'))) %>%
  select(-c('Group1','Group2','DPH'))
format_16S_L6

raw_ITS_L6 <- read.csv('~/Desktop/01maturity/material/02micro/ITS/level-6_for_bar.csv')
format_ITS_L6 <- raw_ITS_L6 %>%
  remove_rownames %>% 
  column_to_rownames(var="X") %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  mutate(across(X__:unidentified.53, ~ . / sum)) %>%
  select(-c('sum')) %>%
  mutate(DPH = rep(c(1,4,7,14,21,28,35,42), each = 40)) %>%
  mutate(Group1 = factor(Group1, levels = c('Duodenum','Jejunum','Ileum','Cecum'))) %>%
  select(-c('Group1','Group2','DPH'))
format_ITS_L6
format_16S_L6_AA <- format_16S_L6
format_ITS_L6_AA <- format_ITS_L6
ab_temp <- read.xlsx('~/Desktop/01maturity/material/01qpcr/01data_for_quantitive.xlsx',sheet = 11)

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
df_pcoa_16S_AA$Segment <- rep(c('Cecum','Duodenum','Ileum','Jejunum'),80)
df_pcoa_16S_AA$Segment <- factor(df_pcoa_16S_AA$Segment, levels = c('Duodenum', 'Jejunum', 'Ileum', 'Cecum'))
f11_5 <- ggplot(df_pcoa_16S_AA, aes(x = PCoA1, y = PCoA2, color = Segment)) +
  geom_point() +
  stat_ellipse() + 
  theme_bw() +
  scale_color_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) +
  xlab('PCoA 1 (23.5%)') +
  ylab('PCoA 2 (16.7%)') +
  ggtitle('Bacteria') +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor = element_line(colour=NA))
f11_5

vegdist_ITS_AA <- vegdist(as.matrix(format_ITS_L6_AA[,1:603]))
pcoa_ITS_AA <- cmdscale(vegdist_ITS_AA, eig = TRUE)
df_pcoa_ITS_AA <- as.data.frame(pcoa_ITS_AA$points)
colnames(df_pcoa_ITS_AA) <- c('PCoA1', 'PCoA2')
PCo1_ITS_AA <- round(eigenvals(pcoa_ITS_AA)[1]/sum(eigenvals(pcoa_ITS_AA)) * 100, 1)
PCo2_ITS_AA <- round(eigenvals(pcoa_ITS_AA)[2]/sum(eigenvals(pcoa_ITS_AA)) * 100, 1)
print(PCo1_ITS_AA)#19.6
print(PCo2_ITS_AA)#14.2
df_pcoa_ITS_AA$Segment <- rep(c('Cecum','Duodenum','Ileum','Jejunum'),80)
df_pcoa_ITS_AA$Segment <- factor(df_pcoa_ITS_AA$Segment, levels = c('Duodenum', 'Jejunum', 'Ileum', 'Cecum'))
f11_6 <- ggplot(df_pcoa_ITS_AA, aes(x = PCoA1, y = PCoA2, color = Segment)) +
  geom_point() +
  stat_ellipse() + 
  theme_bw() +
  scale_color_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) +
  xlab('PCoA 1 (19.6%)') +
  ylab('PCoA 2 (14.2%)') +
  ggtitle('Fungi') +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor = element_line(colour=NA))
f11_6
f11_5 + f11_6 + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

df_pcoa_16S_AA$DPH <- df_pcoa_ITS_AA$DPH <- rep(c(1,4,7,14,21,28,35,42), each = 40)
df_pcoa_16S_AA
f11_7 <- ggplot(df_pcoa_16S_AA, aes(x = as.factor(DPH), y = PCoA1,fill = Segment)) + 
  geom_boxplot(aes(color = 'white'),outlier.shape = NA) + 
  facet_wrap(.~Segment,nrow = 4) + 
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor = element_line(colour=NA),
        legend.position = 'none') + 
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) + 
  scale_color_manual(values = c('black')) +
  xlab('DPH') + 
  ggtitle('Bacteria')
f11_8 <- ggplot(df_pcoa_ITS_AA, aes(x = as.factor(DPH), y = PCoA1,fill = Segment)) + 
  geom_boxplot(aes(color = 'white'),outlier.shape = NA) + 
  facet_wrap(.~Segment,nrow = 4) + 
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),panel.grid.minor = element_line(colour=NA),
        legend.position = 'none') + 
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) + 
  scale_color_manual(values = c('black')) +
  xlab('DPH') + 
  ggtitle('Fungi')

pairwiseAdonis::pairwise.adonis(format_16S_L6_AA,rep(c('Cecum','Duodenum','Ileum','Jejunum'),80))
l_o <- '
AB
CD
CD
CD'

pairwiseAdonis::pairwise.adonis(format_ITS_L6_AA,rep(c('Cecum','Duodenum','Ileum','Jejunum'),80))

pairwise.wilcox.test(df_pcoa_16S_AA$PCoA1[df_pcoa_16S_AA$Segment=='Duodenum'],df_pcoa_16S_AA$DPH[df_pcoa_16S_AA$Segment=='Duodenum'],p.adjust.method = 'BH')

pairwise.wilcox.test(df_pcoa_16S_AA$PCoA1[df_pcoa_16S_AA$Segment=='Jejunum'],df_pcoa_16S_AA$DPH[df_pcoa_16S_AA$Segment=='Jejunum'],p.adjust.method = 'BH')

pairwise.wilcox.test(df_pcoa_16S_AA$PCoA1[df_pcoa_16S_AA$Segment=='Ileum'],df_pcoa_16S_AA$DPH[df_pcoa_16S_AA$Segment=='Ileum'],p.adjust.method = 'BH')

pairwise.wilcox.test(df_pcoa_16S_AA$PCoA1[df_pcoa_16S_AA$Segment=='Cecum'],df_pcoa_16S_AA$DPH[df_pcoa_16S_AA$Segment=='Cecum'],p.adjust.method = 'BH')

pairwise.wilcox.test(df_pcoa_ITS_AA$PCoA1[df_pcoa_ITS_AA$Segment=='Duodenum'],df_pcoa_ITS_AA$DPH[df_pcoa_ITS_AA$Segment=='Duodenum'],p.adjust.method = 'BH')

pairwise.wilcox.test(df_pcoa_ITS_AA$PCoA1[df_pcoa_ITS_AA$Segment=='Jejunum'],df_pcoa_ITS_AA$DPH[df_pcoa_ITS_AA$Segment=='Jejunum'],p.adjust.method = 'BH')

pairwise.wilcox.test(df_pcoa_ITS_AA$PCoA1[df_pcoa_ITS_AA$Segment=='Ileum'],df_pcoa_ITS_AA$DPH[df_pcoa_ITS_AA$Segment=='Ileum'],p.adjust.method = 'BH')

pairwise.wilcox.test(df_pcoa_ITS_AA$PCoA1[df_pcoa_ITS_AA$Segment=='Cecum'],df_pcoa_ITS_AA$DPH[df_pcoa_ITS_AA$Segment=='Cecum'],p.adjust.method = 'BH')

ggplot(ab_bf, aes(x = DPH, y = value,  color = Segment, group = Segment)) +
  geom_smooth(size = 1) +
  geom_jitter(width = 0.1,alpha = .5,size = 0.5) +
  scale_color_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) +
  theme_bw() +
  ylab('Log10(copies per gram)') +
  xlab("DPH") +
  scale_x_continuous(breaks = c(1,7,14,21,28,35,42)) +
  theme(legend.position = "none",panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  facet_wrap(.~variable)

f10_1 + f10_2 +  plot_layout(guides = 'collect')

f10_3 <- ggplot(alpha_all, aes(x = DPH, y = alpha_bacteria,  color = Segment, group = Segment)) +
  geom_smooth(size = 1) +
  geom_jitter(width = 0.1,alpha = .5,size = 0.5) +
  scale_color_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) +
  theme_bw() +
  ylab('Shannon index') +
  xlab("DPH") +
  scale_x_continuous(breaks = c(1,7,14,21,28,35,42)) +
  theme(legend.position = "none",panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_color_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) +
  ggtitle('Bacteria')
f10_4 <- ggplot(alpha_all, aes(x = DPH, y = alpha_fungi,  color = Segment, group = Segment)) +
  geom_smooth(size = 1) +
  geom_jitter(width = 0.1,alpha = .5,size = 0.5) +
  scale_color_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) +
  theme_bw() +
  ylab('Shannon index') +
  xlab("DPH") +
  scale_x_continuous(breaks = c(1,7,14,21,28,35,42)) +
  theme(legend.position = "none",panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  scale_color_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) +
  ggtitle('Fungi') 

f10_3 + f10_4

f8_1 + f8_2 + f8_3 + f8_4 + plot_layout(guides = 'collect',nrow = 1)
f8_1 + theme(legend.position = 'none') +
  f8_2 + theme(legend.position = 'none') +
  f8_3 + theme(legend.position = 'none') +
  f8_4+ theme(legend.position = 'none') +
  plot_layout(guides = 'collect', nrow = 1)

f9_1 + f9_2 + f9_3 + f9_4 + plot_layout(guides = 'collect',nrow = 1)
f9_1 + theme(legend.position = 'none') +
  f9_2 + theme(legend.position = 'none') +
  f9_3 + theme(legend.position = 'none') +
  f9_4+ theme(legend.position = 'none') +
  plot_layout(guides = 'collect', nrow = 1)

ggplot(ab_bf, aes(x = DPH, y = value,  color = Segment, group = Segment)) +
  geom_smooth(size = 1) +
  geom_jitter(width = 0.1,alpha = .5,size = 0.5) +
  scale_color_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) +
  theme_bw() +
  ylab('Log10(copies per gram)') +
  xlab("DPH") +
  scale_x_continuous(breaks = c(1,7,14,21,28,35,42)) +
  theme(legend.position = "none",panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  facet_wrap(.~variable)

ggplot(ab_bf, aes(x = DPH, y = value,  color = Segment, group = Segment)) +
  geom_smooth(size = 1) +
  geom_jitter(width = 0.1,alpha = .5,size = 0.5) +
  scale_color_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) +
  theme_bw() +
  ylab('Log10(copies per gram)') +
  xlab("DPH") +
  scale_x_continuous(breaks = c(1,7,14,21,28,35,42)) +
  theme(legend.position = "none",panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
  facet_wrap(.~variable)

f10_3 + f10_4
