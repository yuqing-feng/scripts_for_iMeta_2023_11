# Part 4
# Find the associated 
rm(list = ls())
library(DirichletMultinomial) #1.34.0
library(dplyr) 
library(reshape2)
library(tibble)
library(openxlsx)
library(ggExtra)
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
  select(-c('Group1','Group2'))%>%
  as.matrix() 
format_ITS_L6

group_store <- data.frame(raw_16S_L6$index,substring(raw_16S_L6$index,1,3),substring(raw_16S_L6$index,7,7))
group_store
colnames(group_store) <- c('ID','DPH','Segment')
data_16S_D <- format_16S_L6[group_store$Segment=='D',] 
data_16S_J <- format_16S_L6[group_store$Segment=='J',]
data_16S_I <- format_16S_L6[group_store$Segment=='I',]
data_16S_C <- format_16S_L6[group_store$Segment=='C',]

data_ITS_D <- format_ITS_L6[group_store$Segment=='D',] 
data_ITS_J <- format_ITS_L6[group_store$Segment=='J',]
data_ITS_I <- format_ITS_L6[group_store$Segment=='I',]
data_ITS_C <- format_ITS_L6[group_store$Segment=='C',]

data_16S_I
data_16S_C

set.seed(1234)
fit_16S_D <- lapply(1:6, dmn, count = data_16S_D, verbose = TRUE)
cluster_16S_D <- data.frame(apply(mixture(fit_16S_D[[which.min(sapply(fit_16S_D, laplace))]]), 1, which.max)) 
set.seed(1234)
fit_16S_J <- lapply(1:6, dmn, count = data_16S_J, verbose = TRUE)
cluster_16S_J <- data.frame(apply(mixture(fit_16S_J[[which.min(sapply(fit_16S_J, laplace))]]), 1, which.max))
set.seed(1234)
fit_16S_I <- lapply(1:6, dmn, count = data_16S_I, verbose = TRUE)
cluster_16S_I <- data.frame(apply(mixture(fit_16S_I[[which.min(sapply(fit_16S_I, laplace))]]), 1, which.max))
set.seed(1234)
fit_16S_C <- lapply(1:6, dmn, count = data_16S_C, verbose = TRUE)
cluster_16S_C <- data.frame(apply(mixture(fit_16S_C[[which.min(sapply(fit_16S_C, laplace))]]), 1, which.max))

set.seed(1234)
fit_ITS_D <- lapply(1:6, dmn, count = data_ITS_D, verbose = TRUE)
cluster_ITS_D <- data.frame(apply(mixture(fit_ITS_D[[which.min(sapply(fit_ITS_D, laplace))]]), 1, which.max)) 
set.seed(1234)
fit_ITS_J <- lapply(1:6, dmn, count = data_ITS_J, verbose = TRUE)
cluster_ITS_J <- data.frame(apply(mixture(fit_ITS_J[[which.min(sapply(fit_ITS_J, laplace))]]), 1, which.max))
set.seed(1234)
fit_ITS_I <- lapply(1:6, dmn, count = data_ITS_I, verbose = TRUE)
cluster_ITS_I <- data.frame(apply(mixture(fit_ITS_I[[which.min(sapply(fit_ITS_I, laplace))]]), 1, which.max))
set.seed(1234)
fit_ITS_C <- lapply(1:6, dmn, count = data_ITS_C, verbose = TRUE)
cluster_ITS_C <- data.frame(apply(mixture(fit_ITS_C[[which.min(sapply(fit_ITS_C, laplace))]]), 1, which.max))

cluster_16S <- cbind(cluster_16S_D,cluster_16S_J,cluster_16S_I,cluster_16S_C)
cluster_ITS <- cbind(cluster_ITS_D,cluster_ITS_J,cluster_ITS_I,cluster_ITS_C)


colnames(cluster_16S) <- colnames(cluster_ITS) <- c('Duodenum','Jejunum','Ileum','Cecum')
cluster_16S$DPH <- rep(c('D1','D4','D7','D14','D21','D28','D35','D42'),each = 10)
# f21_1 <- cluster_16S %>%
#   select(c('Duodenum','DPH')) %>%
#   table() %>%
#   melt() %>%
#   mutate(DPH = factor(DPH, levels = c('D1','D4','D7','D14','D21','D28','D35','D42')))%>%
#   ggplot(aes(x = DPH,y= value, fill = as.factor(Duodenum))) +
#   geom_bar(stat = 'identity') +
#   scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087')) +
#   guides(fill = guide_legend(title = 'Duodenum')) +
#   theme_bw() +
#   ylab('Number') + 
#   theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) +
#   ggtitle('Duodenum (Bacteria)') +
#   scale_x_discrete(label = c(1,4,7,14,21,28,35,42))
# f21_1
# f21_2 <- cluster_16S %>%
#   select(c('Jejunum','DPH')) %>%
#   table() %>%
#   melt() %>%
#   mutate(DPH = factor(DPH, levels = c('D1','D4','D7','D14','D21','D28','D35','D42')))%>%
#   ggplot(aes(x = DPH,y= value, fill = as.factor(Jejunum))) +
#   geom_bar(stat = 'identity') +
#   scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087')) +
#   guides(fill = guide_legend(title = 'Jejunum')) +
#   theme_bw() +
#   ylab('Number') + 
#   theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA))+
#   ggtitle('Jejunum (Bacteria)') +
#   scale_x_discrete(label = c(1,4,7,14,21,28,35,42))
# f21_2
# f21_3 <- cluster_16S %>%
#   select(c('Ileum','DPH')) %>%
#   table() %>%
#   melt() %>%
#   mutate(DPH = factor(DPH, levels = c('D1','D4','D7','D14','D21','D28','D35','D42')))%>%
#   ggplot(aes(x = DPH,y= value, fill = as.factor(Ileum))) +
#   geom_bar(stat = 'identity') +
#   scale_fill_manual(values = c('#E64B35','#00A087')) +
#   guides(fill = guide_legend(title = 'Ileum')) +
#   theme_bw() +
#   ylab('Number') + 
#   theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA))+
#   ggtitle('Ileum (Bacteria)') +
#   scale_x_discrete(label = c(1,4,7,14,21,28,35,42))
# f21_3
# f21_4 <- cluster_16S %>%
#   select(c('Cecum','DPH')) %>%
#   table() %>%
#   melt() %>%
#   mutate(DPH = factor(DPH, levels = c('D1','D4','D7','D14','D21','D28','D35','D42')))%>%
#   ggplot(aes(x = DPH,y= value, fill = as.factor(Cecum))) +
#   geom_bar(stat = 'identity') +
#   scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087')) +
#   guides(fill = guide_legend(title = 'Cecum')) +
#   theme_bw() +
#   ylab('Number') + 
#   theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA))+
#   ggtitle('Cecum (Bacteria)') +
#   scale_x_discrete(label = c(1,4,7,14,21,28,35,42))
# f21_4
# cluster_ITS$DPH <- rep(c('D1','D4','D7','D14','D21','D28','D35','D42'),each = 10)
# f21_5 <- cluster_ITS %>%
#   select(c('Duodenum','DPH')) %>%
#   table() %>%
#   melt() %>%
#   mutate(DPH = factor(DPH, levels = c('D1','D4','D7','D14','D21','D28','D35','D42')))%>%
#   ggplot(aes(x = DPH,y= value, fill = as.factor(Duodenum))) +
#   geom_bar(stat = 'identity') +
#   scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087')) +
#   guides(fill = guide_legend(title = 'Duodenum')) +
#   theme_bw() +
#   ylab('Number') + 
#   theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA))+
#   ggtitle('Duodenum (Fungi)') +
#   scale_x_discrete(label = c(1,4,7,14,21,28,35,42))
# f21_5
# f21_6 <- cluster_ITS %>%
#   select(c('Jejunum','DPH')) %>%
#   table() %>%
#   melt() %>%
#   mutate(DPH = factor(DPH, levels = c('D1','D4','D7','D14','D21','D28','D35','D42')))%>%
#   ggplot(aes(x = DPH,y= value, fill = as.factor(Jejunum))) +
#   geom_bar(stat = 'identity') +
#   scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087')) +
#   guides(fill = guide_legend(title = 'Jejunum')) +
#   theme_bw() +
#   ylab('Number') + 
#   theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA))+
#   ggtitle('Jejunum (Fungi)') +
#   scale_x_discrete(label = c(1,4,7,14,21,28,35,42))
# f21_6
# f21_7 <- cluster_ITS %>%
#   select(c('Ileum','DPH')) %>%
#   table() %>%
#   melt() %>%
#   mutate(DPH = factor(DPH, levels = c('D1','D4','D7','D14','D21','D28','D35','D42')))%>%
#   ggplot(aes(x = DPH,y= value, fill = as.factor(Ileum))) +
#   geom_bar(stat = 'identity') +
#   scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087')) +
#   guides(fill = guide_legend(title = 'Ileum')) +
#   theme_bw() +
#   ylab('Number') + 
#   theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA))+
#   ggtitle('Ileum (Fungi)') +
#   scale_x_discrete(label = c(1,4,7,14,21,28,35,42))
# f21_7
# f21_8 <- cluster_ITS %>%
#   select(c('Cecum','DPH')) %>%
#   table() %>%
#   melt() %>%
#   mutate(DPH = factor(DPH, levels = c('D1','D4','D7','D14','D21','D28','D35','D42')))%>%
#   ggplot(aes(x = DPH,y= value, fill = as.factor(Cecum))) +
#   geom_bar(stat = 'identity') +
#   scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087')) +
#   guides(fill = guide_legend(title = 'Cecum')) +
#   theme_bw() +
#   ylab('Number') + 
#   theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA))+
#   ggtitle('Cecum (Fungi)') +
#   scale_x_discrete(label = c(1,4,7,14,21,28,35,42))
# f21_1 / f21_2 / f21_3 / f21_4 | f21_5 / f21_6 / f21_7 / f21_8

format_16S_L6
dim(format_16S_L6)
binary_16S <- format_16S_L6
binary_ITS <- format_ITS_L6
for(i in 1:nrow(format_16S_L6)){
  for(j in 1:ncol(format_16S_L6)){
    if(format_16S_L6[i,j]!=0){
      binary_16S[i,j] <- 1
    }
  }
}
for(i in 1:nrow(format_ITS_L6)){
  for(j in 1:ncol(format_ITS_L6)){
    if(format_ITS_L6[i,j]!=0){
      binary_ITS[i,j] <- 1
    }
  }
}

ids_bacteria <- read.xlsx('taxa_bacteria.xlsx',sheet = 1)
ids_fungi <- read.xlsx('taxa_fungi.xlsx', sheet = 2)
head(ids_bacteria)
colnames(binary_16S) <- ids_bacteria$Genus
colnames(binary_ITS) <- ids_fungi$Genus
binary_16S_D <- cbind(binary_16S,group_store) %>%
  filter(Segment=='D') %>%
  select(-c(Segment,ID)) %>%
  group_by(DPH) %>% 
  summarise(across(everything(), sum)) %>%
  column_to_rownames(var='DPH')

binary_16S_J <- cbind(binary_16S, group_store) %>%
  filter(Segment=='J') %>%
  select(-c(Segment,ID)) %>%
  group_by(DPH) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames(var='DPH')

binary_16S_I <- cbind(binary_16S, group_store) %>%
  filter(Segment=='I') %>%
  select(-c(Segment,ID)) %>%
  group_by(DPH) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames(var='DPH')

binary_16S_C <- cbind(binary_16S, group_store) %>%
  filter(Segment =='C') %>%
  select(-c(Segment,ID)) %>%
  group_by(DPH) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames(var='DPH')

binary_ITS_D <- cbind(binary_ITS,group_store) %>%
  filter(Segment=='D') %>%
  select(-c(Segment,ID)) %>%
  group_by(DPH) %>% 
  summarise(across(everything(), sum)) %>%
  column_to_rownames(var='DPH')

binary_ITS_J <- cbind(binary_ITS, group_store) %>%
  filter(Segment=='J') %>%
  select(-c(Segment,ID)) %>%
  group_by(DPH) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames(var='DPH')

binary_ITS_I <- cbind(binary_ITS, group_store) %>%
  filter(Segment=='I') %>%
  select(-c(Segment,ID)) %>%
  group_by(DPH) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames(var='DPH')

binary_ITS_C <- cbind(binary_ITS, group_store) %>%
  filter(Segment =='C') %>%
  select(-c(Segment,ID)) %>%
  group_by(DPH) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames(var='DPH')

# pheatmap::pheatmap(binary_16S_I,cluster_rows = F)
con_2_2_pre <- function(data){
  data <- data[,colSums(data)!=0]
  data1 <- data
  for(i in 1:nrow(data)){
    for(j in 1:ncol(data)){
      if(data[i,j]>0 & data[i,j]<9){
        data1[i,j] <- 'B'
      }
    }
  }
  data1[data==9] <- 'C'
  data1[data==10] <- 'C'
  data1[data==0] <- 'A'
  return(data1)
}
bin_16S_C2 <- con_2_2_pre(binary_16S_C)
bin_16S_D2 <- con_2_2_pre(binary_16S_D)
bin_16S_I2 <- con_2_2_pre(binary_16S_I)
bin_16S_J2 <- con_2_2_pre(binary_16S_J)
bin_ITS_C2 <- con_2_2_pre(binary_ITS_C) 
bin_ITS_D2 <- con_2_2_pre(binary_ITS_D)
bin_ITS_I2 <- con_2_2_pre(binary_ITS_I)
bin_ITS_J2 <- con_2_2_pre(binary_ITS_J)



count_16S_C1 <- t(bin_16S_C2) %>% as.data.frame() %>% group_by_all() %>% summarise(COUNT = n()) %>% arrange(desc(COUNT))
count_16S_D1 <- t(bin_16S_D2) %>% as.data.frame() %>% group_by_all() %>% summarise(COUNT = n()) %>% arrange(desc(COUNT))
count_16S_I1 <- t(bin_16S_I2) %>% as.data.frame() %>% group_by_all() %>% summarise(COUNT = n()) %>% arrange(desc(COUNT))
count_16S_J1 <- t(bin_16S_J2) %>% as.data.frame() %>% group_by_all() %>% summarise(COUNT = n()) %>% arrange(desc(COUNT))
count_ITS_C1 <- t(bin_ITS_C2) %>% as.data.frame() %>% group_by_all() %>% summarise(COUNT = n()) %>% arrange(desc(COUNT))
count_ITS_D1 <- t(bin_ITS_D2) %>% as.data.frame() %>% group_by_all() %>% summarise(COUNT = n()) %>% arrange(desc(COUNT))
count_ITS_I1 <- t(bin_ITS_I2) %>% as.data.frame() %>% group_by_all() %>% summarise(COUNT = n()) %>% arrange(desc(COUNT))
count_ITS_J1 <- t(bin_ITS_J2) %>% as.data.frame() %>% group_by_all() %>% summarise(COUNT = n()) %>% arrange(desc(COUNT))


count_16S_D2 <- data.frame(sort(table(apply(bin_16S_D2, 2, paste0, collapse="")),decreasing = T))
count_16S_J2 <- data.frame(sort(table(apply(bin_16S_J2, 2, paste0, collapse="")),decreasing = T))
count_16S_I2 <- data.frame(sort(table(apply(bin_16S_I2, 2, paste0, collapse="")),decreasing = T))
count_16S_C2 <- data.frame(sort(table(apply(bin_16S_C2, 2, paste0, collapse="")),decreasing = T))
count_ITS_D2 <- data.frame(sort(table(apply(bin_ITS_D2, 2, paste0, collapse="")),decreasing = T))
count_ITS_J2 <- data.frame(sort(table(apply(bin_ITS_J2, 2, paste0, collapse="")),decreasing = T))
count_ITS_I2 <- data.frame(sort(table(apply(bin_ITS_I2, 2, paste0, collapse="")),decreasing = T))
count_ITS_C2 <- data.frame(sort(table(apply(bin_ITS_C2, 2, paste0, collapse="")),decreasing = T))

count_16S_C3 <- data.frame(count_16S_C1,count_16S_C2)
count_16S_D3 <- data.frame(count_16S_D1,count_16S_D2)
count_16S_I3 <- data.frame(count_16S_I1,count_16S_I2)
count_16S_J3 <- data.frame(count_16S_J1,count_16S_J2)
count_ITS_C3 <- data.frame(count_ITS_C1,count_ITS_C2)
count_ITS_D3 <- data.frame(count_ITS_D1,count_ITS_D2)
count_ITS_I3 <- data.frame(count_ITS_I1,count_ITS_I2)
count_ITS_J3 <- data.frame(count_ITS_J1,count_ITS_J2)

conv_2_num <- function(data){
  data0 <- data[,1:8]
  data1 <- matrix(,nrow = nrow(data),ncol = 8)
  print(dim(data0))
  rownames(data1) <- data[,10]
  for(i in 1:nrow(data0)){
    for(j in 1:8){
      if(data0[i,j]=='A'){data1[i,j] <- 0}
      else if(data0[i,j]=='B'){data1[i,j] <- 1}
      else if(data0[i,j]=='C'){data1[i,j] <- 2}
    }
  }
  print(nrow(data0))
  ID <- rep(rownames(data1),8)
  Num <- c(data1[,1],data1[,2],data1[,3],data1[,4],data1[,5],data1[,6],data1[,7],data1[,8])
  DPH <- rep(c('1','4','7','14','21','28','35','42'),each = nrow(data1))
  data2 <- data.frame(ID,DPH,Num)
  return(data2)
}
count_16S_C4 <- conv_2_num(count_16S_C3)
count_16S_D4 <- conv_2_num(count_16S_D3)
count_16S_I4 <- conv_2_num(count_16S_I3)
count_16S_J4 <- conv_2_num(count_16S_J3)
count_ITS_C4 <- conv_2_num(count_ITS_C3)
count_ITS_D4 <- conv_2_num(count_ITS_D3)
count_ITS_I4 <- conv_2_num(count_ITS_I3)
count_ITS_J4 <- conv_2_num(count_ITS_J3)

count_16S_C4$ID <- factor(count_16S_C4$ID, levels = rev(count_16S_C3$Var1))
count_16S_C3$Var1 <- factor(count_16S_C3$Var1, levels = rev(count_16S_C3$Var1))
count_16S_D4$ID <- factor(count_16S_D4$ID, levels = rev(count_16S_D3$Var1))
count_16S_D3$Var1 <- factor(count_16S_D3$Var1, levels = rev(count_16S_D3$Var1))
count_16S_I4$ID <- factor(count_16S_I4$ID, levels = rev(count_16S_I3$Var1))
count_16S_I3$Var1 <- factor(count_16S_I3$Var1, levels = rev(count_16S_I3$Var1))
count_16S_J4$ID <- factor(count_16S_J4$ID, levels = rev(count_16S_J3$Var1))
count_16S_J3$Var1 <- factor(count_16S_J3$Var1, levels = rev(count_16S_J3$Var1))

count_ITS_C4$ID <- factor(count_ITS_C4$ID, levels = rev(count_ITS_C3$Var1))
count_ITS_C3$Var1 <- factor(count_ITS_C3$Var1, levels = rev(count_ITS_C3$Var1))
count_ITS_D4$ID <- factor(count_ITS_D4$ID, levels = rev(count_ITS_D3$Var1))
count_ITS_D3$Var1 <- factor(count_ITS_D3$Var1, levels = rev(count_ITS_D3$Var1))
count_ITS_I4$ID <- factor(count_ITS_I4$ID, levels = rev(count_ITS_I3$Var1))
count_ITS_I3$Var1 <- factor(count_ITS_I3$Var1, levels = rev(count_ITS_I3$Var1))
count_ITS_J4$ID <- factor(count_ITS_J4$ID, levels = rev(count_ITS_J3$Var1))
count_ITS_J3$Var1 <- factor(count_ITS_J3$Var1, levels = rev(count_ITS_J3$Var1))

count_16S_C4$DPH <- factor(count_16S_C4$DPH, levels = c(1,4,7,14,21,28,35,42))
f_C_b1 <- ggplot(count_16S_C4,aes(x = as.factor(DPH),y = ID, color = as.factor(Num))) + 
  geom_point() + 
  scale_color_manual(values = c('white','#3C54881f','#3C5488ff')) +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        panel.border = element_blank(),axis.text.y = element_blank(),axis.ticks.length.y = unit(0,'cm'),
        axis.ticks.length.x = unit(0,'cm'),
        legend.position = 'none') + 
  xlab('DPH') + ylab('Pattern')
f_C_b2 <- ggplot(count_16S_C3,aes(x = Var1,y = Freq)) +
  geom_bar(stat = 'identity',fill = '#3C5488') +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        axis.text.y = element_blank(),axis.ticks.length.y = unit(0,'cm'),
        legend.position = 'none') +
  ylab('Counts') + xlab('')

count_16S_D4$DPH <- factor(count_16S_D4$DPH, levels = c(1,4,7,14,21,28,35,42))
f_D_b1 <- ggplot(count_16S_D4,aes(x = as.factor(DPH),y = ID, color = as.factor(Num))) + 
  geom_point() + 
  scale_color_manual(values = c('white','#E64B351f','#E64B35ff')) +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        panel.border = element_blank(),axis.text.y = element_blank(),axis.ticks.length.y = unit(0,'cm'),
        axis.ticks.length.x = unit(0,'cm'),
        legend.position = 'none') + 
  xlab('DPH') + ylab('Pattern')
f_D_b2 <- ggplot(count_16S_D3,aes(x = Var1,y = Freq)) +
  geom_bar(stat = 'identity',fill = '#E64B35') +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        axis.text.y = element_blank(),axis.ticks.length.y = unit(0,'cm'),
        legend.position = 'none') +
  ylab('Counts') + xlab('')
f_D_b1 + f_D_b2

count_16S_I4$DPH <- factor(count_16S_I4$DPH, levels = c(1,4,7,14,21,28,35,42))
f_I_b1 <- ggplot(count_16S_I4,aes(x = as.factor(DPH),y = ID, color = as.factor(Num))) + 
  geom_point() + 
  scale_color_manual(values = c('white','#00A0871f','#00A087ff')) +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        panel.border = element_blank(),axis.text.y = element_blank(),axis.ticks.length.y = unit(0,'cm'),
        axis.ticks.length.x = unit(0,'cm'),
        legend.position = 'none') + 
  xlab('DPH') + ylab('Pattern')
f_I_b2 <- ggplot(count_16S_I3,aes(x = Var1,y = Freq)) +
  geom_bar(stat = 'identity',fill = '#00A087') +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        axis.text.y = element_blank(),axis.ticks.length.y = unit(0,'cm'),
        legend.position = 'none') +
  ylab('Counts') + xlab('')


count_16S_J4$DPH <- factor(count_16S_J4$DPH, levels = c(1,4,7,14,21,28,35,42))
f_J_b1 <- ggplot(count_16S_J4,aes(x = as.factor(DPH),y = ID, color = as.factor(Num))) + 
  geom_point() + 
  scale_color_manual(values = c('white','#4DBBD51f','#4DBBD5ff')) +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        panel.border = element_blank(),axis.text.y = element_blank(),axis.ticks.length.y = unit(0,'cm'),
        axis.ticks.length.x = unit(0,'cm'),
        legend.position = 'none') + 
  xlab('DPH') + ylab('Pattern')
f_J_b2 <- ggplot(count_16S_J3,aes(x = Var1,y = Freq)) +
  geom_bar(stat = 'identity',fill = '#4DBBD5') +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        axis.text.y = element_blank(),axis.ticks.length.y = unit(0,'cm'),
        legend.position = 'none') +
  ylab('Counts') + xlab('')
svg('~/Desktop/figure materials/2nd/figure S7.svg',width = 10.5,height = 12)
f_D_b1 + f_D_b2 +
  f_J_b1 + f_J_b2 + 
  f_I_b1 + f_I_b2 + 
  f_C_b1 + f_C_b2 + plot_layout(nrow = 1)
dev.off()
count_ITS_C4$DPH <- factor(count_ITS_C4$DPH, levels = c(1,4,7,14,21,28,35,42))
f_C_f1 <- ggplot(count_ITS_C4,aes(x = as.factor(DPH),y = ID, color = as.factor(Num))) + 
  geom_point() + 
  scale_color_manual(values = c('white','#3C54881f','#3C5488ff')) +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        panel.border = element_blank(),axis.text.y = element_blank(),axis.ticks.length.y = unit(0,'cm'),
        axis.ticks.length.x = unit(0,'cm'),
        legend.position = 'none') + 
  xlab('DPH') + ylab('Pattern')
f_C_f2 <- ggplot(count_ITS_C3,aes(x = Var1,y = Freq)) +
  geom_bar(stat = 'identity',fill = '#3C5488') +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        axis.text.y = element_blank(),axis.ticks.length.y = unit(0,'cm'),
        legend.position = 'none') +
  ylab('Counts') + xlab('')

count_ITS_D4$DPH <- factor(count_ITS_D4$DPH, levels = c(1,4,7,14,21,28,35,42))
f_D_f1 <- ggplot(count_ITS_D4,aes(x = as.factor(DPH),y = ID, color = as.factor(Num))) + 
  geom_point() + 
  scale_color_manual(values = c('white','#E64B351f','#E64B35ff')) +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        panel.border = element_blank(),axis.text.y = element_blank(),axis.ticks.length.y = unit(0,'cm'),
        axis.ticks.length.x = unit(0,'cm'),
        legend.position = 'none') + 
  xlab('DPH') + ylab('Pattern')
f_D_f2 <- ggplot(count_ITS_D3,aes(x = Var1,y = Freq)) +
  geom_bar(stat = 'identity',fill = '#E64B35') +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        axis.text.y = element_blank(),axis.ticks.length.y = unit(0,'cm'),
        legend.position = 'none') +
  ylab('Counts') + xlab('')
f_D_f1 + f_D_f2

count_ITS_I4$DPH <- factor(count_ITS_I4$DPH, levels = c(1,4,7,14,21,28,35,42))
f_I_f1 <- ggplot(count_ITS_I4,aes(x = as.factor(DPH),y = ID, color = as.factor(Num))) + 
  geom_point() + 
  scale_color_manual(values = c('white','#00A0871f','#00A087ff')) +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        panel.border = element_blank(),axis.text.y = element_blank(),axis.ticks.length.y = unit(0,'cm'),
        axis.ticks.length.x = unit(0,'cm'),
        legend.position = 'none') + 
  xlab('DPH') + ylab('Pattern')
f_I_f2 <- ggplot(count_ITS_I3,aes(x = Var1,y = Freq)) +
  geom_bar(stat = 'identity',fill = '#00A087') +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        axis.text.y = element_blank(),axis.ticks.length.y = unit(0,'cm'),
        legend.position = 'none') +
  ylab('Counts') + xlab('')


count_ITS_J4$DPH <- factor(count_ITS_J4$DPH, levels = c(1,4,7,14,21,28,35,42))
f_J_f1 <- ggplot(count_ITS_J4,aes(x = as.factor(DPH),y = ID, color = as.factor(Num))) + 
  geom_point() + 
  scale_color_manual(values = c('white','#4DBBD51f','#4DBBD5ff')) +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        panel.border = element_blank(),axis.text.y = element_blank(),axis.ticks.length.y = unit(0,'cm'),
        axis.ticks.length.x = unit(0,'cm'),
        legend.position = 'none') + 
  xlab('DPH') + ylab('Pattern')
f_J_f2 <- ggplot(count_ITS_J3,aes(x = Var1,y = Freq)) +
  geom_bar(stat = 'identity',fill = '#4DBBD5') +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        axis.text.y = element_blank(),axis.ticks.length.y = unit(0,'cm'),
        legend.position = 'none') +
  ylab('Counts') + xlab('')
svg('~/Desktop/figure materials/2nd/figure S8.svg',width = 10.5,height = 12)
f_D_f1 + f_D_f2 +
  f_J_f1 + f_J_f2 + 
  f_I_f1 + f_I_f2 + 
  f_C_f1 + f_C_f2 + plot_layout(nrow = 1)
dev.off()

count_16S_C5 <- count_16S_C3[count_16S_C3$D01=='C'|count_16S_C3$D42=='C',]
count_16S_C5
count_16S_D5 <- count_16S_D3[count_16S_D3$D01=='C'|count_16S_D3$D42=='C',]
count_16S_D5
count_16S_I5 <- count_16S_I3[count_16S_I3$D01=='C'|count_16S_I3$D42=='C',]
count_16S_I5
count_16S_J5 <- count_16S_J3[count_16S_J3$D01=='C'|count_16S_J3$D42=='C',]
count_16S_J5

count_ITS_C5 <- count_ITS_C3[count_ITS_C3$D01=='C'|count_ITS_C3$D42=='C',]
count_ITS_C5
count_ITS_D5 <- count_ITS_D3[count_ITS_D3$D01=='C'|count_ITS_D3$D42=='C',]
count_ITS_D5
count_ITS_I5 <- count_ITS_I3[count_ITS_I3$D01=='C'|count_ITS_I3$D42=='C',]
count_ITS_I5
count_ITS_J5 <- count_ITS_J3[count_ITS_J3$D01=='C'|count_ITS_J3$D42=='C',]
count_ITS_J5
sort(table(c(count_16S_C5$Var1,count_16S_D5$Var1,count_16S_I5$Var1,count_16S_J5$Var1,
             count_ITS_C5$Var1,count_ITS_D5$Var1,count_ITS_I5$Var1,count_ITS_J5$Var1)))

names_16S_C <- data.frame(apply(format(bin_16S_C2), 2, paste, collapse=""))
names_16S_D <- data.frame(apply(format(bin_16S_D2), 2, paste, collapse=''))
names_16S_I <- data.frame(apply(format(bin_16S_I2), 2, paste, collapse=''))
names_16S_J <- data.frame(apply(format(bin_16S_J2), 2, paste, collapse=''))
names_ITS_C <- data.frame(apply(format(bin_ITS_C2), 2, paste, collapse=""))
names_ITS_D <- data.frame(apply(format(bin_ITS_D2), 2, paste, collapse=''))
names_ITS_I <- data.frame(apply(format(bin_ITS_I2), 2, paste, collapse=''))
names_ITS_J <- data.frame(apply(format(bin_ITS_J2), 2, paste, collapse=''))
rownames(names_16S_C)[names_16S_C$apply.format.bin_16S_C2...2..paste..collapse......=='CCCCCCCC']
rownames(names_16S_D)[names_16S_D$apply.format.bin_16S_D2...2..paste..collapse......=='CCCCCCCC']
rownames(names_16S_I)[names_16S_I$apply.format.bin_16S_I2...2..paste..collapse......=='CCCCCCCC']
rownames(names_16S_J)[names_16S_J$apply.format.bin_16S_J2...2..paste..collapse......=='CCCCCCCC']
rownames(names_ITS_C)[names_ITS_C$apply.format.bin_ITS_C2...2..paste..collapse......=='CCCCCCCC']
rownames(names_ITS_D)[names_ITS_D$apply.format.bin_ITS_D2...2..paste..collapse......=='CCCCCCCC']
rownames(names_ITS_I)[names_ITS_I$apply.format.bin_ITS_I2...2..paste..collapse......=='CCCCCCCC']
rownames(names_ITS_J)[names_ITS_J$apply.format.bin_ITS_J2...2..paste..collapse......=='CCCCCCCC']
emerge <- c('ABBBBCCC','ABBBBBCC','AABBBCCC','AABBBBCC', 'AAAABBBC', 'ABBBCCCC', 'ABBCCCCC', 'AABBCCCC','AABCCCCC',
            'ABCCCCCC','AACCCCCC', 'ACCCCCCC', 'BBBBCCCC', 'BBBBBBBC','BCCCCCCC')
dismiss <- c('CBBBBBAA','CBBBBBBA','CAAAAAAA','CCCBBBBB','CCBBBBBB','CBBBBBBB','CCCCCCCB')



rownames(names_16S_C)[names_16S_C$apply.format.bin_16S_C2...2..paste..collapse......%in%emerge]
rownames(names_16S_D)[names_16S_D$apply.format.bin_16S_D2...2..paste..collapse......%in%emerge]
rownames(names_16S_I)[names_16S_I$apply.format.bin_16S_I2...2..paste..collapse......%in%emerge]
rownames(names_16S_J)[names_16S_J$apply.format.bin_16S_J2...2..paste..collapse......%in%emerge]
rownames(names_ITS_C)[names_ITS_C$apply.format.bin_ITS_C2...2..paste..collapse......%in%emerge]
rownames(names_ITS_D)[names_ITS_D$apply.format.bin_ITS_D2...2..paste..collapse......%in%emerge]
rownames(names_ITS_I)[names_ITS_I$apply.format.bin_ITS_I2...2..paste..collapse......%in%emerge]
rownames(names_ITS_J)[names_ITS_J$apply.format.bin_ITS_J2...2..paste..collapse......%in%emerge]

rownames(names_16S_C)[names_16S_C$apply.format.bin_16S_C2...2..paste..collapse......%in%dismiss]
rownames(names_16S_D)[names_16S_D$apply.format.bin_16S_D2...2..paste..collapse......%in%dismiss]
rownames(names_16S_I)[names_16S_I$apply.format.bin_16S_I2...2..paste..collapse......%in%dismiss]
rownames(names_16S_J)[names_16S_J$apply.format.bin_16S_J2...2..paste..collapse......%in%dismiss]
rownames(names_ITS_C)[names_ITS_C$apply.format.bin_ITS_C2...2..paste..collapse......%in%dismiss]
rownames(names_ITS_D)[names_ITS_D$apply.format.bin_ITS_D2...2..paste..collapse......%in%dismiss]
rownames(names_ITS_I)[names_ITS_I$apply.format.bin_ITS_I2...2..paste..collapse......%in%dismiss]
rownames(names_ITS_J)[names_ITS_J$apply.format.bin_ITS_J2...2..paste..collapse......%in%dismiss]

for(i in 1:15){print(length(rownames(names_16S_D)[names_16S_D$apply.format.bin_16S_D2...2..paste..collapse......%in%emerge[i]]))}
for(i in 1:15){print(length(rownames(names_16S_J)[names_16S_J$apply.format.bin_16S_J2...2..paste..collapse......%in%emerge[i]]))}
for(i in 1:15){print(length(rownames(names_16S_I)[names_16S_I$apply.format.bin_16S_I2...2..paste..collapse......%in%emerge[i]]))}
for(i in 1:15){print(length(rownames(names_16S_C)[names_16S_C$apply.format.bin_16S_C2...2..paste..collapse......%in%emerge[i]]))}
for(i in 1:7){print(length(rownames(names_16S_D)[names_16S_D$apply.format.bin_16S_D2...2..paste..collapse......%in%dismiss[i]]))}
for(i in 1:7){print(length(rownames(names_16S_J)[names_16S_J$apply.format.bin_16S_J2...2..paste..collapse......%in%dismiss[i]]))}
for(i in 1:7){print(length(rownames(names_16S_I)[names_16S_I$apply.format.bin_16S_I2...2..paste..collapse......%in%dismiss[i]]))}
for(i in 1:7){print(length(rownames(names_16S_C)[names_16S_C$apply.format.bin_16S_C2...2..paste..collapse......%in%dismiss[i]]))}

for(i in 1:15){
  print(length(rownames(names_16S_D)[names_16S_D$apply.format.bin_16S_D2...2..paste..collapse......%in%emerge[i]]))
  print(emerge[i])
}
for(i in 1:7){
  print(length(rownames(names_16S_D)[names_16S_D$apply.format.bin_16S_D2...2..paste..collapse......%in%dismiss[i]]))
  print(dismiss[i])}

for(i in 1:15){print(length(rownames(names_ITS_D)[names_ITS_D$apply.format.bin_ITS_D2...2..paste..collapse......%in%emerge[i]]))}
for(i in 1:15){print(length(rownames(names_ITS_J)[names_ITS_J$apply.format.bin_ITS_J2...2..paste..collapse......%in%emerge[i]]))}
for(i in 1:15){print(length(rownames(names_ITS_I)[names_ITS_I$apply.format.bin_ITS_I2...2..paste..collapse......%in%emerge[i]]))}
for(i in 1:15){print(length(rownames(names_ITS_C)[names_ITS_C$apply.format.bin_ITS_C2...2..paste..collapse......%in%emerge[i]]))}
for(i in 1:7){print(length(rownames(names_ITS_D)[names_ITS_D$apply.format.bin_ITS_D2...2..paste..collapse......%in%dismiss[i]]))}
for(i in 1:7){print(length(rownames(names_ITS_J)[names_ITS_J$apply.format.bin_ITS_J2...2..paste..collapse......%in%dismiss[i]]))}
for(i in 1:7){print(length(rownames(names_ITS_I)[names_ITS_I$apply.format.bin_ITS_I2...2..paste..collapse......%in%dismiss[i]]))}
for(i in 1:7){print(length(rownames(names_ITS_C)[names_ITS_C$apply.format.bin_ITS_C2...2..paste..collapse......%in%dismiss[i]]))}

print(length(rownames(names_16S_D)[names_16S_D$apply.format.bin_16S_D2...2..paste..collapse......=='CCCCCCCC']))
print(length(rownames(names_16S_J)[names_16S_J$apply.format.bin_16S_J2...2..paste..collapse......=='CCCCCCCC']))
print(length(rownames(names_16S_I)[names_16S_I$apply.format.bin_16S_I2...2..paste..collapse......=='CCCCCCCC']))
print(length(rownames(names_16S_C)[names_16S_C$apply.format.bin_16S_C2...2..paste..collapse......=='CCCCCCCC']))

print(length(rownames(names_ITS_D)[names_ITS_D$apply.format.bin_ITS_D2...2..paste..collapse......=='CCCCCCCC']))
print(length(rownames(names_ITS_J)[names_ITS_J$apply.format.bin_ITS_J2...2..paste..collapse......=='CCCCCCCC']))
print(length(rownames(names_ITS_I)[names_ITS_I$apply.format.bin_ITS_I2...2..paste..collapse......=='CCCCCCCC']))
print(length(rownames(names_ITS_C)[names_ITS_C$apply.format.bin_ITS_C2...2..paste..collapse......=='CCCCCCCC']))

num_bac_d <- c(5,0,0,0,0,0,0,0,0,0,0,0,5,0,0,9,0,0,0,0,0,0,0)
num_fun_d <- c(5,0,0,0,0,0,0,0,0,0,0,0,1,0,0,7,0,0,0,0,1,2,1)
num_bac_j <- c(2,0,0,0,0,0,0,0,0,0,0,0,13,0,0,4,0,0,0,0,0,0,1)
num_fun_j <- c(11,0,0,0,0,0,0,0,0,0,0,0,9,0,0,18,1,1,3,0,0,5,1)
num_bac_i <- c(1,0,0,0,0,0,0,0,0,0,0,0,2,0,0,3,0,0,0,0,0,0,1)
num_fun_i <- c(8,0,0,0,0,0,0,0,0,0,0,0,1,1,1,7,0,0,1,1,1,2,1)
num_bac_c <- c(2,1,1,1,1,1,2,4,6,8,14,1,21,0,0,3,0,0,0,0,0,0,0)
num_fun_c <- c(6,0,0,0,0,0,0,0,0,0,0,0,3,0,0,8,0,1,1,1,1,2,1)

numbers_pattern <- data.frame(ID = rep(c('CCCCCCCC',emerge,dismiss),8),
                              Number = c(num_bac_d,num_bac_j,num_bac_i,num_bac_c,num_fun_d,num_fun_j,num_fun_i,num_fun_c),
                              Segment = c(rep(c('Duodenum','Jejunum','Ileum','Cecum'),each = 23),
                                          rep(c('Duodenum','Jejunum','Ileum','Cecum'),each = 23)),
                              Kingdom = rep(c('Bacteria','Fungi'),each = 92))
numbers_pattern

numbers_pattern$ID <- factor(numbers_pattern$ID, levels = c('CCCCCCCC',emerge,dismiss))

corrs <- data.frame(ID = rep(c('CCCCCCCC',emerge,dismiss),8),
                    Pattern = c('Present',rep('Colonize',15),rep('Passenger',7)))
numbers_pattern2 <- merge(numbers_pattern,corrs, by.x = 'ID',by.y = 'ID',all.x=T)
numbers_pattern2$Segment <- factor(numbers_pattern2$Segment,levels = c('Duodenum','Jejunum','Ileum','Cecum'))

f41 <- ggplot(numbers_pattern2,aes(y = ID,x = Segment, fill = Number)) + geom_tile() +
  scale_fill_gradientn(colors = rev(colorRampPalette(c("#E64B35", "white"))(10))[c(1,4,5,7:10)],breaks = c(0,1,5,10,15,21), limits=c(0,21)) +
  facet_grid(Pattern~Kingdom,scales="free_y",space='free_y') +
  theme_bw() + ylab('') + xlab('') +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        axis.text.x = element_text(angle = 90,vjust = 1,hjust = 1),axis.text.y = element_blank(),
        axis.ticks.length.y = unit(0,'cm'),legend.position = 'bottom')
f41
arrays <- data.frame(Cat = c(3,3,3,3,3,3,3,3,3,2,2,2,2,2,1,1,3,2,2,2,2,2,2,1,3,1,1,1,1,1,1,1,
                             3,3,3,2,2,2,2,2,3,3,2,2,2,2,2,2,3,2,2,2,2,2,2,2,3,3,3,3,3,3,3,2,
                             1,2,2,2,2,3,3,3,1,2,2,2,2,2,3,3,1,1,2,2,2,3,3,3,1,1,2,2,2,2,3,3,
                             1,1,1,1,2,2,2,3,1,2,2,2,3,3,3,3,1,2,2,3,3,3,3,3,1,1,2,2,3,3,3,3,
                             1,1,2,3,3,3,3,3,1,2,3,3,3,3,3,3,1,1,3,3,3,3,3,3,1,3,3,3,3,3,3,3,
                             2,2,2,2,3,3,3,3,2,2,2,2,2,2,2,3,2,3,3,3,3,3,3,3),
                     Pattern = rep(c('CCCCCCCC','CBBBBBAA','CBBBBBBA','CAAAAAAA','CCCBBBBB','CCBBBBBB','CBBBBBBB',
                                 'CCCCCCCB','ABBBBCCC','ABBBBBCC','AABBBCCC','AABBBBCC','AAAABBBC','ABBBCCCC','ABBCCCCC',
                                 'AABBCCCC','AABCCCCC','ABCCCCCC','AACCCCCC','ACCCCCCC','BBBBCCCC','BBBBBBBC',
                                 'BCCCCCCC'),each = 8),
                     Position = rep(c('P1','P2','P3','P4','P5','P6','P7','P8'),23))
arrays1 <- merge(arrays, corrs,by.x = 'Pattern', by.y = 'ID',all.x = T)
head(arrays1)
arrays$Pattern <- factor(arrays$Pattern,levels = c('CCCCCCCC','CBBBBBAA','CBBBBBBA','CAAAAAAA','CCCBBBBB','CCBBBBBB','CBBBBBBB',
                                                   'CCCCCCCB','ABBBBCCC','ABBBBBCC','AABBBCCC','AABBBBCC','AAAABBBC','ABBBCCCC','ABBCCCCC',
                                                   'AABBCCCC','AABCCCCC','ABCCCCCC','AACCCCCC','ACCCCCCC','BBBBCCCC','BBBBBBBC',
                                                   'BCCCCCCC'))
head(array)
f42 <- ggplot(arrays1,aes(x = Position,y = Pattern,color = as.factor(Cat))) + geom_point(size = 5) +
  facet_grid(`Pattern.y`~.,scales="free_y",space='free_y') + 
  theme_bw() +
  theme(legend.position = 'bottom',axis.text.y = element_blank(),axis.ticks.length.y = unit(0,'cm'),
        panel.border = element_blank(),panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA)) + 
  ylab('') + xlab('DPH') +
  scale_x_discrete(label = c(1,4,7,14,21,28,35,42)) +
  guides(fill = guide_legend(title = 'Occurrence')) +
  scale_color_manual(values = c('#ededed','#3C54881f','#3C5488ff')) 
layout1 <- 'AABBB'
svg('~/Desktop/figure materials/2nd/figure 4A.svg',width = 6,height = 6.8)
f42 + f41 + plot_layout(design = layout1)
dev.off()
# output the number of each patterns
aggregate(numbers_pattern2$Number[numbers_pattern2$Kingdom=='Bacteria'],list(numbers_pattern2$ID[numbers_pattern2$Kingdom=='Bacteria']),sum)
aggregate(numbers_pattern2$Number[numbers_pattern2$Kingdom=='Fungi'],list(numbers_pattern2$ID[numbers_pattern2$Kingdom=='Fungi']),sum)

library(UpSetR)

emerge <- c('ABBBBCCC','ABBBBBCC','AABBBCCC','AABBBBCC', 'AAAABBBC', 'ABBBCCCC', 'ABBCCCCC', 'AABBCCCC','AABCCCCC',
            'ABCCCCCC','AACCCCCC', 'ACCCCCCC', 'BBBBCCCC', 'BBBBBBBC','BCCCCCCC')
dismiss <- c('CBBBBBAA','CBBBBBBA','CAAAAAAA','CCCBBBBB','CCBBBBBB','CBBBBBBB','CCCCCCCB')

Cecum_I <- c(rownames(names_16S_C)[names_16S_C$apply.format.bin_16S_C2...2..paste..collapse......%in%emerge],
             rownames(names_ITS_C)[names_ITS_C$apply.format.bin_ITS_C2...2..paste..collapse......%in%emerge])
Duodenum_I <- c(rownames(names_16S_D)[names_16S_D$apply.format.bin_16S_D2...2..paste..collapse......%in%emerge],
                rownames(names_ITS_D)[names_ITS_D$apply.format.bin_ITS_D2...2..paste..collapse......%in%emerge])
Ileum_I <- c(rownames(names_16S_I)[names_16S_I$apply.format.bin_16S_I2...2..paste..collapse......%in%emerge],
               rownames(names_ITS_I)[names_ITS_I$apply.format.bin_ITS_I2...2..paste..collapse......%in%emerge])
Jejunum_I <- c(rownames(names_16S_J)[names_16S_J$apply.format.bin_16S_J2...2..paste..collapse......%in%emerge],
               rownames(names_ITS_J)[names_ITS_J$apply.format.bin_ITS_J2...2..paste..collapse......%in%emerge])

Cecum_D <- c(rownames(names_16S_C)[names_16S_C$apply.format.bin_16S_C2...2..paste..collapse......%in%dismiss],
             rownames(names_ITS_C)[names_ITS_C$apply.format.bin_ITS_C2...2..paste..collapse......%in%dismiss])
Duodenum_D <- c(rownames(names_16S_D)[names_16S_D$apply.format.bin_16S_D2...2..paste..collapse......%in%dismiss],
                rownames(names_ITS_D)[names_ITS_D$apply.format.bin_ITS_D2...2..paste..collapse......%in%dismiss])
Ileum_D <- c(rownames(names_16S_I)[names_16S_I$apply.format.bin_16S_I2...2..paste..collapse......%in%dismiss],
             rownames(names_ITS_I)[names_ITS_I$apply.format.bin_ITS_I2...2..paste..collapse......%in%dismiss])
Jejunum_D <- c(rownames(names_16S_J)[names_16S_J$apply.format.bin_16S_J2...2..paste..collapse......%in%dismiss],
               rownames(names_ITS_J)[names_ITS_J$apply.format.bin_ITS_J2...2..paste..collapse......%in%dismiss])

Cecum_P <- c(rownames(names_16S_C)[names_16S_C$apply.format.bin_16S_C2...2..paste..collapse......=='CCCCCCCC'],
          rownames(names_ITS_C)[names_ITS_C$apply.format.bin_ITS_C2...2..paste..collapse......=='CCCCCCCC'])
Duodenum_P <- c(rownames(names_16S_D)[names_16S_D$apply.format.bin_16S_D2...2..paste..collapse......=='CCCCCCCC'],
             rownames(names_ITS_D)[names_ITS_D$apply.format.bin_ITS_D2...2..paste..collapse......=='CCCCCCCC'])
Ileum_P <- c(rownames(names_16S_I)[names_16S_I$apply.format.bin_16S_I2...2..paste..collapse......=='CCCCCCCC'],
          rownames(names_ITS_I)[names_ITS_I$apply.format.bin_ITS_I2...2..paste..collapse......=='CCCCCCCC'])
Jejunum_P <- c(rownames(names_16S_J)[names_16S_J$apply.format.bin_16S_J2...2..paste..collapse......=='CCCCCCCC'],
            rownames(names_ITS_J)[names_ITS_J$apply.format.bin_ITS_J2...2..paste..collapse......=='CCCCCCCC'])
list_P <- list(Cecum = Cecum_P,
               Duodenum = Duodenum_P,
               Ileum = Ileum_P,
               Jejunum = Jejunum_P)
list_I <- list(Cecum = Cecum_I,
               Duodenum = Duodenum_I,
               Ileum = Ileum_I,
               Jejunum = Jejunum_I)
list_D <- list(Cecum = Cecum_D,
               Duodenum = Duodenum_D,
               Ileum = Ileum_D,
               Jejunum = Jejunum_D)
upset(fromList(list_P),order.by = 'freq',point.size = 2)
f43 <- upset(fromList(list_P),order.by = 'freq',point.size = 3,
      main.bar.color = c('gray','#8491B4',rep('gray',8)),
      mb.ratio = c(0.7,0.3),
      shade.color = c('gray','gray'),
      matrix.color= c('#F39B7F'),
      sets.bar.color=c("#E64B35","#4DBBD5","#00A087","#3C5488"))

f44 <- upset(fromList(list_I),order.by = 'freq',point.size = 3,
      main.bar.color = c('gray','gray','gray','gray','#8491B4',rep('gray',5)),
      mb.ratio = c(0.7,0.3),
      shade.color = c('gray','gray'),
      matrix.color= c('#F39B7F'),
      sets.bar.color=c("#3C5488","#E64B35","#4DBBD5","#00A087"))

f45 <- upset(fromList(list_D),order.by = 'freq',point.size = 3,
      main.bar.color = c(rep('gray',8)),
      mb.ratio = c(0.7,0.3),
      shade.color = c('gray','gray'),
      matrix.color= c('#F39B7F'),
      sets.bar.color=c("#E64B35","#00A087","#3C5488","#4DBBD5"))


# 3.5 inch * 6.5 inch

library(dplyr)
library(tibble)
raw_16S_L6 <- read.csv('level-6_bacteria.csv')
format_16S_L6 <- raw_16S_L6 %>%
  remove_rownames %>% 
  column_to_rownames(var="index") %>%
  select(-c('Group1','Group2'))%>%
  as.matrix()
format_16S_L6

raw_ITS_L6 <- read.csv('level-6_fungi.csv')
format_ITS_L6 <- raw_ITS_L6 %>%
  remove_rownames %>% 
  column_to_rownames(var="X") %>%
  select(-c('Group1','Group2')) %>%
  as.matrix()
format_ITS_L6
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
ab_temp <- read.xlsx('01data_for_quantitive.xlsx',sheet = 11)
for(i in 1:nrow(ab_temp)){
  format_16S_L6_AA[i,1:529] <- format_16S_L6[i,1:529] * 10^(ab_temp$`log10(16S)`[i])
  format_ITS_L6_AA[i,1:603] <- format_ITS_L6[i,1:603] * 10^(ab_temp$`log10(ITS)`[i])
}
set_all <- unique(c(Cecum_P,Duodenum_P,Jejunum_P,Ileum_P))
taxa1 <- read.xlsx('taxa1.xlsx',sheet = 1)
taxa2 <- read.xlsx('taxa2.xlsx',sheet = 2)
colnames(format_16S_L6_AA)[1:529] <- taxa1$Genus
colnames(format_ITS_L6_AA)[1:603] <- taxa2$Genus
format_L6_AA <- cbind(format_16S_L6_AA[,1:529],format_ITS_L6_AA)
format_L6_AA2 <- format_L6_AA[,colnames(format_L6_AA)%in%set_all]
dim(format_L6_AA2)
format_L6_AA2$Segment <- c(rep(c('Cecum','Duodenum','Ileum','Jejunum'),80))
format_L6_AA3 <- melt(format_L6_AA2)
format_L6_AA3$DPH <- rep(c(rep(c(1,4,7,14,21,28,35,42),each = 40)),18)
format_L6_AA3$Segment <- factor(format_L6_AA3$Segment,levels = c('Duodenum','Jejunum','Ileum','Cecum'))
svg('~/Desktop/figure materials/2nd/figure S9.svg',width = 7,height = 6)
ggplot(format_L6_AA3, aes(x = as.numeric(DPH), y = log10(value + 1),fill = Segment)) + 
  geom_smooth(aes(y = log10(value + 1), colour=Segment)) + 
  facet_wrap(.~variable) +
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) +
  scale_color_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) +
  theme_bw() +
  ylab('Log10(abundance)') + 
  xlab("Time point") + 
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        legend.position = c(0.8,0.1)) +
  scale_x_continuous(breaks = c(1,7,14,21,28,35,42))
dev.off()
rownames(format_L6_AA)
format_L6_AA$PHAGE <- c(rep('P1',40),rep('P2',40*3),rep('P3',40*4))
format_L6_AA_D <- format_L6_AA %>%
  filter(Group1 == 'Duodenum') %>%
  select(-c('Group1'))
format_L6_AA_J <- format_L6_AA %>%
  filter(Group1 == 'Jejunum') %>%
  select(-c('Group1'))
format_L6_AA_I <- format_L6_AA %>%
  filter(Group1 == 'Ileum') %>%
  select(-c('Group1'))
format_L6_AA_C <- format_L6_AA %>%
  filter(Group1 == 'Cecum') %>%
  select(-c('Group1'))

sig_D <- c()
sig_J <- c()
sig_I <- c()
sig_C <- c()
for(i in 1:1132){
  sig_D[i] <- kruskal.test(format_L6_AA_D[,i]~format_L6_AA_D$PHAGE)$p.value
  sig_J[i] <- kruskal.test(format_L6_AA_J[,i]~format_L6_AA_J$PHAGE)$p.value
  sig_I[i] <- kruskal.test(format_L6_AA_I[,i]~format_L6_AA_I$PHAGE)$p.value
  sig_C[i] <- kruskal.test(format_L6_AA_C[,i]~format_L6_AA_C$PHAGE)$p.value
}
Q_D <- p.adjust(sig_D,method = 'BH')
Q_J <- p.adjust(sig_J,method = 'BH')
Q_I <- p.adjust(sig_I,method = 'BH')
Q_C <- p.adjust(sig_C,method = 'BH')
table(Q_D<0.05)
table(Q_J<0.05)
table(Q_I<0.05)
table(Q_C<0.05)
Qvalue <- data.frame(Q_D,Q_J,Q_I,Q_C,colnames(format_L6_AA_D)[1:1132])
Qvalue[is.na(Qvalue)] <- 1
colnames(Qvalue) <- c('Duodenum','Jejunum','Ileum','Cecum','ID')

Q_common <- Qvalue %>%
  filter(Q_D < 0.05 &Q_J < 0.05& Q_I < 0.05 &Q_C < 0.05)
Q_common
View(Q_common)

list_diff <- list(Cecum = Qvalue$ID[Qvalue$Cecum<0.05],
               Duodenum = Qvalue$ID[Qvalue$Duodenum<0.05],
               Ileum = Qvalue$ID[Qvalue$Ileum<0.05],
               Jejunum = Qvalue$ID[Qvalue$Jejunum<0.05])
svg('~/Desktop/figure materials/2nd/figure S10A.svg',width = 10.5,height = 4)
upset(fromList(list_diff),order.by = 'freq',point.size = 3,
      main.bar.color = c('gray','#8491B4',rep('gray',13)),
      mb.ratio = c(0.7,0.3),
      shade.color = c('gray','gray'),
      matrix.color = c('#F39B7F'),
      sets.bar.color = c("#E64B35","#4DBBD5","#00A087","#3C5488"))
dev.off()
Qvalue$Kingdom <- c(rep('Bacteria',529),rep('Fungi',603))

D_b <- Qvalue %>% filter(Kingdom == 'Bacteria') %>% group_by(Duodenum) %>% slice(which.min(Duodenum)) %>%  select(ID) 
D_f <- Qvalue %>% filter(Kingdom == 'Fungi') %>% group_by(Duodenum) %>% slice(which.min(Duodenum)) %>% select(ID) 
J_b <- Qvalue %>% filter(Kingdom == 'Bacteria') %>% group_by(Jejunum) %>% slice(which.min(Jejunum)) %>% select(ID) 
J_f <- Qvalue %>% filter(Kingdom == 'Fungi') %>% group_by(Jejunum) %>% slice(which.min(Jejunum)) %>% select(ID) 
I_b <- Qvalue %>% filter(Kingdom == 'Bacteria') %>% group_by(Ileum) %>% slice(which.min(Ileum)) %>%  select(ID) 
I_f <- Qvalue %>% filter(Kingdom == 'Fungi') %>% group_by(Ileum) %>% slice(which.min(Ileum)) %>% select(ID) 
C_b <- Qvalue %>% filter(Kingdom == 'Bacteria') %>% group_by(Cecum) %>% slice(which.min(Cecum)) %>% select(ID) 
C_f <- Qvalue %>% filter(Kingdom == 'Fungi') %>% group_by(Cecum) %>% slice(which.min(Cecum)) %>% select(ID) 


L6_D <- format_L6_AA_D[,colnames(format_L6_AA_D)%in%c(D_b$ID[1:20],D_f$ID[1:20])]
L6_J <- format_L6_AA_J[,colnames(format_L6_AA_J)%in%c(J_b$ID[1:20],J_f$ID[1:20])]
L6_I <- format_L6_AA_I[,colnames(format_L6_AA_I)%in%c(I_b$ID[1:20],I_f$ID[1:20])]
L6_C <- format_L6_AA_C[,colnames(format_L6_AA_C)%in%c(C_b$ID[1:20],C_f$ID[1:20])]

median_2_long_data2 <- function(data){
  data0 <- scale(data)
  DPH <- rep(c('D01','D04','D07','D14','D21','D28','D35','D42'),each = 10)
  data2 <- aggregate(data0,list(DPH),mean) 
  rownames(data2) <- data2$Group.1
  data2 <- data2[,-1]
  data3 <- t(data2)
  data4 <- melt(data3)
  data4$Kingdom <- rep(rep(c('Bacteria','Fungi'),each = 20),8)
  data4$Phase <- c(rep('Phase1',40),rep('Phase2',120),rep('Phase3',160))
  colnames(data4) <- c("Genus","DPH","Zscore","Kingdom","Phase")
  return(data4)
}
L6_D2 <- median_2_long_data2(L6_D)

levels(L6_D2$Genus)
L6_D3 <- dcast(L6_D2[,1:3],DPH ~ Genus) 
head(L6_D3)
rownames(L6_D3) <- L6_D3$DPH
L6_D3 <- L6_D3[,-1]
head(L6_D3)
pheatmap::pheatmap(t(L6_D3)[21:40,],cluster_cols = F)

L6_D2$Genus <- factor(L6_D2$Genus,
                      levels = c('Luteimonas',"Erysipelatoclostridium","Macrococcus","Sellimonas" ,
                                 "Ruminiclostridium 5","Eisenbergiella", "[Ruminococcus] torques group",
                                 "Prevotellaceae UCG-001","Pseudomonas","f__Muribaculaceae_3","Parabacteroides",
                                 "Cutibacterium","Akkermansia","Lactobacillus","Prevotella 9","Phascolarctobacterium",
                                 "Helicobacter","Staphylococcus","Sphingobium","Acinetobacter",
                                 "Rhodotorula","o__Saccharomycetales","Nigrograna","Lecanicillium","Naganishia",
                                 "Wallemia","p__Basidiomycota", "Nectriopsis","f__Nectriaceae","Diutina","Trichosporon",
                                 "Tetracladium","Gibellulopsis","Stachybotrys","Oliveonia","Podospora","Schizothecium","Chaetomium","Lasiobolidium", "f__Stachybotryaceae_1"))
f51 <- ggplot(L6_D2,aes(x = DPH,y = Genus, fill = Zscore)) + 
  geom_tile() +
  facet_grid(Kingdom~Phase,scales = c('free'),space = 'free_x') +
  scale_fill_gradient2(high = '#E64B35',low = '#3C5488',mid = 'white') +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        text = element_text(size=6)) +
  ggtitle('Duodenum') +
  ylab('')

L6_J2 <- median_2_long_data2(L6_J)

levels(L6_J2$Genus)

L6_J3 <- dcast(L6_J2[,1:3],DPH ~ Genus) 
head(L6_J3)
rownames(L6_J3) <- L6_J3$DPH
L6_J3 <- L6_J3[,-1]
head(L6_J3)
pheatmap::pheatmap(t(L6_J3)[21:40,],cluster_cols = F)


L6_J2$Genus <- factor(L6_J2$Genus,
                      levels = c("Clostridium sensu stricto 1","Streptococcus","Escherichia-Shigella",
                                 "Candidatus Arthromitus","Romboutsia","Ruminococcaceae UCG-014",
                                 "f__Barnesiellaceae","Sellimonas","Alistipes","Phascolarctobacterium",
                                 "CHKCI001","Barnesiella","Anaerostipes","[Eubacterium] hallii group",
                                 "Butyricicoccus","Ruminiclostridium 9","Bacteroides","Helicobacter",
                                 "Lactobacillus","f__Clostridiales vadinBB60 group_5",
                                 "Rhodocollybia","Zygoascus","Clitocybe","Baeospora","Mycena","Resinicium","Sporisorium",
                                 "Plectosphaerella","Gibellulopsis","c__Sordariomycetes_1","Botryotrichum",
                                 "Leohumicola","c__Dothideomycetes","Fusariella","Tetracladium","p__Ascomycota",
                                 "Scutellinia","Heydenia","Paraphoma","Phlebia"))
f52 <- ggplot(L6_J2,aes(x = DPH,y = Genus, fill = Zscore)) + 
  geom_tile() +
  facet_grid(Kingdom~Phase,scales = c('free'),space = 'free_x') +
  scale_fill_gradient2(high = '#E64B35',low = '#3C5488',mid = 'white') +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        text = element_text(size=6)) +
  ggtitle('Jejunum') +
  ylab('')


L6_I2 <- median_2_long_data2(L6_I)
levels(L6_I2$Genus)
L6_I3 <- dcast(L6_I2[,1:3],DPH ~ Genus) 
head(L6_I3)
rownames(L6_I3) <- L6_I3$DPH
L6_I3 <- L6_I3[,-1]
head(L6_I3)
pheatmap::pheatmap(t(L6_I3)[21:40,],cluster_cols = F)

L6_I2$Genus <- factor(L6_I2$Genus,
                      levels = c("Klebsiella","Escherichia-Shigella","Aerococcus","Corynebacterium 1","Candidatus Arthromitus","Macrococcus",
                                 "f__Peptostreptococcaceae","Erysipelatoclostridium","Bacteroides","Alistipes","Anaerostipes","Blautia",
                                 "[Ruminococcus] torques group","Barnesiella","Ruminiclostridium 9","f__Clostridiales vadinBB60 group_5",
                                 "f__Barnesiellaceae","Ruminococcaceae UCG-014","Helicobacter","Lactobacillus",
                                 "Mrakia","Robbauera","Rhodocollybia",  "Sporobolomyces","Baeospora","Sporisorium","Pseudogymnoascus",
                                 "p__Ascomycota","Mortierella", "Gibellulopsis","Tetracladium","Coprinellus","Botryotrichum", "Scutellinia",
                                 "Alternaria", "Candida","Trichosporon","Aspergillus","Schizothecium","Mycena"))
f53 <- ggplot(L6_I2,aes(x = DPH,y = Genus, fill = Zscore)) + 
  geom_tile() +
  facet_grid(Kingdom~Phase,scales = c('free'),space = 'free_x') +
  scale_fill_gradient2(high = '#E64B35',low = '#3C5488',mid = 'white') +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        text = element_text(size=6)) +
  ggtitle('Ileum') +
  ylab('')

L6_C2 <- median_2_long_data2(L6_C)
levels(L6_C2$Genus)
L6_C3 <- dcast(L6_C2[,1:3],DPH ~ Genus) 
head(L6_C3)
rownames(L6_C3) <- L6_I3$DPH
L6_C3 <- L6_C3[,-1]
head(L6_C3)
pheatmap::pheatmap(t(L6_C3)[1:20,],cluster_cols = F)
pheatmap::pheatmap(t(L6_C3)[21:40,],cluster_cols = F)

L6_C2$Genus <- factor(L6_C2$Genus,
                      levels = c("Escherichia-Shigella", "Oscillibacter","Eisenbergiella","Ruminiclostridium 9","Sellimonas","[Ruminococcus] torques group",
                                 "Butyricicoccus","Lactobacillus","Negativibacillus","Phascolarctobacterium","Coprobacter","Parabacteroides",
                                 "Barnesiella","Bacteroides","o__Rhodospirillales","Helicobacter","Faecalibacterium","Rikenella","f__Barnesiellaceae",
                                 "Christensenellaceae R-7 group","Lecanicillium","Pleospora","Rhodocollybia","Sarocladium","Aspergillus",
                                 "Podospora","Pseudogymnoascus","Acremonium","Chaetomium","Fusarium","Gibellulopsis","Leohumicola","f__Nectriaceae",
                                 "Tetracladium","Schizothecium","Botryotrichum","Plectosphaerella","c__Sordariomycetes_1","Cadophora","Talaromyces"))
f54 <- ggplot(L6_C2,aes(x = DPH,y = Genus, fill = Zscore)) + 
  geom_tile() +
  facet_grid(Kingdom~Phase,scales = c('free'),space = 'free_x') +
  scale_fill_gradient2(high = '#E64B35',low = '#3C5488',mid = 'white') +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        text = element_text(size=6)) +
  ggtitle('Cecum') +
  ylab('')
svg('~/Desktop/figure materials/2nd/figure4C.svg',width = 10.5,height = 5)
f51 + f52 + f53 + f54  + plot_layout(guides = 'collect',nrow = 1)
dev.off()

Q_common
L6_59 <- format_L6_AA[,colnames(format_L6_AA)%in%Q_common$ID]
dim(L6_59)
L6_59$Groups <- rep(c('Cecum','Duodenum','Ileum','Jejunum'),80)
del_str <- function(col){
  gsub('_1$','',col)
}
L6_59_C <- L6_59 %>%
  filter(Groups == 'Cecum') %>%
  select(-c("Groups")) %>%
  scale() %>%
  as.data.frame() %>%
  mutate(DPH = rep(c('D01','D04','D07','D14','D21','D28','D35','D42'),each = 10)) %>%
  group_by(DPH) %>%
  summarise(across(everything(), list(mean))) %>%
  melt() %>%
  rename(Genus0 = variable, Zscore = value) %>%
  mutate(Segment = 'Cecum') %>%
  mutate(Genus = del_str(Genus0)) %>%
  select(-c("Genus0"))
head(L6_59_C)

L6_59_D <- L6_59 %>%
  filter(Groups == 'Duodenum') %>%
  select(-c("Groups")) %>%
  scale() %>%
  as.data.frame() %>%
  mutate(DPH = rep(c('D01','D04','D07','D14','D21','D28','D35','D42'),each = 10)) %>%
  group_by(DPH) %>%
  summarise(across(everything(), list(mean))) %>%
  melt() %>%
  rename(Genus0 = variable, Zscore = value) %>%
  mutate(Segment = 'Duodenum') %>%
  mutate(Genus = del_str(Genus0)) %>%
  select(-c('Genus0'))
L6_59_I <- L6_59 %>%
  filter(Groups == 'Ileum') %>%
  select(-c("Groups")) %>%
  scale() %>%
  as.data.frame() %>%
  mutate(DPH = rep(c('D01','D04','D07','D14','D21','D28','D35','D42'),each = 10)) %>%
  group_by(DPH) %>%
  summarise(across(everything(), list(mean))) %>%
  melt() %>%
  rename(Genus0 = variable, Zscore = value) %>%
  mutate(Segment = 'Ileum') %>%
  mutate(Genus = del_str(Genus0)) %>%
  select(-c('Genus0'))
L6_59_J <- L6_59 %>%
  filter(Groups == 'Jejunum') %>%
  select(-c('Groups')) %>%
  scale() %>%
  as.data.frame() %>%
  mutate(DPH = rep(c('D01','D04','D07','D14','D21','D28','D35','D42'),each = 10)) %>%
  group_by(DPH)%>%
  summarise(across(everything(), list(mean))) %>%
  melt() %>%
  rename(Genus0 = variable, Zscore = value) %>%
  mutate(Segment = 'Jejunum') %>%
  mutate(Genus = del_str(Genus0)) %>%
  select(-c('Genus0'))
head(L6_59_J)
head(L6_59_D)

L6_59_plot <- rbind(L6_59_D, L6_59_J, L6_59_I, L6_59_C)
L6_59_plot$Segment <- factor(L6_59_plot$Segment, levels = c('Duodenum','Jejunum','Ileum','Cecum'))
L6_59_plot2 <- data.frame(dph_seg = paste0(L6_59_plot$DPH,L6_59_plot$Segment),
                         genus = L6_59_plot$Genus,
                         zscore = L6_59_plot$Zscore) %>%
  dcast(., genus ~ dph_seg, value.var = "zscore") %>%
  column_to_rownames(var = 'genus')
pheatmap::pheatmap(L6_59_plot2,cluster_cols = F)
L6_59_plot$Genus <- factor(L6_59_plot$Genus,levels = c("Robbauera","Eutypella","Zygoascus","Escherichia-Shigella","Streptococcus",
                                                 "Sellimonas","[Ruminococcus] torques group","Butyricicoccus","Erysipelatoclostridium","f__Ruminococcaceae_1","Ruminiclostridium 5",
                                                 "f__Peptostreptococcaceae","Lactobacillus","Romboutsia","Blautia","CHKCI001","Candidatus Arthromitus","Rhodotorula","Xerochrysium",
                                                 "Dactylonectria","Heydenia","f__Nectriaceae","p__Ascomycota","Aphanoascus","Fusicolla","Metarhizium","Acremonium","Humicola","c__Sordariomycetes_1",
                                                 "Trichosporon","Paraphoma","Botryotrichum","Pseudogymnoascus","Tetracladium","Mortierella","Scutellinia","Gibellulopsis","Plectosphaerella","Fusariella","Gibberella","Bifidobacterium","Trichoderma", 
                                                 "Sarocladium","Talaromyces","Diutina","Nectriopsis","Stachybotrys","Leohumicola","Podospora",
                                                 "Ruminococcaceae UCG-014","Fournierella","Ruminiclostridium","Faecalibacterium","Alistipes","Christensenellaceae R-7 group",  "Helicobacter","Barnesiella","Bacteroides","Phascolarctobacterium"))
ggplot(L6_59_plot,aes(x = DPH, y = Genus, fill = Zscore)) + 
  geom_tile() +
  facet_grid(.~Segment, scales = c('free'), space = 'free_x') + 
  scale_fill_gradient2(high = '#E64B35', low = '#3C5488', mid = 'white') + 
  theme_bw() + 
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),axis.text.x = element_text(angle = 90))+
  scale_x_discrete(labels = c(1,4,7,14,21,28,35,42))