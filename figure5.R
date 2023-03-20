# bacterial co-occurrence network
# four individual networks
# 
# partial network
library(dplyr)
library(tibble)
library(ggplot2)
library(openxlsx)
# library(ppcor)
raw_16S_L6 <- read.csv('level-6_no_archaea.csv')
format_16S_L6 <- raw_16S_L6 %>%
  remove_rownames %>% 
  column_to_rownames(var="index") %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  mutate(across(Granulicella:uncultured.bacterium.40, ~ . / sum)) %>%
  select(-c('sum')) %>%
  mutate(Group1 = factor(Group1, levels = c('Duodenum','Jejunum','Ileum','Cecum'))) %>%
  select(-c('Group2'))

format_16S_L6
dim(format_16S_L6)
raw_ITS_L6 <- read.csv('level-6_for_bar.csv')
format_ITS_L6 <- raw_ITS_L6 %>%
  remove_rownames %>% 
  column_to_rownames(var="X") %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  mutate(across(X__:unidentified.53, ~ . / sum)) %>%
  select(-c('sum')) %>%
  mutate(Group1 = factor(Group1, levels = c('Duodenum','Jejunum','Ileum','Cecum'))) %>%
  select(-c('Group2'))
format_ITS_L6
ids_bacteria <- read.xlsx('taxa1.xlsx',sheet = 1)
ids_fungi <- read.xlsx('taxa2.xlsx', sheet = 2)
colnames(format_16S_L6)[1:529] <- ids_bacteria$Genus
colnames(format_ITS_L6)[1:603] <- ids_fungi$Genus
ab_temp <- read.xlsx('01data_for_quantitive.xlsx',sheet = 11)
ab_bac_all <- format_16S_L6
ab_fun_all <- format_ITS_L6
for(i in 1:nrow(ab_temp)){
  ab_bac_all[i,1:529] <- format_16S_L6[i,1:529] * 10^(ab_temp$`log10(16S)`[i])
  ab_fun_all[i,1:603] <- format_ITS_L6[i,1:603] * 10^(ab_temp$`log10(ITS)`[i])
}

ab_L6_all <- cbind(ab_bac_all[1:529], ab_fun_all)

ab_L6_C <- ab_L6_all %>% filter(Group1 == 'Cecum') %>% select(-c('Group1'))
ab_L6_D <- ab_L6_all %>% filter(Group1 == 'Duodenum') %>% select(-c('Group1'))
ab_L6_I <- ab_L6_all %>% filter(Group1 == 'Ileum') %>% select(-c('Group1'))
ab_L6_J <- ab_L6_all %>% filter(Group1 == 'Jejunum') %>% select(-c('Group1'))

filt_30_percent <- function(data, percent = 0.3){
  num_zero <- c()
  for(i in 1:ncol(data)){
    num_zero[i] <- sum(data[,i]==0)
  }
  cut_off <- nrow(data)*(1-percent)
  data2 <- data[,num_zero<cut_off]
  print(dim(data2))
  return(data2)
}

ab_L6_C2 <- filt_30_percent(ab_L6_C) # 274
ab_L6_D2 <- filt_30_percent(ab_L6_D) # 182
ab_L6_I2 <- filt_30_percent(ab_L6_I) # 124
ab_L6_J2 <- filt_30_percent(ab_L6_J) # 143

pair_ppcor <- function(data, partial){
  require(ppcor)
  P1 <- c()
  C1 <- c()
  combinations <- t(combn(colnames(data),2))
  for(i in 1:nrow(combinations)){
    storage <- pcor.test(data[,combinations[i,1]],data[,combinations[i,2]],partial, method = 'spearman')
    P1[i] <- storage$p.value
    C1[i] <- storage$estimate
  }
  raw0 <- data.frame(combinations,P1,C1)
  raw0$Q1 <- p.adjust(raw0$P1,method = 'BH')
  raw1 <- raw0[raw0$Q1<0.05,]
  print(dim(raw1))
  return(raw1)
}


par_Cor_C <- pair_ppcor(ab_L6_C2,c(rep(c(1,4,7,14,21,28,35,42),each = 10))) # 19626
par_Cor_D <- pair_ppcor(ab_L6_D2,c(rep(c(1,4,7,14,21,28,35,42),each = 10))) # 5338
par_Cor_I <- pair_ppcor(ab_L6_I2,c(rep(c(1,4,7,14,21,28,35,42),each = 10))) # 2540
par_Cor_J <- pair_ppcor(ab_L6_J2,c(rep(c(1,4,7,14,21,28,35,42),each = 10))) # 5197

par_cor_all_dupli <- rbind(par_Cor_C[,1:2],par_Cor_D[,1:2],par_Cor_I[,1:2],par_Cor_J[,1:2])
par_cor_all <-  par_cor_all_dupli[!duplicated(t(apply(par_cor_all_dupli[1:2], 1, sort))),]

#Part 1 Stats

library(igraph)
graph_C <- graph.data.frame(d = par_Cor_C[,1:2], directed = FALSE)
graph_D <- graph.data.frame(d = par_Cor_D[,1:2], directed = FALSE)
graph_I <- graph.data.frame(d = par_Cor_I[,1:2], directed = FALSE)
graph_J <- graph.data.frame(d = par_Cor_J[,1:2], directed = FALSE)
graph_all <- graph.data.frame(d = par_cor_all, directed = FALSE)

graph_stat <- function(graph,Segment){
  Vcount <- vcount(graph)
  Ecount <- ecount(graph) # edge number
  Cluster <- transitivity(graph) # clustering coefficient
  Separation <- mean_distance(graph) # average separation
  Betweenness <- mean(centr_betw(graph)$res) # average betweenness
  Stats <- data.frame(Segment, Vcount, Ecount, Cluster, Betweenness, Separation)
  return(Stats)
}

Stat_graph_C <- graph_stat(graph_C, 'Cecum')
Stat_graph_D <- graph_stat(graph_D, 'Duodenum')
Stat_graph_I <- graph_stat(graph_I, 'Ileum')
Stat_graph_J <- graph_stat(graph_J, 'Jejunum')
require(reshape2)
statis_networks <- rbind(Stat_graph_C,Stat_graph_D,Stat_graph_I,Stat_graph_J)
statis_networks$Segment <- factor(statis_networks$Segment, levels = c('Duodenum','Jejunum','Ileum','Cecum'))

# save.image('~/Desktop/temp.RData')
plot_g1 <- ggplot(statis_networks,aes(x = Segment, y = Vcount, fill = Segment)) + 
  geom_bar(stat = 'identity') + 
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) +
  theme_bw() + 
  theme(legend.position = 'none',panel.grid.minor = element_line(color = NA),axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 6),
        axis.text.y = element_text(size = 6), axis.title = element_text(size = 8)) +
  ylab('Number of vertices')
plot_g1
plot_g2 <- ggplot(statis_networks, aes(x = Segment, y = Ecount, fill = Segment)) +
  geom_bar(stat = 'identity') + 
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) +
  theme_bw() + 
  theme(legend.position = 'none',panel.grid.minor = element_line(color = NA),axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 6),
        axis.text.y = element_text(size = 6), axis.title = element_text(size = 8)) +
  ylab('Number of edges')
plot_g2
plot_g3 <- ggplot(statis_networks, aes(x = Segment, y = Cluster, fill = Segment)) +
  geom_bar(stat = 'identity') + 
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) +
  theme_bw() + 
  theme(legend.position = 'none',panel.grid.minor = element_line(color = NA),axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 6),
        axis.text.y = element_text(size = 6), axis.title = element_text(size = 8)) +
  ylab('Clustering coefficient')
plot_g3

plot_g4 <- ggplot(statis_networks, aes(x = Segment, y = Betweenness, fill = Segment)) + 
  geom_bar(stat = 'identity') + 
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) +
  theme_bw() + 
  theme(legend.position = 'none',panel.grid.minor = element_line(color = NA),axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 6),
        axis.text.y = element_text(size = 6), axis.title = element_text(size = 8)) +
  ylab('Average betweenness')
plot_g4

plot_g5 <- ggplot(statis_networks, aes(x = Segment, y = Separation, fill = Segment)) + 
  geom_bar(stat = 'identity') + 
  scale_fill_manual(values = c('#E64B35','#4DBBD5','#00A087','#3C5488')) +
  theme_bw() + 
  theme(legend.position = 'none',panel.grid.minor = element_line(color = NA),axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 6),
        axis.text.y = element_text(size = 6), axis.title = element_text(size = 8)) +
  ylab('Average separation')
plot_g5
require(patchwork)


id_all <- rbind(ids_bacteria,ids_fungi)
add_name_for_igraph <- function(graph, id_all){
  vertex_id <- unlist(vertex_attr(graph)$name)
  corrs <- c()
  for(i in 1:vcount(graph)){
    corrs[i] <- id_all$Phylum[vertex_id[i]==id_all$Genus]
  }
  V(graph)$Phylum <- corrs
  return(graph)
}
graph_C2 <- add_name_for_igraph(graph_C, id_all)
graph_D2 <- add_name_for_igraph(graph_D, id_all)
graph_I2 <- add_name_for_igraph(graph_I, id_all)
graph_J2 <- add_name_for_igraph(graph_J, id_all)
overrep <- function(p1 = "Firmicutes", p2 = "Firmicutes",g1 = g.link){
  vertex.phylum <- V(g1)$Phylum
  V(g1)$name <- vertex.phylum
  
  if (p1 != p2){
    v1 = table(vertex.phylum)[p1]
    v2 = table(vertex.phylum)[p2]
    m <- v1*v2*(vcount(g1)-1)/vcount(g1)
    n <- vcount(g1)*(vcount(g1)-1) - m
    k <- ecount(g1)
    x <- length(E(g1)[grep(p1, vertex.phylum) 
                      %--% grep(p2, vertex.phylum)])
    p = phyper(x, m, n, k)}
  else {
    v1 = table(vertex.phylum)[p1]
    m <- v1*c(v1-1)*(vcount(g1)-1)/vcount(g1)
    n <- vcount(g1)*(vcount(g1)-1) - m
    k <- ecount(g1)
    x <- length(E(g1)[grep(p1, vertex.phylum) 
                      %--% grep(p2, vertex.phylum)])
    p = phyper(x, m, n, k)
  }
}
# Cecum 
p_over_C <- matrix(nrow = length(table(V(graph_C2)$Phylum)),ncol = length(table(V(graph_D2)$Phylum)))
for(i in 1:length(table(V(graph_C2)$Phylum))){
  for(j in 1:length(table(V(graph_C2)$Phylum))){
    temp_C <- V(graph_C2)$Phylum
    po_C <- overrep(p1 = names(sort(table(temp_C),decreasing = TRUE))[i],
                    p2 = names(sort(table(temp_C),decreasing = TRUE))[j],
                    graph_D2)
    p_over_C[i,j] <- po_C
  }
  print(i)
}
Q_over_C <- matrix(p.adjust(p_over_C,method = "fdr"),nrow = length(table(V(graph_C2)$Phylum)))
vertex.phylum_C <- V(graph_C2)$Phylum
row.names(Q_over_C) <- names(sort(table(vertex.phylum_C),decreasing = TRUE))
colnames(Q_over_C) <- names(sort(table(vertex.phylum_C),decreasing = TRUE))
p_over_C <- matrix(nrow = length(table(V(graph_C2)$Phylum)),ncol = length(table(V(graph_C2)$Phylum)))
for(i in 1:length(table(V(graph_C2)$Phylum))){
  for(j in 1:length(table(V(graph_C2)$Phylum))){
    temp_C <- V(graph_C2)$Phylum
    po_C <- overrep(p1 = names(sort(table(temp_C),decreasing = TRUE))[i], 
                    p2 = names(sort(table(temp_C),decreasing = TRUE))[j],
                    graph_C2)
    p_over_C[i,j] <- po_C
  }
  print(i)
}
Q_over_C <- matrix(p.adjust(p_over_C,method = "fdr"),nrow = length(table(V(graph_C2)$Phylum)))
vertex.phylum_C <- V(graph_C2)$Phylum
row.names(Q_over_C) <- names(sort(table(vertex.phylum_C),decreasing = TRUE))
colnames(Q_over_C) <- names(sort(table(vertex.phylum_C),decreasing = TRUE))

# Duodenum
p_over_D <- matrix(nrow = length(table(V(graph_D2)$Phylum)),ncol = length(table(V(graph_D2)$Phylum)))
for(i in 1:length(table(V(graph_D2)$Phylum))){
  for(j in 1:length(table(V(graph_D2)$Phylum))){
    temp_D <- V(graph_D2)$Phylum
    po_D <- overrep(p1 = names(sort(table(temp_D),decreasing = TRUE))[i], 
                    p2 = names(sort(table(temp_D),decreasing = TRUE))[j],
                    graph_D2)
    p_over_D[i,j] <- po_D
  }
  print(i)
}
Q_over_D <- matrix(p.adjust(p_over_D,method = "fdr"),nrow = length(table(V(graph_D2)$Phylum)))
vertex.phylum_D <- V(graph_D2)$Phylum
row.names(Q_over_D) <- names(sort(table(vertex.phylum_D),decreasing = TRUE))
colnames(Q_over_D) <- names(sort(table(vertex.phylum_D),decreasing = TRUE))
p_over_D <- matrix(nrow = length(table(V(graph_D2)$Phylum)),ncol = length(table(V(graph_D2)$Phylum)))
for(i in 1:length(table(V(graph_D2)$Phylum))){
  for(j in 1:length(table(V(graph_D2)$Phylum))){
    temp_D <- V(graph_D2)$Phylum
    po_D <- overrep(p1 = names(sort(table(temp_D),decreasing = TRUE))[i], 
                    p2 = names(sort(table(temp_D),decreasing = TRUE))[j],
                    graph_D2)
    p_over_D[i,j] <- po_D
  }
  print(i)
}
Q_over_D <- matrix(p.adjust(p_over_D,method = "fdr"),nrow = length(table(V(graph_D2)$Phylum)))
vertex.phylum_D <- V(graph_D2)$Phylum
row.names(Q_over_D) <- names(sort(table(vertex.phylum_D),decreasing = TRUE))
colnames(Q_over_D) <- names(sort(table(vertex.phylum_D),decreasing = TRUE))



# Ileum
p_over_I <- matrix(nrow = length(table(V(graph_I2)$Phylum)),ncol = length(table(V(graph_I2)$Phylum)))
for(i in 1:length(table(V(graph_I2)$Phylum))){
  for(j in 1:length(table(V(graph_I2)$Phylum))){
    temp_I <- V(graph_I2)$Phylum
    po_I <- overrep(p1 = names(sort(table(temp_I),decreasing = TRUE))[i], 
                    p2 = names(sort(table(temp_I),decreasing = TRUE))[j],
                    graph_I2)
    p_over_I[i,j] <- po_I
  }
  print(i)
}
Q_over_I <- matrix(p.adjust(p_over_I,method = "fdr"),nrow = length(table(V(graph_I2)$Phylum)))
vertex.phylum_I <- V(graph_I2)$Phylum
row.names(Q_over_I) <- names(sort(table(vertex.phylum_I),decreasing = TRUE))
colnames(Q_over_I) <- names(sort(table(vertex.phylum_I),decreasing = TRUE))
p_over_I <- matrix(nrow = length(table(V(graph_I2)$Phylum)),ncol = length(table(V(graph_I2)$Phylum)))
for(i in 1:length(table(V(graph_I2)$Phylum))){
  for(j in 1:length(table(V(graph_I2)$Phylum))){
    temp_I <- V(graph_I2)$Phylum
    po_I <- overrep(p1 = names(sort(table(temp_I),decreasing = TRUE))[i], 
                    p2 = names(sort(table(temp_I),decreasing = TRUE))[j],
                    graph_I2)
    p_over_I[i,j] <- po_I
  }
  print(i)
}
Q_over_I <- matrix(p.adjust(p_over_I,method = "fdr"),nrow = length(table(V(graph_I2)$Phylum)))
vertex.phylum_I <- V(graph_I2)$Phylum
row.names(Q_over_I) <- names(sort(table(vertex.phylum_I),decreasing = TRUE))
colnames(Q_over_I) <- names(sort(table(vertex.phylum_I),decreasing = TRUE))

# Jejunum
p_over_J <- matrix(nrow = length(table(V(graph_J2)$Phylum)),ncol = length(table(V(graph_J2)$Phylum)))
for(i in 1:length(table(V(graph_J2)$Phylum))){
  for(j in 1:length(table(V(graph_J2)$Phylum))){
    temp_J <- V(graph_J2)$Phylum
    po_J <- overrep(p1 = names(sort(table(temp_J),decreasing = TRUE))[i], 
                    p2 = names(sort(table(temp_J),decreasing = TRUE))[j],
                    graph_J2)
    p_over_J[i,j] <- po_J
  }
  print(i)
}
Q_over_J <- matrix(p.adjust(p_over_J,method = "fdr"),nrow = length(table(V(graph_J2)$Phylum)))
vertex.phylum_J <- V(graph_J2)$Phylum
row.names(Q_over_J) <- names(sort(table(vertex.phylum_J),decreasing = TRUE))
colnames(Q_over_J) <- names(sort(table(vertex.phylum_J),decreasing = TRUE))
p_over_J <- matrix(nrow = length(table(V(graph_J2)$Phylum)),ncol = length(table(V(graph_J2)$Phylum)))
for(i in 1:length(table(V(graph_J2)$Phylum))){
  for(j in 1:length(table(V(graph_J2)$Phylum))){
    temp_J <- V(graph_J2)$Phylum
    po_J <- overrep(p1 = names(sort(table(temp_J),decreasing = TRUE))[i], 
                    p2 = names(sort(table(temp_J),decreasing = TRUE))[j],
                    graph_J2)
    p_over_J[i,j] <- po_J
  }
  print(i)
}
Q_over_J <- matrix(p.adjust(p_over_J,method = "fdr"),nrow = length(table(V(graph_J2)$Phylum)))
vertex.phylum_J <- V(graph_J2)$Phylum
row.names(Q_over_J) <- names(sort(table(vertex.phylum_J),decreasing = TRUE))
colnames(Q_over_J) <- names(sort(table(vertex.phylum_J),decreasing = TRUE))

require(corrplot)
par(mfrow=c(1,4))
D_phylum <- rownames(Q_over_D)
corrplot(corr = 1-Q_over_D[D_phylum,D_phylum],col = c('#E64B35'),tl.col = "black",tl.cex=.6,cl.pos = "n",method = "square",
         insig = "blank",sig.level = 0.05, p.mat = Q_over_D[D_phylum,D_phylum],type = "lower")
J_phylum <- rownames(Q_over_J)
corrplot(corr = 1-Q_over_J[J_phylum,J_phylum],col = c('#E64B35'),tl.col = "black",tl.cex=.6,cl.pos = "n",method = "square",
         insig = "blank",sig.level = 0.05, p.mat = Q_over_J[J_phylum,J_phylum],type = "lower")
I_phylum <- rownames(Q_over_I)
corrplot(corr = 1-Q_over_I[I_phylum,I_phylum],col = c('#E64B35'),tl.col = "black",tl.cex=.6,cl.pos = "n",method = "square",
         insig = "blank",sig.level = 0.05, p.mat = Q_over_I[I_phylum,I_phylum],type = "lower")
C_phylum <- rownames(Q_over_C)
corrplot(corr = 1-Q_over_C[C_phylum,C_phylum],col = c('#E64B35'),tl.col = "black",tl.cex=.6,cl.pos = "n",method = "square",
         insig = "blank",sig.level = 0.05, p.mat = Q_over_C[C_phylum,C_phylum],type = "lower")
par(mfrow=c(1,1))

# generalist and specialist
E(graph_C) #19626
E(graph_D) # 5338
E(graph_I) # 2540
E(graph_J) # 5197
E(graph_all)
# 342 24830
edge_list <- list(get.edgelist(graph_D),get.edgelist(graph_J),get.edgelist(graph_I),get.edgelist(graph_C))
edge_list2 <- lapply(edge_list, FUN = function(x){paste(x[,1],x[,2],sep = '@')})
edgelist.unique <- unique(unlist(edge_list2))
edgelist_df <- matrix(0, nrow = 24830, ncol = 4)
for(i in 1:4){
  edgelist_df[edgelist.unique %in% edge_list2[[i]],i] <- 1
}
rownames(edgelist_df) <- edgelist.unique

dim(edgelist_df[(edgelist_df[,1]==1)&(rowSums(edgelist_df)==1),]) # 2555
dim(edgelist_df[(edgelist_df[,2]==1)&(rowSums(edgelist_df)==1),]) # 1835
dim(edgelist_df[(edgelist_df[,3]==1)&(rowSums(edgelist_df)==1),]) # 183
dim(edgelist_df[(edgelist_df[,4]==1)&(rowSums(edgelist_df)==1),]) # 15315

co_exist_edges <- dim(edgelist_df[rowSums(edgelist_df)==4,])
unique_edges <- c()
for(i in 1:4){
  unique_edges[i] <- dim(edgelist_df[(edgelist_df[,i]==1)&(rowSums(edgelist_df)==1),])[1]
}
total_edge <- c(ecount(graph_D),ecount(graph_J),ecount(graph_I),ecount(graph_C))
total_edge

shared_edges <- data.frame(total_edge,unique_edges,
                           Segment = c('Duodenum','Jejunum','Ileum','Cecum'))
shared_edges$co_exist <- co_exist_edges[1]
shared_edges
shared_edges$Generalist <- (shared_edges$total_edge - shared_edges$unique_edges)/shared_edges$total_edge * 100
shared_edges$Specialist <- shared_edges$unique_edges/shared_edges$total_edge * 100
shared_edges2 <- melt(shared_edges[,c(3,5,6)])
shared_edges2
shared_edges2$Segment <- factor(shared_edges2$Segment, levels = c('Duodenum','Jejunum','Ileum','Cecum'))
shared_edges2$variable <- factor(shared_edges2$variable, levels = c('Specialist','Generalist'))
figure5c <- ggplot(shared_edges2, aes(x = Segment, y = value, fill = variable)) +
  geom_bar(stat = 'identity') + 
  guides(fill = guide_legend(title = 'Edge type')) + 
  theme_bw() + 
  scale_fill_manual(values = c('#E64B35','#3C5488')) +
  ylab('Percentage (%)') +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        legend.position = 'bottom',axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 6),axis.text.y = element_text(size = 6),
        legend.key.size = unit(0.5, 'cm'))
figure5c

# core taxa
core_taxa <- unique(c(names(sort(degree(graph_D),decreasing = T)[1:10]),
                      names(sort(degree(graph_J),decreasing = T)[1:10]),
                      names(sort(degree(graph_I),decreasing = T)[1:11]),
                      names(sort(degree(graph_C),decreasing = T)[1:10])))
core_taxa_sub <- list(names(sort(degree(graph_D),decreasing = T)[1:10]),
                      names(sort(degree(graph_J),decreasing = T)[1:10]),
                      names(sort(degree(graph_I),decreasing = T)[1:11]),
                      names(sort(degree(graph_C),decreasing = T)[1:10]))
core_taxa_degree <- matrix(,nrow = length(core_taxa),ncol = 4)
for(i in 1:4){
  core_taxa_degree[,i] <- core_taxa%in%core_taxa_sub[[i]]
}
core_taxa_degree[core_taxa_degree==TRUE] <- 1
core_taxa_degree[core_taxa_degree==FALSE] <- 0
rownames(core_taxa_degree) <- core_taxa
colnames(core_taxa_degree) <- c('Duodenum','Jejunum','Ileum','Cecum')
core_taxa_degree2 <- melt(core_taxa_degree)
core_taxa_degree2$Var2 <- factor(core_taxa_degree2$Var2,levels = c('Cecum','Ileum','Jejunum','Duodenum'))
core_taxa_degree2$Var1 <- factor(core_taxa_degree2$Var1,levels = aggregate(core_taxa_degree2$value,list(core_taxa_degree2$Var1),sum)$Group.1[order(aggregate(core_taxa_degree2$value,list(core_taxa_degree2$Var1),sum)$x,decreasing = T)])
ab_L6_C_degree <- data.frame(degree(graph_C))
ab_L6_D_degree <- data.frame(degree(graph_D))
ab_L6_I_degree <- data.frame(degree(graph_I))
ab_L6_J_degree <- data.frame(degree(graph_J))
head(core_taxa_degree2)
head(core_taxa_degree2)
core_taxa_degree2$M1 <- paste(core_taxa_degree2$Var2,core_taxa_degree2$Var1,sep = '@')
colnames(ab_L6_C_degree) <- 'Degree'
colnames(ab_L6_D_degree) <- 'Degree'
colnames(ab_L6_I_degree) <- 'Degree'
colnames(ab_L6_J_degree) <- 'Degree'
ab_L6_C_degree$M2 <- paste('Cecum',rownames(ab_L6_C_degree),sep = '@')
ab_L6_D_degree$M2 <- paste('Duodenum',rownames(ab_L6_D_degree),sep = '@')
ab_L6_I_degree$M2 <- paste('Ileum',rownames(ab_L6_I_degree),sep = '@')
ab_L6_J_degree$M2 <- paste('Jejunum',rownames(ab_L6_J_degree),sep = '@')
head(ab_L6_C_degree)
ab_L6_degree <- rbind(ab_L6_C_degree,ab_L6_D_degree,ab_L6_I_degree,ab_L6_J_degree)
core_taxa_degree3 <- merge(core_taxa_degree2,ab_L6_degree,by.x = 'M1',by.y = 'M2',all.x = T)
head(core_taxa_degree3)
#svg('~/Desktop/figure materials/2nd/figure5E.svg',width = 7,height = 2.5)
core_taxa_degree4 <- core_taxa_degree3[core_taxa_degree3$value==1,]
figure5d <- ggplot(core_taxa_degree4,aes(x = Var1,y = Var2,color = Degree)) +
  geom_point() +
  theme_bw() +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 6),
        axis.text.y = element_text(size = 6),
        legend.key.size = unit(0.5, 'cm')) +
  ylab('Segment') + xlab('Genus') +
  scale_color_continuous(low = '#E64B351f',high = '#E64B35')
figure5d

# svg('~/Desktop/01maturity/figure materials/3rd_svg/figure5CD.svg',width = 6.69, height = 2.5)
# figure5c + figure5d + plot_layout(widths = c(1,3))
# dev.off()

generate_ratio_pn <- function(par_cor){
  genus_list <- unique(c(par_cor[,1],par_cor[,2]))
  par_temp <- data.frame(genus = c(par_cor[,1],par_cor[,2]),
                         cor = c(par_cor[,4],par_cor[,4]))
  RA <- c()
  temp_P <- c()
  temp_N <- c()
  for(i in 1:length(genus_list)){
    temp_P[i] <- sum(par_temp$cor[par_temp$genus==genus_list[i]]>0)
    temp_N[i] <- sum(par_temp$cor[par_temp$genus==genus_list[i]]<0)
    if(temp_P[i] == 0){temp_P[i] = temp_P[i] + 1}
    else if(temp_N[i] == 0){temp_N[i] = temp_N[i] + 1}
    RA[i] <- temp_N[i]/temp_P[i]
  }
  return(data.frame(genus_list,RA))
}
C_PN <- generate_ratio_pn(par_Cor_C) # N/P
D_PN <- generate_ratio_pn(par_Cor_D)
I_PN <- generate_ratio_pn(par_Cor_I)
J_PN <- generate_ratio_pn(par_Cor_J)
head(C_PN)
degree_C <- data.frame(degree(graph_C))
degree_D <- data.frame(degree(graph_D))
degree_I <- data.frame(degree(graph_I))
degree_J <- data.frame(degree(graph_J))
colnames(degree_C) <- colnames(degree_D) <- colnames(degree_I) <- colnames(degree_J) <- 'Degree'
C_PN_degree <- merge(C_PN, degree_C, by.x = 'genus_list', by.y = 'row.names', all.x = T)
D_PN_degree <- merge(D_PN, degree_D, by.x = 'genus_list', by.y = 'row.names', all.x = T)
I_PN_degree <- merge(I_PN, degree_I, by.x = 'genus_list', by.y = 'row.names', all.x = T)
J_PN_degree <- merge(J_PN, degree_J, by.x = 'genus_list', by.y = 'row.names', all.x = T)
head(D_PN_degree)
D_PN_degree$Category[D_PN_degree$RA>1] <- "Neg"
D_PN_degree$Category[D_PN_degree$RA<1] <- 'Pos'
J_PN_degree$Category[J_PN_degree$RA>1] <- "Neg"
J_PN_degree$Category[J_PN_degree$RA<1] <- 'Pos'
I_PN_degree$Category[I_PN_degree$RA>1] <- "Neg"
I_PN_degree$Category[I_PN_degree$RA<1] <- 'Pos'
C_PN_degree$Category[C_PN_degree$RA>1] <- "Neg"
C_PN_degree$Category[C_PN_degree$RA<1] <- 'Pos'

D_PN_degree <- D_PN_degree[D_PN_degree$genus_list!='Enterobacter',]

figure5e1 <- ggplot(D_PN_degree, aes(x = RA, y = Degree, color = Category)) + 
  geom_point(size = 0.5) +
  scale_x_log10() +
  scale_color_manual(values = c('#E64B35','#3C5488'))+
  theme_bw() +
  xlab('Negative/Positive') +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        legend.position = 'none') +
  geom_vline( xintercept = c(1),color = "gray") +
  ggtitle('Duodenum')

figure5e2 <- ggplot(J_PN_degree, aes(x = RA, y = Degree, color = Category)) + 
  geom_point(size = 0.5) +
  scale_x_log10() +
  scale_color_manual(values = c('#E64B35','#3C5488'))+
  theme_bw() +
  xlab('Negative/Positive') +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        legend.position = 'none') +
  geom_vline( xintercept = c(1),color = "gray") +
  ggtitle('Jejunum')

I_PN_degree <- I_PN_degree[-87,]
figure5e3 <- ggplot(I_PN_degree, aes(x = RA, y = Degree, color = Category)) + 
  geom_point(size = 0.5) + 
  scale_x_log10() +
  scale_color_manual(values = c('#E64B35','#3C5488'))+
  theme_bw() +
  xlab('Negative/Positive') +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        legend.position = 'none') +
  geom_vline( xintercept = c(1),color = "gray") +
  ggtitle('Ileum')

figure5e4 <- ggplot(C_PN_degree, aes(x = RA, y = Degree, color = Category)) + 
  geom_point(size = 0.5) + 
  scale_x_log10() + 
  scale_color_manual(values = c('#E64B35','#3C5488')) +
  theme_bw() +
  xlab('Negative/Positive') +
  theme(panel.grid.major = element_line(color = NA), panel.grid.minor = element_line(color = NA),
        legend.position = 'none') +
  geom_vline( xintercept = c(1),color = "gray") +
  ggtitle('Cecum') 


 svg('~/Desktop/01maturity/figure materials/3rd_svg/figure5e.svg',width = 6.69,height = 4)
 figure5e1 + figure5e2 + figure5e3 + figure5e4 + plot_layout(design = 'AB
                                                             CD')
dev.off()

# output sub-networks for gephi
# format requirement:
# 1) correlation with pos and neg;
# 2) add phylum to the nodes
prepare_edge_for_gephi <- function(data){
  data1 <- data[,c(1,2,4)]
  data1$PN[data1$C1>0] <- 'P'
  data1$PN[data1$C1<0] <- 'N'
  data2 <- data1[,c(1,2,4)]
  return(data2)
}
par_Cor_D2 <- prepare_edge_for_gephi(par_Cor_D)
par_Cor_J2 <- prepare_edge_for_gephi(par_Cor_J)
par_Cor_I2 <- prepare_edge_for_gephi(par_Cor_I)
par_Cor_C2 <- prepare_edge_for_gephi(par_Cor_C)
prepare_node_for_gephi <- function(data,id_all){
  ids_genus <- unique(c(data[,1],data[,2]))
  print(length(ids_genus))
  id2 <- id_all[id_all$Genus%in%ids_genus,c(1,2,3)]
  return(id2)
}
par_Cor_D3 <- prepare_node_for_gephi(par_Cor_D,id_all)
par_Cor_J3 <- prepare_node_for_gephi(par_Cor_J,id_all)
par_Cor_I3 <- prepare_node_for_gephi(par_Cor_I,id_all)
par_Cor_C3 <- prepare_node_for_gephi(par_Cor_C,id_all)
require(openxlsx)
write.xlsx(par_Cor_D2,'~/Desktop/01maturity/material/06network/gephi_4th_0521/edge_D.xlsx')
write.xlsx(par_Cor_J2,'~/Desktop/01maturity/material/06network/gephi_4th_0521/edge_J.xlsx')
write.xlsx(par_Cor_I2,'~/Desktop/01maturity/material/06network/gephi_4th_0521/edge_I.xlsx')
write.xlsx(par_Cor_C2,'~/Desktop/01maturity/material/06network/gephi_4th_0521/edge_C.xlsx')
write.xlsx(par_Cor_D3,'~/Desktop/01maturity/material/06network/gephi_4th_0521/node_D.xlsx')
write.xlsx(par_Cor_J3,'~/Desktop/01maturity/material/06network/gephi_4th_0521/node_J.xlsx')
write.xlsx(par_Cor_I3,'~/Desktop/01maturity/material/06network/gephi_4th_0521/node_I.xlsx')
write.xlsx(par_Cor_C3,'~/Desktop/01maturity/material/06network/gephi_4th_0521/node_C.xlsx')

# First round reversion
# Date: Mar 3rd, 2023
# Author: Yuqing Feng
# Aim: Analysis the associations between weight and the microbiobes
ls(pattern = 'ab_L6_[DJIC]2')
weight <- read.xlsx('weight.xlsx',sheet = 1, rowNames = T)
head(weight)
head(ab_L6_C2)
# sig_glm[i] <- summary(glm(m_all_metlin_log_scale[,i]~DPH))$coefficients[2,4]
glm_L6_weight <- function(L6, weight){
  if(nrow(L6) == nrow(weight)){
    print('The format is OK')
  }
  sig_glm <- c()
  cor_glm <- c()
  for(i in 1:ncol(L6)){
    sig_glm[i] <- summary(glm(L6[,i] ~ weight$Date))$coefficients[2,4]
    cor_glm[i] <- summary(glm(L6[,i] ~ weight$Date))$coefficients[2,1]
  }
  print(i)
  taxa <- colnames(L6)
  print(sum(sig_glm<0.05))
  sig_glm2 <- p.adjust(sig_glm, method = 'BH')
  print(sum(sig_glm2<0.05))
  color <- c()
  color[(sig_glm2 < 0.05)&(cor_glm>0)] <- 'A'
  color[(sig_glm2 < 0.05)&(cor_glm<0)] <- 'B'
  color[sig_glm2 >= 0.05] <- 'C'
  taxa[sig_glm2 >= 0.05] <- ''
  pva <- -log10(sig_glm)
  est <- sign(cor_glm) * log10(abs(cor_glm))
  return(data.frame(pva, est, color, taxa))
}
glm_D <- glm_L6_weight(ab_L6_D2, weight)
glm_J <- glm_L6_weight(ab_L6_J2, weight)
glm_I <- glm_L6_weight(ab_L6_I2, weight)
glm_C <- glm_L6_weight(ab_L6_C2, weight)
# library(scales)
# library(ggsci)

s1 <- ggplot(glm_D, aes(x = est, y = pva, color = color)) + 
  geom_point(size = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = 'none') +
  geom_text(data = glm_D, mapping = aes(x = est, y = pva + 0.2, label = taxa),size = 1) +
  scale_color_manual(values = c('#d3d3d3')) +
  ylab('-log10(pvalue)') +
  xlab('log10(estimate)') +
  ggtitle('Duodenum')+ ylim(c(0,9)) +xlim(c(-10,10))
s1  

s2 <- ggplot(glm_J, aes(x = est, y = pva, color = color)) + 
  geom_point(size = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = 'none') +
  geom_text(data = glm_J, mapping = aes(x = est , y = pva + 0.2, label = taxa),size = 1) +
  scale_color_manual(values = c('#d3d3d3')) +
  ylab('-log10(pvalue)') +
  xlab('log10(estimate)') +
  ggtitle('Jejunum')+ ylim(c(0,9)) +xlim(c(-10,10))
s2  
s3 <- ggplot(glm_I, aes(x = est, y = pva, color = color)) + 
  geom_point(size = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = 'none') +
  geom_text(data = glm_I, mapping = aes(x = est , y = pva + 0.2, label = taxa),size = 1) +
  scale_color_manual(values = c('#E64B35','#d3d3d3')) +
  ylab('-log10(pvalue)') +
  xlab('log10(estimate)') +
  ggtitle('Ileum')+ ylim(c(0,9)) +xlim(c(-10,10))
s3  
s4 <- ggplot(glm_C, aes(x = est, y = pva, color = color)) + 
  geom_point(size = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = 'none') +
  geom_text(data = glm_C, mapping = aes(x = est, y = pva + 0.2, label = taxa),size = 1) +
  scale_color_manual(values = c('#E64B35','#3C5488','#d3d3d3')) +
  ylab('-log10(pvalue)') +
  xlab('log10(estimate)') +
  ggtitle('Cecum') + ylim(c(0,9)) +xlim(c(-10,10))
s4  
# library(patchwork)
s1 + s2 + s3 + s4 + plot_annotation(tag_levels = "a")
