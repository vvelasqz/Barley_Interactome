# make conditional subnetwork and analyse results, clustering, GO, overlap and difference network, wilcoxon for essentiality

calc_corr <- function(interaction_table, expression_table, met="dist", log=F, genotypes= c("wt","mla6","rar3","dm","bln1")){
  library(energy)
  library(reshape2)
  expression2 <- expression_table
  nodes<- unique(c(interaction_table$source, interaction_table$target))
  expression2<- expression2[rownames(expression2) %in% nodes, ]
  if(log==T){
    expression2<- log2(expression2+1)
  }
  expression2<- expression2[,grep(paste(genotypes, collapse = "|"), colnames(expression2))]
  net_x <- as.data.frame(t(subset(expression2, rownames(expression2) %in% interaction_table$source)))
  net_y <- as.data.frame(t(subset(expression2, rownames(expression2) %in% interaction_table$target)))
  
  # calculate the correlation more efficiently
  if(met=="dist"){
    interaction_table$weight<- apply(interaction_table, 1, 
                                     FUN=function(z){ifelse(z[1] %in% rownames(expression2) & z[2] %in% rownames(expression2),
                                                            dcor(x = net_x[,z[1]], y = net_y[,z[2]]), NA)})
  }else{
    interaction_table$weight<- apply(interaction_table, 1, 
                                     FUN=function(z){ifelse(z[1] %in% rownames(expression2) & z[2] %in% rownames(expression2),
                                                            cor(x= net_x[,z[1]], y=net_y[,z[2]], method=met,use="pairwise"), NA)})
  }
  return(interaction_table)
}

interactome_R2_HC_MLA<- read.csv(file='d_weighted_interactome_R2_HC2_MLA.csv', stringsAsFactors = F)
expression_table<- read.csv("hv_R2_genes_norm_counts_de_tax_sp_filtered.csv", row.names = 1)

cor_data_R<- calc_corr(interaction_table=interactome_R2_HC_MLA, expression_table=expression_table, genotypes=c("wt","bln1"))
cor_data_S<- calc_corr(interaction_table=interactome_R2_HC_MLA, expression_table=expression_table, genotypes=c("mla6","rar3","dm"))
cor_interactomes<- merge(interactome_R2_HC_MLA, cor_data_R, by=c("source", "target"), all=T)
cor_interactomes<- merge(cor_interactomes, cor_data_S, by=c("source", "target"), all=T)
colnames(cor_interactomes)[3:5]<- c("weight_all", "weight_R", "weight_S")

perm_corR<- read.csv("permuted_wt_bln1_dist_tritex.csv")
perm_corR<- perm_corR[perm_corR$freq!=0,]
perm_corR$p_val<- sapply(FUN=function(x)sum(perm_corR$freq[1:x])/sum(perm_corR$freq), 1:nrow(perm_corR))
perm_corS<- read.csv("permuted_mla6_rar3_dm_dist_tritex.csv")
perm_corS<- perm_corS[perm_corS$freq!=0,]
perm_corS$p_val<- sapply(FUN=function(x)sum(perm_corS$freq[1:x])/sum(perm_corS$freq), 1:nrow(perm_corS))

#subset HvInt based on a threshold
HvInt_R<- cor_data_R[cor_data_R$weight> perm_corR[perm_corR$p_val==min(perm_corR$p_val[perm_corR$p_val> 0.95]), "cor"] & !is.na(cor_data_R$weight), ]
nodes_R<- unique(c(HvInt_R$source, HvInt_R$target))

HvInt_S<- cor_data_S[cor_data_S$weight> perm_corS[perm_corS$p_val==min(perm_corS$p_val[perm_corS$p_val> 0.95]), "cor"] & !is.na(cor_data_S$weight), ]
nodes_S<- unique(c(HvInt_S$source, HvInt_S$target))

#difference network R-S
HvInt_R<- HvInt_R[!duplicated(t(apply(HvInt_R, 1, sort))), ]
HvInt_S<- HvInt_S[!duplicated(t(apply(HvInt_S, 1, sort))), ]

HvInt_R$link<- paste0(HvInt_R$source, "_", HvInt_R$target)
HvInt_S$link<- paste0(HvInt_S$source, "_", HvInt_S$target)

HvInt_diff<- HvInt_R[!HvInt_R$link %in% HvInt_S$link, ]
nodes_diff<- unique(c(HvInt_diff$source, HvInt_diff$target))

HvInt_diff_S_R<- HvInt_S[!HvInt_S$link %in% HvInt_R$link, ]
nodes_diff_S_R<- unique(c(HvInt_diff_S_R$source, HvInt_diff_S_R$target))
#number of common links R and S
HvInt_R_S_common<- intersect(HvInt_R[,1:2], HvInt_S[,1:2])
nodes_R_S_common<- unique(c(HvInt_R_S_common$source, HvInt_R_S_common$target))
barley_net_diff <- graph_from_data_frame(HvInt_diff, directed = FALSE, vertices = NULL)
barley_net_diff<-igraph::simplify(barley_net_diff)

barley_net_diff_S_R <- graph_from_data_frame(HvInt_diff_S_R, directed = FALSE, vertices = NULL)
barley_net_diff_S_R<-igraph::simplify(barley_net_diff_S_R)

components(barley_net_diff)$csize
components_barley_diff<- decompose(barley_net_diff, mode = "weak")
barley_net_d1 <- components_barley_diff[[1]]
clusters_d1<- cluster_walktrap(barley_net_d1, steps = 100)
clusters_net_df_d<-data.frame(node=clusters_d1$names, cluster=clusters_d1$membership)
table(clusters_net_df_d$cluster)

#wilcoxon tests
wilcox_test_net<- function(vec1, vec2){
  library(effsize)
  #table trts by bins
  wilcox_test <- data.frame(matrix(ncol = 4, nrow = 1)) 
  colnames(wilcox_test)<- c("V", "pvalue", "Size_effect_VD", "Effect")
  c<- wilcox.test(vec1, vec2, paired=F)
  # #vargha implementation
  d<-VD.A(d = c(vec1, vec2), f = c(rep("R-S", length(vec1)),rep("R=S", length(vec2))))
  wilcox_test[1,]<- c(c["statistic"], c["p.value"], d["estimate"], d["magnitude"])
  return(wilcox_test)
}

nodes_stats_R<- read.csv(file = "nodes_stats_R.csv", row.names = 1)

HvInt_R<- read.csv(file = "HvInt_R.csv", row.names = 1)
HvInt_R<- HvInt_R[!duplicated(t(apply(HvInt_R, 1, sort))), ]
HvInt_R$link<- paste0(HvInt_R$source, "_", HvInt_R$target)
HvInt_diff<- read.csv(file = "HvInt_diff.csv", row.names = 1)
HvInt_R$interaction<- ifelse(HvInt_R$link %in% HvInt_diff$link, "R", "R_S")

#for R-S
w_diffslc<- wilcox_test_net(nodes_stats_R[nodes_stats_R$diff==T, "diffslc_09"], nodes_stats_R[nodes_stats_R$diff==F, "diffslc_09"])
w_degree<- wilcox_test_net(nodes_stats_R[nodes_stats_R$diff==T, "degree"], nodes_stats_R[nodes_stats_R$diff==F, "degree"])
w_betwe<- wilcox_test_net(nodes_stats_R[nodes_stats_R$diff==T, "betweenness"], nodes_stats_R[nodes_stats_R$diff==F, "betweenness"])
w_cor<- wilcox_test_net(HvInt_R[HvInt_R$interaction=="R_S", "weight"], HvInt_R[HvInt_R$interaction=="R", "weight"])

w_table<- rbind(w_diffslc, w_degree, w_betwe, w_cor)
w_table$Property<- c("DiffSLC", "Degree", "Betweenness", "Distance Correlation")

