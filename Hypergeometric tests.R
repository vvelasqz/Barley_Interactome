#HG test for enrichment of top prots in HvInt clusters 
HG_test<- function(cluster_matrix, top_table, cluster_table){
  library(effsize)
  #table trts by bins
  hypergeom_test <- data.frame(matrix(ncol = 1+2*nrow(top_table), nrow = nrow(cluster_matrix))) 
  hypergeom_test[,1]<- rownames(cluster_matrix)
  colnames(hypergeom_test)<- c("cluster", unlist(sapply(top_table$cluster_diffslc, function(x) c(paste0(x,"_HG_pvalue"), paste0(x,"_HG_padj")))))
  
  for(m in top_table$cluster_diffslc){
    
    #total cluster diffscl, total cluster diffscl in cluster, sum not cluster diffscl in cluster, genes in cluster
    hypergeom_test[,paste0(m,"_HG_pvalue")]<- phyper(cluster_matrix[,m], sum(cluster_matrix[,m]), 
                                                     sum(cluster_table$Freq - cluster_matrix[,m]), cluster_table$Freq, lower.tail = FALSE)
    
    hypergeom_test[, paste0(m,"_HG_padj")]<- p.adjust(hypergeom_test[,paste0(m,"_HG_pvalue")], method = "BH")
    
  }
  hypergeom_test<- hypergeom_test[rowSums(is.na(hypergeom_test[,-1])) != ncol(hypergeom_test[,-1]), ]
  hypergeom_test<- merge(hypergeom_test, cluster_table, by.x="cluster",by.y="Var1", all.x=T)
  
  return(hypergeom_test)
}
clusters_net_df<- read.csv("clusters_HvInt_walktrap.csv", row.names = 1)
clusters_net_df_table<- data.frame(table(clusters_net_df$cluster))
node_stats<- read.csv("node_stats.csv", row.names = 1)
node_stats$rank<- rank(-node_stats$diffslc_09)
node_stats$cluster_diffslc<- ifelse(node_stats$rank %in% 1:100, "Top 100", ifelse(node_stats$rank %in% 101:200, "Top 101-200", ifelse(node_stats$rank %in% 201:500, "Top 201-500",  ifelse(node_stats$rank %in% 501:1000, "Top 501-1000", "More_1000"))))
node_stats$cluster_diffslc<-factor(node_stats$cluster_diffslc, levels = c("Top 100", "Top 101-200", "Top 201-500", "Top 501-1000", "More_1000"))
node_stats<- merge(node_stats, clusters_net_df, by="node", all.x=T)
node_stats_1000<- node_stats[!node_stats$cluster_diffslc %in% "More_1000",]

top_table<- data.frame(table(node_stats_1000$cluster_diffslc))[-5,]
colnames(top_table)<- c("cluster_diffslc", "total")
cluster_table<- data.frame(table(clusters_net_df$cluster))
library(reshape2)
nodes_1000_matrix<- acast(node_stats_1000[, 10:11], cluster~cluster_diffslc)
cluster_table<- cluster_table[cluster_table$Var1 %in% rownames(nodes_1000_matrix),]
HG_cluster<- HG_test(cluster_matrix=nodes_1000_matrix, top_table, cluster_table)


#eqtl enrichment
HG_test_eqtl_cluster<- function(cluster_table, eqtl_table){
  library(effsize)
  #table trts by bins
  trts<- unique(eqtl_table$timepoint)
  hypergeom_test <- data.frame(matrix(ncol = 1+2*length(trts), nrow = length(unique(cluster_table$cluster)))) 
  hypergeom_test[,1]<- unique(cluster_table$cluster)
  colnames(hypergeom_test)<- c("cluster", unlist(sapply(trts, function(x) c(paste0(x,"_HG_pvalue"), paste0(x,"_HG_padj")))))
  
  for(m in 1:length(trts)){
    DE_gene_list<- eqtl_table[eqtl_table$timepoint==trts[m], "gene"]
    cluster_table_filtered<- cluster_table[cluster_table$node %in% DE_gene_list, c("cluster","node")]
    DE_genes_TF<-data.frame(table(cluster_table_filtered[,"cluster"]))
    cluster_table_filtered<- merge(cluster_table_filtered, DE_genes_TF, by.x="cluster", by.y="Var1", all.x=T)
    DE_TF<- data.frame(table(cluster_table[,"cluster"]))
    DE_TF<- DE_TF[DE_TF$Var1 %in% DE_genes_TF$Var1,]
    #total DE in bin, total DE trt, sum not DE in bins listed in trt, genes in bin
    hypergeom_test[match(DE_genes_TF$Var1, unlist(hypergeom_test$cluster)), 
                   paste0(trts[m],"_HG_pvalue")]<- phyper(DE_genes_TF$Freq, length(DE_gene_list), 
                                                          sum(DE_TF$Freq - DE_genes_TF$Freq), 
                                                          DE_TF$Freq, lower.tail = FALSE)
    
    hypergeom_test[unlist(hypergeom_test$cluster) %in% DE_genes_TF$Var1,paste0(trts[m],"_HG_padj")]<- p.adjust(hypergeom_test[unlist(hypergeom_test$cluster) %in% DE_genes_TF$Var1,paste0(trts[m],"_HG_pvalue")], method = "BH")
  }
  hypergeom_test<- hypergeom_test[rowSums(is.na(hypergeom_test[,-1])) != ncol(hypergeom_test[,-1]), ]
  return(hypergeom_test)
}

clusters_mla_df<- read.csv("clusters_mla_noMLA_df.csv", row.names = 1)
Mla_associations<- read.csv("Mla_associations.csv", row.names = 1)
node_stats_mla<- read.csv("node_stats_MLA.csv", row.names = 1)
Mla_associations<- Mla_associations[Mla_associations$gene %in% unique(node_stats_mla$node), ]
HG_mla_eqtl<-HG_test_eqtl_cluster(cluster_table=clusters_mla_df, eqtl_table=Mla_associations)

#look for DE genes in modules
HG_test_cluster<- function(trts, go_cluster_table, DE_table, p_adj, name){
  library(effsize)
  #table trts by bins
  go_cluster_table$ID_group<- paste0(go_cluster_table$group, "/", go_cluster_table$ID)
  cluster_table<- data.frame(ID=unlist(apply(go_cluster_table, 1, FUN=function(x)rep(x[1], length(strsplit(x[8],"/")[[1]])))), 
                             ID_cluster=unlist(apply(go_cluster_table, 1, FUN=function(x)rep(x[14], length(strsplit(x[8],"/")[[1]])))), 
                             Target= unlist(apply(go_cluster_table, 1, FUN=function(x)strsplit(x[8],"/")[[1]])))
  hypergeom_test <- data.frame(matrix(ncol = 1+2*length(trts), nrow = length(unique(cluster_table$ID_cluster)))) 
  hypergeom_test[,1]<- unique(cluster_table$ID_cluster)
  colnames(hypergeom_test)<- c("ID_cluster", unlist(sapply(trts, function(x) c(paste0(x,"_HG_pvalue"), paste0(x,"_HG_padj")))))
  
  DE_genes <- data.frame(matrix(ncol = 1+length(trts), nrow = length(unique(cluster_table$ID)))) 
  DE_genes[,1]<- unique(cluster_table$ID)
  colnames(DE_genes)<- c("ID", unlist(sapply(trts, function(x) paste0(x,"_DE_genes"))))
  
  for(m in 1:length(trts)){
    DE_gene_list<- DE_table[DE_table[,paste0(trts[m],"_padj")]<= p_adj ,"gene"]
    DE_gene_list<- DE_gene_list[grep("HORVU", DE_gene_list)]
    cluster_table_filtered<- cluster_table[cluster_table$Target %in% DE_gene_list, c("ID", "ID_cluster","Target")]
    if(nrow(cluster_table_filtered)>0){
      DE_genes[match(unique(cluster_table_filtered$ID), unlist(DE_genes$ID)), paste0(trts[m],"_DE_genes")]<- sapply(FUN=function(x)paste(cluster_table_filtered[cluster_table_filtered$ID==x,"Target"], collapse = ";"), unique(cluster_table_filtered$ID))
      #plot_targets(cluster_table_filtered, trt=trts[m], name)
      
      DE_genes_TF<-data.frame(table(cluster_table_filtered[,"ID_cluster"]))
      cluster_table_filtered<- merge(cluster_table_filtered, DE_genes_TF, by.x="ID_cluster", by.y="Var1", all.x=T)
      DE_TF<- data.frame(table(cluster_table[,"ID_cluster"]))
      DE_TF<- DE_TF[DE_TF$Var1 %in% DE_genes_TF$Var1,]
      #total DE in bin, total DE trt, sum not DE in bins listed in trt, genes in bin
      hypergeom_test[match(DE_genes_TF$Var1, unlist(hypergeom_test$ID_cluster)), 
                     paste0(trts[m],"_HG_pvalue")]<- phyper(DE_genes_TF$Freq, length(DE_gene_list), 
                                                            sum(DE_TF$Freq - DE_genes_TF$Freq), 
                                                            DE_TF$Freq, lower.tail = FALSE)
      
      hypergeom_test[unlist(hypergeom_test$ID_cluster) %in% DE_genes_TF$Var1,paste0(trts[m],"_HG_padj")]<- p.adjust(hypergeom_test[unlist(hypergeom_test$ID_cluster) %in% DE_genes_TF$Var1,paste0(trts[m],"_HG_pvalue")], method = "BH")
      
    }}
  hypergeom_test<- hypergeom_test[rowSums(is.na(hypergeom_test[,-1])) != ncol(hypergeom_test[,-1]), ]
  if(length(trts)==1){
    DE_genes<- DE_genes[!is.na(DE_genes[,-1]), ]
  }else{
    DE_genes<- DE_genes[rowSums(is.na(DE_genes[,-1])) != ncol(DE_genes[,-1]), ]
  }
  return(list(hypergeom_test, DE_genes))
}

FC_padj_table<-read.csv("hv_R2_genes_deseq2_results_pairwise_genotype_logfc_padj_tax_sp.csv", stringsAsFactors = F)
go_result_16<- read.csv("GO_results_16.csv", row.names = 1)
DE_trts<- paste0("t16_wt_mla6")
HG_test_GO_16<-HG_test_cluster(trts=DE_trts, go_cluster_table=go_result_16, DE_table=FC_padj_table, p_adj=0.001, name="GO_16_DE_genes")

