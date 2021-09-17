#now go term analysis of top ranked nodes with diffslc
library(clusterProfiler)
library(GO.db)
library(AnnotationDbi)

GO_terms<- read.delim("annotation_HvR2.txt", sep = ";", stringsAsFactors = F, header = F, fill=T)
colnames(GO_terms)<- c("gene", "description", "go_terms")
s <- strsplit(GO_terms$go_terms, split = ",")
GO_terms<- data.frame(gene = rep(GO_terms$gene, sapply(s, length)), description = rep(GO_terms$description, sapply(s, length)),GO = unlist(s))
rm(s)

term2gene=GO_terms[, c("GO", "gene")]
term2name=GO_terms[, c("GO", "description")]
uniqueGO<- read.csv("uniqueGO.csv")

# Now clustering for HvInt and GO analysis 
interactome_R2_HC<- read.csv(file='d_weighted_interactome_R2_HC2.csv')
interactome_R2_HC$weight[is.na(interactome_R2_HC$weight)]<- min(interactome_R2_HC$weight, na.rm = T)

library(igraph)
barley_net <- graph_from_data_frame(interactome_R2_HC, directed = FALSE, vertices = NULL)
barley_net<-igraph::simplify(barley_net)
components(barley_net)$csize
components_barley<- decompose(barley_net, mode = "weak")
barley_net1 <- components_barley[[1]]
clusters_1<- cluster_walktrap(barley_net1, steps = 100)
clusters_net_df<-data.frame(node=clusters_1$names, cluster=clusters_1$membership)
table(clusters_net_df$cluster)

go_result<- data.frame()
pdf(paste0("GO_HvInt_clusters.pdf"), width = 12, height = 12, fonts = "ArialMT", pointsize = 18)
for (cluster in unique(clusters_net_df$cluster)) {
  gene_list_cluster<- clusters_net_df[clusters_net_df$cluster==cluster, "node"]
  x_horvu <- enricher(gene_list_cluster, TERM2GENE=term2gene, TERM2NAME=term2name)
  
  if(!is.null(x_horvu)){
    print(barplot(x_horvu, title = paste0(cluster, " Hv Int cluster enriched proteins")))
    print(dotplot(x_horvu, title = paste0(cluster, " Hv Int cluster enriched proteins")))
    x_horvu<- merge(as.data.frame(x_horvu), uniqueGO, by.x="ID", by.y="GO")
    x_horvu$group<- cluster
    go_result<- rbind(go_result, x_horvu)
  }
}
dev.off()

#write.csv(go_result, "GO_HvInt.csv")



