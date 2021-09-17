# This function calculates net statistics

calculate_net_stats<- function(weighted_network, annotation){
  barley_net<- weighted_network
  
  getECC <- function(graph, e) {
    ve <- ends(graph, e)
    numerator <- length(intersect(neighbors(graph, v = ve[1], mode = 1), neighbors(graph, v = ve[2], mode = 1))) + 1
    denominator <- min(degree(graph = graph, v = ve[1], loops = F), degree(graph = graph, v = ve[2], loops = F))
    ECC <- numerator / denominator
    ECC
  }
  
  E(barley_net)$ecc <- sapply(X = E(barley_net), simplify = T, FUN = function(x){getECC(barley_net, x)})
  
  Edges_barley_net <- data.frame(
    Source=as_edgelist(barley_net)[,1], 
    Target=as_edgelist(barley_net)[,2],
    weight= E(barley_net)$weight,
    ecc= E(barley_net)$ecc)
  
  rankingByLambda <- function(graph, lambda, compOne, compTwo) {
    # check if given components/attributes actually exist
    if(!(compOne %in% edge_attr_names(graph))) {
      stop(paste(compOne,"isn't one of the edge attributes of graph"))
    }
    if(!(compTwo %in% edge_attr_names(graph))) {
      stop(paste(compTwo, "isn't one of the edge attributes of graph"))
    }
    # return the score
    sapply(X = V(graph), simplify = TRUE, FUN = function(v){
      incidentEdges <- E(graph)[unlist(graph[[v,,edges=TRUE]])]
      sum((lambda * edge_attr(graph = graph, name = compOne, index = incidentEdges)) +
            ((1-lambda) * edge_attr(graph = graph, name = compTwo, index = incidentEdges)))
    })
  }
  
  V(barley_net)$bdc <- rankingByLambda(graph = barley_net, lambda = 0.8, compOne = "ecc", compTwo = "weight")
  V(barley_net)$degree <- degree(graph = barley_net)
  V(barley_net)$closeness <- closeness(graph = barley_net)
  V(barley_net)$betweenness <- betweenness(graph = barley_net)
  V(barley_net)$eigcent <- evcent(graph = barley_net)$vector
  #(barley_net)$sgc <- subgraph.centrality(graph = barley_net, diag = FALSE)
  
  graphRankingStats_barley_net <- data.frame(
    degree = V(barley_net)$degree,
    betweenness = V(barley_net)$betweenness,
    closeness = V(barley_net)$closeness,
    eigcent = V(barley_net)$eigcent,
    bdc = V(barley_net)$bdc,
    row.names = V(barley_net)$name)
  
  om <- 0.9
  graphRankingStats_barley_net$diffslc_09 <- (om*V(barley_net)$eigcent) + ((1-om)*graphRankingStats_barley_net$bdc)
  graphRankingStats_barley_net$node<- rownames(graphRankingStats_barley_net)
  graphRankingStats_barley_net<- merge(graphRankingStats_barley_net, annotation, by.x = "node", by.y="gene", all.x = T)
  return(list(Edges_barley_net, graphRankingStats_barley_net))
}

#example calculating the statistics using HvInt
library(igraph)
annotation<- read.table("annotation_HvR2.txt", sep=';', stringsAsFactors = F, quote = "\"")[,1:2]
colnames(annotation)<-c("gene","description")
interactome_R2_HC<- read.csv(file='d_weighted_interactome_R2_HC2.csv')
interactome_R2_HC$weight[is.na(interactome_R2_HC$weight)]<- min(interactome_R2_HC$weight, na.rm = T)/10
barley_net <- graph_from_data_frame(interactome_R2_HC, directed = FALSE, vertices = NULL)
barley_net<-igraph::simplify(barley_net)
net_stats<- calculate_net_stats(barley_net, annotation)
#write.csv(net_stats[[2]], file = "node_stats.csv")
#write.csv(net_stats[[1]], file = "edges_stats.csv")
