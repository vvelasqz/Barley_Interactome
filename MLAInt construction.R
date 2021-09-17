#find shortest paths among MLA interactors (not passing throug MLA), use interactome_R2_HC or rm MLA

library(igraph)
interactome_R2_HC<- read.csv(file='d_weighted_interactome_R2_HC2.csv')
interactome_R2_HC$weight[is.na(interactome_R2_HC$weight)]<- min(interactome_R2_HC$weight, na.rm = T)/10
barley_net0 <- graph_from_data_frame(interactome_R2_HC, directed = FALSE, vertices = NULL)
barley_net0<-igraph::simplify(barley_net0)
components(barley_net0)$csize
components_barley<- decompose(barley_net0, mode = "weak")
barley_net0 <- components_barley[[1]]

MLA_interactors_all<- read.csv("MLA_targets.csv")[,1]
MLA_interactors_all<- unique(c(MLA_interactors_all, interactome_R2_HC[interactome_R2_HC$source=="HORVU.MOREX.r2.1HG0009580", "target"]))
MLA_interactors<- MLA_interactors_all[MLA_interactors_all %in% V(barley_net0)$name]

"HORVU.MOREX.r2.1HG0009580" %in% V(barley_net0)$name
barley_net_mla<- barley_net0

MLA_numbers<- which(V(barley_net_mla)$name %in% MLA_interactors)
for(i in 1:length(MLA_numbers)){
  a<- all_shortest_paths(barley_net_mla, from=MLA_numbers[i], to = V(barley_net_mla)[MLA_numbers[-1:-i]])
  if(i==1){
    result<- unlist(a$res)
  }else{
    result<- c(result, unlist(a$res))
  }
}
names_result<- unique(names(result))
result<- unique(result)

#get subnetwork short paths and 2nd neighbors
#mla "HORVU.MOREX.r2.1HG0009580"
neighbors(barley_net_mla, names_result)
induced_subgraph(barley_net_mla,names_result)
MLA_net<- induced_subgraph(barley_net_mla, c(names(neighbors(barley_net_mla, neighbors(barley_net_mla, names_result))), names(neighbors(barley_net_mla, names_result)), names_result))

Edges_MLA_net <- data.frame(Source=as_edgelist(MLA_net)[,1], Target=as_edgelist(MLA_net)[,2])
MLA_interactions<- data.frame(Source= "HORVU.MOREX.r2.1HG0009580", Target=MLA_interactors_all)
Edges_MLA_net <-rbind(Edges_MLA_net, MLA_interactions)
Edges_barley_net<- read.csv("edges_stats_MLA.csv", row.names = 1)
Edges_MLA_interactors<- Edges_barley_net[Edges_barley_net$Source %in% MLA_interactors_all| Edges_barley_net$Target %in% MLA_interactors_all, c("Source", "Target")]
Edges_MLA_net <-rbind(Edges_MLA_net, Edges_MLA_interactors)
Edges_MLA_net<- Edges_MLA_net[!duplicated(t(apply(Edges_MLA_net, 1, sort))), ]
Edges_MLA_net<- merge(Edges_MLA_net, interactome_R2_HC_MLA, by.x=c("Source", "Target"), by.y=c("source", "target"), all.x = T)
Edges_MLA_net<- merge(Edges_MLA_net, interactome_R2_HC_MLA, by.x=c("Source", "Target"), by.y=c("target", "source"), all.x = T)
Edges_MLA_net$weight<- ifelse(is.na(Edges_MLA_net$weight.x), Edges_MLA_net$weight.y, Edges_MLA_net$weight.x)
Edges_MLA_net<- Edges_MLA_net[, c(1,2,5)]
Edges_MLA_net$weight[is.na(Edges_MLA_net$weight)]<- min(interactome_R2_HC_MLA$weight, na.rm = T)/10

#write.csv(Edges_MLA_net, file = "MLAInt.csv")
