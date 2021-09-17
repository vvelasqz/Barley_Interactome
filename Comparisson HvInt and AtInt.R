#Comparisson HvInt and AtInt

interactome_R2_HC<- read.csv(file='d_weighted_interactome_R2_HC2.csv')
arab_interactome<- read.csv(file='arab_interactome.csv', row.names = F)

library(igraph)
barley_net <- graph_from_data_frame(interactome_R2_HC, directed = FALSE, vertices = NULL)
barley_net<-igraph::simplify(barley_net)
arab_net<- graph_from_data_frame(arab_interactome, directed = FALSE, vertices = NULL)
arab_net<-igraph::simplify(arab_net)

#fit power laws
fit_power_law(degree(barley_net), xmin = NULL, start = 2, force.continuous = FALSE,implementation = "plfit")
fit_power_law(degree(arab_net), xmin = NULL, start = 2, force.continuous = FALSE,implementation = "plfit")

#calculate random network to check small-world 
simulated_random_barley<- erdos.renyi.game(vcount(barley_net), ecount(barley_net), type =  "gnm")
simulated_random_arab<- erdos.renyi.game(vcount(arab_net), ecount(arab_net), type =  "gnm")

paste0("Diameter of Arabidopsis network: ",diameter(arab_net, directed = F))
paste0("Radius of Arabidopsis network: ",radius(arab_net))

paste0("Global clustering coefficient of Arabidopsis network: ",transitivity(arab_net))
paste0("Average shortest path length of Arabidopsis network: ",mean_distance(arab_net, directed = F))

paste0("Global clustering coefficient of simulated Arabidopsis random network: ",transitivity(simulated_random_arab))
paste0("Average shortest path length of simulated Arabidopsis random network: ",mean_distance(simulated_random_arab, directed = F))

paste0("Global clustering coefficient ratio arab_net/simulated random network: ",transitivity(arab_net)/transitivity(simulated_random_arab))
paste0("Average shortest path length ratio arab_net/simulated random network: ", mean_distance(arab_net, directed = F)/mean_distance(simulated_random_arab, directed = F))
