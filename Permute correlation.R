
#this file is for calculation of correlations for networks

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
expression_table<- read.csv("hv_R2_genes_norm_counts_de_tax_sp.csv", row.names = 1)

permute_cor <- function(interaction_table, expression_table, m="dist", genotypes= c("wt","mla6","rar3","dm","bln1")){  # or spearman, pearson
  expression_table<- expression_table[rowSums(expression_table)>100,grep(paste(genotypes, collapse = "|"), colnames(expression_table))]
  b <- data.frame(gene=sample(rownames(expression_table), nrow(expression_table), replace = FALSE, prob = NULL))
  for(i in seq(1,ncol(expression_table),3)){
    b <- cbind(b, expression_table[match(b$gene, rownames(expression_table)),c(i:(i+2))])
    #b <- merge(b, expression_table[,c(1,i:(i+2))], by='gene', all.x = TRUE)
    b$gene <- sample(rownames(expression_table), nrow(expression_table), replace = FALSE, prob = NULL)
  }
  rownames(b)<- b$gene
  b$gene <- NULL
  cor_data_b<- calc_corr(interaction_table, expression_table=b, m, log=F, genotypes)
  return(data.frame(table(round(cor_data_b$weight, digits = 4))))
}

genotype_list=list(c("wt","bln1"))
#genotype_list=list(c("rar3","mla6","dm"))

for (gen in 1:length(genotype_list)){
  result<- data.frame(cor=seq(0,1, 0.0001), freq=0)
  for (j in 1:10000){
    b <- permute_cor(interactome_R2_HC_MLA, expression_table, m="dist", genotypes= genotype_list[[gen]])
    result$freq[match(b$Var1, result$cor)] <- result$freq[match(b$Var1, result$cor)] + b$Freq
    print(j)
  }
  write.csv(result, file = paste0("permuted_",paste(genotype_list[[gen]], collapse = "_"),"_dist_tritex.csv"), row.names = F)
}

