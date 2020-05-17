
## Script to create edgelists using adjacency matrices from C-script as input 
rm(list=ls())
require(igraph)
require(tidyverse)

## Make edge list from adjacency matrix. 
## Input: 
## - c_webfile: file containig the adjacency matrix (web.out)
## Output: 
## - edgelist: edgelist, i.e who eats whom (web.csv) 
make_edgelist <- function(cwebfile ){
  Amat <- read.table(cwebfile) %>% as.matrix ## make web into a matrix 
  colnames(Amat) <- rownames(Amat) <- paste0("species", sprintf("%03d", 1:nrow(Amat))) 
  g <- graph_from_adjacency_matrix(Amat, mode = "directed") ## generate graph
  if(!is_dag(g)) print("Web is cyclic.") ## test if web contains cycles
  if(is_dag(g)){
    g2 <- graph_from_adjacency_matrix(Amat, mode = "directed") ## regenerate graph 
    web <- g2 %>% get.edgelist %>% as_tibble %>% 
      rename(consumer=V1, resource=V2)
  }
  return(web)
}

webs <- Sys.glob("../data/c_webs/A_*.csv") ## list of web files
outputdir <- Sys.glob("../data/c_webs/") ## list of web files
w <- webs[1]

for(w in webs){
  temp <- w %>% str_split(.,"A_",simplify=TRUE)
  outfile <- paste0(temp[1], "web", temp[,2])
  make_edgelist(w) %>% write_csv(outfile)
}


