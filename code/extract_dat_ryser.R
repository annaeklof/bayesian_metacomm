
##----------------------- Packages --------------------##

require(igraph)
require(tidyverse)

##----------------------- Functions --------------------##

## Make edge list from adjacency matrix. 
## Input: 
## - webfile: file containig the weighted adjacency matrix (web.out)
## Output: 
## - edgelist: edgelist, i.e who eats whom (web.csv) 
make_edgelist <- function(webfile){
  A <- read.table(webfile) %>% as.matrix ## make web into a matrix 
  Amat <- ifelse(A > 0, 1, 0) ## make into binary matrix (0,1)
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

'%ni%' <- function(x,y)!('%in%'(x,y))

numextract <- function(string){ 
  temp <- gsub("[^[:digit:].]", "", string)
  as.numeric(str_extract(temp, "\\-*\\d+\\.*\\d*"))
}

##----------------------- Script --------------------##

ryser <- list.dirs("../data/ryser")[-1]

webfiles <- list.files(path=ryser, pattern = "^web(.*).out$", all.files = TRUE, 
                         full.names = TRUE, no.. = TRUE, include.dirs = TRUE, recursive=TRUE) 

paramfiles <- list.files(path=ryser, pattern = "^params(.*).out$", all.files = TRUE, 
                           full.names = TRUE, no.. = TRUE, include.dirs = TRUE, recursive=TRUE) 

for(k in 1:length(paramfiles)){
  
  webfile <- webfiles[k] 
  name1 <- dirname(webfile) %>% numextract()
  name2 <- basename(webfile) %>% numextract()
  
  web <- webfiles[k] %>% read.table() %>% as.matrix() %>% as_tibble() %>% 
    write_csv(paste0("../data/ryser_webs/Amat_",name1,"_",name2,".csv"))
  
  edgelist <- make_edgelist(webfile) %>% 
    write_csv(paste0("../data/ryser_webs/web_",name1,"_",name2,".csv"))
  
  S <- web %>% nrow
  S_c <- ((web %>% rowSums > 0) * 1) %>% sum
  S_b <- ((web %>% rowSums == 0) * 1) %>% sum
  # S_b + S_c == S

  params <- paramfiles[k] %>% read_csv
  name3 <- dirname(paramfiles[k]) %>% numextract()
  name4 <- basename(paramfiles[k]) %>% numextract()
  
  name1 <- basename(webfile) %>% str_split(".out",simplify = TRUE) %>%
    as_tibble() %>% pull(V1) 
  
  coords <-  params %>% filter(variable == "patch location XY") %>% #extract coords
    pull(value) #vector with XY alternating

  x <- coords[seq(1, length(coords), by = 2)] #extract X coordinates
  y <- coords[seq(2, length(coords), by = 2)] #extract Y coordinates
  Z <- length(x)
  # plot(x,y, xlim=c(0,1), ylim=c(0,1))
  landscape <- cbind(x,y) %>% as_tibble %>% #generate landscape matrix
    write_csv(paste0("../data/ryser_landscapes/landscape_",name3,"_",name4,".csv"))
  
  # distmat <- dist(landscape) %>% as.matrix #distance matrix containing the euclidean distances between all patche
  
  D_0 <-  params %>% filter(variable == "RGG dispersal distance D_0") %>% #extract D_0
    pull(value)  %>% #scaling factor for landscape size
    as_tibble() %>% rename(D_0 = value) %>% 
    write_csv(paste0("../data/ryser_landscapes/D0_",name3,"_",name4,".csv"))
  
  SW <-  params %>% filter(variable == "SW dispersal matrix") %>% #extract success-weighted dispersal matrix
   pull(value) 
  
  con <- edgelist %>% pull(consumer) %>% unique
  res <- edgelist %>% pull(resource) %>% unique
  species_list <- c(con,res) %>% unique %>% sort
  
  s <- 0
  tmp <- SW[(1+s*Z*Z):(Z*Z+s*Z*Z)] %>% as_tibble() %>% rename(!!(species_list[s+1]) := value)
  
  for(s in 1:(S-1))
    tmp <- tmp %>% add_column(!!(species_list[s+1]) := SW[(1+s*Z*Z):(Z*Z+s*Z*Z)])
  
  # s <- 0
  # SW <- tmp$species001
  # SW_mat <- matrix(SW[(1+s*Z*Z):(Z*Z+s*Z*Z)], ncol=Z, nrow=Z)
  
  tmp %>% write_csv(paste0("../data/ryser_landscapes/swmat_",name3,"_",name4,".csv"))
  
}


