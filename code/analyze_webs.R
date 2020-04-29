rm(list=ls())
require(igraph)
require(tidyverse)
source("functions.R")
   
   
   ## web information: 
   webfiles <- Sys.glob("../data/webs/*.csv") ## list of web files
   
   i <- 1
   for(webfile in webfiles){
     webname <- str_split(basename(webfile), (".csv"), simplify=TRUE)[1] %>% 
                        str_split(.," ",simplify = TRUE) %>% paste(.,collapse="_")
     Amat <- read_csv(webfile, col_types=cols()) %>% create_adj_matrix ## create adjacency matrix from links
     g <- graph_from_adjacency_matrix(Amat, mode = "directed") ## generate graph
     if(!is_dag(g)) print("Web is not acyclic.") 
     if(is_dag(g)){
      S = (length(Amat[1,]))
      Sb =  length(which(rowSums(Amat) == 0)) 
      Sc =  length(which(rowSums(Amat) != 0)) 
      if(Sb+Sc!=S) print("Species numbers do not add up.")
      tib <- tibble(web = webname, 
                 S = S,
                 Sb = Sb,  
                 Sc = Sc,
                 connectance = calc_NetInd(Amat)$connectance,  
                 generality = calc_NetInd(Amat)$generality, 
                 sd.generality = calc_NetInd(Amat)$sd.generality, 
                 vulnerability  = calc_NetInd(Amat)$vulnerability,
                 sd.vulnerability = calc_NetInd(Amat)$sd.vulnerability, 
                 modularity = calc_NetInd(Amat)$modularity,  
                 TL.max = calc_NetInd(Amat)$TL.max, 
                 TL.mean = calc_NetInd(Amat)$TL.mean, 
                 TL.sd = calc_NetInd(Amat)$TL.sd, 
                 OI.max = calc_NetInd(Amat)$OI.max, 
                 OI.mean = calc_NetInd(Amat)$OI.mean, 
                 OI.sd = calc_NetInd(Amat)$OI.sd
                 )
   
      if (i == 1) out <- tib else out <- rbind(out, tib) 
      i <- i+1
     }
   }  
   
   out %>% write_csv("../data/webinfo.csv") ## write web info to file 
    
## species distribution across trophic levels: 
dat <- readRDS( "../data/summaries/modified_summary.rds")

species_per_ITL <- list()
species_per_TL05 <- list()

for(webs in unique(dat$web)){
   for(i in unique(dat$ITL))
      species_per_ITL[[webs]][[i]] <- dat %>% filter(web==webs, ITL==i) %>% pull(species) %>% unique %>% length
   
   for(j in 1:length(unique(dat$TL05)))
      species_per_TL05[[webs]][j] <- dat %>% filter(web==webs, TL05==unique(dat$TL05)[j]) %>% pull(species) %>% unique %>% length
}

t <- as_tibble(species_per_ITL) %>% mutate(TL = unique(dat$ITL)) ## make into tibble
write_csv(t, "../data/species_per_ITL.csv") ## save
xtable(t) 

t <- as_tibble(species_per_TL05) %>% mutate(TL = unique(dat$TL05)) ## make into tibble
write_csv(t, "../data/species_per_TL05.csv") ## save
xtable(t) 
View(t)

## analyze input parameters xi and pi: 
range_pi <- list()
range_xi <- list()

for(webs in unique(dat$web))
{
   range_pi[[webs]] <- dat %>%
      filter(web==webs) %>% pull(pi) %>% unique %>% range %>% round(.,4)
   range_xi[[webs]] <- dat %>%
      filter(web==webs) %>% pull(xi) %>% unique %>% range %>% round(.,4)
}

t_pi <- range_pi %>% as_tibble %>% write_csv("../data/range_pi.csv") ## make into tibble and save
t_xi <- range_xi %>% as_tibble %>% write_csv("../data/range_xi.csv") ## make into tibble and save


