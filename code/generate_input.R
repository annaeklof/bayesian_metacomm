rm(list=ls())
require(NetIndices)
require(igraph)
require(tidyverse)

source("functions.R")

clargs <- commandArgs(trailingOnly=TRUE)
if (length(clargs) > 0) { ## command-line arguments
  outputdir <- clargs[1] ## output directory 
  webdir <- clargs[2] ## web directory 
  webs <- clargs[3:length(clargs)] ## webs
} else { ## sample input parameters, if no command line arguments are given
  outputdir <-  "../data/species_input/"
  webdir <-  "../data/webs/"
  webs <- Sys.glob("../data/webs/web50350.csv")
}

webs <- c("../data/webs/web503501.csv", 
          "../data/webs/web503502.csv", 
          "../data/webs/web503503.csv", 
          "../data/webs/web503504.csv", 
          "../data/webs/web503505.csv")

webs <- sub('\\.csv$', '', basename(webs)) 

for (w in webs) {
  webfile <- paste0(webdir, w, ".csv")
  print(webfile)
  kshape <- "Exponential"
  
  ## create species_input file 
  A <- read_csv(webfile) %>% create_adj_matrix ## create adjacency matrix from links
  
  n <- 1 ## offset for file ending
  
  ## Eklöf et al. 2013: pi = 0.2 * TL/mean(TL) 
  pi0 <- 0.2
  xivec <- seq(0.01, 0.1, l=5)
  
  ## xi same for all species 
  for(xi0 in xivec){
    dat <- tibble(species = rownames(A)) %>% 
      add_column(pi = 0, xi = 0, kernel = kshape, TL = TrophInd(t(A))$TL, 
                 OI = TrophInd(t(A))$OI, FT = get_feeding_type(webfile))
    
    for(sp in 1:length(dat$species)){
      dat$pi[sp] <- pi0 * (dat$TL[sp]/mean(dat$TL))
      dat$xi[sp] <- xi0
    }
    dat$kernel <- kshape
    
    outfile <- paste0(outputdir, w, sprintf("_%03d.csv", n))
    dat %>% write_csv(outfile)     
    n <- n + 1
  }
  
  ## xi increases with trophic level 
  dat <- tibble(species = rownames(A)) %>% 
    add_column(pi = 0, xi = 0, kernel = kshape, TL = TrophInd(t(A))$TL, 
               OI = TrophInd(t(A))$OI, FT =get_feeding_type(webfile))
  
  for(sp in 1:length(dat$species)){
    dat$pi[sp] <- pi0 * (dat$TL[sp]/mean(dat$TL))
    dat$xi[sp] <- mean(xivec) * (dat$TL[sp]/mean(dat$TL))
  }
  dat$kernel <- kshape
  
  outfile <- paste0(outputdir, w, sprintf("_%03d.csv", n))
  dat %>% write_csv(outfile)     
  n <- n + 1
  
  ## pi same for all species (average over TL) 
  dat <- tibble(species = rownames(A)) %>% 
    add_column(pi = 0, xi = 0, kernel = kshape, TL = TrophInd(t(A))$TL, 
               OI = TrophInd(t(A))$OI, FT = get_feeding_type(webfile))
  
  for(sp in 1:length(dat$species))
    dat$xi[sp] <- mean(xivec) * (dat$TL[sp]/mean(dat$TL))
  
  dat$pi <- pi0
  dat$kernel <- kshape
  
  outfile <- paste0(outputdir, w, sprintf("_%03d.csv", n))
  dat %>% write_csv(outfile)     
  n <- n + 1
  
  ## xi and pi same for all species (average over TL) 
  dat <- tibble(species = rownames(A)) %>% 
    add_column(pi = 0, xi = 0, kernel = kshape, TL = TrophInd(t(A))$TL, 
               OI = TrophInd(t(A))$OI, FT = get_feeding_type(webfile))
  
  dat$pi <- pi0 
  dat$xi <- mean(xivec)
  dat$kernel <- kshape
  
  outfile <- paste0(outputdir, w, sprintf("_%03d.csv", n))
  dat %>% write_csv(outfile)     
  n <- n + 1
  
}



