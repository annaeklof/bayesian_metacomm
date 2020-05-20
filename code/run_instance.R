
require(Rcpp, quietly=TRUE)
sourceCpp("c_functions.cpp")
require(NetIndices, quietly=TRUE)
require(igraph, quietly=TRUE)
require(tidyverse, quietly=TRUE)
source("functions.R")

clargs <- commandArgs(trailingOnly=TRUE)
if (length(clargs) > 0) { ## command-line arguments
  web <- clargs[1] ## name of web file (w/ path & extension)
  landsc <- clargs[2] ## name of landscape file (w/ path & extension)
  spinput <- clargs[3] ## name of species input file (w/ path & extension)
  spf <- clargs[4] ## focal species for patch removal ("basal" or "top")
  outputdir <- clargs[5] ## ouptut directory
  alpha <- as.numeric(clargs[6]) ## 1st parameter of beta distr (Bayesian network)
  beta <- as.numeric(clargs[7]) ## 2nd parameter of beta distr (Bayesian network)
} else { ## sample input parameters, if no command line arguments are given
  web <- "../data/webs/web50350.csv"
  landsc <- "../data/landscapes/gridlike_2.csv"
  spinput <- "../data/species_input/web50350_008.csv"
  spf <- "basal"
  outputdir <- "../data/results/"
  alpha <- 1
  beta <- 1
}

seed <- 548261 ## seed for random numbers 
rem <- 10 ## number of patches removed per patch loss iteration
nreps <- 1000 ## number of iterations for Bayesian network
scenarios <- c("random", "best", "worst") ## set of scenarios

## 1: simulate metacommunity in the full landscape (pre patch loss)  
A <- create_adj_matrix(read_csv(web, col_types=cols())) ## create adjacency matrix from links
landscape <- as.matrix(read_csv(landsc, col_types=cols())) ## read landscape file
N <- length(landscape[,1]) ## number of patches in the landscape

## combine species and landscape information
spdat0 <- generate_spdat(web, spinput, landsc) %>%
  add_column(delta=0, lambda=0, p=0.5, pvalues=0)

spdat0 <- sim_metacomm(A, spdat0, nreps, alpha, beta) %>%
  add_column(prem=0,
             web=web,
             spinput=spinput,
             landsc=landsc,
             spf=spf,
             seed=seed,
             alpha=alpha,
             beta=beta,
             scenario="all")

spdat <- spdat0

flag <- FALSE
## if one of the following conditions is TRUE we stop here else we continue and simulate the habitat loss scenarios: 
if ((all(spdat0 %>% filter(FT!="basal") %>% pull(lambda) %>% unique < 1)) || ## stop if no consumer can persist 
    (max(spdat0$p)<(1e-10)) || ## stop if persistence too low 
    ((spf!="basal") & (spdat0 %>% filter(species==rownames(A)[nrow(A)]) %>% pull(lambda) %>% unique<1))){ ## stop if top species is spf and extinct 
    flag <- TRUE
}

focal_sp <- ifelse(spf=="basal", rownames(A)[1], rownames(A)[nrow(A)]) ## spf name

## print to log file: 
print(web)
print(landsc)
## baseline extinction probability pi: uniform or increases with trophic level (TLB)
print(paste("pi =", ifelse(length(unique(read_csv(spinput, col_types = cols())$pi))==1, unique(read_csv(spinput, col_types = cols())$pi), "TLB")))
## dispersal distance xi: uniform or increases with trophic level (TLB)
print(paste("xi =", ifelse(length(unique(read_csv(spinput, col_types = cols())$xi))==1, unique(read_csv(spinput, col_types = cols())$xi), "TLB")))
print(paste("spf patch loss:", spf, focal_sp)) 
## if one of these is TRUE we stop at the full landscape and do not simulate habitat loss: 
print(paste("all consumers extinct:", all(spdat0 %>% filter(FT!="basal") %>% pull(lambda) %>% unique < 1))) 
print(paste("top species extinct:", spdat0 %>% filter(species==rownames(A)[nrow(A)]) %>% pull(lambda) %>% unique<1)) ## only relevant if spf=="top"
print(paste("persistence too low:", max(spdat0$p)<(1e-10)))
print(paste("stopped after full landscape:", flag))
print(paste0("scenarios running for ",spf, " ", focal_sp, ": ", !flag))

if(!flag){
## 2: simualte metacommunities with patch loss
  for (scen in scenarios) {
    N0 <- N ## number of patches in landscape initially
    j <- 0
    rem0 <- rem
    pselec <- integer(0) ## initialize pselec
    crashed <- FALSE ## test if populations are still extant
    spdat1 <- spdat0 %>% mutate(scenario=scen)
    while (!crashed) {
      j <- j + rem0
      if(scen %in% "best"){
        ## remove patches starting with lowest relative patch values
        pselec <- spdat1 %>% ## vector of patches selected for removal
          filter(species==focal_sp) %>%
          top_n(-rem0, pvalues) %>% ## negative to select bottom rem0 rows
          pull(patch)
      }
      if(scen %in% "worst"){
        ## remove patches starting with highest relative patch values
        pselec <- spdat1 %>%  ## vector of patches selected for removal
          filter(species==focal_sp) %>%  
          top_n(rem0, pvalues) %>%
          pull(patch) 
      }
      if(scen=="random"){ ## remove patches randomly
        pselec <- sample(unique(spdat1$patch), rem0) ## vector of patches selected for removal
      }
      
      pselec <- ifelse(length(pselec)<=rep(rem,length(pselec)), pselec, pselec[1:rem]) ## make sure pselec has a max. length of rem
      spdat_rem_patches <- spdat1 %>%
        filter(patch %ni% pselec) %>% ## remove selected patches
        mutate(prem=j)
      spdat1 <- sim_metacomm(A, spdat_rem_patches, nreps, alpha, beta)
      spdat <- rbind(spdat, spdat1)
      N0 <- length(unique(spdat1$patch)) ## number of patches in landscape
      if (N0<=rem0) rem0 <- 1
      if ((spf!= "basal") & (spdat1 %>% filter(species==rownames(A)[nrow(A)]) %>% pull(lambda) %>% unique<1)){ ## stop if top species is spf and extinct
          crashed <- TRUE
          print(paste("spf extinct:", spf, focal_sp, scen, N0)) ## print to log file 
         
      } 
      if((all((spdat1 %>% filter(FT!="basal") %>% pull(lambda) %>% unique)<1)) || ## stop if all consumers have crashed
          (max(spdat1$p)<=(1e-10)) || (N0<=2)) { ## stop if persistence too low or less than 2 patches
        crashed <- TRUE
        print(paste("metacommunity crashed:", spf, focal_sp, scen, N0)) ## print to log file
      }
    }
  }
}

outfile <- paste0(outputdir, str_replace(basename(spinput), ".csv", ""), "_",
                  str_replace(basename(landsc), ".csv", ""), "_", spf, "_",
                  alpha, "_", beta, ".rds")
saveRDS(spdat, outfile)
