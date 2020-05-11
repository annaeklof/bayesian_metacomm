
require(Rcpp, quietly=TRUE)
sourceCpp("c_functions.cpp")
require(NetIndices, quietly=TRUE)
require(igraph, quietly=TRUE)
require(tidyverse, quietly=TRUE)
source("functions.R")

##spdat_ref <- read_csv("../data/results/test_R_alpha5_beta5.csv")

web <- "../data/webs/test.csv"
landsc <- "../data/landscapes/landscape_001.csv"
spinput <- "../data/species_input/test.csv"
spf <- "basal"
outputdir <- "../data/results/"
alpha <- 5
beta <- 5

seed <- 548261 ## seed for random numbers 
nreps <- 10000 ## number of iterations for Bayesian network

A <- create_adj_matrix(read_csv(web, col_types=cols())) ## create adjacency matrix from links
landscape <- as.matrix(read_csv(landsc, col_types=cols())) ## read landscape file
N <- length(landscape[,1]) ## number of patches in the landscape

## combine species and landscape information
spdat <- generate_spdat(web, spinput, landsc) %>%
  add_column(delta=0, lambda=0, p=0.5, pvalues=0)

system.time(
  spdat <- sim_metacomm(A, spdat, nreps, alpha, beta) %>%
    add_column(prem=0,
               web=web,
               spinput=spinput,
               landsc=landsc,
               spf=spf,
               seed=seed,
               alpha=alpha,
               beta=beta,
               scenario="all")
)

##write_csv(spdat, "../data/results/test_R_alpha5_beta5.csv")

par(mfrow=c(2,2))
plot(spdat_ref$p, spdat$p)
plot(spdat_ref$pvalues, spdat$pvalues)
plot(spdat_ref$delta, spdat$delta)
plot(spdat_ref$lambda, spdat$lambda)
