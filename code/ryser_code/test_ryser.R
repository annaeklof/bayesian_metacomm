
require(Rcpp)
sourceCpp("c_functions.cpp")
require(NetIndices)
require(igraph)
require(tidyverse)
source("functions.R")

web <- "../../data/ryser_data/webs/web_1_1.csv"
landsc <- "../../data/ryser_data/landscapes/lsc_70.csv"
spinput <- "../../data/ryser_data/spinput/test_input.csv"
spf <- "basal"
outputdir <- "../../data/results/"
alpha <- 1
beta <- 1
nreps <- 1000 ## number of iterations for Bayesian network

A <- create_adj_matrix(read_csv(web, col_types=cols())) ## adj matrix from links
landscape <- as.matrix(read_csv(landsc, col_types=cols())) ## read landscape file
N <- length(landscape[,1]) ## number of patches in landscape

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
               seed=0,
               alpha=alpha,
               beta=beta,
               scenario="all")
)

##write_csv(spdat, "../../data/results/test_ryser.csv")

