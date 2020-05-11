
require(Rcpp, quietly=TRUE)
sourceCpp("bayes_marginal.cpp")
require(NetIndices, quietly=TRUE)
require(igraph, quietly=TRUE)
require(tidyverse, quietly=TRUE)
source("functions_old.R")

##spdat_ref <- read_csv("../data/results/test_serengeti.csv")

web <- "../data/webs/serengeti.csv"
landsc <- "../data/landscapes/landscape_001.csv"
spinput <- "../data/species_input/serengeti_006.csv"
spf <- "basal"
outputdir <- "../data/results/"
alpha <- 1
beta <- 1

seed <- 548261 ## seed for random numbers 
nreps <- 1000 ## number of iterations for Bayesian network

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

##write_csv(spdat, "../data/results/test_serengeti.csv")

par(mfrow=c(2,2))
plot(spdat_ref$p, spdat$p)
plot(spdat_ref$pvalues, spdat$pvalues)
plot(spdat_ref$delta, spdat$delta)
plot(spdat_ref$lambda, spdat$lambda)

bind_rows(
  spdat_ref %>%
    select(patch, species, TL, p) %>%
    mutate(data="ref"),
  spdat %>%
    select(patch, species, TL, p) %>%
    mutate(data="sim")
) %>%
  spread(data, p) %>%
  ggplot +
  aes(x=ref, y=sim, colour=TL) +
  geom_point() +
  theme_bw()
