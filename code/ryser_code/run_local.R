
require(Rcpp)
sourceCpp("c_functions.cpp")
require(NetIndices)
require(igraph)
require(tidyverse)
source("functions.R")

clargs <- commandArgs(trailingOnly=TRUE)
if (length(clargs) > 0) { ## command-line arguments
  outfile <- clargs[1]
} else { ## sample input parameters, if no command line arguments are given
  outfile <- ""
}

nreps <- 1000 ## number of iterations for Bayesian network
tab <- expand_grid(
  web="../../data/ryser_data/webs/web_1_1.csv",
  landsc=Sys.glob("../../data/ryser_data/landscapes/*.csv"),
  spinput=Sys.glob("../../data/ryser_data/spinput/*.csv")[1],
  spf=c("basal"),
  outputdir="../../data/ryser_data/results/",
  alpha=c(1, 5),
  beta=c(1, 5),
  xiscale=10^(seq(-2, 0.1, l=11))/0.055/5)

heatmap <- tibble(xi=rep(0, nrow(tab)), N=0, alpha=0, beta=0, tl=0)

for (r in 1:nrow(tab)) {
  web <- tab$web[r]
  landsc <- tab$landsc[r]
  spinput <- tab$spinput[r]
  spf <- tab$spf[r]
  outputdir <- tab$outputdir[r]
  alpha <- tab$alpha[r]
  beta <- tab$beta[r]
  A <- create_adj_matrix(read_csv(web, col_types=cols())) ## adj matrix from links
  landscape <- as.matrix(read_csv(landsc, col_types=cols())) ## landscape file
  spdat0 <- generate_spdat(web, spinput, landsc) %>%
    add_column(delta=0, p=0.5) %>%
    mutate(xi=xi*tab$xiscale[r])
  spdat <- sim_metacomm(A, spdat0, nreps, alpha, beta) %>%
    add_column(prem=0,
               web=web,
               spinput=spinput,
               landsc=landsc,
               spf=spf,
               seed=0,
               alpha=alpha,
               beta=beta,
               scenario="all")
  heatmap$xi[r] <- tab$xiscale[r]
  heatmap$N[r] <- nrow(landscape)
  heatmap$alpha[r] <- tab$alpha[r]
  heatmap$beta[r] <- tab$beta[r]
  heatmap$tl[r] <- spdat %>% filter(p>1e-10) %>% pull(TL) %>% max
  cat(paste(r, "/", nrow(tab), "finished \n"))
}

if (outfile!="") write_csv(heatmap, outfile)

