
require(XRJulia) ## package allowing for Julia interfacing
require(igraph)
require(tidyverse) ## for efficient data manipulation and plotting

## Julia packages for generating random numbers and using Beta distributions
juliaUsing("Random, Distributions")

## Monte Carlo evaluation of Bayesian network, written in Julia
bayes_ntw_MC_jul <- juliaEval("
function bayesianNetworkMonteCarlo(A, pBaseline, a, b, reps)
    sA = sum(A, dims = 2)
    nRows = length(sA)
    lastBasal = maximum(findall(sA .== minimum(sA)))[1][1]
    marginals = zeros(nRows)
    for rep in 1:1:reps
        extant = ones(nRows)
        for i in 1:1:lastBasal
            extant[i] = rand() < pBaseline[i] ? 0 : 1
        end 
        for i in (lastBasal + 1):1:size(A, 1)
            frac = 1 - (extant' * A[i,:])[1] / sum(A[i,:])
            pExt = pBaseline[i] + (1 - pBaseline[i]) * cdf(Beta(a, b), frac)
            extant[i] = rand() < pExt ? 0 : 1
        end
        marginals = marginals + extant
    end
    return(vec(marginals) / reps)
end
")

## turn the above Julia code to an R function
bayes_jul <- JuliaFunction(bayes_ntw_MC_jul)
bayes_ntw_MC <- function(A, Pb, a, b, N) juliaGet(bayes_jul(A, Pb, a, b, N))

## Create adjacency matrix from a data frame of edges.
## Input:
## - edgelist: data frame with two columns of taxon names, where consumers
##             are in column 1 and resources in column 2
## Output:
## - The adjacency matrix, with A[i,j]=1 if species i eats j and 0 otherwise.
##   Species are topologically sorted to make the matrix upper triangular.
##   The rows and columns of the matrix are labeled by taxon names.
create_adj_matrix <- function(edgelist) {
  ## Create adjacency matrix A; links point from resource to consumer, so
  ## the two columns of the edge list data frame are flipped
  A <- graph_from_data_frame(d=edgelist[,2:1]) %>% ## create igraph graph
    as_adjacency_matrix %>% ## convert to igraph adjacency matrix
    as.matrix %>% ## convert to a regular matrix
    t ## transpose result (so A[i,j] is 1 if i eats j and 0 otherwise)
  ## Find a sorting of A's rows and columns to make A lower triangular
  o <- graph_from_adjacency_matrix(A) %>% topo_sort(mode="in")
  ## Return sorted, lower triangular adjacency matrix A
  return(A[o,o])
}

## for benchmarking purposees, here is the equivalent R function
bayesian_network_monte_carlo <- function(A, Pb, a=1, b=1, reps=1) {
    ind_basal <- max(which(rowSums(A)==0))
    marginals <- rep(0, nrow(A))
    for (i in 1:reps) {
        extant <- rep(1, nrow(A))
        for (i in 1:ind_basal) if (runif(1)<Pb[i]) extant[i] <- 0
        for (i in (ind_basal+1):nrow(A)) {
            frac <- 1 - as.numeric(A[i,] %*% extant / sum(A[i,]))
            Pext <- Pb[i] + (1 - Pb[i]) * pbeta(frac, a, b)
            if (runif(1)<Pext) extant[i] <- 0
        }
        marginals <- marginals + extant
    }
    return(marginals/reps)
}

serengeti <- read_csv("../data/serengeti/serengeti_edgelist.csv") %>%
    create_adj_matrix ## load Serengeti food web edgelist, and create adj matrix
Pb <- rep(0.1, nrow(serengeti)) ## baseline extinction probabilities

a_par <- 1 ## alpha parameter of Beta CFD
b_par <- 1 ## beta parameter of Beta CFD
N <- 10000 ## number of iterations of Monte Carlo simulation

## benchmark no. 1: Julia code
t0_jl <- Sys.time() ## record time at start of evaluation of Julia function
bayes_ntw_MC(serengeti, Pb, a_par, b_par, N) ## Bayes eval
t1_jl <- Sys.time() ## record time at the end of evaluation

## benchmark no. 2: R code
t0_R <- Sys.time() ## record time at start of evaluation of R function
bayesian_network_monte_carlo(serengeti, Pb, a_par, b_par, N) %>%
    print ## Bayes eval
t1_R <- Sys.time() ## record time at the end of evaluation

## print results of benchmarks on screen
t1_jl-t0_jl ## display time it took to run the Julia function
t1_R-t0_R ## display time it took to run the R function

