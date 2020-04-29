
require(XRJulia) ## package allowing for Julia interfacing
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
bayes_ntw_MC <- JuliaFunction(bayes_ntw_MC_jul)

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

A <- matrix(c(0,1,1,0,0,1,0,0,0), 3, 3) ## adj. matrix of simple omnivory module
Pb <- c(0.1, 0.7, 0.4) ## baseline extinction probabilities

## benchmark no. 1: Julia code
t0_jl <- Sys.time() ## record time at start of evaluation of Julia function
bayes_ntw_MC(A, Pb, 1, 1, 10000) %>% ## one million evals of BN
    juliaGet %>% ## convert vector from Julia to R format
    print ## display marginal persistence probabilities on screen
t1_jl <- Sys.time() ## record time at the end of evaluation

## benchmark no. 2: R code
t0_R <- Sys.time() ## record time at start of evaluation of R function
bayesian_network_monte_carlo(A, Pb, 1, 1, 10000) %>% ## one million evals of BN
    print ## display marginal persistence probabilities on screen
t1_R <- Sys.time() ## record time at the end of evaluation

## print results of benchmarks on screen
t1_jl-t0_jl ## display time it took to run the Julia function
t1_R-t0_R ## display time it took to run the R function

