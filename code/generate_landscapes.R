
rm(list=ls())

require(igraph)
require(tidyverse)

source("functions.R")

## Randomly draw patch coordinates. 
## Input:
## - N: number of patches in the landscape, 
## - dim: dimension of the landscape (1, 2 or 3). 
## Output: 
## - matrix with N rows and dim columns. 
generate_coordinates <- function(N, dim) {
    return(matrix(runif(N*dim, 0, 1), nrow=N, ncol=dim))
}

## Generate random landscapes
N <- 300
dim <- 2
k <- 1
seed <- 54321

for(N in N){
  dir <- paste0("../data/landscapes/")
  if(file.exists(dir)){
    print("Directory already exists.")
  }else{
    dir.create(file.path(dir)) ## create direcotry for ouput files
  } 
  for(k in 1:50){
    set.seed(seed*k)
    filename <- sprintf("landscape_%03d.csv", k)
    landscape <- generate_coordinates(N, dim) %>% as_tibble %>% rename(x=V1, y=V2) %>% 
                  write_csv(paste0(dir,filename))
  }
}