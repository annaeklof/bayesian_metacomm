
library(igraph)
library(NetIndices)
library(tidyverse)
source("functions.R")

A <- read_csv("../data/webs/web5351.csv") %>% create_adj_matrix
lay<-matrix(nrow=nrow(A),ncol=2) # create a matrix with one column as runif, the other as trophic level
lay[,1]<-1:nrow(A)
lay[,2]<-TrophInd(t(A))$TL
par(mar=c(.1,.1,.1,.1))
plot.igraph(A %>% graph_from_adjacency_matrix, 
            layout=lay,vertex.label=NA,vertex.size=2,edge.arrow.size=.5,edge.width=.5)
