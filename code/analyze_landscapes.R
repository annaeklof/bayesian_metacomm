rm(list=ls())
require(igraph)
require(tidyverse)
require(ggplot2)
require(foreach)

source("functions.R")

####------------------- landscape characteristics ---------------------####
landscapefiles <- Sys.glob("../data/landscapes/*.csv") ## list of summary files

list <- lapply(landscapefiles, read_csv) 
matlist <- foreach(i = list) %do% as.matrix(i)
mean <- foreach(i = matlist) %do% mean(i)
sd <- foreach(i = matlist) %do% sd(i)
range <- foreach(i = matlist) %do% range(i)
quantile <- foreach(i = matlist) %do% quantile(i)
dist <- foreach(i = matlist) %do% dist(i) 

mean(unlist(mean))
range(unlist(mean))
mean(unlist(range))
range(unlist(range))
mean(unlist(dist))
range(unlist(dist))

i <- 1
for(w in landscapefiles){
  landsc <- (str_split(basename(w), (".csv")) %>% unlist)[1]
  tib <- tibble(landsc = landsc,
                N = length(list[[1]]$x),
                mean_dist = mean[[i]], 
                sd_dist = sd[[i]], 
                min_dist = range[[i]][1],
                max_dist = range[[i]][2])
    if (i == 1) out <- tib else out <- rbind(out, tib) 
    i <- i+1
}  

out %>% write_csv("../data/landscinfo.csv")

####------------------- plot landscapes and connectivity ---------------------####
webfile = "../data/webs/rc.csv"
landscapefile = "../data/landscapes/landscape_001.csv"

par(pty="s")
par(mfrow = c(1,1))
read_csv(landscapefile) %>% as.matrix %>% plot

coords <- read_csv(landscapefile) %>% as.matrix
N <- length(coords[,1]) ## number of patches

A <- create_adj_matrix(read_csv(webfile))
rownames(A)
xi <- seq(0.01, 1, by = 0.01)[1:length(rownames(A))]
colors <- rev(rainbow(length(rownames(A))))
kshape <- "Exponential"
spp <- 0

for(sp in rownames(A)){
  if(spp == 0) par(mfrow=c(1,length(rownames(A))))
  spp <- spp+1
  M <- get_matrix(coords = coords, ## create sp's landscape matrix 
                kernel = kshape, ## for plotting Rectanuglar (i.e. 0/1 links) needed
                xi = xi[spp])
  M[which(M<10^-3)] <- 0
  plot(coords, pch = 1, xaxt = "n", yaxt = "n", 
       xlab = "", ylab = "", main = sp)
  for(i in 1:N){
    for(j in 1:N){
      if(M[i,j]>0){
        x1=coords[i,1]
        y1=coords[i,2]
        x2=coords[j,1]
        y2=coords[j,2]
        segments(x1,y1,x2,y2,col=colors[spp])
      }
    }
  }
}


