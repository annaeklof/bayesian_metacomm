#Load libraries
require(tidyverse)
require(igraph)

#Load R functions 
source("functions.R")
source("generate_summary_function.R")

#Helper functions
my_round_any = function(x, accuracy, f=round){
  f(x/ accuracy) * accuracy
}

numextract <- function(string){ 
  temp <- gsub("[^[:digit:].]", "", string)
  as.numeric(str_extract(temp, "\\-*\\d+\\.*\\d*"))
}


infiles <- Sys.glob("../data/results/web*_basal_*.rds")
outfile <- paste0("../data/summaries/summary_sym_basal.rds")

out <- generate_summary(infiles=infiles)
saveRDS(out, outfile)

  
