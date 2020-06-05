#Load libraries
require(tidyverse)
require(igraph)

#Load R functions 
source("functions.R")
source("generate_summary_function_weighted.R")

#Helper functions
my_round_any = function(x, accuracy, f=round){
  f(x/ accuracy) * accuracy
}

numextract <- function(string){ 
  temp <- gsub("[^[:digit:].]", "", string)
  as.numeric(str_extract(temp, "\\-*\\d+\\.*\\d*"))
}



for(i in 1:5){
	infiles <- Sys.glob(paste0("/data/idiv_brose/johanna/bayesian_metacomm_round2/data/results/web50350",i,"_*.rds"))
	# infiles <- Sys.glob(paste0("../data/results/web50350",i,"_*.rds"))
	outfile <- paste0("/data/idiv_brose/johanna/bayesian_metacomm_round2/data/summaries/summary_web50350",i,".rds")

	out <- generate_summary(infiles=infiles)
	saveRDS(out, outfile)
}
  
