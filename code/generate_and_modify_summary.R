## Generate summary file: 
## Summarizes over all patches for each simulation run. 
## Input: 
## - resultfiles 
## Output: 
## - outfile 

require(igraph)
require(tidyverse)
source("functions.R")
source("generate_summary_function.R")

resultfiles <- Sys.glob("../data/results/serengeti_009_*.rds")
outfile <- "../data/summaries/summary_serengeti_009_new.rds"

out <- function_to_generate_and_modify_summary(resultfiles)

saveRDS(out, outfile)

