## Generate summary file: 
## Summarizes over all patches for each simulation run. 
## Input: 
## - resultfiles 
## Output: 
## - outfile 

require(tidyverse)
source("functions.R")

resultfiles <- Sys.glob("../data/results/*.rds")
outfile <- "../data/summaries/summary_2.rds"

i <- 1
for (resultfile in resultfiles) {
  sol <- readRDS(resultfile)
  stop_for_problems(sol)
  temp <- sol %>%
    group_by(web, landsc, spinput, seed, spf, scenario, prem, alpha, beta,
             species, TL, OI, FT) %>%
    summarise(pi=unique(pi), xi=unique(xi), 
              sump=sum(p), meanp=mean(p), minp=min(p), maxp=max(p),
              lambda=mean(lambda), meand=mean(delta), mind=min(delta),
              maxd=max(delta)) %>%
    ungroup
  if (i==1) out <- temp else out <- rbind(out, temp)
  i <- i+1
}

saveRDS(out, outfile)

