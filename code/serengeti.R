
rm(list=ls())
require(tidyverse)
require(igraph)
source("functions.R")

#Generate additional group based species input file for the serengeti web
serengeti_files <- Sys.glob(paste0(outputdir, "serengeti_*.csv")) 
groups <- read_csv("../data/serengeti/serengeti_groups.csv")

for(serengeti_file in serengeti_files)
  serengeti_file %>% read_csv %>% left_join(groups) %>% #add group as additional column to all serengeti species input files 
  write_csv(serengeti_file)

dat <- read_csv(paste0(outputdir, "serengeti_006.csv")) ## 'pi' and 'xi' TLB
for(sp in 1:length(dat$species)){
  dat$pi[sp] <- 0.2*(15-dat$group[sp])/mean(dat$group) #group based 'baseline extinction probability 'pi'
  dat$xi[sp] <- 0.055*(15-dat$group[sp])/mean(dat$group) #group based 'diserpersal distance 'xi'
}
dat %>% write_csv(paste0(outputdir, sprintf("serengeti_%03d.csv", n)))

