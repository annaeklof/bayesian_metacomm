
rm(list=ls())
require(tidyverse)

dat <- "../data/summaries/modified_summary.rds" %>% readRDS ## import the data

dat <- dat %>% mutate(extinct=case_when(
  lambda>=1 ~ 0,
  TRUE ~ 1))    


listlist <- list()
for(scens in unique(dat$scenario))
  for(webs in unique(dat$web))
    listlist[[scens]][[webs]] <- dat %>% filter(web == webs, scenario == scens) %>% pull(extinct) %>% sum



listlist$`best-case scenario`/listlist$`worst-case scenario`
listlist$`best-case scenario`/listlist$`random scenario`
listlist$`worst-case scenario`/listlist$`random scenario`


boxplot(dat$extinct ~ dat$scenario)


all <- best <- worst <- random <- list()
i <- 1
for(webs in unique(dat$web)){
  all[[webs]] <- dat %>% filter(web == webs, scenario == "all") %>% pull(extinct) %>% sum
  best[[webs]] <- dat %>% filter(web == webs, scenario == "best-case scenario") %>% pull(extinct) %>% sum
  worst[[webs]] <- dat %>% filter(web == webs, scenario == "worst-case scenario") %>% pull(extinct) %>% sum
  random[[webs]] <- dat %>% filter(web == webs, scenario == "random scenario") %>% pull(extinct) %>% sum
}



