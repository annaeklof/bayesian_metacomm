
rm(list=ls())
require(tidyverse)

logs <- Sys.glob("../logs/webs_*.log")
#logs <- Sys.glob("../logs/soil*.log")

topspfext <- consext <- lowp <- basalspf <- topspf <- webs <- landscs <- pis <- xis <- c()
ls <- list() ## fill list with simualtions that crashed before simulating habitat loss 

for(i in logs){
  log <- read_log(i, skip=24)
  ls[[i]] <- log$X2
  if(("stopped after full landscape: TRUE" %in% ls[[i]]) && (any(str_detect(ls[[i]], "scenarios running for basal")))) basalspf[i] <- i
  if(("stopped after full landscape: TRUE" %in% ls[[i]]) && (any(str_detect(ls[[i]], "scenarios running for top")))) topspf[i] <- i
  if(("stopped after full landscape: TRUE" %in% ls[[i]]) && (any(str_detect(ls[[i]], "top species extinct: TRUE")))) topspfext[i] <- i
  if(("stopped after full landscape: TRUE" %in% ls[[i]]) && (any(str_detect(ls[[i]], "all consumers extinct: TRUE")))) consext[i] <- i
  if(("stopped after full landscape: TRUE" %in% ls[[i]]) && (any(str_detect(ls[[i]], "persistence too low: TRUE")))) lowp[i] <- i
  if("stopped after full landscape: TRUE" %in% ls[[i]]) webs[i] <- ls[[i]][1]
  if("stopped after full landscape: TRUE" %in% ls[[i]]) landscs[i] <- ls[[i]][2]
  if("stopped after full landscape: TRUE" %in% ls[[i]]) pis[i] <- ls[[i]][3]
  if("stopped after full landscape: TRUE" %in% ls[[i]]) xis[i] <- ls[[i]][4]
}

## list of simulations thata all webs: 
length(ls) ## total of 720 simulations 160 spf basal, 160 spf top; 
length(logs) ## total of 720 simulations 160 spf basal, 160 spf top; 
length(basalspf) ## 0
length(topspf) ## 84
length(topspfext) ## all of these stopped because top is focal species for patch loss (spf) and already extinct
length(consext) ## 0 
length(lowp) ## 0
  
sum(webs %>% str_detect("largeweb")) ## 39 400-species webs
sum(webs %>% str_detect("serengeti")) ## 5 serengeti

levels(as.factor(webs))
levels(as.factor(landscs))
levels(as.factor(pis))
levels(as.factor(xis)) 




