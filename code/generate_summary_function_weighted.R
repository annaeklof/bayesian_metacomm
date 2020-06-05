#Load libraries
require(tidyverse)
require(igraph)

#Load R functions 
source("functions.R")


#Helper functions
my_round_any = function(x, accuracy, f=round){
  f(x/ accuracy) * accuracy
}

numextract <- function(string){ 
  temp <- gsub("[^[:digit:].]", "", string)
  as.numeric(str_extract(temp, "\\-*\\d+\\.*\\d*"))
}


## Summarize over all patches for each simulation run and modify the generated summary. 
## Input: - infiles: simulation output files (w/ path & extension)
## Output: -out: summarized output
generate_summary <- function(infiles){
    i <- 1
    for (infile in infiles) {
      print(infile)
      sol <- infile %>% readRDS
      problems(sol)
      #stop_for_problems(sol)
      sol <- sol %>% mutate(weighted = str_detect(infile, pattern="weighted"))
      
      temp <- sol %>%
        group_by(web, landsc, spinput, seed, spf, scenario, prem, alpha, beta,
                 species, TL, OI, FT, weighted) %>%
        summarise(pi=unique(pi), xi=unique(xi), 
                  sump=sum(p), meanp=mean(p), minp=min(p), maxp=max(p),
                  lambda=mean(lambda), meand=mean(delta), mind=min(delta),
                  maxd=max(delta)) %>%
        ungroup
      if (i==1) dat <- temp else dat <- rbind(dat, temp)
      i <- i+1
    }
    
    #Modify summary
    out <- dat %>% mutate(
      topsp = dat %>% pull(species) %>% tail(n = 1),
      web = case_when(
        numextract(dat$web) %>% unique == "503501" ~ "Model food web with 350 consumer and 50 basal species (1)",
        numextract(dat$web) %>% unique == "503502" ~ "Model food web with 350 consumer and 50 basal species (2)",
        numextract(dat$web) %>% unique == "503503" ~ "Model food web with 350 consumer and 50 basal species (3)",
        numextract(dat$web) %>% unique == "503504" ~ "Model food web with 350 consumer and 50 basal species (4)",
        numextract(dat$web) %>% unique == "503505" ~ "Model food web with 350 consumer and 50 basal species (5)"),
      landsc = numextract(str_split(basename(dat$landsc), "_", simplify = TRUE)[, 2]), 
      ITL = round(dat$TL), 
      TL05 = my_round_any(dat$TL, 0.5), 
      params = case_when(
        numextract(str_split(basename(dat$spinput),"_", simplify = TRUE)[, 2]) == 1 ~ "pi[i]*' TLB,  '*xi[i]*' = 0.01'",
        numextract(str_split(basename(dat$spinput),"_", simplify = TRUE)[, 2]) == 2 ~ "pi[i]*' TLB,  '*xi[i]*' = 0.0325'",
        numextract(str_split(basename(dat$spinput),"_", simplify = TRUE)[, 2]) == 3 ~ "pi[i]*' TLB,  '*xi[i]*' = 0.055'",
        numextract(str_split(basename(dat$spinput),"_", simplify = TRUE)[, 2]) == 4 ~ "pi[i]*' TLB,  '*xi[i]*' = 0.0775'",
        numextract(str_split(basename(dat$spinput),"_", simplify = TRUE)[, 2]) == 5 ~ "pi[i]*' TLB,  '*xi[i]*' = 0.1'",
        numextract(str_split(basename(dat$spinput),"_", simplify = TRUE)[, 2]) == 6 ~ "pi[i]*' and '*xi[i]*' TLB'",
        numextract(str_split(basename(dat$spinput),"_", simplify = TRUE)[, 2]) == 7 ~ "pi[i]*' = 0.2,  '*xi[i]*' TLB'",
        numextract(str_split(basename(dat$spinput),"_", simplify = TRUE)[, 2]) == 8 ~ "pi[i]*' = 0.2,  '*xi[i]*' = 0.055'"
        ),
      scenario = case_when(
        dat$scenario == "best" ~ "best-case scenario",
        dat$scenario == "worst" ~ "worst-case scenario",
        dat$scenario == "random" ~ "random scenario",
        TRUE ~ "pre-patch loss"
      )
    )
    return(out)
  }
