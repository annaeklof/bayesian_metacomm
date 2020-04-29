
rm(list=ls())

## Modify summary for plotting: 
require(tidyverse)
require(igraph)

source("functions.R")

## Helper function: 
my_round_any = function(x, accuracy, f=round){
  f(x/ accuracy) * accuracy
}



data<- data %>% mutate(params = recode(params,
                                       "xi and pi TLB" = "pi[i] and xi[i] TLB",
                                       "xi TLB, pi = 0.2" = "pi[i] = 0.2 xi[i] TLB", 
                                       "xi = 0.055, pi = 0.2" = "pi[i] = 0.2 xi[i] = 0.05", 
                                       "xi = 0.055, pi = TLB" = "pi[i] TLB xi[i] = 0.05"
))


modify_summary <- function(infiles){
  dat <- infiles %>% map_df(readRDS) 
  out <- dat %>% mutate(topsp = case_when(
    web == ("../data/webs/serengeti.csv") ~ read_csv("../data/webs/serengeti.csv") %>% 
      create_adj_matrix %>% rownames %>% tail(n = 1),
    web %>% str_detect("largeweb") ~ read_csv("../data/webs/largeweb1003001.csv") %>% 
      create_adj_matrix %>% rownames %>% tail(n = 1),
    TRUE ~ read_csv("../data/webs/web10301.csv") %>% create_adj_matrix %>% rownames %>% tail(n = 1)),
    web = recode(web,  
      "../data/webs/web30101.csv" = "Model food web with 30 consumer and 10 basal species",
      "../data/webs/web3551.csv"  = "Model food web with 35 consumer and 5 basal species",
      "../data/webs/web25151.csv"  = "Model food web with 25 consumer and 15 basal species",
      "../data/webs/web20101.csv" = "Model food web with 20 consumer and 20 basal species",
      "../data/webs/web25151.csv" = "Model food web with 300 consumer and 100 basal species",
      "../data/webs/largeweb350501.csv" = "Model food web with 350 consumer and 50 basal species",
      "../data/webs/largeweb2501501.csv" = "Model food web with 250 consumer and 150 basal species",
      "../data/webs/largeweb2002001.csv" = "Model food web with 200 consumer and 200 basal species",
      "../data/webs/serengeti.csv" = "Serengeti food web with 32 consumer and 129 basal species"
    ), 
    landsc = numextract(str_split(basename(landsc), "_", simplify = TRUE)[, 2]), 
    ITL = round(TL), 
    TL05 = my_round_any(TL, 0.5), 
    params = recode(
      numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]) == 1 ~ "pi_i TLB; xi_i = 0.01",
      numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]) == 2 ~ "pi_i TLB; xi_i = 0.0325",
      numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]) == 3 ~ "pi_i TLB; xi_i = 0.055",
      numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]) == 4 ~ "pi_i TLB; xi_i = 0.0775",
      numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]) == 5 ~ "pi_i TLB; xi_i = 0.1",
      numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]) == 6 ~ "pi_i and xi_i TLB",
      numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]) == 7 ~ "pi_i = 0.2; xi_i TLB",
      numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]) == 8 ~ "pi_i = 0.2; xi_i = 0.055",
      numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]) == 9 ~ "pi_i and xi_i SGB",
      TRUE ~ as.character(numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]))
    ), 
    pi0 = case_when(
      numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]) == 1 ~ "TLB",
      numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]) == 2 ~ "TLB",
      numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]) == 3 ~ "TLB",
      numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]) == 4 ~ "TLB",
      numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]) == 5 ~ "TLB",
      numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]) == 6 ~ "TLB",
      numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]) == 7 ~ "0.2",
      numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]) == 8 ~ "0.2",
      numextract(str_split(basename(spinput),"_", simplify = TRUE)[, 2]) == 9 ~ "SGB",
    ), 
    scenario = case_when(
      scenario == "best" ~ "best-case scenario",
      scenario == "worst" ~ "worst-case scenario",
      scenario == "random" ~ "random scenario",
      TRUE ~ "pre-patch loss"
    )
  )
  return(out)
}

library(tidyverse)
dat <- "../data/summaries/summary_serengeti_linear.rds" %>% read_rds
dat2 <- "../data/summaries/summary_serengeti_non-linear.rds" %>% read_rds

dat21 <- "../data/summaries/modified_summary_serengeti_groups.rds" %>% read_rds

dat3 <- bind_rows(dat, dat2, dat21)
saveRDS(dat3, "../data/summaries/modified_summary_serengeti_newnew.rds" )

dat4 <- generate_summary("../data/summaries/modified_summary_serengeti_new.rds")



dat2 <- dat %>% mutate(params = recode(params,
                                       "xi and pi TLB" = "pi[i] and xi[i] TLB",
                                       "xi TLB, pi = 0.2" = "pi[i] = 0.2 xi[i] TLB", 
                                       "xi = 0.055, pi = 0.2" = "pi[i] = 0.2 xi[i] = 0.05", 
                                       "xi = 0.055, pi = TLB" = "pi[i] TLB xi[i] = 0.05",
                                       "xi = 0.0325, pi = TLB" = "pi[i] TLB xi[i] = 0.0325",
                                       "xi = 0.0775, pi = TLB" = "pi[i] TLB xi[i] = 0.0775",
                                       "xi = 0.01, pi = TLB" = "pi[i] TLB xi[i] = 0.01",
                                       "xi = 0.1, pi = TLB" = "pi[i] TLB xi[i] = 0.1", 
                                       "pi_i and xi_i SGB" = "pi[i] and xi[i] SGB", 
))


# saveRDS(dat2,"../data/summaries/summary_serengeti_groups.rds")

## Modifiy summaries 40-species webs: 
infiles = "../data/summaries/summary_webs_linear.rds" %>% Sys.glob
outfile = "../data/summaries/modified_summary_webs_testing.rds" 
modify_summary(infiles) %>% saveRDS(outfile)

## Modifiy summaries 400-species webs (seperately for different fforms): 
infiles = "../data/summaries/summary_largewebs_*.rds" %>% Sys.glob
for(infile in infiles){
  outfile = paste0("../data/summaries/modified_",basename(infile))
  modify_summary(infile) %>% saveRDS(outfile)
}

## Modifiy summaries Serengeti web: 
infiles = "../data/summaries/summary_serengeti_*.rds" %>% Sys.glob
outfile = "../data/summaries/modified_summary_serengeti.rds" 
out <- modify_summary(infiles) 
groups <- read_csv("../data/serengeti/serengeti_groups.csv") ## get spatial groups
out <- out %>% left_join(groups, by="species") ## add them to summary 
out %>% saveRDS(outfile)

