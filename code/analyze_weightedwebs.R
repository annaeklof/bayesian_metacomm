
#Load libraries 
require(tidyverse)
require(cowplot)
require(magick)
require(RColorBrewer)

#Load R functions
source("plot_functional_forms.R") #Functional forms 
source("plot_patchloss.R") 

figuredir <- "../figures/"
if(file.exists(figuredir)){
  print("Directory already exists.")
}else{
  dir.create(file.path(figuredir)) #create direcotry for figures, if it does not exist yet
} 

#Generate .png plots for the different consumer responses to resource loss (functional forms) 
infiles <- "../data/summaries_round2/summary_web50350*.rds" %>% Sys.glob() #summary files
infile <- "../data/summaries_round2/summary_web503501.rds" #summary files

for(infile in infiles){
  dat_full <- infile %>% read_rds 
  # pparams <- unique(dat_full$params)
  # spinput <- unique(dat_full$spinput) 
  # nums <- str_split(spinput,pattern="_",simplify=TRUE)[,3] %>% numextract()
  # dat <- dat_full %>% filter(params == pparams) ## selected pi_i and xi_i combinations
  
  for(WEIGHTED in c(TRUE,FALSE)){
    dat <- dat_full %>% filter(weighted == WEIGHTED) 
    dat$web <- ifelse(WEIGHTED==TRUE, paste(dat$web, "Weighted"),  paste(dat$web, "Not weighted"))
    
    pas <- unique(dat$alpha)
    pbs <- unique(dat$beta)
    # pas <- 1
    # pbs <- 1
  
    for(pa in pas){ #loop over functional forms
      for(pb in pbs){
      
      if(pa==1 && pb==1) ptitles="A"
      if(pa==5 && pb==1) ptitles="B"
      if(pa==1 && pb==5) ptitles="C"
      if(pa==5 && pb==5) ptitles="D"
      
      ptitles <- ifelse(WEIGHTED==TRUE,paste(ptitles,dat$web), paste(ptitles, dat$web))
      
      # plot_functional_form(figuredir, pa, pb)
      
      #Add prem=0 to scenarios
      temp <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "best-case scenario")
      temp2 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "worst-case scenario")
      temp3 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "random scenario")
      dat_plot <- dat %>% bind_rows(temp) %>% bind_rows(temp2) %>% bind_rows(temp3) %>% 
        filter(scenario != "pre-patch loss")
      
      #Define input for plotting functions
      pweb <- unique(dat_plot$web) #web
      pscenarios <-  unique(dat_plot$scenario) #scenarios
      pspfs <- "basal" ## which focal species for patch loss
      pbasals <- FALSE ## if TRUE basal species are plotted if FALSE only consumer species
      
      #Metapopulation capacity ~ patches removed
      pparams <- unique(dat_plot$params)#[c(3,6:8)] ## selected pi_i and xi_i combinations

      #Create named vectors for facet labels with greek letters
      var_names <- pparams
      names(var_names) <- pparams
        
      #Plot name
      lambda <- paste0(figuredir, "lambda_", pweb %>% str_split(" ", simplify=TRUE) %>% as_tibble %>% 
                   paste(.,collapse="_") %>% unname, "_", pspfs, "_", pa, "_", pb,".pdf")
    
      #Generate the plot
      p <- plot_patchloss(data = dat_plot, webs = pweb, paramss = pparams, scenarioss = pscenarios, spfs = pspfs,
                          a = pa, b = pb, basals = pbasals, ptitle = ptitles, labnames = var_names)
      
      #Add functional form inset
      # g <- ggdraw(p) + draw_image(fform_png, x = 0.835, y = 1.01, hjust = 0, vjust = 1, width = 0.2, height = 0.2)
      # g
      
      #Save it
      ggsave(filename = lambda, width=14, height=8, plot=p)
    } 
    }
  }
}


#Analyze webs
require(igraph)
require(tidyverse)
require(NetIndices)
require(xtable)
source("functions.R")

#web information: 
webfiles <- Sys.glob("../data/webs/web50350*.csv")[1:5] ## list of web files

i <- 1
for(webfile in webfiles){
  webname <- str_split(basename(webfile), (".csv"), simplify=TRUE)[1] %>% 
    str_split(.," ",simplify = TRUE) %>% paste(.,collapse="_")
  Amat <- read_csv(webfile, col_types=cols()) %>% create_adj_matrix ## create adjacency matrix from links
  g <- graph_from_adjacency_matrix(Amat, mode = "directed") ## generate graph
  if(!is_dag(g)) print("Web is not acyclic.") 
  if(is_dag(g)){
    S = (length(Amat[1,]))
    Sb =  length(which(rowSums(Amat) == 0)) 
    Sc =  length(which(rowSums(Amat) != 0)) 
    if(Sb+Sc!=S) print("Species numbers do not add up.")
    tib <- tibble(web = webname, 
                  S = S,
                  Sb = Sb,  
                  Sc = Sc,
                  connectance = calc_NetInd(Amat)$connectance,  
                  generality = calc_NetInd(Amat)$generality, 
                  sd.generality = calc_NetInd(Amat)$sd.generality, 
                  vulnerability  = calc_NetInd(Amat)$vulnerability,
                  sd.vulnerability = calc_NetInd(Amat)$sd.vulnerability, 
                  modularity = calc_NetInd(Amat)$modularity,  
                  TL.max = calc_NetInd(Amat)$TL.max, 
                  TL.mean = calc_NetInd(Amat)$TL.mean, 
                  TL.sd = calc_NetInd(Amat)$TL.sd, 
                  OI.max = calc_NetInd(Amat)$OI.max, 
                  OI.mean = calc_NetInd(Amat)$OI.mean, 
                  OI.sd = calc_NetInd(Amat)$OI.sd
    )
    
    if (i == 1) out <- tib else out <- rbind(out, tib) 
    i <- i+1
  }
}  

out %>% write_csv("../data/notweighted_webinfo.csv") ## write web info to file 

#species distribution across trophic levels: 
dats <- Sys.glob( "../data/summaries_round2/summary_web50350*.rds")

for(m in 2:5){
  dat <- dats[m] %>% readRDS

  for(i in unique(dat$ITL))
    tmp <- dat %>% group_by(ITL) %>% summarize(no.of.species = species %>% unique %>% length)
  if(m==1) 
    tib <- tmp
  else 
    tib <- right_join(tib, tmp, by="ITL")
}

oldnames <- colnames(tib)[-1]
newnames <- paste0("web50350",1:5)
tib <- tib %>% rename_at(vars(oldnames), ~ newnames)

write_csv(tib, "../data/no.of.species.TL.weightedwebs.csv") ## save
# xtable(tib)


## analyze input parameters xi and pi: 
range_pi <- list()
range_xi <- list()

for(webs in unique(dat$web))
{
  range_pi[[webs]] <- dat %>%
    filter(web==webs) %>% pull(pi) %>% unique %>% range %>% round(.,4)
  range_xi[[webs]] <- dat %>%
    filter(web==webs) %>% pull(xi) %>% unique %>% range %>% round(.,4)
}

t_pi <- range_pi %>% as_tibble %>% write_csv("../data/range_pi.csv") ## make into tibble and save
t_xi <- range_xi %>% as_tibble %>% write_csv("../data/range_xi.csv") ## make into tibble and save

