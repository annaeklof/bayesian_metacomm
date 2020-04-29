
#Load libraries 
require(tidyverse)
require(cowplot)
require(magick)
require(RColorBrewer)


#Load R functions
source("plot_functional_forms.R") #generate .png plots for the different consumer responses to resource loss (functional forms) 
source("plot_patchloss.R") #Figure 2
source("plot_patchloss_groups.R") #Figure 3

figuredir <- "../figures/"

#Figure 2: 
infiles <- "../data/summaries/modified_summary_largewebs_*.rds" %>% Sys.glob() #summary files

for(infile in infiles){
  dat <- infile %>% readRDS %>% filter(web == "Model food web with 300 consumer and 100 basal species") 
  pas <- unique(dat$alpha)
  pbs <- unique(dat$beta)
  
  for(pa in pas){ #loop over functional forms
    for(pb in pbs){
      
      if(pa==5 && pb==1) ptitles="A"
      if(pa==5 && pb==5) ptitles="B"
      if(pa==1 && pb==1) ptitles="C"
      if(pa==1 && pb==5) ptitles="D"
      
      fform_png <-  paste0("/homes/jh57masa/github/bayesian_metacomm/figures/fform_",pa,"_",pb,".png")
      
      #Add prem=0 to scenarios
      temp <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "best-case scenario")
      temp2 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "worst-case scenario")
      temp3 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "random scenario")
      dat_plot <- dat %>% bind_rows(temp) %>% bind_rows(temp2) %>% bind_rows(temp3)
      
      #Define input for plotting functions
      (pweb <- unique(dat_plot$web)) ## which web
      pscenarios <- c("best-case scenario", "worst-case scenario", "random scenario") ## which scenarios
      pspfs <- "basal" ## which focal species for patch loss
      pbasals <- FALSE ## if TRUE basal s pecies are plotted if FALSE only consumers
      
      #Metapopulation capacity ~ patches removed
      pparams <- unique(dat$params)[c(3,6:8)] ## selected pi_i and xi_i combinations
      
      var_names <- c("xi = 0.055, pi TLB"  = "pi[i]*' TLB,  '*xi[i]*' = 0.055'", 
                     "xi and pi TLB" = "pi[i]*' and '*xi[i]*' TLB'",
                     "xi TLB, pi = 0.2" = "pi[i]*' = 0.2,  '*xi[i]*' TLB'" , 
                     "xi = 0.055, pi = 0.2" = "pi[i]*' = 0.2,  '*xi[i]*' = 0.055'" )
      
      #Plot name
      lambda <- paste0(figuredir, "lambda_", pweb %>% str_split(" ", simplify=TRUE) %>% as_tibble %>% 
                         paste(.,collapse="_") %>% unname, "_", pspfs, "_", pa, "_", pb,".pdf")
      
      #Generate the plot
      p <- plot_patchloss(data = dat_plot, webs = pweb, paramss = pparams, scenarioss = pscenarios, spfs = pspfs,
                          a = pa, b = pb, basals = pbasals, ptitle=ptitles, labnames = var_names)
      
      #Add functional form 
      g <- ggdraw(p) + draw_image(fform_png, x = 0.835, y = 1.01, hjust = 0, vjust = 1, width = 0.2, height = 0.2)
      
      #Save it
      ggsave(filename = lambda, width=12, height=8, plot=g)
    } 
  }
}


#Figure 3: 
infile <- "../data/summaries/modified_summary_serengeti_groups.rds" #summary file
dat <- infile %>% readRDS

pas <- unique(dat$alpha) 
pbs <- unique(dat$beta)

for(pa in pas){ #loop over functional forms
  for(pb in pbs){
    
    if(pa==5 && pb==1) ptitles="A"
    if(pa==5 && pb==5) ptitles="B"
    if(pa==1 && pb==1) ptitles="C"
    if(pa==1 && pb==5) ptitles="D"
    
    fform_png <- paste0("/homes/jh57masa/github/bayesian_metacomm/figures/fform_",pa,"_",pb,".png")
    
    temp <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "best-case scenario")
    temp2 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "worst-case scenario")
    temp3 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "random scenario")
    dat_plot <- dat %>% bind_rows(temp) %>% bind_rows(temp2) %>% bind_rows(temp3)
    
    ## Define input
    (pweb <- unique(dat_plot$web))
    pscenarios <- c("best-case scenario", "worst-case scenario", "random scenario")
    pspfs <- "basal"
    pbasals <- TRUE
    
    #Metapopulation capacity ~ patches removed
    pparams <- unique(dat_plot$params) ## pi and xi

    #Plot name  
    lambda <- paste0(figuredir, "lambda_", pweb %>% str_split(" ", simplify=TRUE) %>% as_tibble %>% 
                       paste(.,collapse="_") %>% unname, "_", pspfs, "_", pa, "_", pb,"_SGB.pdf")
    
    #Generate the plot
    p <- plot_patchloss_groups(data = dat_plot, webs = pweb, paramss = pparams, scenarioss = pscenarios, spfs = pspfs,
                             a = pa, b = pb, basals = pbasals, ptitle = ptitles)
    
    #Add functional form
    g <- ggdraw(p) + draw_image(fform_png, x = 0.795, y = 1.01, hjust = 0, vjust = 1, width = 0.2, height = 0.2)

    #Save it
    ggsave(filename = lambda, width=8, height=7, plot=g) 
  }
}

