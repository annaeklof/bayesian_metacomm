
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

infile <- infiles[1]


for(infile in infiles){
  dat_full <- infile %>% read_rds 
  dat_full <- dat_full %>% filter(params == unique(dat_full$params)[c(3,6:8)]) ## selected pi_i and xi_i combinations
  
  for(w in c(TRUE,FALSE)){
    dat <- dat_full %>% filter(weighted == w)
  
    pas <- unique(dat$alpha)
    pbs <- unique(dat$beta)
  
    for(pa in pas){ #loop over functional forms
      for(pb in pbs){
      
      if(pa==1 && pb==1) ptitles="A"
      if(pa==5 && pb==1) ptitles="B"
      if(pa==1 && pb==5) ptitles="C"
      if(pa==5 && pb==5) ptitles="D"
      
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
      # ggsave(filename = lambda, width=14, height=8, plot=p)
    } 
    }
  }
}

