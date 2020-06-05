
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

for(infile in infiles[1:2]){
  dat_full <- infile %>% read_rds 
  pparams <- unique(dat_full$params)[c(3,6:8)]
  # spinput <- unique(dat_full$spinput)[c(3,6:8)] 
  # nums <- str_split(spinput,pattern="_",simplify=TRUE)[,3] %>% numextract()
  dat <- dat_full %>% filter(params == pparams) ## selected pi_i and xi_i combinations
  
  for(WEIGHTED in c(TRUE,FALSE)){
    dat <- dat_full %>% filter(weighted == WEIGHTED) 
    dat$web <- ifelse(WEIGHTED==TRUE, paste(dat$web, "Weighted"),  paste(dat$web, "Not weighted"))
    
    # pas <- unique(dat$alpha)
    # pbs <- unique(dat$beta)
    pas <- 1
    pbs <- 1
  
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

