#Load libraries 
require(tidyverse) #for data handling
require(cowplot) #for inlay plots
require(magick) #for inlay plots
require(RColorBrewer) #for mixing palettes
theme_set(theme_bw()) #set black and white theme

#Load R functions
source("plot_functional_forms.R") #generate .png plots for the different consumer responses to resource loss (functional forms) 
source("plot_persistence.R")
source("plot_patchloss.R")

figuredir <- "../SI/"
if(file.exists(figuredir)){
  print("Directory already exists.")
}else{
  dir.create(file.path(figuredir)) #create direcotry for figures, if it does not exist yet
} 

#Figure S1
# dpalette <- brewer.pal(9, "Blues")[c(4, 5, 6, 7, 8)]
# 
# expand.grid(piR=c(seq(0, 0.99, l=5001), seq(0.99, 1, l=15001)),
#             lR=c(0.5, 1, 2, 3, 5)) %>%
#   as_tibble %>%
#   mutate(lMR=-lR*log(1-piR)) %>%
#   filter(lMR>=0, lMR<=5) %>%
#   ggplot +
#   aes(x=piR, y=lMR, colour=factor(lR), size=factor(lR), linetype=factor(lR)) +
#   geom_path() +
#   scale_x_continuous(name=expression(paste(pi[R])), labels=abbreviate) +
#   scale_y_continuous(name=expression(paste(lambda[M[R]]))) +
#   scale_colour_manual(name=expression(paste(lambda[R])), values=dpalette) +
#   scale_size_manual(name=expression(paste(lambda[R])),
#                     values=c(0.5, 0.6, 0.5, 0.5, 0.5)) +
#   scale_linetype_manual(name=expression(paste(lambda[R])),
#                         values=c(2, 1, 2, 2, 2))
# ## ggsave("../figures/R_persist.pdf", width=3, height=2.4)
# 
# expand.grid(s=seq(0, 1, l=20002), lMC=c(0.5, 1, 2, 3, 4)) %>%
#   as_tibble %>%
#   mutate(lR=s/(s-exp(-lMC))) %>%
#   ggplot +
#   aes(x=s, y=lR, colour=factor(lMC)) +
#   geom_path(size=0.5) +
#   scale_x_continuous(name=expression(paste((1-pi[C])(1-pi[R]))),
#                      labels=abbreviate) +
#   scale_y_continuous(name=expression(paste(lambda[R])), limits=c(1, 5)) +
#   scale_colour_manual(values=dpalette, name=expression(paste(lambda[M[C]])))
# ## ggsave("../figures/C_persist.pdf", width=3, height=2.4)

#Figures S5, S7, S9, S11
#Generate .png plots for the different consumer responses to resource loss (functional forms) 
infiles <- "../data/summaries_round2/summary_largewebs_*.rds" %>% Sys.glob() #summary files
infile <- "../data/summaries_round2/summary_serengeti.rds"  #summary files

for(infile in infiles){
  dat0 <- infile %>% read_rds 
  webs <- dat0$web %>% unique
  
  for(pweb in webs){
    dat <- dat0 %>% filter(web == pweb)
  
  pas <- unique(dat$alpha)
  pbs <- unique(dat$beta)
  
  for(pa in pas){ #loop over functional forms
    for(pb in pbs){
      
      if(pa==1 && pb==1) ptitles="A"
      if(pa==5 && pb==1) ptitles="B"
      if(pa==1 && pb==5) ptitles="C"
      if(pa==5 && pb==5) ptitles="D"
      
      #Add prem=0 to scenarios
      temp <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "best-case scenario")
      temp2 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "worst-case scenario")
      temp3 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "random scenario")
      dat_plot <- dat %>% bind_rows(temp) %>% bind_rows(temp2) %>% bind_rows(temp3) %>% 
        filter(scenario != "pre-patch loss")
      
      #Define input for plotting functions
      # pweb <- unique(dat_plot$web) #web
      pscenarios <-  unique(dat_plot$scenario) #scenarios
      pspfs <- "basal" #unique(dat_plot$spf) ## which focal species for patch loss
      print(pspfs)
      pbasals <- ifelse(pweb=="Serengeti food web with 32 consumer and 129 basal species", 
                        TRUE, FALSE) ## if TRUE basal species are plotted if FALSE only consumer species

      #Metapopulation capacity ~ patches removed
      pparams <- c("pi[i]*' = 0.2,  '*xi[i]*' = 0.055'", "pi[i]*' TLB,  '*xi[i]*' = 0.055'", 
        "pi[i]*' and '*xi[i]*' TLB'", "pi[i]*' = 0.2,  '*xi[i]*' TLB'")
      
      #Create named vectors for facet labels with greek letters
      var_names <- pparams
      names(var_names) <- pparams
      
      #Plot name
      lambda <- paste0(figuredir, "lambda_", pweb %>% str_split(" ", simplify=TRUE) %>% as_tibble %>%
                         paste(.,collapse="_") %>% unname, "_", pspfs, "_", pa, "_", pb,".pdf")

      #Generate the plot
      p <- plot_patchloss(data = dat_plot, webs = pweb, paramss = pparams, scenarioss = pscenarios, spfs = pspfs,
                          a = pa, b = pb, basals = pbasals, ptitle=ptitles, labnames = var_names)

      #Save it
      ggsave(filename = lambda, width=14, height=8, plot=p)
  
      #Average landscape persistence ~ dispersal distance 
      # pparams <- unique(dat$params)[c(1:5)] #pi TLB
      vec <- which(unique(dat$params) %>% str_detect(pattern="TLB,") == TRUE)
      pparams <- unique(dat$params)[vec]

      #Plot name
      sump <- paste0(figuredir, "sump_", pweb %>% str_split(" ", simplify=TRUE) %>% as_tibble %>% 
                         paste(.,collapse="_") %>% unname, "_", pspfs, "_", pa, "_", pb,".pdf")
      
      p2 <- plot_persistence(data = dat_plot, webs = pweb, paramss = pparams, scenarioss = pscenarios, spfs = pspfs,
                          a = pa, b = pb, basals = pbasals, ptitle=ptitles)

      #Save it
      ggsave(filename = sump, width=14, height=8, plot=p2)
    } 
  }
 }
}

#Figures S13
#Generate .pdf plots for the different model food webs with the top species as focal species for patch loss 
infile <- "../data/summaries_round1/summary_largewebs_1_1.rds" #summary files
# infile <- "../data/summaries_round1/summary_serengeti.rds" #summary files

  dat0 <- infile %>% read_rds 
  webs <- dat0$web %>% unique
  
  for(pweb in webs){
    dat <- dat0 %>% filter(web == pweb)
    
      pa <- 1 
      pb <- 1 
        
        if(pweb=="Model food web with 350 consumer and 50 basal species") ptitles=paste("A", pweb) 
        if(pweb=="Model food web with 300 consumer and 100 basal species") ptitles=paste("B", pweb) 
        if(pweb=="Model food web with 250 consumer and 150 basal species") ptitles=paste("C", pweb) 
        if(pweb=="Model food web with 200 consumer and 200 basal species") ptitles=paste("D", pweb) 
        if(pweb=="Serengeti food web with 32 consumer and 129 basal species") ptitles=paste("E", pweb) 

        #Add prem=0 to scenarios
        temp <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "best-case scenario")
        temp2 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "worst-case scenario")
        temp3 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "random scenario")
        dat_plot <- dat %>% bind_rows(temp) %>% bind_rows(temp2) %>% bind_rows(temp3) %>% 
          filter(scenario != "pre-patch loss")
        
        #Define input for plotting functions
        # pweb <- unique(dat_plot$web) #web
        pscenarios <-  unique(dat_plot$scenario) #scenarios
        pspfs <- "top" #unique(dat_plot$spf) ## which focal species for patch loss
        pbasals <- ifelse(pweb=="Serengeti food web with 32 consumer and 129 basal species", TRUE, FALSE) ## if TRUE basal species are plotted if FALSE only consumer species
        
        #Metapopulation capacity ~ patches removed
        pparams <- c("pi[i]*' = 0.2,  '*xi[i]*' = 0.055'", "pi[i]*' TLB,  '*xi[i]*' = 0.055'", 
                     "pi[i]*' and '*xi[i]*' TLB'", "pi[i]*' = 0.2,  '*xi[i]*' TLB'")
        
        #Create named vectors for facet labels with greek letters
        var_names <- pparams
        names(var_names) <- pparams
        
        #Plot name
        lambda <- paste0(figuredir, "lambda_", pweb %>% str_split(" ", simplify=TRUE) %>% as_tibble %>%
                           paste(.,collapse="_") %>% unname, "_", pspfs, "_", pa, "_", pb,".pdf")

        #Generate the plot
        p <- plot_patchloss(data = dat_plot, webs = pweb, paramss = pparams, scenarioss = pscenarios, spfs = pspfs,
                            a = pa, b = pb, basals = pbasals, ptitle=ptitles, labnames = var_names)

        #Save it
        ggsave(filename = lambda, width=14, height=8, plot=p)
      } 

  
  #Figures S14
  #Generate .pdf plot for a weigthted food web 
  infile <- "../data/summaries_round2/summary_web503501.rds" #summary files
  # infile <- "../data/summaries_round1/summary_serengeti.rds" #summary files
  
  dat0 <- infile %>% read_rds 
  pweb <- dat0$web %>% unique
  
  dat <- dat0 %>% filter(weighted == TRUE)

  ptitles=paste("Weighted model food web with 350 consumer and 50 basal species") 

  #Add prem=0 to scenarios
  temp <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "best-case scenario")
  temp2 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "worst-case scenario")
  temp3 <- dat %>% filter(scenario == "pre-patch loss", alpha == pa, beta == pb) %>% mutate(scenario = "random scenario")
  dat_plot <- dat %>% bind_rows(temp) %>% bind_rows(temp2) %>% bind_rows(temp3) %>% 
    filter(scenario != "pre-patch loss")
    
  #Define input for plotting functions
  # pweb <- unique(dat_plot$web) #web
  pscenarios <-  unique(dat_plot$scenario) #scenarios
  pspfs <- "basal" #unique(dat_plot$spf) ## which focal species for patch loss
  pbasals <- FALSE ## if TRUE basal species are plotted if FALSE only consumer species
    
  #Metapopulation capacity ~ patches removed
  pparams <- c("pi[i]*' = 0.2,  '*xi[i]*' = 0.055'", "pi[i]*' TLB,  '*xi[i]*' = 0.055'", 
                 "pi[i]*' and '*xi[i]*' TLB'", "pi[i]*' = 0.2,  '*xi[i]*' TLB'")
    
  #Create named vectors for facet labels with greek letters
  var_names <- pparams
  names(var_names) <- pparams
    
  #Plot name
  lambda <- paste0(figuredir, "lambda_", ptitles %>% str_split(" ", simplify=TRUE) %>% as_tibble %>%
                      paste(.,collapse="_") %>% unname, "_", pspfs, "_", pa, "_", pb,".pdf")
    
  #Generate the plot
  p <- plot_patchloss(data = dat_plot, webs = pweb, paramss = pparams, scenarioss = pscenarios, spfs = pspfs,
                      a = pa, b = pb, basals = pbasals, ptitle=ptitles, labnames = var_names)
    
  #Save it
  ggsave(filename = lambda, width=14, height=8, plot=p)
 
  
