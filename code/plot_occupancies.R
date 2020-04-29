    
require(RColorBrewer)
require(tidyverse)
source("functions.R")
   
## Define theme and colour palette:
theme_set(theme_bw()) ## set black and white theme
getPalette  <- colorRampPalette(brewer.pal(11,"RdYlBu"))
cpalette <- getPalette(11)[c(10, 2)]

solfiles <- Sys.glob("../data/results/largeweb100300*_1_1.rds")

for(q in solfiles){
    name <- basename(q) %>% str_split(".rds",simplify=TRUE) %>% max
    sol <- readRDS(q)
    sol <- bind_rows(sol[1,], sol) ## add fake point with p=0 to scale
    sol$p[1] <- 0 ## alpha such that white=zero

    sol <- sol %>% filter(species %in% c(min(species), max(species)))
    
    sol_plot <- sol %>% 
        rename(`patches removed`=prem) %>%
        mutate(species=case_when(
            species==min(species) ~ "a basal species",
            TRUE ~ "a top predator")) %>%  
        mutate(scenario=case_when(
            scenario=="best" ~ "A) best-case scenario",
            scenario=="worst" ~ "B) worst-case scenario",
            scenario=="random" ~ "C) random scenario",
            TRUE ~ "all")) %>% 
        mutate(params=case_when(
            numextract(str_split(basename(spinput),
                             "_", simplify = TRUE)[,2])==1 ~ "xi = 0.01, pi TLB",
            numextract(str_split(basename(spinput),
                             "_", simplify = TRUE)[,2])==2 ~ "xi = 0.0325, pi TLB",
            numextract(str_split(basename(spinput),
                             "_", simplify = TRUE)[,2])==3 ~ "xi = 0.055, pi TLB",
            numextract(str_split(basename(spinput),
                             "_", simplify = TRUE)[,2])==4 ~ "xi = 0.0775, pi TLB",
            numextract(str_split(basename(spinput),
                             "_", simplify = TRUE)[,2])==5 ~ "xi = 0.1, pi TLB",
            numextract(str_split(basename(spinput),
                             "_", simplify = TRUE)[,2])==6 ~ "xi and pi TLB",
            numextract(str_split(basename(spinput),
                             "_", simplify = TRUE)[,2])==7 ~ "xi TLB, pi = 0.2",
            numextract(str_split(basename(spinput),
                             "_", simplify = TRUE)[,2])==8 ~ "xi = 0.055, pi = 0.2",
            TRUE ~ as.character(numextract(str_split(basename(spinput),
                                                 "_", simplify = TRUE)[,2]))))
    
    for(i in c("A) best-case scenario", "B) worst-case scenario", "C) random scenario")){
        sol_plot %>% 
            filter(scenario %in% c("all", i), `patches removed`%in% c(0,100,200,290)) %>%
            mutate(p=ifelse(p<1e-10,0,p)) %>%
            mutate(p=factor(p)) %>%
            ggplot +
            aes(x=x, y=y, alpha=factor(p), colour=species, fill=species) +
            geom_point(shape=21) +
            geom_point(alpha=0.3, shape=1, colour="black") +
            facet_grid(factor(species, levels=c("a basal species", "a top predator"))~`patches removed`, 
                       labeller=labeller(.rows = label_value, .cols = label_both)) +
            scale_fill_manual(values=cpalette, guide=FALSE) +
            scale_colour_manual(values=cpalette, guide=FALSE) +
            ## use discrete here as p was converted into a factor
            scale_alpha_discrete(range=c(0, 1), guide=FALSE) +
            scale_x_continuous(name="", label=abbreviate, limits=c(0, 1)) +
            scale_y_continuous(name="", label=abbreviate, limits=c(0, 1)) +
            labs(title = i) + #subtitle = sol_plot$params %>% unique) +
            theme(plot.title = element_text(size=15), 
                  #plot.subtitle = element_text(size=15), 
                  text = element_text(size=15))
        ggsave(gsub(" ", "_", paste0("../figures/",name,"_",i,".pdf"),fixed=TRUE), width=8.5, height=4.5)
    }
}
