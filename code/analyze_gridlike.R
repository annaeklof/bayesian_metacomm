
require(tidyverse)
require(RColorBrewer)

getPalette  <- colorRampPalette(brewer.pal(11,"RdYlBu"))
path <- "../data/summaries/summaries_gridlike/summary_gridlike_" ## summary files
psizes <- c(0, exp(seq(log(0.01), log(0.5), l=19))) ## perturbation sizes

gl <- tibble() ## gridlike, extracted & summarized
for (i in 1:20) {
  gl <- readRDS(paste0(path, i, ".rds")) %>% ## read ith gridlike summary file
    select(scenario, spinput, prem, alpha, beta, TL, lambda) %>%
    filter(scenario!="pre-patch loss") %>% ## drop pre-patch loss scenarios
    mutate(TL=as.integer(round(TL))) %>% ## integer trophic levels
    filter(TL>1) %>% ## only trophic levels 2 and up
    group_by(scenario, spinput, alpha, beta, prem, TL) %>% ## for each combo:
    summarise(mlambda=mean(lambda)) %>% ## get mean metapop capacity
    ungroup %>%
    filter(mlambda<=1) %>% ## keep only those mlambdas below the threshold
    group_by(scenario, spinput, alpha, beta, TL) %>% ## for each param combo & TL:
    filter(mlambda==max(mlambda)) %>% ## find maximum mlambda
    ungroup %>%
    group_by(scenario, spinput, alpha, beta, TL) %>% ## for each param combo & TL:
    filter(prem==min(prem)) %>% ## minimum no of patches to remove to cause ext.
    ungroup %>%
    mutate(psize=psizes[i]) %>% ## add size of perturbation of landscape from grid
    arrange(psize, scenario, spinput, alpha, beta, TL) %>% ## reorder rows
    select(psize, scenario, spinput, alpha, beta, TL, prem) %>% ## mlambda dropped!
    bind_rows(gl) ## add to extracted data table
  cat(paste0(path, i, ".rds processed","\n"))
}
gl <- gl %>% ## remove bloat from species input file names
  mutate(spinput=gsub("../data/species_input/web50350_00", "", spinput)) %>%
  mutate(spinput=gsub(".csv", "", spinput)) %>% ## remove everything except number
  mutate(spinput=case_when( ## mutate numbers into parameter descriptions
    spinput=="3" ~ "pi[i]~TLB~~xi[i]==0.055",
    spinput=="6" ~ "pi[i]~and~xi[i]~TLB",
    spinput=="7" ~ "pi[i]==0.2~~xi[i]~TLB",
    spinput=="8" ~ "pi[i]==0.2~~xi[i]==0.055",
    TRUE ~ spinput)
  )

gl %>%
  mutate(TL=factor(TL)) %>%
  filter(alpha==1, beta==1) %>%
  ggplot +
  aes(x=psize, y=prem, colour=TL) +
  geom_line(linetype="dashed") +
  geom_point() +
  facet_grid(scenario~spinput,
             labeller=labeller(.rows=label_value, .cols=label_parsed)) +
  scale_colour_manual(values=colorRampPalette(brewer.pal(11, "RdYlBu"))(7),
                      name="Trophic levels") +
  scale_x_continuous(name="Perturbation size to perfect grid structure") +
  scale_y_continuous(name="Patches to remove for eliminating trophic level") +
  theme_bw()
##ggsave("../figures/gridlike.pdf", width=8, height=5.33)
