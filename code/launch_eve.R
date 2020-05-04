
require(tidyverse)

landscs <- Sys.glob("../data/landscapes/gridlike*.csv")
spinputs <- Sys.glob("../data/species_input/web50350*.csv")[c(3, 6, 7, 8)]
spfs <- "basal" ## spfs <- c("basal", "top")
outputdir <- "../data/results/"
alpha <- c(1, 5)
beta <- c(5, 1)  

#note: expand_grid does not work on the cluster (probably because of an older R version)
tab <- expand.grid(landsc=landscs, spinput=spinputs, spf=spfs,
                   outputdir=outputdir, alpha=alpha, beta=beta, stringsAsFactors = FALSE) %>% 
    mutate(web=paste0("../data/webs/",
                      str_split_fixed(basename(spinput), "_", 2)[,1], ".csv")) %>%
    select(web, landsc, spinput, spf, outputdir, alpha, beta) ##%>%
    ##filter(alpha!=beta) ## only non-linear fforms of prey loss

for (r in 1:nrow(tab)) {
  jobname <- sprintf("webs_%03d.sh",r)
  f <- file(jobname)
  writeLines(c(
    paste("#!/bin/bash"),
    paste("#$ -S /bin/bash"),
    #paste("#$ -M johanna.haeussler@idiv.de"),
    paste("#$ -N", jobname),
    paste("#$ -o /work/$USER/$JOB_NAME-$JOB_ID.log"),
    paste("#$ -wd /work/$USER"),
    paste("#$ -l h_rt=36:00:00"),
    paste("#$ -l h_vmem=8G"),
    paste("#$ -j y"),
    paste("#$ -cwd"),
    paste("#$ -binding linear:1"),        	
    paste("#"),
    paste("# Run single task in foreground"),
    paste("module load R/3.5.1-1"),
    paste("Rscript run_instance.R", paste(tab[r,], collapse=" ")),
    paste("#"),
    paste("# Script ends here")
  ), f)
  close(f)
  system(paste("qsub", jobname))
  Sys.sleep(0.1)
}
