
require(tidyverse)

for(i in 2:5){
  spinputs <- Sys.glob(paste0("../data/species_input/web50350",i,"*.csv"))[c(3, 6, 7, 8)]
  landscs <- Sys.glob("../data/landscapes/landscape_*.csv")
  spfs <- "basal" ## spfs <- c("basal", "top")
  outputdir <- "../data/results/"
  alpha <- c(1, 5)
  beta <- c(5, 1)  
  
  #note: expand_grid does not work on the cluster (probably because of an older R version)
  tab <- expand.grid(landsc=landscs, spinput=spinputs, spf=spfs,
                     outputdir=outputdir, alpha=alpha, beta=beta, stringsAsFactors = FALSE) %>% 
      mutate(web=paste0("../data/webs/",
                        str_split_fixed(basename(spinput), "_", 2)[,1], ".csv")) %>%
    dplyr::select(web, landsc, spinput, spf, outputdir, alpha, beta) 
  
  for (r in 1:nrow(tab)) {
    jobname <- sprintf(paste0("weighted_", i, "_%03d.sh"),r)
    f <- file(jobname)
    writeLines(c(
      paste("#!/bin/bash"),
      paste("#$ -S /bin/bash"),
      #paste("#$ -M johanna.haeussler@idiv.de"),
      paste("#$ -N", jobname),
      paste("#$ -o /work/$USER/$JOB_NAME-$JOB_ID.log"),
      paste("#$ -wd /work/$USER"),
      paste("#$ -l h_rt=240:00:00"),
      paste("#$ -l h_vmem=8G"),
      paste("#$ -j y"),
      paste("#$ -cwd"),
      paste("#$ -binding linear:1"),        	
      paste("#"),
      paste("# Run single task in foreground"),
      paste("module load R/3.5.1-1"),
      paste("Rscript run_instance_weighted.R", paste(tab[r,], collapse=" ")),
      paste("#"),
      paste("# Script ends here")
    ), f)
    close(f)
    system(paste("qsub", jobname))
    Sys.sleep(0.1)
  }
  
  for (r in 1:nrow(tab)) {
    jobname <- sprintf(paste0("unweighted_", i, "_%03d.sh"),r)
    f <- file(jobname)
    writeLines(c(
      paste("#!/bin/bash"),
      paste("#$ -S /bin/bash"),
      #paste("#$ -M johanna.haeussler@idiv.de"),
      paste("#$ -N", jobname),
      paste("#$ -o /work/$USER/$JOB_NAME-$JOB_ID.log"),
      paste("#$ -wd /work/$USER"),
      paste("#$ -l h_rt=240:00:00"),
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
}


