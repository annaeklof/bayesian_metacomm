
require(tidyverse)

landscs <- Sys.glob("../data/landscapes/landscape*.csv")[1:5]
spinputs <- Sys.glob("../data/species_input/web*.csv")
spfs <- c("basal", "top")
outputdir <- "../data/results/"
alpha <- 1
beta <- 1

tab <- expand.grid(landsc=landscs, spinput=spinputs, spf=spfs,
                   outputdir=outputdir, alpha=alpha, beta=beta,
                   stringsAsFactors=FALSE) %>%
    as_tibble %>%
    mutate(web=paste0("../data/webs/",
                      str_split_fixed(basename(spinput), "_", 2)[,1], ".csv")) %>%
    select(web, landsc, spinput, spf, outputdir, alpha, beta)

tab <- tab %>% filter(web=="../data/webs/web5351.csv",
                      spf=="top",
                      spinput %in% c("../data/species_input/web5351_001.csv",
                                     "../data/species_input/web5351_006.csv",
                                     "../data/species_input/web5351_007.csv"))

for (r in 1:nrow(tab)) {
    jobname <- sprintf("SBN_%03d.sh", r)
    f <- file(jobname)
    writeLines(c(
        paste("#!/bin/bash"),
        paste("#SBATCH -J", jobname),
        paste("#SBATCH -A liu-2018-29"),
        paste("#SBATCH -t 05:00:00"),
        paste("#SBATCH --mem=16000"),
        paste("#SBATCH -n 1"),
        paste("#"),
        paste("# Run single task in foreground"),
        paste("module add R/3.5.1-nsc1-gcc-2018a"),
        paste("cd", Sys.getenv("PWD")),
        paste("Rscript run_instance.R", paste(tab[r,], collapse=" ")),
        paste("#"),
        paste("# Script ends here")
    ), f)
    close(f)
    system(paste("sbatch", jobname))
    Sys.sleep(0.1)
}

