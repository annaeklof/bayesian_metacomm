
require(tidyverse)

landscs <- Sys.glob("../../data/landscapes/landscape*.csv")[1]
spinputs <- Sys.glob("../../data/species_input/linkremoval/web*.csv")
spfs <- "basal" ## spfs <- c("basal", "top")
outputdir <- "../../data/linkrem_results/"
alpha <- c(1, 5)
beta <- c(1, 5)

tab <- expand_grid(landsc=landscs, spinput=spinputs, spf=spfs,
                   outputdir=outputdir, alpha=alpha, beta=beta) %>%
  mutate(web=paste0("../../data/webs/",
                    str_split_fixed(basename(spinput), "_", 2)[,1], ".csv")) %>%
  dplyr::select(web, landsc, spinput, spf, outputdir, alpha, beta)

for (r in 1:nrow(tab)) {
  jobname <- sprintf("BN%02d.sh",r)
  f <- file(jobname)
  writeLines(c(
    paste("#!/bin/bash"),
    paste("#SBATCH -J", jobname),
    paste("#SBATCH -A liu-2018-29"),
    paste("#SBATCH -t 24:00:00"),
    paste("#SBATCH --mem=16000"),
    paste("#SBATCH -n 1"),
    paste("#"),
    paste("# Run single task in foreground"),
    paste("module add R/3.6.0-nsc1-gcc-7.3.0"),
    paste("cd", Sys.getenv("PWD")),
    paste("Rscript run_instance.R", paste(tab[r,], collapse=" ")),
    paste("#"),
    paste("# Script ends here")
  ), f)
  close(f)
  system(paste("sbatch", jobname))
  Sys.sleep(0.1)
}
