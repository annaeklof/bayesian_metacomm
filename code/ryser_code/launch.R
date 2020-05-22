
require(tidyverse)

tab <- expand_grid(alpha=c(1, 5), beta=c(1, 5), pB=c(0.2, 0.3, 0.4))
tab <- tab %>% mutate(outfile=paste0(sprintf("out_%03d.csv", 1:nrow(tab))))

for (r in 1:nrow(tab)) {
    jobname <- sprintf("SBN_%03d.sh", r)
    f <- file(jobname)
    writeLines(c(
        paste("#!/bin/bash"),
        paste("#SBATCH -J", jobname),
        paste("#SBATCH -A liu-2017-00089-2"), ## liu-2018-29
        paste("#SBATCH -t 10:00:00"),
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
