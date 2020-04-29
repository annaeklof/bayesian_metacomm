
require(tidyverse)

jobname <- sprintf("summarize.sh")
f <- file(jobname)
writeLines(c(
    paste("#!/bin/bash"),
    paste("#$ -S /bin/bash"),
    #paste("#$ -M johanna.haeussler@idiv.de"),
    paste("#$ -N", jobname),
    paste("#$ -o /work/$USER/$JOB_NAME-$JOB_ID.log"),
    paste("#$ -wd /work/$USER"),
    paste("#$ -l h_rt=12:00:00"),
    paste("#$ -l h_vmem=8G"),
    paste("#$ -j y"),
    paste("#$ -cwd"),
    paste("#$ -binding linear:1"),        	
    paste("#"),
    paste("# Run single task in foreground"),
    paste("module load R/3.5.1-1"),
    paste("Rscript generate_summary.R"), 
    paste("#"),
    paste("# Script ends here")
  ), f)
  close(f)
  system(paste("qsub", jobname))
  Sys.sleep(0.1)
