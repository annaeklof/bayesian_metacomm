#!/bin/bash
#$ -S /bin/bash
#$ -N summarize_serengeit_009.sh
#$ -o /work/$USER/$JOB_NAME-$JOB_ID.log
#$ -wd /work/$USER
#$ -l h_rt=12:00:00
#$ -l h_vmem=8G
#$ -j y
#$ -cwd
#$ -binding linear:1
#
# Run single task in foreground
module load R/3.5.1-1
Rscript generate_and_modify_summary.R
#
# Script ends here
