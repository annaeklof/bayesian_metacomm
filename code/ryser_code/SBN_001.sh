#!/bin/bash
#SBATCH -J SBN_001.sh
#SBATCH -A liu-2017-00089-2
#SBATCH -t 10:00:00
#SBATCH --mem=16000
#SBATCH -n 1
#
# Run single task in foreground
module add R/3.6.0-nsc1-gcc-7.3.0
cd /home/x_gyoba/bayesian_metacomm_round2/code/ryser_code/
Rscript run_local.R ../../data/ryser_data/results/tldat.csv
#
# Script ends here
