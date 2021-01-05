#!/bin/bash
#SBATCH --job-name=bth-ShanEntr
#SBATCH --partition=serial
#SBATCH --time=1:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/bhoover/ATMS/GDAS_ensemble/Shannon_entropy/bth-ShanEntr_output.txt
source /etc/bashrc
/data/users/bhoover/ATMS/GDAS_ensemble/Shannon_entropy/run_Compute_Shannon_Entropy.sh  # job runs here
