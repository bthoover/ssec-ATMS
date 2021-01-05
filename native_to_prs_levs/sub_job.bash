#!/bin/bash
#SBATCH --job-name=compute_ensemble_stats
#SBATCH --partition=serial
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/data/users/bhoover/ATMS/GDAS_ensemble/native_to_prs_levs/compute_ensemble_stats_output.txt

inp01=${1}
inp02=${2}
inp03=${3}
inp04=${4}
inp05=${5}
inp06=${6}
inp07=${7}
inp08=${8}
inp09=${9}
inp10=${10}
inp11=${11}
inp12=${12}
inp13=${13}
inp14=${14}
inp15=${15}
inp16=${16}
inp17=${17}
inp18=${18}
inp19=${19}
inp20=${20}
inp21=${21}
inp22=${22}

source /etc/bashrc
module purge

  # job runs here
sbatch /data/users/bhoover/ATMS/GDAS_ensemble/native_to_prs_levs/exec_compute_ensemble_stats.bash ${inp01} ${inp02} ${inp03} ${inp04} ${inp05} ${inp06} ${inp07} ${inp08} ${inp09} ${inp10} ${inp11} ${inp12} ${inp13} ${inp14} ${inp15} ${inp16} ${inp17} ${inp18} ${inp19} ${inp20} ${inp21} ${inp22} 
