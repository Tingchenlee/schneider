#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1:00:00
#SBATCH --mem=20GB
#SBATCH --job-name=schneider_O2try
#SBATCH --output=logs/CSTR_try/schneider.%a.log
#SBATCH --error=logs/CSTR_try/schneider.%a.slurm.log
#SBATCH --partition=short

#an array for the job.
#SBATCH --array=1-199


####################################################
source activate rmg_env
python -u CSTR_script_try.py