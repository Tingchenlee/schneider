#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1:00:00
#SBATCH --mem=20GB
#SBATCH --job-name=schneider_Rh111
#SBATCH --output=logs/schneider.%a.log
#SBATCH --error=logs/schneider.%a.slurm.log
#SBATCH --partition=short
#SBATCH --mail-user=lee.ting@northeastern.edu
#SBATCH --mail-type=FAIL,END

#an array for the job.
#SBATCH --array=1-99


####################################################
source activate rmg_env
python -u CSTR_script_Rh111.py