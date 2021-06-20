#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1:00:00
#SBATCH --mem=20GB
#SBATCH --job-name=schneider_no_gasrxn_PFR_temp
#SBATCH --output=logs/PFR/schneider.%a.log
#SBATCH --error=logs/PFR/schneider.%a.slurm.log
#SBATCH --partition=short
#SBATCH --mail-user=lee.ting@northeastern.edu
#SBATCH --mail-type=FAIL,END

#an array for the job.
#SBATCH --array=1-22


####################################################
source activate rmg_env
python -u PFR_script_temp.py