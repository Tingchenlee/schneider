#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1:00:00
#SBATCH --mem=20GB
#SBATCH --job-name=schneider_Pt111_time
#SBATCH --output=logs/schneider.%a.log
#SBATCH --error=logs/schneider.%a.slurm.log
#SBATCH --exclude=c5003
#SBATCH --partition=west
#SBATCH --mail-user=lee.ting@northeastern.edu
#SBATCH --mail-type=FAIL,END

#an array for the job.
#SBATCH --array=1-199


####################################################
source activate rmg_env
python -u CSTR_script_time.py