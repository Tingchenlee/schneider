#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1:00:00
#SBATCH --mem=20GB
#SBATCH --job-name=schneider_Pt211_time
#SBATCH --output=logs/schneider.%a.log
#SBATCH --error=logs/schneider.%a.slurm.log
#SBATCH --exclude=c5003
#SBATCH --partition=west


#an array for the job.
#SBATCH --array=1-99


####################################################
source activate rmg_env
python -u CSTR_script_time.py