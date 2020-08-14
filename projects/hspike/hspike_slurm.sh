#!/bin/bash
#SBATCH --job-name=Hspike
#SBATCH --partition=normal
#SBATCH --time=99:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --chdir=.
#SBATCH --output=/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/slurm_output/output-%j_%a-%x.txt
#SBATCH --error=/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/slurm_output/error-%j_%a-%x.txt
#SBATCH --array=1-7

module load MATLAB/R2019b
matlab -nodesktop -softwareopengl -nosplash -nodisplay -r "hspike_cluster($SLURM_ARRAY_TASK_ID);"
sleep 5;
