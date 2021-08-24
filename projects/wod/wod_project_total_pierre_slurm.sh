#!/bin/bash

#SBATCH --job-name=wod_project_antoine

#SBATCH --partition=normal,bigmem

#SBATCH --time=99:99:99

#SBATCH --mem=12G

#SBATCH --cpus-per-task=1

#SBATCH --chdir=.

#SBATCH --output=/network/lustre/iss01/charpier/analyses/wod/Pierre/slurm-output/output-%j_%a-%x.txt

#SBATCH --error=/network/lustre/iss01/charpier/analyses/wod/Pierre/slurm-output/error-%j_%a-%x.txt

#SBATCH --array=0

module load MATLAB/R2019b

matlab -nodesktop -softwareopengl -nosplash -nodisplay -r "wod_project_total_pierre($SLURM_ARRAY_TASK_ID);"

sleep 5;

