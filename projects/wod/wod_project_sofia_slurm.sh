#!/bin/bash

#SBATCH --job-name=wod

#SBATCH --partition=bigmem

#SBATCH --time=99:99:99

#SBATCH --mem=120G

#SBATCH --cpus-per-task=2

#SBATCH --chdir=.

#SBATCH --output=/network/lustre/iss01/charpier/analyses/wod/Sofia/slurm-output/-%j_%a-%x.txt

#SBATCH --error=/network/lustre/iss01/charpier/analyses/wod/Sofia/slurm-output/error-%j_%a-%x.txt

#SBATCH --array=1-16

module load MATLAB/R2019b

matlab -nodesktop -softwareopengl -nosplash -nodisplay -r "wod_project_sofia($SLURM_ARRAY_TASK_ID);"

sleep 5;

