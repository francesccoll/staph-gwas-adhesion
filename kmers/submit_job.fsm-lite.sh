#!/usr/bin/env bash
# The name to show in queue lists for this job:
#SBATCH -J fsm_file

# Number of desired cpus:
#SBATCH --cpus-per-task=1

# Amount of RAM needed for this job:
#SBATCH --mem=258gb

# The time the job will be running:
#SBATCH --time=48:00:00

# To use GPUs you have to request them:
##SBATCH --gres=gpu:1
#SBATCH --constraint=cal

# Set output and error files
#SBATCH --error=fsm_file.%J.err
#SBATCH --output=fsm_file.%J.out

# the program to execute with its parameters:
#time ./mi_programa argumentos${SLURM_ARRAYID}.jpg > out_$SLURM_ARRAYID.out
hostname

# $SLURM_CPUS_PER_TASK has the same value requested in --cpus-per-task 

fsm-lite -l fsm_file_list.txt -m 30 -M 100 -v -t mzIE_phen_samples.fsm_kmers > mzIE_phen_samples.fsm_kmers.30-100bp.txt