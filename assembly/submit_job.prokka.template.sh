#!/usr/bin/env bash
# The name to show in queue lists for this job:
#SBATCH -J template_job

# Number of desired cpus:
#SBATCH --cpus-per-task=8

# Amount of RAM needed for this job:
#SBATCH --mem=32gb

# The time the job will be running:
#SBATCH --time=48:00:00

# To use GPUs you have to request them:
##SBATCH --gres=gpu:1
#SBATCH --constraint=cal

# Set output and error files
#SBATCH --error=template_job.%J.err
#SBATCH --output=template_job.%J.out

# the program to execute with its parameters:
#time ./mi_programa argumentos${SLURM_ARRAYID}.jpg > out_$SLURM_ARRAYID.out
hostname

# $SLURM_CPUS_PER_TASK has the same value requested in --cpus-per-task 
