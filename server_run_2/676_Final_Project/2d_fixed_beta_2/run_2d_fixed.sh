#!/bin/sh
#SBATCH --job-name 2d_fixed_beta
#SBATCH --partition=normal

#SBATCH --mail-type=BEGIN,END,FAIL         # ALL,NONE,BEGIN,END,FAIL,REQUEUE,..
#SBATCH --mail-user=tsong8@gmu.edu     # Put your GMU email address here

#SBATCH --output=/home/tsong8/676_Final_Project/2d_fixed_beta/%x-%N-%j.out  # Output file
#SBATCH --error=/home/tsong8/676_Final_Project/2d_fixed_beta/%x-%N-%j.err   # Error file

#SBATCH --mem=32G         # Total memory needed for your job (suffixes: K,M,G,T)
#SBATCH --time=3-00:01   # Total time needed for your job: Days-Hours:Minutes

## Load the relevant modules needed for the job
module load r/3.6.3

## Start the job
Rscript "2d_fixed_beta.R"