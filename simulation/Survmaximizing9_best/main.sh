#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --output=log/estimator_rep-%a.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --array=1-50  # Number of Rscript submitting

module load r/4.1.0

# Define and create directories for output
mkdir -p "log"
mkdir -p "data"

# Run the R script with specified parameters
Rscript "simul_Survmaximizing.R"
