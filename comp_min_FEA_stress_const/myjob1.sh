#!/bin/bash 

#SBATCH -p  shared
#SBATCH -c 1
#SBATCH --mem=24G
#SBATCH -t 00-02:00:00 
#SBATCH --mail-type=END # Choose from BEGIN, END, FAIL, ALL

module load matlab/R2023a 
matlab -nodisplay -r main


