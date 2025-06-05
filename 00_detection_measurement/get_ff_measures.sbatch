#!/usr/bin/bash

#SBATCH --job-name=ffsummary
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=500GB
#SBATCH --mail-user=vrtilek.1@osu.edu
#SBATCH --mail-type=FAIL
#SBATCH --account=PAS1986

module load R/4.0.2-gnu9.1
Rscript get_ff_measures3.R $args