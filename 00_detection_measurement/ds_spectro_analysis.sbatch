#!/usr/bin/bash

#SBATCH --job-name=peakf
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=24GB
#SBATCH --mail-user=vrtilek.1@osu.edu
#SBATCH --mail-type=FAIL
#SBATCH --account=PAS1986

module load gcc/12.3.0
module load R/4.4.0
Rscript ds_spectro_analysis_Marcelofix.R $args