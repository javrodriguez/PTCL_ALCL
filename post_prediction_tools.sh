#!/bin/bash
#SBATCH -J post_prediction_tools
#SBATCH --mem=20gb 
#SBATCH --time=24:00:00 
#SBATCH -c 8
#SBATCH --output=logs-post_prediction_tools/%J.out
#SBATCH --error=logs-post_prediction_tools/%J.err

module load python/cpu/3.8.11
module load r/4.1.2

cp ../post_prediction_tools.R ./
Rscript post_prediction_tools.R
