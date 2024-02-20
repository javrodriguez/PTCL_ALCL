#!/bin/bash
#SBATCH -J corigami-train
#SBATCH -c16
#SBATCH -p gpu4_long,gpu8_long,gpu4_short
#SBATCH --gres=gpu:4
#SBATCH --mem=180G
#SBATCH --time=20-00:00:00
#SBATCH --output=logs-train/%j.log

source /gpfs/home/rodrij92/home_abl/miniconda3/etc/profile.d/conda.sh
conda activate corigami_dno3

corigami-train --data-root ./data/ --save_path checkpoints_dno3 --celltype dno3 --assembly hg38 --num-gpu 4 --batch-size 4 --num-workers 16
(base) [rodrij92@bigpurple-ln2 C.Origami]$ sbatch train.sh 
