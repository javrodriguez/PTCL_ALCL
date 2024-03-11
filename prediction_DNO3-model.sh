#!/bin/bash
#SBATCH -J predBatch
#SBATCH --mem=1gb 
#SBATCH --time=4:00:00 
#SBATCH --output=logs-prediction-DNO3-model/%J.log
#SBATCH --error=logs-prediction-DNO3-model/%J.out

sample_sheet=/gpfs/data/abl/home/rodrij92/DNO3_TRAIN_C.ORIGAMI/C.Origami/sample_atac_ctcf_dno3.txt
model=/gpfs/data/abl/home/rodrij92/DNO3_TRAIN_C.ORIGAMI/C.Origami/checkpoints_dno3/models/epoch=35-step=10728.ckpt
target_points=/gpfs/data/abl/home/rodrij92/DNO3_TRAIN_C.ORIGAMI/C.Origami/hg38_startPos_w2097152.txt
seq=/gpfs/data/abl/home/rodrij92/DNO3_TRAIN_C.ORIGAMI/C.Origami/data/hg38/dna_sequence/
outdir="prediction_DNO3-model"

cp $sample_sheet sample_sheet.txt
cp $target_points target_points.txt

n_samples=`cat sample_sheet.txt| wc -l`

sbatch --array=1-$n_samples prediction_sample_point.sh ${model} ${seq} ${outdir}
