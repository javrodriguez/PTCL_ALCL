#!/bin/bash
#SBATCH -J prediction
#SBATCH --mem=10gb 
#SBATCH --time=4:00:00 
#SBATCH --output=logs-prediction/%J.log
#SBATCH --error=logs-prediction/%J.out

outdir="prediction_ITK_locus"
celltype="DNO3"
model=/gpfs/data/abl/home/rodrij92/DNO3_TRAIN_C.ORIGAMI/C.Origami/checkpoints_dno3/models/epoch=35-step=10728.ckpt
seq=/gpfs/data/abl/home/rodrij92/DNO3_TRAIN_C.ORIGAMI/C.Origami/data/hg38/dna_sequence/
ctcf=/gpfs/data/abl/home/rodrij92/DNO3_TRAIN_C.ORIGAMI/C.Origami/data/hg38/dno3/genomic_features/ctcf_log2fc.bw
atac=/gpfs/data/abl/home/rodrij92/DNO3_TRAIN_C.ORIGAMI/C.Origami/data/hg38/dno3/genomic_features/atac.bw

chr="chr5"
start=156094357

module purge
module add default-environment

source /gpfs/home/rodrij92/home_abl/miniconda3/etc/profile.d/conda.sh
conda activate corigami_dno3

corigami-predict --out ${outdir} --celltype ${celltype} --chr ${chr} --start ${start} --model ${model} --seq ${seq} --ctcf ${ctcf} --atac ${atac}
