#!/bin/bash
#SBATCH -J screen
#SBATCH --mem=10gb 
#SBATCH --time=1:00:00 
#SBATCH --output=logs-ALCL_screen/%J.out
#SBATCH --error=logs-ALCL_screen/%J.err

model=/gpfs/data/abl/home/rodrij92/DNO3_TRAIN_C.ORIGAMI/C.Origami/checkpoints_dno3/models/epoch=35-step=10728.ckpt
seq=/gpfs/data/abl/home/rodrij92/DNO3_TRAIN_C.ORIGAMI/C.Origami/data/hg38/dna_sequence/
ctcf=/gpfs/data/abl/home/rodrij92/DNO3_TRAIN_C.ORIGAMI/C.Origami/data/hg38/dno3/genomic_features/ctcf_log2fc.bw
atac=/gpfs/data/abl/home/rodrij92/DNO3_TRAIN_C.ORIGAMI/C.Origami/data/hg38/dno3/genomic_features/atac.bw
outdir="ALCL_screen"

chr=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$1}" ALCL_screen_hg38.bed)
start=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$2}" ALCL_screen_hg38.bed)
end=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$3}" ALCL_screen_hg38.bed)
peak_id=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$4}" ALCL_screen_hg38.bed)
peak_length=$((end-start))

source /gpfs/home/rodrij92/home_abl/miniconda3/etc/profile.d/conda.sh
conda activate corigami_dno3

sleep 5
corigami-screen --out ${outdir} --celltype ${peak_id} --chr ${chr} --model ${model} --seq ${seq} --ctcf ${ctcf} --atac ${atac} --screen-start ${start} --screen-end ${end} --perturb-width ${peak_length} --step-size ${peak_length} --save-bedgraph --padding zero --save-frames
