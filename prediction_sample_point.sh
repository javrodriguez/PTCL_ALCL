#!/bin/bash
#SBATCH -J SamplePoint
#SBATCH --mem=4gb 
#SBATCH --time=24:00:00 
#SBATCH --output=logs-prediction-sample_point/%J.log
#SBATCH --error=logs-prediction-sample_point/%J.out

model=$1
seq=$2
outdir=$3

sample=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$1}" sample_sheet.txt)
atac=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$2}" sample_sheet.txt)
ctcf=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$3}" sample_sheet.txt)

source /gpfs/home/rodrij92/home_abl/miniconda3/etc/profile.d/conda.sh
conda activate corigami_analysis

declare -a chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX")
#for chr in "${chromosomes[@]}"; do read -a points <<< `grep -w "$chr" target_points.txt | cut -f2| tr '\n' ' '`; for point in "${points[@]}"; do corigami-predict --out ${outdir} --celltype ${sample} --chr "${chr}" --start "${point}" --model ${model} --seq ${seq} --ctcf ${ctcf} --atac ${atac}; done; done

for chr in "${chromosomes[@]}"; do read -a points <<< `grep -w "$chr" target_points.txt | cut -f2| tr '\n' ' '`; for point in "${points[@]}"; do if [ ! -f ${outdir}/${sample}/prediction/npy/"$chr"_"${point}".npy ]; then corigami-predict --out ${outdir} --celltype ${sample} --chr "${chr}" --start "${point}" --model ${model} --seq ${seq} --ctcf ${ctcf} --atac ${atac}; else echo "File exists. Skipping."; fi; done; done
