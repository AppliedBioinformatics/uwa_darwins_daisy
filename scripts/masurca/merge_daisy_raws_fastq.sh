#!/bin/bash
#SBATCH --job-name=merge_fastq
#SBATCH --cpus-per-task=3
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=merge_fastq_%A_%a.out
#SBATCH --error=merge_fastq_%A_%a.err
#SBATCH --array=1-1  # Replace N with the number of samples in samples.txt
#SBATCH --partition=work

# Merge fastq.gz files by sample identifier and read direction
# Usage: bash merge_fastq.sh

shopt -s nullglob

# Get the sample name for this array job
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)

for read in 1 2; do
    files=( ${sample}_*_${read}.fq.gz )
    if [ ${#files[@]} -gt 0 ]; then
        out="../merged_raw_fastqs/${sample}_R${read}.fastq.gz"
        echo "Merging ${#files[@]} files for $sample read $read into $out"
        cat "${files[@]}" > "$out"
    else
        echo "No files found for $sample read $read"
    fi
done
