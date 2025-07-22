#!/bin/bash
#SBATCH --job-name=extract
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --partition=work
#SBATCH --output=out
#SBATCH --error=err
#SBATCH --mail-user=jack.bruton@uwa.edu.au
#SBATCH --mem=4G
#SBATCH --array=1-34

# Generic SLURM commands
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --clusters=setonix
#SBATCH --account=pawsey0149
#SBATCH --mail-type=ALL
#SBATCH --export=NONE

# Run this script using sbatch [slurm].sh

# Set these variables for conda paths.
source /software/projects/pawsey0149/jbruton/miniconda3/etc/profile.d/conda.sh
conda activate /software/projects/pawsey0149/jbruton/miniconda3/envs/biogenerics

echo "========================================="
echo "SLURM_JOB_ID = $SLURM_JOB_ID"
echo "SLURM_NODELIST = $SLURM_NODELIST"

if [ ! -z $SLURM_ARRAY_TASK_ID ]; then
	echo "SLURM_ARRAY_TASK_ID = $SLURM_ARRAY_TASK_ID"
fi
echo "========================================="

# Script to run (srun -m command recommended by Pawsey to pack threads)
time srun -m block:block:block

mkdir -p extracted_unmapped
mkdir -p logs

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)
echo "Processing sample: $sample"

# Define raw reads
raw_r1="../merged_raw_fastqs/${sample}_R1.fastq.gz"
raw_r2="../merged_raw_fastqs/${sample}_R2.fastq.gz"

# Find unmapped files (assume only one file per sample)
unmapped_r1_file="../samtools_unmapped/${sample}_merged_unmapped_R1.fastq.gz"
unmapped_r2_file="../samtools_unmapped/${sample}_merged_unmapped_R2.fastq.gz"

if [[ -f "$unmapped_r1_file" && -f "$unmapped_r2_file" && -f "$raw_r1" && -f "$raw_r2" ]]; then

    # Extract read IDs without /1 or /2 suffixes
    zcat "$unmapped_r1_file" | awk 'NR % 4 == 1' | sed 's/^@//; s/\/[12]$//' > "${sample}_unmapped_ids.txt"

    # Extract unmapped reads from raw files using seqtk
    seqtk subseq "$raw_r1" "${sample}_unmapped_ids.txt" | gzip > "extracted_unmapped/${sample}_R1.unmapped.fastq.gz"
    seqtk subseq "$raw_r2" "${sample}_unmapped_ids.txt" | gzip > "extracted_unmapped/${sample}_R2.unmapped.fastq.gz"

    echo "Done processing sample: $sample"

    # Cleanup temporary files
    rm "${sample}_unmapped_ids.txt"

else
    echo "Missing file(s) for $sample"
fi
