#!/bin/bash
#SBATCH --job-name=diamond_annot
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=diamond_%j.out
#SBATCH --error=diamond_%j.err
#SBATCH --partition=work

source /software/projects/pawsey0149/jbruton/miniconda3/etc/profile.d/conda.sh
conda activate /software/projects/pawsey0149/jbruton/miniconda3/envs/biogenerics

diamond blastp -q protein_clean.fa -d uniprot_sprot.dmnd -o diamond_results.tsv -p 64 --sensitive -e 1e-5 -k 1 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids stitle sscinames sskingdoms skingdoms sphylums
