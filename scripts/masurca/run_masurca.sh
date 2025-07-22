#!/bin/bash
#SBATCH --job-name=masurca_pooled
#SBATCH --output=masurca_pooled_%j.out
#SBATCH --error=masurca_pooled_%j.err
#SBATCH --cpus-per-task=128
#SBATCH --mem=230G
#SBATCH --time=24:00:00
#SBATCH --partition=work

# Set these variables
SAMPLES_DIR="/scratch/pawsey0149/jbruton/daisy_pangenome/untrimmed_unmapped_17072025"  # <-- CHANGE THIS
RESULTS_DIR="/scratch/pawsey0149/jbruton/daisy_pangenome/masurca_build"  # <-- CHANGE THIS
WORK_DIR="/scratch/pawsey0149/jbruton/daisy_pangenome/masurca_build/_masurca_tmp"  # <-- CHANGE THIS

# MaSuRCA Configuration Parameters # < --- UPDATE THESE AS NEEDED
INSERT_SIZE=231
INSERT_STDEV=33.2
GRAPH_KMER_SIZE="auto"
USE_LINKING_MATES=1
LIMIT_JUMP_COVERAGE=300
CA_PARAMETERS="cgwErrorRate=0.15"
KMER_COUNT_THRESHOLD=1
JF_SIZE=200000000
SOAP_ASSEMBLY=0
EXTEND_JUMP_READS=0
LHE_COVERAGE=25
CLOSE_GAPS=1
FLYE_ASSEMBLY=1

# Make directories.
mkdir -p ${RESULTS_DIR}
mkdir -p ${WORK_DIR}
mkdir -p logs

# Set these variables for conda paths.
source /software/projects/pawsey0149/jbruton/miniconda3/etc/profile.d/conda.sh
conda activate /software/projects/pawsey0149/jbruton/miniconda3/envs/masurca

# Generate sample_ids.txt if it does not exist
if [ ! -f sample_ids.txt ]; then
  ls $SAMPLES_DIR/*_R1.unmapped.fastq.gz | sed 's#.*/\(.*\)_R1.unmapped.fastq.gz#\1#' > sample_ids.txt
fi

echo "Starting pooled MaSuRCA assembly"
echo "Job ID: ${SLURM_JOB_ID}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

# Create working directory for pooled assembly
POOLED_WORK_DIR="${WORK_DIR}/pooled_assembly"
mkdir -p ${POOLED_WORK_DIR}
cd ${POOLED_WORK_DIR}

# Find all sample files
PE1_FILES=($(find ${SAMPLES_DIR} -name "*_R1.unmapped.fastq.gz" | sort))
PE2_FILES=($(find ${SAMPLES_DIR} -name "*_R2.unmapped.fastq.gz" | sort))

echo "Found ${#PE1_FILES[@]} R1 files and ${#PE2_FILES[@]} R2 files"

# Check that we have equal numbers of R1 and R2 files
if [[ ${#PE1_FILES[@]} -ne ${#PE2_FILES[@]} ]]; then
    echo "ERROR: Unequal number of R1 and R2 files found"
    echo "R1 files: ${#PE1_FILES[@]}"
    echo "R2 files: ${#PE2_FILES[@]}"
    exit 1
fi

if [[ ${#PE1_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No input files found in ${SAMPLES_DIR}"
    exit 1
fi

# Start creating MaSuRCA configuration file
cat > config.txt << EOF
# MaSuRCA configuration file for pooled assembly
DATA
EOF

# Add PE entries for each sample pair, won't break until 676 unique ID's are used. 
LETTERS=(A B C D E F G H I J K L M N O P Q R S T U V W X Y Z)

generate_library_id() {
    local num=$1
    local first_idx=$(( num / 26 ))
    local second_idx=$(( num % 26 ))
    
    if [[ $first_idx -ge 26 ]]; then
        echo "ERROR: Too many samples (>676)"
        exit 1
    fi
    
    echo "${LETTERS[$first_idx]}${LETTERS[$second_idx]}"
}

LIBRARY_COUNTER=0
for i in "${!PE1_FILES[@]}"; do
    PE1="${PE1_FILES[$i]}"
    PE2="${PE2_FILES[$i]}"
    
    # Generate proper 2-character ID: AA, AB, AC, ..., AZ, BA, BB, ..., ZZ
    LIBRARY_ID=$(generate_library_id $LIBRARY_COUNTER)
    
    # Extract sample name for logging
    SAMPLE=$(basename ${PE1} | sed 's/_merged_unmapped_R1.fastq.gz//')
    
    echo "Adding sample ${SAMPLE} as library ${LIBRARY_ID}"
    echo "  R1: ${PE1}"
    echo "  R2: ${PE2}"
    
    # Add PE entry to config
    echo "PE= ${LIBRARY_ID} ${INSERT_SIZE} ${INSERT_STDEV} ${PE1} ${PE2}" >> config.txt
    
    # Increment counter
    ((LIBRARY_COUNTER++))
done

# Finish the configuration file
cat >> config.txt << EOF
END

PARAMETERS
USE_GRID=0
GRID_ENGINE=SLURM
GRID_QUEUE=work
GRID_BATCH_SIZE=500000000
EXTEND_JUMP_READS=${EXTEND_JUMP_READS}
GRAPH_KMER_SIZE = ${GRAPH_KMER_SIZE}
USE_LINKING_MATES = ${USE_LINKING_MATES}
LHE_COVERAGE=${LHE_COVERAGE}
LIMIT_JUMP_COVERAGE = ${LIMIT_JUMP_COVERAGE}
CA_PARAMETERS = ${CA_PARAMETERS}
CLOSE_GAPS=${CLOSE_GAPS}
KMER_COUNT_THRESHOLD = ${KMER_COUNT_THRESHOLD}
NUM_THREADS = ${SLURM_CPUS_PER_TASK}
JF_SIZE = ${JF_SIZE}
SOAP_ASSEMBLY = ${SOAP_ASSEMBLY}
FLYE_ASSEMBLY=${FLYE_ASSEMBLY}
END
EOF

echo "Created configuration file with ${#PE1_FILES[@]} sample pairs"
echo "Configuration file contents:"
cat config.txt

# Run Masurca.
echo "Starting MaSuRCA for pooled assembly..."
masurca config.txt

if [[ ! -f "assemble.sh" ]]; then
    echo "Error: assemble.sh was not created"
    exit 1
fi

chmod +x assemble.sh
./assemble.sh

# Check if assembly completed successfully
if [[ -f "CA/final.genome.scf.fasta" ]]; then
    echo "Pooled assembly completed successfully"
    
    # Copy results to output directory
    cp CA/final.genome.scf.fasta ${RESULTS_DIR}/pooled_assembly.fasta
    cp runCA.out ${RESULTS_DIR}/pooled_runCA.log
    cp config.txt ${RESULTS_DIR}/pooled_config.txt
    
    # Optional: compress and copy additional outputs
    # tar -czf ${RESULTS_DIR}/pooled_full_output.tar.gz *
    
    echo "Results copied to Results directory"
    echo "Final assembly: ${RESULTS_DIR}/pooled_assembly.fasta"
else
    echo "Error: Pooled assembly failed"
    exit 1
fi

# Clean up temporary files (optional)
# cd ${WORK_DIR}
# rm -rf ${POOLED_WORK_DIR}
echo "Pooled Masurca assembly job completed."
