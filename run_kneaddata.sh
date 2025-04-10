#!/bin/bash
# Source configuration
source ~/projects/metagenomics/scripts/config.sh

# Set KneadData specific parameters
KNEADDATA_DB="/home/cclare/projects/metagenomics/databases/kneaddata/"

# Create a separate directory for KneadData outputs
KNEADDATA_DIR="${PROCESSED_DATA_DIR}/kneaddata"
mkdir -p "${KNEADDATA_DIR}"

for sample in "${ALL_SAMPLES[@]}"; do
    echo "Processing KneadData for sample: ${sample}"
    # Create sample-specific directories
    mkdir -p "${KNEADDATA_DIR}/${sample}"
    mkdir -p "${LOGS_DIR}/${sample}"
    
    # Check if KneadData output already exists
    if ([ -f "${KNEADDATA_DIR}/${sample}/${sample}_kneaddata_paired_1.fastq" ] && 
        [ -f "${KNEADDATA_DIR}/${sample}/${sample}_kneaddata_paired_2.fastq" ]) || 
       ([ -f "${KNEADDATA_DIR}/${sample}/${sample}_kneaddata_paired_1.fastq.gz" ] && 
        [ -f "${KNEADDATA_DIR}/${sample}/${sample}_kneaddata_paired_2.fastq.gz" ]); then
	echo "KneadData already completed for sample: ${sample}, skipping..."
        continue
    fi
    
    # Run KneadData
    echo "Running KneadData for sample: ${sample} with 8 threads..."
    kneaddata --input1 ${RAW_DATA_DIR}/${sample}${R1_SUFFIX} \
              --input2 ${RAW_DATA_DIR}/${sample}${R2_SUFFIX} \
              --output ${KNEADDATA_DIR}/${sample} \
              --output-prefix ${sample}_kneaddata \
              --reference-db ${KNEADDATA_DB} \
              -t 8 \
	      --bypass-trf \
              --trimmomatic ${CONDA_PREFIX}/share/trimmomatic \
              --trimmomatic-options "ILLUMINACLIP:${CONDA_PREFIX}/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36" \
              --remove-intermediate-output \
              &> ${LOGS_DIR}/${sample}/kneaddata.log
    
    # Compress output files to save space
    echo "Compressing output files for sample: ${sample}"
    find "${PROCESSED_DATA_DIR}/kneaddata/${sample}" -name "*.fastq" -exec gzip {} \;
    
    echo "KneadData completed for sample: ${sample}"
done

# Run KneadData read count table to summarize results
echo "Generating KneadData read count summary..."
mkdir -p "${KNEADDATA_DIR}/stats"
kneaddata_read_count_table --input ${KNEADDATA_DIR} \
                          --output ${KNEADDATA_DIR}/stats/kneaddata_read_counts.txt \
                          &> ${LOGS_DIR}/kneaddata_stats.log

echo "Quality control processing with KneadData completed for all samples!"
