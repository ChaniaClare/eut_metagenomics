#!/bin/bash

# Source configuration
source /Users/chaniaclare/Library/CloudStorage/OneDrive-UniversityCollegeLondon/Documents/PhD/Barnes-PhD/metagenomics/projects/eut_metagenomics/scripts/config.sh


for sample in "${ALL_SAMPLES[@]}"; do
    echo "Processing QC for sample: ${sample}"

    # Create sample-specific directories
    mkdir -p "${PROCESSED_DATA_DIR}/${sample}"
    mkdir -p "${LOGS_DIR}/${sample}"
    
    # Check if trimmed output already exists
    if [ -f "${PROCESSED_DATA_DIR}/${sample}/${sample}_trimmed_R1.fastq.gz" ] && [ -f "${PROCESSED_DATA_DIR}/${sample}/${sample}_trimmed_R2.fastq.gz" ]; then
        echo "QC already completed for sample: ${sample}, skipping..."
        continue
    fi

    
    # Run FastQC
    echo "Running FastQC..."
    fastqc ${RAW_DATA_DIR}/${sample}${R1_SUFFIX} ${RAW_DATA_DIR}/${sample}${R2_SUFFIX} \
        -o ${PROCESSED_DATA_DIR}/${sample} \
        &> ${LOGS_DIR}/${sample}/fastqc.log 
        
    # Run Trimmomatic
    echo "Running Trimmomatic..."
    trimmomatic PE \
        ${RAW_DATA_DIR}/${sample}${R1_SUFFIX} \
        ${RAW_DATA_DIR}/${sample}${R2_SUFFIX} \
        ${PROCESSED_DATA_DIR}/${sample}/${sample}_trimmed_R1.fastq.gz \
        ${PROCESSED_DATA_DIR}/${sample}/${sample}_unpaired_R1.fastq.gz \
        ${PROCESSED_DATA_DIR}/${sample}/${sample}_trimmed_R2.fastq.gz \
        ${PROCESSED_DATA_DIR}/${sample}/${sample}_unpaired_R2.fastq.gz \
        ILLUMINACLIP:${CONDA_PREFIX}/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
        &> ${LOGS_DIR}/${sample}/trimmomatic.log
        
    echo "QC completed for sample: ${sample}"
done
    
# Run MultiQC to summarize all FastQC results
echo "Running MultiQC to summarize results..."
mkdir -p "${PROCESSED_DATA_DIR}/multiqc"
multiqc ${PROCESSED_DATA_DIR} -o ${PROCESSED_DATA_DIR}/multiqc &> ${LOGS_DIR}/multiqc.log
        
echo "Quality control processing completed for all samples!"
