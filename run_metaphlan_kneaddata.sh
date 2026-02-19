#!/bin/bash

# Define fixed paths directly
METAPHLAN_DB="/home/cclare/miniconda3/envs/metagenomics_env/lib/python3.9/site-packages/metaphlan/metaphlan_databases"

METAPHLAN_CMD="${CONDA_PREFIX}/bin/metaphlan"
METAPHLAN_INDEX="mpa_vJun23_CHOCOPhlAnSGB_202403"

# Source config for other paths
source ~/projects/metagenomics/scripts/config.sh

# Define KneadData directory
KNEADDATA_DIR="${PROCESSED_DATA_DIR}/kneaddata"

# Count total samples
TOTAL_SAMPLES=${#ALL_SAMPLES[@]}
PROCESSED=0
SUCCESSFUL=0
FAILED=0

echo "=================================="
echo "Starting MetaPhlAn analysis on ${TOTAL_SAMPLES} samples using KneadData output"
echo "Using database: ${METAPHLAN_DB}"
echo "Using index: ${METAPHLAN_INDEX}"
echo "Start time: $(date)"
echo "=================================="



# Process each sample
for sample in "${ALL_SAMPLES[@]}"; do
    # Create output directories
    mkdir -p "${TAXONOMY_DIR}/${sample}"
    mkdir -p "${LOGS_DIR}/${sample}"
    
    # Update counter 
    PROCESSED=$((PROCESSED + 1))
    
    echo "--------------------------------"
    echo "Processing sample ${sample} [${PROCESSED}/${TOTAL_SAMPLES}]"
    echo "Start time: $(date)"
    
   # Check if output files already exist
   if [[ -f "${TAXONOMY_DIR}/${sample}/${sample}_profile.txt" ]] && \
      [[ -f "${TAXONOMY_DIR}/${sample}/${sample}.bowtie2.bz2" ]] && \
      [[ -f "${TAXONOMY_DIR}/${sample}_profile.txt" ]]; then
       echo "Output files already exist for sample ${sample}, skipping..."
       SKIPPED=$((SKIPPED + 1))
       echo "Progress: ${PROCESSED}/${TOTAL_SAMPLES} (Successful: ${SUCCESSFUL}, Failed: ${FAILED}, Skipped: ${SKIPPED})"
       echo "--------------------------------"
       continue
   fi


    # Clear any existing files
    rm -f "${TAXONOMY_DIR}/${sample}/${sample}.bowtie2.bz2"
    rm -f "${TAXONOMY_DIR}/${sample}/${sample}_profile.txt"
    
    # Define KneadData output files (check for both paired and trimmed formats)
    # First format (paired)
    R1_KNEADDATA_PAIRED_GZ="${KNEADDATA_DIR}/${sample}/${sample}_kneaddata_paired_1.fastq.gz"
    R2_KNEADDATA_PAIRED_GZ="${KNEADDATA_DIR}/${sample}/${sample}_kneaddata_paired_2.fastq.gz"
    R1_KNEADDATA_PAIRED="${KNEADDATA_DIR}/${sample}/${sample}_kneaddata_paired_1.fastq"
    R2_KNEADDATA_PAIRED="${KNEADDATA_DIR}/${sample}/${sample}_kneaddata_paired_2.fastq"
    
    # Second format (trimmed)
    R1_KNEADDATA_TRIMMED_GZ="${KNEADDATA_DIR}/${sample}/${sample}_kneaddata.trimmed.1.fastq.gz"
    R2_KNEADDATA_TRIMMED_GZ="${KNEADDATA_DIR}/${sample}/${sample}_kneaddata.trimmed.2.fastq.gz"
    R1_KNEADDATA_TRIMMED="${KNEADDATA_DIR}/${sample}/${sample}_kneaddata.trimmed.1.fastq"
    R2_KNEADDATA_TRIMMED="${KNEADDATA_DIR}/${sample}/${sample}_kneaddata.trimmed.2.fastq"
    
    R1_KNEADDATA_UNMATCHED_GZ="${KNEADDATA_DIR}/${sample}/${sample}_kneaddata_unmatched_1.fastq.gz"
    R2_KNEADDATA_UNMATCHED_GZ="${KNEADDATA_DIR}/${sample}/${sample}_kneaddata_unmatched_2.fastq.gz"
    R1_KNEADDATA_UNMATCHED="${KNEADDATA_DIR}/${sample}/${sample}_kneaddata_unmatched_1.fastq"
    R2_KNEADDATA_UNMATCHED="${KNEADDATA_DIR}/${sample}/${sample}_kneaddata_unmatched_2.fastq"

# Determine which files to use
if [[ -f "$R1_KNEADDATA_PAIRED_GZ" && -f "$R2_KNEADDATA_PAIRED_GZ" ]]; then
    R1="$R1_KNEADDATA_PAIRED_GZ"
    R2="$R2_KNEADDATA_PAIRED_GZ"
    echo "Using compressed KneadData paired output files"
    
    # Check file sizes
    R1_SIZE=$(stat -c %s "$R1")
    R2_SIZE=$(stat -c %s "$R2")
    
    # If paired files are too small, check unmatched files
    if [[ "$R1_SIZE" -lt 1000 || "$R2_SIZE" -lt 1000 ]]; then
        echo "WARNING: Paired files are too small (${R1_SIZE} and ${R2_SIZE} bytes)"
        
        # Try unmatched files
        if [[ -f "$R1_KNEADDATA_UNMATCHED_GZ" && -f "$R2_KNEADDATA_UNMATCHED_GZ" ]]; then
            R1="$R1_KNEADDATA_UNMATCHED_GZ"
            R2="$R2_KNEADDATA_UNMATCHED_GZ"
            echo "Falling back to compressed unmatched files"
        elif [[ -f "$R1_KNEADDATA_UNMATCHED" && -f "$R2_KNEADDATA_UNMATCHED" ]]; then
            R1="$R1_KNEADDATA_UNMATCHED"
            R2="$R2_KNEADDATA_UNMATCHED"
            echo "Falling back to uncompressed unmatched files"
        fi
    fi
elif [[ -f "$R1_KNEADDATA_PAIRED" && -f "$R2_KNEADDATA_PAIRED" ]]; then
    R1="$R1_KNEADDATA_PAIRED"
    R2="$R2_KNEADDATA_PAIRED"
    echo "Using uncompressed KneadData paired output files"
elif [[ -f "$R1_KNEADDATA_UNMATCHED_GZ" && -f "$R2_KNEADDATA_UNMATCHED_GZ" ]]; then
    R1="$R1_KNEADDATA_UNMATCHED_GZ"
    R2="$R2_KNEADDATA_UNMATCHED_GZ"
    echo "Using compressed KneadData unmatched output files"
elif [[ -f "$R1_KNEADDATA_UNMATCHED" && -f "$R2_KNEADDATA_UNMATCHED" ]]; then
    R1="$R1_KNEADDATA_UNMATCHED"
    R2="$R2_KNEADDATA_UNMATCHED"
    echo "Using uncompressed KneadData unmatched output files"
elif [[ -f "$R1_KNEADDATA_TRIMMED_GZ" && -f "$R2_KNEADDATA_TRIMMED_GZ" ]]; then
    R1="$R1_KNEADDATA_TRIMMED_GZ"
    R2="$R2_KNEADDATA_TRIMMED_GZ"
    echo "Using compressed KneadData trimmed output files"
elif [[ -f "$R1_KNEADDATA_TRIMMED" && -f "$R2_KNEADDATA_TRIMMED" ]]; then
    R1="$R1_KNEADDATA_TRIMMED"
    R2="$R2_KNEADDATA_TRIMMED"
    echo "Using uncompressed KneadData trimmed output files"
else
    echo "ERROR: KneadData output files missing for ${sample}"
    echo "Checked for: $R1_KNEADDATA_PAIRED_GZ, $R2_KNEADDATA_PAIRED_GZ, $R1_KNEADDATA_PAIRED, $R2_KNEADDATA_PAIRED, $R1_KNEADDATA_UNMATCHED_GZ, $R2_KNEADDATA_UNMATCHED_GZ, $R1_KNEADDATA_UNMATCHED, $R2_KNEADDATA_UNMATCHED, $R1_KNEADDATA_TRIMMED_GZ, $R2_KNEADDATA_TRIMMED_GZ, $R1_KNEADDATA_TRIMMED, $R2_KNEADDATA_TRIMMED"
    
    # Don't use fallback files, just skip this sample
    echo "ERROR: KneadData files not found for ${sample}, skipping this sample"
    FAILED=$((FAILED + 1))
    echo "Progress: ${PROCESSED}/${TOTAL_SAMPLES} (Successful: ${SUCCESSFUL}, Failed: ${FAILED})"
    echo "--------------------------------"
    continue
fi

echo "Input files:"
echo "R1: $R1"
echo "R2: $R2"
echo "Running MetaPhlAn..."



    
    # Create log file
    LOG_FILE="${LOGS_DIR}/${sample}/metaphlan_kneaddata.log"
    echo "Starting MetaPhlAn for ${sample} with KneadData output at $(date)" > "$LOG_FILE"
    
    # Create temporary uncompressed files
    echo "Decompressing files for MetaPhlAn..."
    R1_TEMP="${KNEADDATA_DIR}/${sample}/${sample}_kneaddata_paired_1_temp.fastq"
    R2_TEMP="${KNEADDATA_DIR}/${sample}/${sample}_kneaddata_paired_2_temp.fastq"

    # Decompress the files
    zcat "${R1}" > "${R1_TEMP}"
    zcat "${R2}" > "${R2_TEMP}"

    # Check if files were created and have content
    echo "Checking decompressed files:"
    ls -la "${R1_TEMP}"
    ls -la "${R2_TEMP}"
    echo "First few lines of R1:"
    head -n 4 "${R1_TEMP}"

    # Run MetaPhlAn with decompressed files
    ${METAPHLAN_CMD} \
        "${R1_TEMP},${R2_TEMP}" \
        --input_type fastq \
        --bowtie2db "${METAPHLAN_DB}" \
        --index "${METAPHLAN_INDEX}" \
        --nproc "${METAPHLAN_THREADS}" \
        --sample_id "${sample}" \
        --add_viruses \
        -t "${METAPHLAN_ANALYSIS_TYPE}" \
        --bowtie2out "${TAXONOMY_DIR}/${sample}/${sample}.bowtie2.bz2" \
        --samout "${SAM_DIR}/${sample}/${sample}.sam.bz2" \
        -o "${TAXONOMY_DIR}/${sample}/${sample}_profile.txt" \
        2>&1 | tee -a "$LOG_FILE"

    # Remove temporary files after processing
    rm "${R1_TEMP}" "${R2_TEMP}"


    METAPHLAN_EXIT=$?
    
    # Check if successful
    if [[ $METAPHLAN_EXIT -eq 0 && -f "${TAXONOMY_DIR}/${sample}/${sample}_profile.txt" && -s "${TAXONOMY_DIR}/${sample}/${sample}_profile.txt" ]]; then
        echo "MetaPhlAn completed successfully for sample: ${sample}"
        # Copy to main directory
        cp "${TAXONOMY_DIR}/${sample}/${sample}_profile.txt" "${TAXONOMY_DIR}/"
        echo "Copied profile to ${TAXONOMY_DIR}/"
        SUCCESSFUL=$((SUCCESSFUL + 1))
    else
        echo "ERROR: MetaPhlAn failed for sample: ${sample} with exit code ${METAPHLAN_EXIT}"
        FAILED=$((FAILED + 1))
    fi
        echo "End time: $(date)"
    echo "Progress: ${PROCESSED}/${TOTAL_SAMPLES} (Successful: ${SUCCESSFUL}, Failed: ${FAILED})"
    echo "--------------------------------"
done

# Merge all profiles into a single table if successful profiles exist
if [[ ${SUCCESSFUL} -gt 0 ]]; then
    echo "Merging taxonomic profiles..."
    /opt/anaconda3/envs/metagenomics_x86/bin/merge_metaphlan_tables.py ${TAXONOMY_DIR}/*_profile.txt > "${TAXONOMY_DIR}/merged_profiles.txt"

    # Convert to more accessible format (relative abundance only)
    grep -v "#" "${TAXONOMY_DIR}/merged_profiles.txt" | sed 's/^.*|//g' > "${TAXONOMY_DIR}/merged_profiles_cleaned.txt"
else
    echo "No successful profiles to merge."
fi

echo "===========================================" 
echo "Taxonomic profiling completed!"
echo "Total samples: ${TOTAL_SAMPLES}"
echo "Successfully processed: ${SUCCESSFUL}"
echo "Failed: ${FAILED}"
echo "End time: $(date)"
echo "==========================================="
