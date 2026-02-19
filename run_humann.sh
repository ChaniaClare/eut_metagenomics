#!/bin/bash

# HUMAnN3 functional profiling + EUT screening pipeline
# Usage: ./run_humann.sh [fresh_start]

source ~/projects/metagenomics/scripts/config.sh

# Check for fresh start option
FRESH_START=false
if [[ "$1" == "fresh_start" ]]; then
    FRESH_START=true
    echo "Fresh start requested - will reprocess all samples"
fi

mkdir -p "${RESULTS_DIR}/humann"
mkdir -p "${RESULTS_DIR}/eut_screening"
mkdir -p "${LOGS_DIR}/humann"
mkdir -p "${LOGS_DIR}/eut_screening"

# Checkpoint file
CHECKPOINT_FILE="${LOGS_DIR}/humann/completed_samples.txt"

# Handle fresh start
if [[ "$FRESH_START" == true ]]; then
    echo "Removing previous results for fresh start..."
    rm -rf "${RESULTS_DIR}/humann"
    rm -rf "${RESULTS_DIR}/eut_screening"
    rm -f "$CHECKPOINT_FILE"
    mkdir -p "${RESULTS_DIR}/humann"
    mkdir -p "${RESULTS_DIR}/eut_screening"
    mkdir -p "${LOGS_DIR}/humann"
    mkdir -p "${LOGS_DIR}/eut_screening"
fi

# Create checkpoint file if it doesn't exist
touch "$CHECKPOINT_FILE"

# Find HUMAnN3 command
if command -v humann3 &> /dev/null; then
    HUMANN_CMD="humann3"
elif command -v humann &> /dev/null; then
    HUMANN_CMD="humann"
elif [[ -f "${CONDA_PREFIX}/bin/humann3" ]]; then
    HUMANN_CMD="${CONDA_PREFIX}/bin/humann3"
elif [[ -f "${CONDA_PREFIX}/bin/humann" ]]; then
    HUMANN_CMD="${CONDA_PREFIX}/bin/humann"
else
    echo "ERROR: HUMAnN3 not found"
    exit 1
fi

HUMANN_THREADS=8
KNEADDATA_DIR="${PROCESSED_DATA_DIR}/kneaddata"
HUMANN_RESULTS_DIR="${RESULTS_DIR}/humann"
METAPHLAN_RESULTS_DIR="/home/cclare/projects/metagenomics/taxonomic_profiles"

# Build lookup arrays for sample group assignment from config.sh arrays
declare -A SAMPLE_GROUP_MAP
declare -A SAMPLE_PATIENT_MAP
declare -A SAMPLE_TIMEPOINT_MAP
declare -A SAMPLE_DESCRIPTION_MAP

for s in "${CONTROL_SAMPLES[@]}"; do
    SAMPLE_GROUP_MAP[$s]="Control"
    SAMPLE_PATIENT_MAP[$s]="NA"
    SAMPLE_TIMEPOINT_MAP[$s]="NA"
    SAMPLE_DESCRIPTION_MAP[$s]="Control sample"
done

for s in "${HEALTHY_SAMPLES[@]}"; do
    SAMPLE_GROUP_MAP[$s]="Healthy"
    SAMPLE_PATIENT_MAP[$s]="NA"
    SAMPLE_TIMEPOINT_MAP[$s]="NA"
    SAMPLE_DESCRIPTION_MAP[$s]="Healthy control sample"
done

for s in "${EXT_HEALTHY_SAMPLES[@]}"; do
    SAMPLE_GROUP_MAP[$s]="ExtHealthy"
    SAMPLE_PATIENT_MAP[$s]="${s}"
    SAMPLE_TIMEPOINT_MAP[$s]="NA"
    SAMPLE_DESCRIPTION_MAP[$s]="External healthy control"
done

# For treated samples, derive patient ID and timepoint from sample name convention
# Convention: PATIENTID_TIMEPOINT_S##  e.g. 10_1_S1 -> patient 10, timepoint 1
for s in "${TREATED_SAMPLES[@]}"; do
    SAMPLE_GROUP_MAP[$s]="Treated"
    patient_id=$(echo "$s" | cut -d'_' -f1)
    timepoint=$(echo "$s" | cut -d'_' -f2)
    SAMPLE_PATIENT_MAP[$s]="${patient_id}"
    SAMPLE_TIMEPOINT_MAP[$s]="${timepoint}"
    if [[ "$timepoint" == "1" ]]; then
        SAMPLE_DESCRIPTION_MAP[$s]="Pre-transplant"
    else
        SAMPLE_DESCRIPTION_MAP[$s]="Post-transplant"
    fi
done

echo "=================================================================="
echo "HUMANN3 FUNCTIONAL PROFILING + EUT SCREENING"
echo "Start time: $(date)"
echo "=================================================================="

# Check databases
HUMANN_DB_DIR="/home/cclare/projects/metagenomics/databases/humann"
CHOCOPHLAN_DB="${HUMANN_DB_DIR}/chocophlan"
UNIREF_DB="${HUMANN_DB_DIR}/uniref"

if [[ ! -d "$CHOCOPHLAN_DB" ]]; then
    echo "ERROR: ChocoPhlAn database not found at $CHOCOPHLAN_DB"
    exit 1
fi

if [[ ! -d "$UNIREF_DB" ]]; then
    echo "ERROR: UniRef database not found at $UNIREF_DB"
    exit 1
fi

# Check for MetaPhlAn results
if [[ -d "$METAPHLAN_RESULTS_DIR" ]]; then
    USE_METAPHLAN_BYPASS=true
else
    USE_METAPHLAN_BYPASS=false
fi

echo "PART 1: RUNNING HUMANN3"

TOTAL_SAMPLES=${#ALL_SAMPLES[@]}
PROCESSED=0
SUCCESSFUL=0
FAILED=0
SKIPPED=0

# Initialise metadata file with header if it does not exist
METADATA_FILE="${RESULTS_DIR}/eut_screening/sample_metadata.tsv"
if [[ ! -f "$METADATA_FILE" ]]; then
    echo -e "Sample\tGroup\tTimepoint\tPatient_ID\tDescription" > "$METADATA_FILE"
    echo "Created new sample metadata file: $METADATA_FILE"
fi

for sample in "${ALL_SAMPLES[@]}"; do
    mkdir -p "${HUMANN_RESULTS_DIR}/${sample}"
    mkdir -p "${LOGS_DIR}/humann/${sample}"

    PROCESSED=$((PROCESSED + 1))
    echo "Processing sample ${sample} [${PROCESSED}/${TOTAL_SAMPLES}]"

    # Find KneadData files - prioritise unmatched files
    R1_KNEADDATA_UNMATCHED_GZ="${KNEADDATA_DIR}/${sample}/${sample}_kneaddata_unmatched_1.fastq.gz"
    R2_KNEADDATA_UNMATCHED_GZ="${KNEADDATA_DIR}/${sample}/${sample}_kneaddata_unmatched_2.fastq.gz"
    R1_KNEADDATA_UNMATCHED="${KNEADDATA_DIR}/${sample}/${sample}_kneaddata_unmatched_1.fastq"
    R2_KNEADDATA_UNMATCHED="${KNEADDATA_DIR}/${sample}/${sample}_kneaddata_unmatched_2.fastq"

    R1_KNEADDATA_PAIRED_GZ="${KNEADDATA_DIR}/${sample}/${sample}_kneaddata_paired_1.fastq.gz"
    R2_KNEADDATA_PAIRED_GZ="${KNEADDATA_DIR}/${sample}/${sample}_kneaddata_paired_2.fastq.gz"
    R1_KNEADDATA_PAIRED="${KNEADDATA_DIR}/${sample}/${sample}_kneaddata_paired_1.fastq"
    R2_KNEADDATA_PAIRED="${KNEADDATA_DIR}/${sample}/${sample}_kneaddata_paired_2.fastq"
    R1_KNEADDATA_TRIMMED_GZ="${KNEADDATA_DIR}/${sample}/${sample}_kneaddata.trimmed.1.fastq.gz"
    R2_KNEADDATA_TRIMMED_GZ="${KNEADDATA_DIR}/${sample}/${sample}_kneaddata.trimmed.2.fastq.gz"
    R1_KNEADDATA_TRIMMED="${KNEADDATA_DIR}/${sample}/${sample}_kneaddata.trimmed.1.fastq"
    R2_KNEADDATA_TRIMMED="${KNEADDATA_DIR}/${sample}/${sample}_kneaddata.trimmed.2.fastq"

    if [[ -f "$R1_KNEADDATA_UNMATCHED_GZ" && -f "$R2_KNEADDATA_UNMATCHED_GZ" ]]; then
        INPUT_FILE="$R1_KNEADDATA_UNMATCHED_GZ"
        INPUT_FILE2="$R2_KNEADDATA_UNMATCHED_GZ"
        echo "Using compressed KneadData unmatched output files"
    elif [[ -f "$R1_KNEADDATA_UNMATCHED" && -f "$R2_KNEADDATA_UNMATCHED" ]]; then
        INPUT_FILE="$R1_KNEADDATA_UNMATCHED"
        INPUT_FILE2="$R2_KNEADDATA_UNMATCHED"
        echo "Using uncompressed KneadData unmatched output files"
    elif [[ -f "$R1_KNEADDATA_PAIRED_GZ" && -f "$R2_KNEADDATA_PAIRED_GZ" ]]; then
        INPUT_FILE="$R1_KNEADDATA_PAIRED_GZ"
        INPUT_FILE2="$R2_KNEADDATA_PAIRED_GZ"
        echo "Using compressed KneadData paired output files"
    elif [[ -f "$R1_KNEADDATA_PAIRED" && -f "$R2_KNEADDATA_PAIRED" ]]; then
        INPUT_FILE="$R1_KNEADDATA_PAIRED"
        INPUT_FILE2="$R2_KNEADDATA_PAIRED"
        echo "Using uncompressed KneadData paired output files"
    elif [[ -f "$R1_KNEADDATA_TRIMMED_GZ" && -f "$R2_KNEADDATA_TRIMMED_GZ" ]]; then
        INPUT_FILE="$R1_KNEADDATA_TRIMMED_GZ"
        INPUT_FILE2="$R2_KNEADDATA_TRIMMED_GZ"
        echo "Using compressed KneadData trimmed output files"
    elif [[ -f "$R1_KNEADDATA_TRIMMED" && -f "$R2_KNEADDATA_TRIMMED" ]]; then
        INPUT_FILE="$R1_KNEADDATA_TRIMMED"
        INPUT_FILE2="$R2_KNEADDATA_TRIMMED"
        echo "Using uncompressed KneadData trimmed output files"
    else
        echo "ERROR: KneadData files missing for ${sample}"
        FAILED=$((FAILED + 1))
        continue
    fi

    CONCAT_INPUT="${HUMANN_RESULTS_DIR}/${sample}/${sample}_concat.fastq"

    if [[ "${INPUT_FILE}" == *.gz ]]; then
        zcat "${INPUT_FILE}" "${INPUT_FILE2}" > "${CONCAT_INPUT}"
    else
        cat "${INPUT_FILE}" "${INPUT_FILE2}" > "${CONCAT_INPUT}"
    fi

    LOG_FILE="${LOGS_DIR}/humann/${sample}/humann3.log"

    HUMANN_CMD_FULL="${HUMANN_CMD} \
        --input \"${CONCAT_INPUT}\" \
        --output \"${HUMANN_RESULTS_DIR}/${sample}\" \
        --nucleotide-database \"$CHOCOPHLAN_DB\" \
        --protein-database \"$UNIREF_DB\" \
        --threads ${HUMANN_THREADS} \
        --remove-temp-output \
        --verbose"

    if [[ "$USE_METAPHLAN_BYPASS" == true ]]; then
        METAPHLAN_FILE="${METAPHLAN_RESULTS_DIR}/${sample}/${sample}_profile.txt"
        if [[ -f "$METAPHLAN_FILE" ]]; then
            HUMANN_CMD_FULL="${HUMANN_CMD_FULL} --taxonomic-profile \"$METAPHLAN_FILE\""
        fi
    fi

    eval $HUMANN_CMD_FULL > "$LOG_FILE" 2>&1
    HUMANN_EXIT=$?

    rm -f "${CONCAT_INPUT}"

    if [[ $HUMANN_EXIT -eq 0 ]] && \
       [[ -f "${HUMANN_RESULTS_DIR}/${sample}/${sample}_concat_genefamilies.tsv" ]]; then

        mv "${HUMANN_RESULTS_DIR}/${sample}/${sample}_concat_genefamilies.tsv" \
           "${HUMANN_RESULTS_DIR}/${sample}/${sample}_genefamilies.tsv"
        mv "${HUMANN_RESULTS_DIR}/${sample}/${sample}_concat_pathcoverage.tsv" \
           "${HUMANN_RESULTS_DIR}/${sample}/${sample}_pathcoverage.tsv"
        mv "${HUMANN_RESULTS_DIR}/${sample}/${sample}_concat_pathabundance.tsv" \
           "${HUMANN_RESULTS_DIR}/${sample}/${sample}_pathabundance.tsv"

        echo "${sample}" >> "$CHECKPOINT_FILE"

        # Append to metadata file if this sample is not already present
        if ! grep -q "^${sample}\t" "$METADATA_FILE"; then
            echo -e "${sample}\t${SAMPLE_GROUP_MAP[$sample]}\t${SAMPLE_TIMEPOINT_MAP[$sample]}\t${SAMPLE_PATIENT_MAP[$sample]}\t${SAMPLE_DESCRIPTION_MAP[$sample]}" >> "$METADATA_FILE"
            echo "Added ${sample} to sample metadata"
        fi

        SUCCESSFUL=$((SUCCESSFUL + 1))
        echo "Sample ${sample} completed successfully"
    else
        echo "ERROR: HUMAnN3 failed for sample: ${sample}"
        FAILED=$((FAILED + 1))
    fi

    echo "Progress: ${PROCESSED}/${TOTAL_SAMPLES} (Successful: ${SUCCESSFUL}, Failed: ${FAILED}, Skipped: ${SKIPPED})"
done

echo ""
echo "PART 2: MERGING RESULTS"

if [[ ${SUCCESSFUL} -gt 0 ]]; then
    MERGE_NEEDED=true
    if [[ "$FRESH_START" == false ]] && [[ -f "${HUMANN_RESULTS_DIR}/merged_genefamilies.tsv" ]]; then
        echo "Merged files exist - checking if update needed..."
        NEWEST_INDIVIDUAL=$(find "${HUMANN_RESULTS_DIR}" -name "*_genefamilies.tsv" -not -name "merged_*" -newer "${HUMANN_RESULTS_DIR}/merged_genefamilies.tsv" | wc -l)

        if [[ ${NEWEST_INDIVIDUAL} -eq 0 ]]; then
            echo "Merged files are up to date - skipping merge step"
            MERGE_NEEDED=false
        else
            echo "New individual results found - re-merging..."
        fi
    fi

    if [[ "$MERGE_NEEDED" == true ]]; then
        echo "Merging HUMAnN3 results..."

        if command -v humann_join_tables &> /dev/null; then
            HUMANN_JOIN_CMD="humann_join_tables"
        elif command -v humann3_join_tables &> /dev/null; then
            HUMANN_JOIN_CMD="humann3_join_tables"
        elif [[ -f "${CONDA_PREFIX}/bin/humann_join_tables" ]]; then
            HUMANN_JOIN_CMD="${CONDA_PREFIX}/bin/humann_join_tables"
        elif [[ -f "${CONDA_PREFIX}/bin/humann3_join_tables" ]]; then
            HUMANN_JOIN_CMD="${CONDA_PREFIX}/bin/humann3_join_tables"
        else
            echo "ERROR: humann_join_tables not found"
            exit 1
        fi

        ${HUMANN_JOIN_CMD} \
            --input "${HUMANN_RESULTS_DIR}" \
            --file_name "genefamilies.tsv" \
            --output "${HUMANN_RESULTS_DIR}/merged_genefamilies.tsv" \
            --search-subdirectories

        ${HUMANN_JOIN_CMD} \
            --input "${HUMANN_RESULTS_DIR}" \
            --file_name "pathcoverage.tsv" \
            --output "${HUMANN_RESULTS_DIR}/merged_pathcoverage.tsv" \
            --search-subdirectories

        ${HUMANN_JOIN_CMD} \
            --input "${HUMANN_RESULTS_DIR}" \
            --file_name "pathabundance.tsv" \
            --output "${HUMANN_RESULTS_DIR}/merged_pathabundance.tsv" \
	    --search-subdirectories

        echo "Merging complete"
    fi
else
    echo "ERROR: No successful samples to merge"
    exit 1
fi

echo ""
echo "PART 3: EUT SCREENING"

EUT_SCREEN_DIR="${RESULTS_DIR}/eut_screening"
LOG_FILE="${LOGS_DIR}/eut_screening/eut_screening.log"

EUT_EXTRACTION_NEEDED=true
if [[ "$FRESH_START" == false ]] && \
   [[ -f "${EUT_SCREEN_DIR}/eut_genes_raw.tsv" ]] && \
   [[ -f "${EUT_SCREEN_DIR}/sample_groups.tsv" ]]; then
    if [[ "${EUT_SCREEN_DIR}/eut_genes_raw.tsv" -nt "${HUMANN_RESULTS_DIR}/merged_genefamilies.tsv" ]]; then
        echo "EUT extraction already completed and up to date - skipping"
        EUT_EXTRACTION_NEEDED=false
    else
        echo "Merged files updated - re-extracting EUT genes..."
    fi
fi

if [[ "$EUT_EXTRACTION_NEEDED" == true ]]; then
    echo "Extracting EUT genes..."

    EUT_GENES=("eutA" "eutB" "eutC" "eutD" "eutE" "eutG" "eutH" "eutJ" "eutK" "eutL" "eutM" "eutN" "eutP" "eutQ" "eutR" "eutS" "eutT")
    EUT_PATTERN=$(IFS='|'; echo "${EUT_GENES[*]}")

    grep -iE "(${EUT_PATTERN})" "${HUMANN_RESULTS_DIR}/merged_genefamilies.tsv" > "${EUT_SCREEN_DIR}/eut_genes_raw.tsv" 2>/dev/null

    EUT_GENE_COUNT=$(wc -l < "${EUT_SCREEN_DIR}/eut_genes_raw.tsv" 2>/dev/null || echo "0")

    if [[ "$EUT_GENE_COUNT" -eq 0 ]]; then
        grep -iE "(ethanolamine|eut)" "${HUMANN_RESULTS_DIR}/merged_genefamilies.tsv" > "${EUT_SCREEN_DIR}/eut_genes_raw.tsv" 2>/dev/null
        EUT_GENE_COUNT=$(wc -l < "${EUT_SCREEN_DIR}/eut_genes_raw.tsv" 2>/dev/null || echo "0")

        if [[ "$EUT_GENE_COUNT" -eq 0 ]]; then
            touch "${EUT_SCREEN_DIR}/eut_genes_raw.tsv"
        fi
    fi

    echo "Found ${EUT_GENE_COUNT} EUT-related gene entries"
fi

echo "EUT extraction completed. Run the separate Python plotting script for detailed results."

echo ""
echo "=================================================================="
echo "ANALYSIS COMPLETED"
echo "End time: $(date)"
echo "=================================================================="
echo ""
echo "PROCESSING SUMMARY:"
echo "Total samples: ${TOTAL_SAMPLES}"
echo "Successfully completed: ${SUCCESSFUL}"
echo "Failed: ${FAILED}"
echo "Skipped (already done): ${SKIPPED}"
echo ""
echo "CHECKPOINT INFO:"
echo "Completed samples logged in: ${CHECKPOINT_FILE}"
echo "Resume: Run script again without 'fresh_start' to continue from checkpoint"
echo "Fresh start: Run with 'fresh_start' argument to reprocess all samples"
echo ""
echo "RESULTS LOCATIONS:"
echo "- HUMAnN3 results: ${HUMANN_RESULTS_DIR}/"
echo "- Merged gene families: ${HUMANN_RESULTS_DIR}/merged_genefamilies.tsv"
echo "- Merged pathways: ${HUMANN_RESULTS_DIR}/merged_pathabundance.tsv"
echo "- EUT screening: ${EUT_SCREEN_DIR}/"
echo "- EUT genes: ${EUT_SCREEN_DIR}/eut_genes_raw.tsv"
echo "- Sample metadata: ${METADATA_FILE}"
echo ""
