#!/bin/bash

# Steps 3-5: Normalisation, annotation, and EUT gene identification

source ~/projects/metagenomics/scripts/config.sh

HUMANN_RESULTS_DIR="${RESULTS_DIR}/humann"
EUT_SCREEN_DIR="${RESULTS_DIR}/eut_screening"

echo "=================================================================="
echo "CORRECT HUMANN EUT GENE ANALYSIS PROTOCOL"
echo "Start time: $(date)"
echo "=================================================================="

cd "$HUMANN_RESULTS_DIR"

# Check if mapping databases are properly installed
echo "Verifying utility mapping database installation..."
if humann_rename_table --help | grep -q "uniref90"; then
    echo "UniRef90 naming is available"
else
    echo "ERROR: UniRef90 naming still not available"
    echo "Please check the utility mapping database installation"
    exit 1
fi

echo ""
echo "STEP 3: NORMALISING GENE FAMILY ABUNDANCES"
echo "Converting from RPK to copies per million (CPM) units..."

humann_renorm_table \
    --input "merged_genefamilies.tsv" \
    --units cpm \
    --output "merged_genefamilies_cpm.tsv"

if [[ ! -f "merged_genefamilies_cpm.tsv" ]]; then
    echo "ERROR: Normalisation failed"
    exit 1
fi

NORM_ROWS=$(wc -l < "merged_genefamilies_cpm.tsv")
NORM_COLS=$(head -1 "merged_genefamilies_cpm.tsv" | tr '\t' '\n' | wc -l)
echo "Normalised file created: ${NORM_ROWS} rows, ${NORM_COLS} columns"

echo ""
echo "STEP 4: ADDING GENE ANNOTATIONS"
echo "Attaching UniRef90 gene names to normalised data..."

humann_rename_table \
    --input "merged_genefamilies_cpm.tsv" \
    --names uniref90 \
    --output "merged_genefamilies_cpm_annotated.tsv"

if [[ ! -f "merged_genefamilies_cpm_annotated.tsv" ]]; then
    echo "ERROR: Gene annotation failed"
    exit 1
fi

ANNOT_ROWS=$(wc -l < "merged_genefamilies_cpm_annotated.tsv")
echo "Annotated file created: ${ANNOT_ROWS} rows"

# Check annotation success rate
NAMED_GENES=$(grep -v "NO_NAME" "merged_genefamilies_cpm_annotated.tsv" | grep -v "^#" | wc -l)
TOTAL_GENES=$(grep -v "^#" "merged_genefamilies_cpm_annotated.tsv" | wc -l)
ANNOTATION_RATE=$(echo "scale=1; $NAMED_GENES * 100 / $TOTAL_GENES" | bc)
echo "Annotation success: ${NAMED_GENES}/${TOTAL_GENES} genes (${ANNOTATION_RATE}%)"

echo ""
echo "STEP 4b: CLEANING COLUMN NAMES"
echo "Stripping _concat_Abundance-RPKs suffix from sample column headers..."

# Strip HUMAnN-appended suffix from column names so they match metadata sample IDs
head -1 "merged_genefamilies_cpm_annotated.tsv" | \
    sed 's/_concat_Abundance-RPKs//g' > /tmp/clean_header.txt
tail -n +2 "merged_genefamilies_cpm_annotated.tsv" | \
    cat /tmp/clean_header.txt - > "merged_genefamilies_cpm_annotated_clean.tsv"
rm -f /tmp/clean_header.txt

echo "Cleaned header columns:"
head -1 "merged_genefamilies_cpm_annotated_clean.tsv" | tr '\t' '\n' | head -5

echo ""
echo "STEP 5: SEARCHING FOR EUT GENES"
echo "Using correct functional terminology for ethanolamine utilisation..."

mkdir -p "$EUT_SCREEN_DIR"

echo "APPROACH 1: HIGHLY SPECIFIC FUNCTIONAL SEARCH"

SPECIFIC_EUT_TERMS=(
    "ethanolamine ammonia-lyase"
    "ethanolamine kinase"
    "phosphoethanolamine cytidylyltransferase"
    "ethanolamine utilization protein"
    "ethanolamine phosphate cytidylyltransferase"
    "ethanolamine ammonia lyase"
    "microcompartment shell protein"
    "carboxysome shell protein"
)

echo "Searching with highly specific terms:"
for term in "${SPECIFIC_EUT_TERMS[@]}"; do
    echo "  - $term"
done

# Write header first, then search
head -1 "merged_genefamilies_cpm_annotated_clean.tsv" > "${EUT_SCREEN_DIR}/eut_specific_functional.tsv"

for term in "${SPECIFIC_EUT_TERMS[@]}"; do
    grep -i "$term" "merged_genefamilies_cpm_annotated_clean.tsv" | \
    grep -v "Coprococcus_eutactus" >> "${EUT_SCREEN_DIR}/eut_specific_functional.tsv" 2>/dev/null
done

# Remove duplicate data rows but preserve header
HEADER_LINE=$(head -1 "${EUT_SCREEN_DIR}/eut_specific_functional.tsv")
tail -n +2 "${EUT_SCREEN_DIR}/eut_specific_functional.tsv" | \
    sort -u | \
    cat <(echo "$HEADER_LINE") - > "${EUT_SCREEN_DIR}/eut_specific_functional_clean.tsv"
mv "${EUT_SCREEN_DIR}/eut_specific_functional_clean.tsv" "${EUT_SCREEN_DIR}/eut_specific_functional.tsv"

SPECIFIC_COUNT=$(tail -n +2 "${EUT_SCREEN_DIR}/eut_specific_functional.tsv" | wc -l)
echo "Found ${SPECIFIC_COUNT} specific EUT pathway genes (excluding C. eutactus)"

echo ""
echo "APPROACH 2: MANUAL EUT GENE IDENTIFICATION"

head -1 "merged_genefamilies_cpm_annotated_clean.tsv" > "${EUT_SCREEN_DIR}/eut_named_proteins.tsv"
grep -E "ethanolamine utilization protein Eut[A-T]|EutM|EutE|EutB|EutC|EutD|EutA" \
    "merged_genefamilies_cpm_annotated_clean.tsv" | \
    grep -v "Coprococcus_eutactus" >> "${EUT_SCREEN_DIR}/eut_named_proteins.tsv" 2>/dev/null

NAMED_COUNT=$(tail -n +2 "${EUT_SCREEN_DIR}/eut_named_proteins.tsv" | wc -l)
echo "Found ${NAMED_COUNT} genes with explicit EUT protein names"

echo ""
echo "ALTERNATIVE APPROACH: REGROUPING TO KEGG ORTHOLOGY"

humann_regroup_table \
    --input "merged_genefamilies_cpm.tsv" \
    --groups uniref90_ko \
    --output "merged_genefamilies_cpm_ko.tsv"

if [[ -f "merged_genefamilies_cpm_ko.tsv" ]]; then
    echo "KEGG KO regrouping successful"

    # Clean KEGG column names too
    head -1 "merged_genefamilies_cpm_ko.tsv" | \
        sed 's/_concat_Abundance-RPKs//g' > /tmp/clean_ko_header.txt
    tail -n +2 "merged_genefamilies_cpm_ko.tsv" | \
        cat /tmp/clean_ko_header.txt - > "merged_genefamilies_cpm_ko_clean.tsv"
    rm -f /tmp/clean_ko_header.txt

    echo ""
    echo "APPROACH 3: KEGG-SPECIFIC EUT GENES"

    EUT_KEGG_KOS=(
        "K04020"  # ethanolamine ammonia-lyase large subunit
        "K03735"  # ethanolamine ammonia-lyase small subunit
        "K15732"  # ethanolamine kinase
        "K04021"  # ethanolamine utilization protein EutE
        "K04022"  # ethanolamine utilization protein EutG
        "K04023"  # ethanolamine utilization protein EutH
        "K04024"  # ethanolamine utilization protein EutJ
        "K04025"  # ethanolamine utilization protein EutK
        "K04026"  # ethanolamine utilization protein EutL
        "K04027"  # ethanolamine utilization protein EutM
        "K04028"  # ethanolamine utilization protein EutN
        "K04029"  # ethanolamine utilization protein EutP
        "K04030"  # ethanolamine utilization protein EutQ
        "K04031"  # ethanolamine utilization protein EutR
        "K04032"  # ethanolamine utilization protein EutS
        "K04033"  # ethanolamine utilization protein EutT
    )

    KEGG_PATTERN=$(IFS='|'; echo "${EUT_KEGG_KOS[*]}")

    head -1 "merged_genefamilies_cpm_ko_clean.tsv" > "${EUT_SCREEN_DIR}/eut_kegg_specific.tsv"
    grep -E "(${KEGG_PATTERN})" "merged_genefamilies_cpm_ko_clean.tsv" | \
        grep -v "Coprococcus_eutactus" >> "${EUT_SCREEN_DIR}/eut_kegg_specific.tsv" 2>/dev/null

    KEGG_SPECIFIC_COUNT=$(tail -n +2 "${EUT_SCREEN_DIR}/eut_kegg_specific.tsv" | wc -l)
    echo "Found ${KEGG_SPECIFIC_COUNT} genes with specific EUT KEGG KO identifiers"
else
    echo "WARNING: KEGG KO regrouping failed"
    KEGG_SPECIFIC_COUNT=0
fi

echo ""
echo "RESULTS COMPARISON"
echo "=================="
echo "Specific functional search: ${SPECIFIC_COUNT} genes"
echo "Named EUT proteins: ${NAMED_COUNT} genes"
echo "KEGG-specific EUT KOs: ${KEGG_SPECIFIC_COUNT} genes"

# Determine best approach
BEST_COUNT=0
BEST_FILE=""
BEST_TYPE=""

if [[ "$SPECIFIC_COUNT" -gt "$BEST_COUNT" ]]; then
    BEST_COUNT="$SPECIFIC_COUNT"
    BEST_FILE="${EUT_SCREEN_DIR}/eut_specific_functional.tsv"
    BEST_TYPE="Specific functional search"
fi

if [[ "$NAMED_COUNT" -gt "$BEST_COUNT" ]]; then
    BEST_COUNT="$NAMED_COUNT"
    BEST_FILE="${EUT_SCREEN_DIR}/eut_named_proteins.tsv"
    BEST_TYPE="Named EUT proteins"
fi

if [[ "$KEGG_SPECIFIC_COUNT" -gt "$BEST_COUNT" ]]; then
    BEST_COUNT="$KEGG_SPECIFIC_COUNT"
    BEST_FILE="${EUT_SCREEN_DIR}/eut_kegg_specific.tsv"
    BEST_TYPE="KEGG-specific KOs"
fi

echo ""
echo "COMPREHENSIVE EUT GENE ANALYSIS"

if [[ "$BEST_COUNT" -gt 0 ]]; then
    echo "Using ${BEST_TYPE} approach (${BEST_COUNT} genes found)"

    echo -e "Gene_ID\tFunction\tTotal_Abundance\tDetected_Samples\tTop_Sample\tTop_Abundance" > "${EUT_SCREEN_DIR}/eut_comprehensive_analysis.tsv"

    echo "Analysing EUT gene abundance patterns..."

    HEADER=$(head -1 "$BEST_FILE")
    SAMPLE_COUNT=$(echo "$HEADER" | tr '\t' '\n' | wc -l)

    while IFS=$'\t' read -r gene_entry rest_of_line; do
        if [[ "$gene_entry" == \#* ]]; then
            continue
        fi

        FUNCTION=$(echo "$gene_entry" | sed 's/.*: //' | sed 's/|.*//')
        if [[ -z "$FUNCTION" || "$FUNCTION" == "$gene_entry" ]]; then
            FUNCTION="Unknown function"
        fi

        TOTAL_ABUNDANCE=$(echo "$rest_of_line" | tr '\t' '\n' | awk '{sum+=$1} END {print sum}')
        DETECTED_SAMPLES=$(echo "$rest_of_line" | tr '\t' '\n' | awk '$1>0{count++} END {print count}')

        MAX_ABUNDANCE=0
        TOP_SAMPLE="None"

        for ((i=2; i<=SAMPLE_COUNT; i++)); do
            SAMPLE_NAME=$(echo "$HEADER" | cut -f$i)
            ABUNDANCE=$(echo "$rest_of_line" | cut -f$((i-1)))

            if (( $(echo "$ABUNDANCE > $MAX_ABUNDANCE" | bc -l) )); then
                MAX_ABUNDANCE="$ABUNDANCE"
                TOP_SAMPLE="$SAMPLE_NAME"
            fi
        done

        echo -e "${gene_entry}\t${FUNCTION}\t${TOTAL_ABUNDANCE}\t${DETECTED_SAMPLES}\t${TOP_SAMPLE}\t${MAX_ABUNDANCE}" >> "${EUT_SCREEN_DIR}/eut_comprehensive_analysis.tsv"

    done < <(tail -n +2 "$BEST_FILE")

    echo ""
    echo "EUT GENE ANALYSIS RESULTS"
    echo "========================"

    echo "EUT genes found by functional description:"
    for term in "${SPECIFIC_EUT_TERMS[@]}"; do
        COUNT=$(grep -Fci "$term" "$BEST_FILE" | head -1 || echo "0")
        if [[ "$COUNT" -gt 0 ]]; then
            echo "  $term: $COUNT genes"
        fi
    done

    echo ""
    echo "Top 10 most abundant EUT genes:"
    tail -n +2 "${EUT_SCREEN_DIR}/eut_comprehensive_analysis.tsv" | sort -k3 -nr | head -10 | \
        awk -F'\t' '{printf "%-60s %10.2f CPM\n", substr($2, 1, 60), $3}'

    echo ""
    echo "EUT gene prevalence across samples:"
    tail -n +2 "${EUT_SCREEN_DIR}/eut_comprehensive_analysis.tsv" | \
        awk -F'\t' '{if($4>0) print $4" samples: "$2}' | sort -nr | head -10

    echo ""
    echo "SAMPLE-WISE EUT ABUNDANCE ANALYSIS"

    echo -e "Sample\tTotal_EUT_CPM\tEUT_Gene_Count\tTop_EUT_Gene" > "${EUT_SCREEN_DIR}/sample_eut_summary.tsv"

    for ((i=2; i<=SAMPLE_COUNT; i++)); do
        SAMPLE_NAME=$(echo "$HEADER" | cut -f$i)

        TOTAL_EUT=$(cut -f$i "$BEST_FILE" | tail -n +2 | awk '{sum+=$1} END {print sum}')
        EUT_COUNT=$(cut -f$i "$BEST_FILE" | tail -n +2 | awk '$1>0{count++} END {print count}')
        TOP_GENE=$(paste <(cut -f1 "$BEST_FILE" | tail -n +2) <(cut -f$i "$BEST_FILE" | tail -n +2) | \
                   sort -k2 -nr | head -1 | cut -f1 | sed 's/.*: //' | sed 's/|.*//')

        echo -e "${SAMPLE_NAME}\t${TOTAL_EUT}\t${EUT_COUNT}\t${TOP_GENE}" >> "${EUT_SCREEN_DIR}/sample_eut_summary.tsv"
    done

    echo ""
    echo "Top 10 samples by total EUT gene abundance:"
    tail -n +2 "${EUT_SCREEN_DIR}/sample_eut_summary.tsv" | sort -k2 -nr | head -10 | \
        awk -F'\t' '{printf "%-25s %10.2f CPM (%2d genes)\n", $1, $2, $3}'

    echo ""
    echo "SAMPLE GROUP ANALYSIS"
    echo "Analysing EUT abundance by experimental groups..."

    declare -A SAMPLE_GROUPS
    for sample in "${CONTROL_SAMPLES[@]}"; do
        SAMPLE_GROUPS["$sample"]="Control"
    done
    for sample in "${HEALTHY_SAMPLES[@]}"; do
        SAMPLE_GROUPS["$sample"]="Healthy"
    done
    for sample in "${EXT_HEALTHY_SAMPLES[@]}"; do
        SAMPLE_GROUPS["$sample"]="ExtHealthy"
    done
    for sample in "${TREATED_SAMPLES[@]}"; do
        SAMPLE_GROUPS["$sample"]="Treated"
    done

    echo -e "Group\tMean_EUT_CPM\tStd_Dev\tSample_Count\tSamples_With_EUT" > "${EUT_SCREEN_DIR}/group_eut_summary.tsv"

    for group in Control Healthy ExtHealthy Treated; do
        GROUP_ABUNDANCES=()
        SAMPLES_WITH_EUT=0
        TOTAL_SAMPLES=0

        while IFS=$'\t' read -r sample total_eut gene_count top_gene; do
            if [[ "${SAMPLE_GROUPS[$sample]}" == "$group" ]]; then
                GROUP_ABUNDANCES+=("$total_eut")
                TOTAL_SAMPLES=$((TOTAL_SAMPLES + 1))
                if (( $(echo "$total_eut > 0" | bc -l) )); then
                    SAMPLES_WITH_EUT=$((SAMPLES_WITH_EUT + 1))
                fi
            fi
        done < <(tail -n +2 "${EUT_SCREEN_DIR}/sample_eut_summary.tsv")

        if [[ ${#GROUP_ABUNDANCES[@]} -gt 0 ]]; then
            MEAN=$(printf '%s\n' "${GROUP_ABUNDANCES[@]}" | awk '{sum+=$1} END {print sum/NR}')
            STD_DEV=$(printf '%s\n' "${GROUP_ABUNDANCES[@]}" | awk -v mean="$MEAN" '{sum+=($1-mean)^2} END {print sqrt(sum/NR)}')
            echo -e "${group}\t${MEAN}\t${STD_DEV}\t${TOTAL_SAMPLES}\t${SAMPLES_WITH_EUT}" >> "${EUT_SCREEN_DIR}/group_eut_summary.tsv"
        fi
    done

    echo ""
    echo "EUT abundance by experimental group:"
    tail -n +2 "${EUT_SCREEN_DIR}/group_eut_summary.tsv" | \
        awk -F'\t' '{printf "%-12s: %8.2f +/- %6.2f CPM (%d/%d samples with EUT genes)\n", $1, $2, $3, $5, $4}'

else
    echo "WARNING: No EUT genes found using either functional annotation or KEGG KO approaches"
fi

echo ""
echo "=================================================================="
echo "EUT GENE ANALYSIS COMPLETE"
echo "End time: $(date)"
echo "=================================================================="
echo ""
echo "OUTPUT FILES CREATED:"
echo "- Normalised gene families: merged_genefamilies_cpm.tsv"
echo "- Annotated gene families: merged_genefamilies_cpm_annotated.tsv"
echo "- Annotated (clean headers): merged_genefamilies_cpm_annotated_clean.tsv"
echo "- KEGG KO regrouped: merged_genefamilies_cpm_ko.tsv"
echo "- KEGG KO (clean headers): merged_genefamilies_cpm_ko_clean.tsv"
echo "- EUT genes (functional): ${EUT_SCREEN_DIR}/eut_specific_functional.tsv"
echo "- EUT genes (named): ${EUT_SCREEN_DIR}/eut_named_proteins.tsv"
echo "- EUT genes (KEGG): ${EUT_SCREEN_DIR}/eut_kegg_specific.tsv"
echo "- Comprehensive analysis: ${EUT_SCREEN_DIR}/eut_comprehensive_analysis.tsv"
echo "- Sample EUT summary: ${EUT_SCREEN_DIR}/sample_eut_summary.tsv"
echo "- Group EUT summary: ${EUT_SCREEN_DIR}/group_eut_summary.tsv"
echo ""
echo "NEXT STEPS:"
echo "1. Review the comprehensive analysis to identify specific EUT genes"
echo "2. Use the sample and group summaries for statistical analysis"
echo "3. Run the Python plotting script for visualisation"
