#!/bin/bash

# Project directories

PROJECT_DIR=~/projects/metagenomics
RAW_DATA_DIR=${PROJECT_DIR}/raw_data
PROCESSED_DATA_DIR=${PROJECT_DIR}/processed_data
TAXONOMY_DIR=${PROJECT_DIR}/taxonomic_profiles
FUNCTIONAL_DIR=${PROJECT_DIR}/functional_analysis
EUT_DIR=${PROJECT_DIR}/eut_operon_analysis
SCRIPTS_DIR=${PROJECT_DIR}/scripts
LOGS_DIR=${PROJECT_DIR}/logs
RESULTS_DIR=${PROJECT_DIR}/results

# MetaPhlAn parameters
METAPHLAN_DB="/home/cclare/miniconda3/envs/metagenomics_env/lib/python3.9/site-packages/metaphlan/metaphlan_databases"
METAPHLAN_INDEX="mpa_vJun23_CHOCOPhlAnSGB_202403"
METAPHLAN_THREADS=8
METAPHLAN_ANALYSIS_TYPE="rel_ab_w_read_stats"

# Sample information
CONTROL_SAMPLES=("Negative_control_S32" "Positive_control_S31")
HEALTHY_SAMPLES=("H3_S25" "H4_S27" "H10_S29" "H11_S30")
EXT_HEALTHY_SAMPLES=("C1")
TREATED_SAMPLES=("10_1_S1" "10_7_S2" "12_1_S4" "12_3_S5" "14_1_S7" "14_5_S8" "27_1_S10" "27_3_S11" "43_1_S13" "43_5_S14" "45_1_S16" "45_6_S17" "46_1_S19" "46_4_S20" "53_1_S22" "53_8_S23")
#ALL_SAMPLES=("${CONTROL_SAMPLES[@]}" "${HEALTHY_SAMPLES[@]}" "${EXT_HEALTHY_SAMPLES[@]}" "${TREATED_SAMPLES[@]}")
ALL_SAMPLES=("${EXT_HEALTHY_SAMPLES[@]}")


# File extensions
R1_SUFFIX="_R1_001.fastq.gz"
R2_SUFFIX="_R2_001.fastq.gz"


