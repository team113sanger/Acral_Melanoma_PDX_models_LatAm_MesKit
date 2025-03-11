#!/bin/bash

set -e

###########################
## PARSE INPUT ARGUMENTS ##
###########################

while getopts m:i:n:j:o: FLAG
do
    case "${FLAG}" in
        # Path to MAF file 1.
        m) MAF1=${OPTARG};;
        # Identifier for MAF 1.
        i) MAF1_ID=${OPTARG};;
        # Path to MAF file 2.
        n) MAF2=${OPTARG};;
        # Identifier for MAF 2.
        j) MAF2_ID=${OPTARG};;
        # Output directory.
        o) OUTDIR=${OPTARG};;
    esac
done

###########################
## CHECK INPUT ARGUMENTS ##
###########################

# Iterate over MAF files...
for ARGUMENT in ${MAF1} ${MAF2}; do

    # Check if MAF files exist.
    if [ ! -f ${ARGUMENT} ]; then
        echo "ERROR: ${ARGUMENT} does not exist. Please \
              review your input."
        exit 1
    fi

    # Check if MAF files have really been provided.
    if [ -z ${ARGUMENT} ]; then
        echo "ERROR: at least one MAF has not been provided. Please \
              review your input."
        exit 1
    fi
done

# Iterate over MAF identifiers...
for ARGUMENT in ${MAF1_ID} ${MAF2_ID}; do

    # Check if MAF identifiers have really been provided.
    if [ -z ${ARGUMENT} ]; then
        echo "ERROR: at least one MAF identifier has not been provided. Please \
              review your input."
        exit 1
    fi
done

# Check if the output directory exists.
if [ ! -d ${OUTDIR} ]; then
    echo "ERROR: ${OUTDIR} does not exist. Please \
          review your input."
    exit 1
fi

# Check if an output directory has been provided.
if [ -z ${OUTDIR} ]; then
    echo "ERROR: no output directory has been provided. Please \
          review your input."
    exit 1
fi

##############################
# GENERATE COMBINED MAF FILE #
##############################

# Determine output file names.
COMBINED_MAF="${OUTDIR}/${MAF1_ID}_${MAF2_ID}_combined_intermediate_file.maf"
FILTERED_COMBINED_MAF="${OUTDIR}/${MAF1_ID}_${MAF2_ID}_combined_filtered_for_meskit.maf"

# Print content from MAF1 to new MAF file.
cat ${MAF1} > ${COMBINED_MAF}

# Print content from MAF2 (except first line, i.e. header) to new MAF file.
cat ${MAF2} | tail -n +2 >> ${COMBINED_MAF}

# Get cases from conbined MAF file for which there are at least 2 samples.
SAMPLES=$(cat ${COMBINED_MAF} | awk '{print $11}' | sort -u | cut -c 1-7 | sort | uniq -d)

# Filter MAF to keep only cases with at least 2 samples.
cat ${COMBINED_MAF} | grep -e 'Tumor_Sample_Barcode' -e "${SAMPLES}" > ${FILTERED_COMBINED_MAF}

# Adjust columns for MesKit.
sed -i 's/VAF_tum/VAF/g' ${FILTERED_COMBINED_MAF}
sed -i 's/t_ref_count/Ref_allele_depth/g' ${FILTERED_COMBINED_MAF}
sed -i 's/t_alt_count/Alt_allele_depth/g' ${FILTERED_COMBINED_MAF}

rm -f ${COMBINED_MAF}

###################################
# GENERATE COMBINED METADATA FILE #
###################################

COMBINED_MDATA="${OUTDIR}/${MAF1_ID}_${MAF2_ID}_metadata_for_meskit.tsv"

# Get tumour sample barcodes from MAF file.
TUMOUR_BARCODES=$(cat ${FILTERED_COMBINED_MAF} | awk '{print $11}' | grep "PD" | sort -u)

# Get patient IDs.
PATIENT_IDS=$(echo "${TUMOUR_BARCODES}" | cut -c 1-7)

# Get tumour IDs.
TUMOUR_IDS=$(echo "${TUMOUR_BARCODES}" | cut -c 8)

# Get tumour sample labels.
TUMOUR_LABELS=$(echo "${TUMOUR_BARCODES}" | cut -c 8-12)

echo -e "Tumor_Sample_Barcode\tTumor_ID\tPatient_ID\tTumor_Sample_Label" > ${COMBINED_MDATA}

paste <(echo "${TUMOUR_BARCODES}") \
    <(echo "${TUMOUR_IDS}") \
    <(echo "${PATIENT_IDS}") \
    <(echo "${TUMOUR_LABELS}") \
    --delimiters '\t' >> ${COMBINED_MDATA}