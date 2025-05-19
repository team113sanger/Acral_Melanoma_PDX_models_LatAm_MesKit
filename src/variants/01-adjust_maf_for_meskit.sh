#!/bin/bash

####################
## DEFINING PATHS ##
####################

DATA_DIR="data/variants"
MAF="${DATA_DIR}/6633_2729_3248-filtered_mutations_matched_allTum_keep.maf"

##############################
# GENERATE ADJUSTED MAF FILE #
##############################

# Determine output file name.
ADJUSTED_MAF="${DATA_DIR}/6633_2729_3248-filtered_mutations_matched_allTum_keep_adjusted_for_meskit.maf"

# Get cases from conbined MAF file for which there are at least 2 samples.
SAMPLES=$(cat "${MAF}" | awk '{print $11}' | sort -u | cut -c 1-7 | sort | uniq -d)

# Filter MAF to keep only cases with at least 2 samples.
cat "${MAF}" | grep -e 'Tumor_Sample_Barcode' -e "${SAMPLES}" > "${ADJUSTED_MAF}"

# Adjust columns for MesKit.
sed -i 's/VAF_tum/VAF/g' "${ADJUSTED_MAF}"
sed -i 's/t_ref_count/Ref_allele_depth/g' "${ADJUSTED_MAF}"
sed -i 's/t_alt_count/Alt_allele_depth/g' "${ADJUSTED_MAF}"

##########################
# GENERATE CLINICAL FILE #
##########################

CLINICAL="${DATA_DIR}/6633_2729_3248-filtered_mutations_matched_allTum_keep_clinical_for_meskit.tsv"

# Get tumour sample barcodes from MAF file.
TUMOUR_BARCODES=$(cat "${ADJUSTED_MAF}" | awk '{print $11}' | grep "PD" | sort -u)

# Get patient IDs.
PATIENT_IDS=$(echo "${TUMOUR_BARCODES}" | cut -c 1-7)

# Get tumour IDs.
TUMOUR_IDS=$(echo "${TUMOUR_BARCODES}" | cut -c 8)

# Get tumour sample labels.
TUMOUR_LABELS=$(echo "${TUMOUR_BARCODES}" | cut -c 8-12)

echo -e "Tumor_Sample_Barcode\tTumor_ID\tPatient_ID\tTumor_Sample_Label" > ${CLINICAL}

paste <(echo "${TUMOUR_BARCODES}") \
    <(echo "${TUMOUR_IDS}") \
    <(echo "${PATIENT_IDS}") \
    <(echo "${TUMOUR_LABELS}") \
    --delimiters '\t' >> ${CLINICAL}
