#!/bin/bash

set -e

###########################
## PARSE INPUT ARGUMENTS ##
###########################

while getopts m:i:n:j:o: FLAG
do
    case "${FLAG}" in
        # Path to manifest file 1.
        m) MANIFEST1=${OPTARG};;
        # Identifier for manifest 1.
        i) MANIFEST1_ID=${OPTARG};;
        # Path to manifest file 2.
        n) MANIFEST2=${OPTARG};;
        # Identifier for manifest 2.
        j) MANIFEST2_ID=${OPTARG};;
        # Output directory.
        o) OUTDIR=${OPTARG};;
    esac
done

###########################
## CHECK INPUT ARGUMENTS ##
###########################

# Iterate over manifest files...
for ARGUMENT in ${MANIFEST1} ${MANIFEST2}; do

    # Check if manifest files exist.
    if [ ! -f ${ARGUMENT} ]; then
        echo "ERROR: ${ARGUMENT} does not exist. Please \
              review your input."
        exit 1
    fi

    # Check if manifest files have really been provided.
    if [ -z ${ARGUMENT} ]; then
        echo "ERROR: at least one manifest has not been provided. Please \
              review your input."
        exit 1
    fi
done

# Iterate over manifest identifiers...
for ARGUMENT in ${MANIFEST1_ID} ${MANIFEST2_ID}; do

    # Check if manifest identifiers have really been provided.
    if [ -z ${ARGUMENT} ]; then
        echo "ERROR: at least one manifest identifier has not been provided. Please \
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

# Create new manifest file.
touch ${OUTDIR}/${MANIFEST1_ID}_${MANIFEST2_ID}_combined_manifest.txt

# Print content from MANIFEST1 to new manifest file.
cat ${MANIFEST1} >> ${OUTDIR}/${MANIFEST1_ID}_${MANIFEST2_ID}_combined_manifest.txt

echo -n "\n" >> ${OUTDIR}/${MANIFEST1_ID}_${MANIFEST2_ID}_combined_manifest.txt

# Print content from MANIFEST2 (except first line, i.e. header) to new manifest file.
cat ${MANIFEST2} | tail -n +2 >> ${OUTDIR}/${MANIFEST1_ID}_${MANIFEST2_ID}_combined_manifest.txt