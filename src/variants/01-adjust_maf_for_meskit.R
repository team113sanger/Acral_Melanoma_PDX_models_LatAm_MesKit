#!/usr/bin/env Rscript

library(here)
library(dplyr)
library(stringr)
library(logger)

logger::log_threshold(logger::INFO)

########################
#### Defining paths ####
########################

data_dir <- here::here("data/")
maf_file <- paste0(data_dir, "variants/6633_2729_3248-filtered_mutations_matched_allTum_keep.maf")

mdata_dir <- here::here("metadata/")
mdata_file <- paste0(mdata_dir, "6633_2729_3248_METADATA_PDX_from_Latin_America_WES.txt")

######################
#### Loading data ####
######################

logger::log_info("Reading MAF...")
maf <- read.csv(maf_file, sep = "\t", check.names = FALSE)

logger::log_info("Loading metadata...")
metadata <- read.csv(mdata_file, sep = "\t", check.names = FALSE)
metadata$patient_id <- stringr::str_extract(metadata$`Case ID`, "^AM[0-9]+[ab]?")

####################################
#### Adjust MAF file for MesKit ####
####################################

# Add sample IDs.
maf$sample_id <- metadata$`Case ID`[match(
    maf$Tumor_Sample_Barcode, metadata$final_sample_name_used
)]

# Add patient IDs.
maf$patient_id <- metadata$patient_id[match(
    maf$Tumor_Sample_Barcode, metadata$final_sample_name_used
)]

logger::log_info("Adjusting MAF...")

adjusted_maf_path <- paste0(
    data_dir,
    "variants/6633_2729_3248-filtered_mutations_matched_allTum_keep_adjusted_for_meskit.maf"
)

#  Get IDs of patients for whom there are at least 2 samples available.
patients_with_at_least_two_samples <- maf |>
    dplyr::group_by(patient_id) |>
    dplyr::summarise(n_samples = n_distinct(sample_id)) |>
    dplyr::filter(n_samples >= 2) |>
    dplyr::select(patient_id) |>
    dplyr::pull()

maf <- maf |>
    # Replace PD IDs with sample IDs.
    dplyr::select(!Tumor_Sample_Barcode) |>
    dplyr::rename(Tumor_Sample_Barcode = sample_id) |>
    # Keep only patients for whom there are at least 2 samples.
    dplyr::filter(patient_id %in% patients_with_at_least_two_samples) |>
    # Remove the column with patient IDs.
    dplyr::select(!patient_id) |>
    # Rename columns.
    dplyr::rename(VAF = VAF_tum) |>
    dplyr::rename(Ref_allele_depth = t_ref_count) |>
    dplyr::rename(Alt_allele_depth = t_alt_count) |>
    # Tidy up sample IDs.
    dplyr::mutate(Tumor_Sample_Barcode = gsub("_gDNA_tumour", "", Tumor_Sample_Barcode)) |>
    # Save to output file.
    readr::write_delim(file = adjusted_maf_path, delim = "\t")


logger::log_info(paste0("Adjusted MAF saved to ", adjusted_maf_path))

###########################################
#### Generate clinical file for MesKit ####
###########################################

logger::log_info("Generating clinical file...")

clinical <- maf |>
    dplyr::select(Tumor_Sample_Barcode) |>
    dplyr::filter(!duplicated(Tumor_Sample_Barcode))

# Create columns.
clinical$Tumor_ID <- clinical$Tumor_Sample_Barcode
clinical$Patient_ID <- stringr::str_extract(clinical$Tumor_Sample_Barcode, "^AM[0-9]+[ab]?")
clinical$Tumor_Sample_Label <- clinical$Tumor_Sample_Barcode

clinical_file_path <- paste0(
    data_dir,
    "variants/6633_2729_3248-filtered_mutations_matched_allTum_keep_clinical_for_meskit.tsv"
)

clinical |>
    dplyr::arrange(Patient_ID, Tumor_Sample_Barcode) |>
    readr::write_delim(file = clinical_file_path, delim = "\t")

logger::log_info(paste0("Clinical file saved to ", clinical_file_path))
