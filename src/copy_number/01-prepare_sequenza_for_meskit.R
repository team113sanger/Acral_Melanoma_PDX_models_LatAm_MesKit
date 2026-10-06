library(here)
library(MesKit)
library(ggpubr)
library(dplyr)
library(tidyr)
library(stringr)

########################
#### Defining paths ####
########################

data_dir <- here::here("data/copy_number/")
segments_files <- fs::dir_ls(paste0(data_dir, "segments/"), recurse = TRUE, glob = "*_segments.txt$")
ploidy_file <- paste0(data_dir, "Adjusted_ploidy_table.tsv")

mdata_dir <- here::here("metadata/")
mdata_file <- paste0(mdata_dir, "6633_2729_3248_METADATA_PDX_from_Latin_America_WES.txt")

#####################################
#### Loading and processing data ####
#####################################

adjust_copynumber <- function(segments, ploidy_table) {
    segments <- segments |>
        dplyr::rename(Sample = Tumor_Sample_Barcode) |>
        dplyr::left_join(
            ploidy_table,
            by = dplyr::join_by(Sample)
        ) |>
        dplyr::mutate(
            CopyNumber_adjusted = dplyr::case_when(
                CopyNumber == 0 ~ 0, # Deletion
                CopyNumber != 0 & CopyNumber < mean_ploidy_adj ~ 1, #  Loss
                CopyNumber == mean_ploidy_adj ~ 2, # Neutral
                CopyNumber > mean_ploidy_adj & CopyNumber < (2 * mean_ploidy_adj) ~ 3, # Gain
                CopyNumber >= (2 * mean_ploidy_adj) ~ 4, #  Amplification
            ),
            .before = CopyNumber
        ) |>
        dplyr::select(!c(CopyNumber, mean_ploidy_adj)) |>
        dplyr::rename(CopyNumber = CopyNumber_adjusted) |>
        dplyr::rename(Tumor_Sample_Barcode = Sample)

    return(segments)
}

load_copynumber <- function(files, metadata, ploidy_table) {
    segments <- readr::read_tsv(files, col_names = TRUE, id = "Tumor_Sample_Barcode") |>
        dplyr::mutate(
            Tumor_Sample_Barcode = stringr::str_remove(basename(Tumor_Sample_Barcode), "_segments.txt"),
            Patient_ID = substring(Tumor_Sample_Barcode, 1, 7)
        ) |>
        dplyr::rename(
            Chromosome = chromosome,
            Start_Position = start.pos,
            End_Position = end.pos,
            CopyNumber = CNt,
            Major_CN = A,
            Minor_CN = B
        ) |>
        dplyr::select(
            Patient_ID,
            Tumor_Sample_Barcode,
            Chromosome,
            Start_Position,
            End_Position,
            CopyNumber,
            Major_CN,
            Minor_CN
        )

    segments <- adjust_copynumber(segments, ploidy_table)

    segments$Patient_ID <- metadata$patient_id[match(
        segments$Tumor_Sample_Barcode, metadata$final_sample_name_used
    )]

    segments$Tumor_Sample_Barcode <- metadata$`Case ID`[match(
        segments$Tumor_Sample_Barcode, metadata$final_sample_name_used
    )]
    segments$Tumor_Sample_Barcode <- gsub("_gDNA_tumour", "", segments$Tumor_Sample_Barcode)

    return(segments)
}

ploidy <- readr::read_table(ploidy_file)

metadata <- read.csv(mdata_file, sep = "\t", check.names = FALSE)
metadata$patient_id <- stringr::str_extract(metadata$`Case ID`, "^AM[0-9]+[ab]?")

#  Load segments.
segments <- load_copynumber(
    files = segments_files,
    metadata = metadata,
    ploidy_table = ploidy
) |> tidyr::drop_na()

write.table(
    segments,
    file = paste0(data_dir, "segments_for_meskit.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)
