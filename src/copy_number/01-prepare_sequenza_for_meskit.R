library(here)
library(MesKit)
library(ggpubr)

########################
#### Defining paths ####
########################

data_dir <- here::here("data/copy_number/")

human_1st_model_files <- fs::dir_ls(paste0(data_dir, "segments/human/Sequenza_1stModels/"), recurse = TRUE, glob = "*_segments.txt$")
human_alt_model_files <- fs::dir_ls(paste0(data_dir, "segments/human/Sequenza_AltModels/"), recurse = TRUE, glob = "*_segments.txt$")

pdx_1st_model_files <- fs::dir_ls(paste0(data_dir, "segments/pdx/Sequenza_1stModels/"), recurse = TRUE, glob = "*_segments.txt$")
pdx_alt_model_files <- fs::dir_ls(paste0(data_dir, "segments/pdx/Sequenza_AltModels/"), recurse = TRUE, glob = "*_segments.txt$")

ploidy_file <- paste0(data_dir, "ploidy_table_jan2025_data.tsv")

#####################################
#### Loading and processing data ####
#####################################

adjust_copynumber <- function(segments) {
    segments <- segments |>
        dplyr::left_join(
            ploidy_table,
            by = dplyr::join_by(Tumor_Sample_Barcode)
        ) |>
        dplyr::mutate(Ploidy = round(Ploidy)) |>
        dplyr::mutate(
            CopyNumber_adjusted = dplyr::case_when(
                CopyNumber == 0 ~ 0, # Deletion
                CopyNumber != 0 & CopyNumber < Ploidy ~ 1, #  Loss
                CopyNumber == Ploidy ~ 2, # Neutral
                CopyNumber > Ploidy & CopyNumber < (2 * Ploidy) ~ 3, # Gain
                CopyNumber >= (2 * Ploidy) ~ 4, #  Amplification
            ),
            .before = CopyNumber
        ) |>
        dplyr::select(!c(CopyNumber, Ploidy)) |>
        dplyr::rename(CopyNumber = CopyNumber_adjusted)

    return(segments)
}

load_copynumber <- function(files, ploidy, source) {
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

    segments <- adjust_copynumber(segments)

    if (source == "pdx") {
        segments <- segments |> dplyr::mutate(Tumor_Sample_Barcode = paste0(Tumor_Sample_Barcode, "_hum"))
    }

    return(segments)
}

ploidy_table <- readr::read_table(ploidy_file)
colnames(ploidy_table) <- c("Tumor_Sample_Barcode", "Ploidy")

#  Load and adjust segments for human data.
segments_human_1st_models <- load_copynumber(
    files = human_1st_model_files,
    ploidy = ploidy_table,
    source = "human"
)
segments_human_alt_models <- load_copynumber(
    files = human_alt_model_files,
    ploidy = ploidy_table,
    source = "human"
)

#  Load and adjust segments for PDX data.
segments_pdx_1st_models <- load_copynumber(
    files = pdx_1st_model_files,
    ploidy = ploidy_table,
    source = "pdx"
)
segments_pdx_alt_models <- load_copynumber(
    files = pdx_alt_model_files,
    ploidy = ploidy_table,
    source = "pdx"
)

all_segments <- dplyr::bind_rows(
    segments_human_1st_models,
    segments_pdx_1st_models,
    segments_human_alt_models,
    segments_pdx_alt_models
) |> tidyr::drop_na()

write.table(
    all_segments,
    file = paste0(data_dir, "combined_Human+PDX_with_1st+AltModels_segments_for_meskit.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)
