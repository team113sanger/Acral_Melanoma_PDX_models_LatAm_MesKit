library(optparse)
library(dplyr)
library(readr)
library(fs)

option_list <- list(
    make_option(c("--human"),
        dest = "human",
        action = "store",
        default = NA,
        type = "character",
        help = "Path to directory with Sequenza outputs for human tumours."
    ),
    make_option(c("--pdx"),
        dest = "pdx",
        action = "store",
        default = NA,
        type = "character",
        help = "Path to directory with Sequenza outputs for PDX samples."
    ),
    make_option(c("--outdir"),
        dest = "outdir",
        action = "store",
        default = NA,
        type = "character",
        help = "Path to output directory."
    ),
    make_option(c("--prefix"),
        dest = "prefix",
        action = "store",
        default = "",
        type = "character",
        help = "Prefix for output segment file."
    ),
    make_option(c("--ploidy"),
        dest = "ploidy_file",
        action = "store",
        default = NA,
        type = "character",
        help = "Path to a file containing each sample's ploidy; if provided, copy number results will be normalised relative to the ploidy of each sample."
    )
)

parser <- OptionParser(
    usage = "prepare_sequenza_for_meskit.R [options] --human path/to/human --pdx path/to/pdx --outdir path/to/outdir --prefix prefix --ploidy path/to/ploidy.tsv",
    option_list = option_list
)
arguments <- parse_args(parser, positional_arguments = 0)

human <- arguments$options$human
pdx <- arguments$options$pdx
outdir <- arguments$options$outdir
prefix <- arguments$options$prefix
ploidy_file <- arguments$options$ploidy_file

if (!is.na(ploidy_file)) {
    ploidy_table <- readr::read_table(ploidy_file)
    colnames(ploidy_table) <- c("Tumor_Sample_Barcode", "Ploidy")
    warning("A ploidy table has been provided. Copy number results will be adjusted relative to the ploidy values!")
} else {
    warning("A ploidy table has not been provided. Proceeding without adjusting copy number results...")
}

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
                CopyNumber < Ploidy ~ 1, #  Loss
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

    if (!is.na(ploidy_file)) {
        segments <- adjust_copynumber(segments)
    }

    if (source == "pdx") {
        segments <- segments |> dplyr::mutate(Tumor_Sample_Barcode = paste0(Tumor_Sample_Barcode, "_hum"))
    }

    return(segments)
}

segments_human <- load_copynumber(
    files = fs::dir_ls(human, recurse = TRUE, glob = "*_segments.txt$"),
    ploidy = ploidy_table,
    source = "human"
)
segments_pdx <- load_copynumber(
    files = fs::dir_ls(pdx, recurse = TRUE, glob = "*_segments.txt$"),
    ploidy = ploidy_table,
    source = "pdx"
)

all_segments <- dplyr::bind_rows(segments_human, segments_pdx) |> tidyr::drop_na()

write.table(
    all_segments,
    file = paste0(outdir, "/", prefix, "_segments_for_meskit.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)
