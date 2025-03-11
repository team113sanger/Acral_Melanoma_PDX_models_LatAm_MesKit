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
    make_option(c("-o", "--outdir"),
        dest = "outdir",
        action = "store",
        default = NA,
        type = "character",
        help = "Path to output directory."
    ),
    make_option(c("--ploidy"),
        dest = "ploidy",
        action = "store",
        default = NA,
        type = "character",
        help = "Path to a file containing each sample's ploidy; if provided, copy number results will be normalised relative to the ploidy of each sample."
    )
)

parser <- OptionParser(
    usage = "prepare_sequenza_for_meskit.R [options] --human path/to/human --pdx path/to/pdx --outdir path/to/outdir --ploidy path/to/ploidy.tsv",
    option_list = option_list
)
arguments <- parse_args(parser, positional_arguments = 0)

human <- arguments$options$human
pdx <- arguments$options$pdx
outdir <- arguments$options$outdir
ploidy <- arguments$options$ploidy

sequenza_files_human <- fs::dir_ls(human, recurse = TRUE, glob = "*_segments.txt$")
sequenza_files_pdx <- fs::dir_ls(pdx, recurse = TRUE, glob = "*_segments.txt$")
if (!is.na(ploidy)) {
    ploidy_table <- readr::read_table(ploidy) |> dplyr::rename(Tumor_Sample_Barcode = sample)
} else {
    warning("A ploidy table has not been provided. Proceeding without normalising copy number results...")
}

seg_human <- readr::read_tsv(sequenza_files_human, col_names = TRUE, id = "Tumor_Sample_Barcode") |>
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

if (!is.na(ploidy)) {
    seg_human <- seg_human |>
        dplyr::left_join(
            ploidy_table,
            by = dplyr::join_by(Tumor_Sample_Barcode)
        ) |>
        dplyr::mutate(
            CopyNumber_normalised = round(CopyNumber / ploidy.mean.cn),
            .before = CopyNumber
        ) |>
        dplyr::select(!c(CopyNumber, ploidy.mean.cn)) |>
        dplyr::rename(CopyNumber = CopyNumber_normalised)
}

seg_pdx <- readr::read_tsv(sequenza_files_pdx, col_names = TRUE, id = "Tumor_Sample_Barcode") |>
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

if (!is.na(ploidy)) {
    seg_pdx <- seg_pdx |>
        dplyr::left_join(
            ploidy_table,
            by = dplyr::join_by(Tumor_Sample_Barcode)
        ) |>
        dplyr::mutate(
            CopyNumber_normalised = round(CopyNumber / ploidy.mean.cn),
            .before = CopyNumber
        ) |>
        dplyr::select(!c(CopyNumber, ploidy.mean.cn)) |>
        dplyr::rename(CopyNumber = CopyNumber_normalised) |>
        dplyr::mutate(
            Tumor_Sample_Barcode = paste0(Tumor_Sample_Barcode, "_hum")
        )
}

seg <- dplyr::bind_rows(seg_human, seg_pdx) |> tidyr::drop_na()

write.table(
    seg,
    file = paste0(outdir, "/segments_for_meskit.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)
