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
    )
)

parser <- OptionParser(
    usage = "prepare_sequenza_for_meskit.R [options] --human path/to/human --pdx path/to/pdx --outdir path/to/outdir",
    option_list = option_list
)
arguments <- parse_args(parser, positional_arguments = 0)

human <- arguments$options$human
pdx <- arguments$options$pdx
outdir <- arguments$options$outdir

sequenza_files_human <- fs::dir_ls(human, recurse = TRUE, glob = "*_segments.txt$")
sequenza_files_pdx <- fs::dir_ls(pdx, recurse = TRUE, glob = "*_segments.txt$")

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
    ) |>
    dplyr::mutate(
        Tumor_Sample_Barcode = paste0(Tumor_Sample_Barcode, "_hum")
    )

seg <- dplyr::bind_rows(seg_human, seg_pdx) |> tidyr::drop_na()

write.table(
    seg, 
    file = paste0(outdir, "/segments_for_meskit.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)
