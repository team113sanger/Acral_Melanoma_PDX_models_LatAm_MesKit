library(optparse)
library(dplyr)
library(readr)

option_list <- list(
    make_option(c("-s", "--segments"),
        dest = "segments",
        action = "store",
        default = NA,
        type = "character",
        help = "Path to segments file."
    ),
    make_option(c("-o", "--outdir"),
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
    usage = "prepare_ascat_for_meskit.R [options] --segments path/to/segments --outdir path/to/outdir",
    option_list = option_list
)
arguments <- parse_args(parser, positional_arguments = 0)

segments <- arguments$options$segments
outdir <- arguments$options$outdir
ploidy_file <- arguments$options$ploidy_file

if (!is.na(ploidy)) {
    ploidy_table <- readr::read_table(ploidy_file, col_names = FALSE)
    colnames(ploidy_table) <- c("Tumor_Sample_Barcode", "Purity", "Ploidy")
    warning("A ploidy table has been provided. Copy number results will be adjusted relative to the ploidy values!")
} else {
    warning("A ploidy table has not been provided. Proceeding without normalising copy number results...")
}

#  Read segments file.
seg <- read.table(segments, header = TRUE, row.names = NULL)
colnames(seg) <- c(
    "Tumor_Sample_Barcode",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Major_CN",
    "Minor_CN",
    "Sex",
    "Purity",
    "Ploidy",
    "CN",
    "Size"
)

seg$Patient_ID <- substring(seg$Tumor_Sample_Barcode, 1, 7)
seg$CopyNumber <- (seg$Major_CN + seg$Minor_CN)
seg <- seg %>% select(
    Patient_ID,
    Tumor_Sample_Barcode,
    Chromosome,
    Start_Position,
    End_Position,
    CopyNumber,
    Major_CN,
    Minor_CN,
    Sex,
    Purity,
    Ploidy,
    CN,
    Size
)

if (!is.na(ploidy_file)) {
    seg <- seg |>
        dplyr::left_join(
            ploidy_table,
            by = dplyr::join_by(Tumor_Sample_Barcode)
        ) |>
        dplyr::mutate(Ploidy = round(Ploidy)) |>
        dplyr::mutate(
            CopyNumber_adjusted = dplyr::case_when(
                CopyNumber == 0 ~ 0, # Deletion
                CopyNumber < Ploidsy ~ 1, #  Loss
                CopyNumber == Ploidy ~ 2, # Neutral
                CopyNumber > Ploidy & CopyNumber < (2 * Ploidy) ~ 3, # Gain
                CopyNumber >= (2 * Ploidy) ~ 4 #  Amplification
            ),
            .before = CopyNumber
        ) |>
        dplyr::select(!c(CopyNumber, Ploidy)) |>
        dplyr::rename(CopyNumber = CopyNumber_adjusted)
}

write.table(
    seg,
    file = paste0(outdir, "/", prefix, "_segments_for_meskit.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)
