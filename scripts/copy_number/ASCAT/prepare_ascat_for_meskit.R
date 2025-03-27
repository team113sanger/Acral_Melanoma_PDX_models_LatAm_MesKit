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
        dest = "use_ploidy",
        action = "store_true",
        default = FALSE,
        type = "logical",
        help = "Whether ploidy values should be used to adjust copy number results [default = FALSE]."
    )
)

parser <- OptionParser(
    usage = "prepare_ascat_for_meskit.R [options] --segments path/to/segments --outdir path/to/outdir",
    option_list = option_list
)
arguments <- parse_args(parser, positional_arguments = 0)

segments <- arguments$options$segments
outdir <- arguments$options$outdir
prefix <- arguments$options$prefix
use_ploidy <- arguments$options$use_ploidy

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

if (use_ploidy == TRUE) {
    warning("Copy number results will be adjusted relative to the ploidy values!")
    seg <- seg |>
        dplyr::mutate(Ploidy = round(Ploidy)) |>
        dplyr::mutate(
            CopyNumber_adjusted = dplyr::case_when(
                CopyNumber == 0 ~ 0, # Deletion
                CopyNumber != 0 & CopyNumber < Ploidy ~ 1, #  Loss
                CopyNumber == Ploidy ~ 2, # Neutral
                CopyNumber > Ploidy & CopyNumber < (2 * Ploidy) ~ 3, # Gain
                CopyNumber >= (2 * Ploidy) ~ 4 #  Amplification
            ),
            .before = CopyNumber
        ) |>
        dplyr::select(!c(CopyNumber)) |>
        dplyr::rename(CopyNumber = CopyNumber_adjusted)
} else {
    warning("--ploidy has not been specified, therefore copy number results will not be adjusted...")
}

seg <- seg |> dplyr::select(
    Patient_ID,
    Tumor_Sample_Barcode,
    Chromosome,
    Start_Position,
    End_Position,
    CopyNumber,
    Major_CN,
    Minor_CN
)

write.table(
    seg,
    file = paste0(outdir, "/", prefix, "_segments_for_meskit.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)
