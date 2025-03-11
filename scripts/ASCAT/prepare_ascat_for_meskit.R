library(optparse)
library(dplyr)

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
    )
)

parser <- OptionParser(
    usage = "prepare_ascat_for_meskit.R [options] --segments path/to/segments --outdir path/to/outdir",
    option_list = option_list
)
arguments <- parse_args(parser, positional_arguments = 0)

segments <- arguments$options$segments
outdir <- arguments$options$outdir

# Read segments file.
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

write.table(
    seg, 
    file = paste0(outdir, "/segments_for_meskit.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)