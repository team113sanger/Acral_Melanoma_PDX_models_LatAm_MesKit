library(MesKit)
library(dplyr)
library(ggpubr)
library(stringr)
library(optparse)

option_list <- list(
    make_option(c("--segments-file"),
        dest = "segments_file",
        action = "store",
        default = NA,
        type = "character",
        help = "Path to a file with Sequenza segments."
    ),
    make_option(c("--manifest"),
        dest = "manifest",
        action = "store",
        default = NA,
        type = "character",
        help = "Path to a manifest file."
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
    usage = "cna_analysis.R [options] --segments-file path/to/segments.tsv --manifest path/to/manifest.txt ---outdir path/to/outdir",
    option_list = option_list
)
arguments <- parse_args(parser, positional_arguments = 0)

segments_file <- arguments$options$segments_file
manifest <- arguments$options$manifest
outdir <- arguments$options$outdir

segments <- readSegment(segFile = segments_file)

#  Load manifest.
manifest <- read.csv(manifest, sep = "\t")

#  Check the last character in the directory path is a slash symbol ("/").
if (substr(outdir, nchar(outdir), nchar(outdir)) != "/") {
    outdir <- paste0(outdir, "/")
}
dir.create(paste0(outdir, "plots/"))

for (patient in segments %>% names()) {
    if (length(unique(segments[[patient]]$Tumor_Sample_Barcode)) == 1) next

    dir.create(paste0(outdir, "plots/", patient))

    segments[[patient]]$Tumor_Sample_Label <- manifest$tumour_sample_ID_in_COSMIC[match(
        segments[[patient]]$Tumor_Sample_Barcode, manifest$Tumour.Sample
    )]

    #  Adjust Tumor_Sample_Label IDs to be shorter and fit better in the plots later on.
    segments[[patient]]$Tumor_Sample_Label <- gsub("_tumor_skin_melanoma", "", segments[[patient]]$Tumor_Sample_Label)
    segments[[patient]]$Tumor_Sample_Label <- gsub("_gDNA", "", segments[[patient]]$Tumor_Sample_Label)
    segments[[patient]]$Tumor_Sample_Label <- gsub("_Lymph_node", "", segments[[patient]]$Tumor_Sample_Label)
    segments[[patient]]$Tumor_Sample_Label <- gsub("_Tumour - local_recurrence", "", segments[[patient]]$Tumor_Sample_Label)
    segments[[patient]]$Tumor_Sample_Label <- gsub("_Tumour - primary", "", segments[[patient]]$Tumor_Sample_Label)

    plot <- plotCNA(
        segments,
        patient.id = patient,
        use.tumorSampleLabel = TRUE,
        sample.bar.height = 0.25,
        chrom.bar.height = 0.1
    )

    ggexport(plot, filename = paste0(outdir, "plots/", patient, "/", patient, ".pdf"), width = 15, height = 4)
}
