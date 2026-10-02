library(here)
library(MesKit)
library(ggpubr)

########################
#### Defining paths ####
########################

data_dir <- here::here("data/copy_number/")
segments_file <- paste0(data_dir, "segments_for_meskit.tsv")

mdata_dir <- here::here("metadata/")
mdata_file <- paste0(mdata_dir, "6633_2729_3248_METADATA_PDX_from_Latin_America_WES.txt")

outdir <- here::here("results/copy_number/")

######################
#### Loading data ####
######################

#  Load segments.
segments <- MesKit::readSegment(segFile = segments_file)

#  Load metadata table.
metadata <- read.csv(mdata_file, sep = "\t", check.names = FALSE)

#################################
#### Generating CNA heatmaps ####
#################################

for (patient in segments %>% names()) {
    # Skip patient if there is only one sample available.
    if (length(unique(segments[[patient]]$Tumor_Sample_Barcode)) == 1) next

    # Create subdirectory for a given patient.
    dir.create(paste0(outdir, patient))

    # Redefine the Tumor_Sample_Label IDs.
    segments[[patient]]$Tumor_Sample_Label <- metadata$`Case ID`[match(
        segments[[patient]]$Tumor_Sample_Barcode, metadata$final_sample_name_used
    )]

    # Create CNA heatmap.
    plot <- MesKit::plotCNA(
        segments,
        patient.id = patient,
        refBuild = "hg38",
        use.tumorSampleLabel = TRUE,
        sample.bar.height = 0.25,
        chrom.bar.height = 0.1
    )

    # Export CNA heatmap.
    ggpubr::ggexport(plot, filename = paste0(outdir, patient, "/", patient, ".pdf"), width = 15, height = 4)
}
