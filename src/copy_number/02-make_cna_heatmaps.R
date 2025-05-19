library(here)
library(MesKit)
library(ggpubr)

########################
#### Defining paths ####
########################

data_dir <- here::here("data/copy_number/")
segments_file <- paste0(data_dir, "combined_Human+PDX_with_1st+AltModels_segments_for_meskit.tsv")

mdata_dir <- here::here("metadata/")
mdata_file <- paste0(mdata_dir, "6633_METADATA_table.tsv")

outdir <- here::here("results/copy_number/")

######################
#### Loading data ####
######################

#  Load segments.
segments <- MesKit::readSegment(segFile = segments_file)

#  Load metadata table.
metadata <- read.csv(mdata_file, sep = "\t")

#################################
#### Generating CNA heatmaps ####
#################################

for (patient in segments %>% names()) {
    # Skip patient if there is only one sample available.
    if (length(unique(segments[[patient]]$Tumor_Sample_Barcode)) == 1) next

    # Create subdirectory for a given patient.
    dir.create(paste0(outdir, patient))

    # Redefine the Tumor_Sample_Label IDs.
    segments[[patient]]$Tumor_Sample_Label <- metadata$sample_ID_in_COSMIC[match(
        segments[[patient]]$Tumor_Sample_Barcode, metadata$Tumor_Sample_Barcode
    )]

    # Adjust Tumor_Sample_Label IDs to be shorter and fit better in the plots later on.
    for (id_substring in c("_tumor_skin_melanoma", "_gDNA", "_Lymph_node", "_Tumour - local_recurrence", "_Tumour - primary")) {
        segments[[patient]]$Tumor_Sample_Label <- gsub(id_substring, "", segments[[patient]]$Tumor_Sample_Label)
    }

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
