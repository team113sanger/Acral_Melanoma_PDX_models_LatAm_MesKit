library(MesKit)
library(dplyr)
library(ggpubr)
library(stringr)

proj_dir <- "/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025"
setwd(proj_dir)

segments <- readSegment(segFile = "segments_for_meskit.tsv")

#  Load manifest.
manifest <- read.csv(
    "/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES/analysis/release_v2/combined_2729_3248/2729_3248_combined_manifest.txt",
    sep = "\t"
)

dir.create("plots/")
setwd("plots/")

for (patient in segments %>% names()) {
    if (length(unique(segments[[patient]]$Tumor_Sample_Barcode)) == 1) next

    dir.create(patient)

    segments[[patient]]$Tumor_Sample_Label <- manifest$tumour_sample_ID_in_COSMIC[match(
        segments[[patient]]$Tumor_Sample_Barcode, manifest$Tumour.Sample
    )]

    #  Adjust Tumor_Sample_Label IDs to be shorter and fit better in the plots later on.
    segments[[patient]]$Tumor_Sample_Label <- gsub("_tumor_skin_melanoma", "", segments[[patient]]$Tumor_Sample_Label)
    segments[[patient]]$Tumor_Sample_Label <- gsub("_gDNA", "", segments[[patient]]$Tumor_Sample_Label)
    segments[[patient]]$Tumor_Sample_Label <- gsub("_Lymph_node", "", segments[[patient]]$Tumor_Sample_Label)
    segments[[patient]]$Tumor_Sample_Label <- gsub("_Tumour - local_recurrence", "", segments[[patient]]$Tumor_Sample_Label)
    segments[[patient]]$Tumor_Sample_Label <- gsub("_Tumour - primary", "", segments[[patient]]$Tumor_Sample_Label)

    # segments[[patient]] <- segments[[patient]] |> tidyr::drop_na()

    plot <- plotCNA(segments, patient.id = patient, use.tumorSampleLabel = TRUE)

    ggexport(plot, filename = paste0(patient, "/", patient, ".pdf"), width = 8, height = 4)
}
