library(MesKit)
library(dplyr)
library(ggpubr)
library(stringr)

#  Set project directory.
proj_dir <- "/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES/analysis/release_v2/combined_2729_3248/"
setwd(proj_dir)

##############
#### Data ####
##############

#  Load manifest.
manifest <- read.csv("2729_3248_combined_manifest.txt", sep = "\t")

#  Load COSMIC's Cancer Gene Census (CGC) file.
cgc <- read.csv(
    "/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES/resources/COSMIC/cancer_gene_census.v97.csv",
    header = TRUE,
    row.names = NULL,
    check.names = FALSE
)

#  Create MAF object.
maf <- readMaf(
    mafFile = "2729_matched_tum_3248_matched_tum_combined_filtered_for_meskit.maf",
    clinicalFile = "2729_matched_tum_3248_matched_tum_metadata_for_meskit.tsv",
    refBuild = "hg38"
)

###################
#### Functions ####
###################

make_readable <- function(maf, patient) {
    maf[[patient]]@data$Tumor_Sample_Label <- manifest$tumour_sample_ID_in_COSMIC[match(
        maf[[patient]]@data$Tumor_Sample_Barcode, manifest$Tumour.Sample
    )]

    maf[[patient]]@data$Tumor_Sample_Label <- gsub("_tumor_skin_melanoma", "", maf[[patient]]@data$Tumor_Sample_Label)
    maf[[patient]]@data$Tumor_Sample_Label <- gsub("_gDNA", "", maf[[patient]]@data$Tumor_Sample_Label)
    maf[[patient]]@data$Tumor_Sample_Label <- gsub("_Lymph_node", "", maf[[patient]]@data$Tumor_Sample_Label)
    maf[[patient]]@data$Tumor_Sample_Label <- gsub("_Tumour - local_recurrence", "", maf[[patient]]@data$Tumor_Sample_Label)
    maf[[patient]]@data$Tumor_Sample_Label <- gsub("_Tumour - primary", "", maf[[patient]]@data$Tumor_Sample_Label)

    return(maf)
}

make_tree <- function(maf, patient) {
    setwd(paste0(proj_dir, substring(patient, 1, 7)))

    tree <- getPhyloTree(maf, patient.id = patient, min.vaf = 0.06)

    plot <- plotPhyloTree(tree, use.tumorSampleLabel = TRUE)

    ggpubr::ggexport(
        plot,
        filename = paste0(patient, "_tree_publication.pdf"),
        width = 8,
        height = 6
    )
}

make_heatmap <- function(maf, patient) {
    setwd(paste0(proj_dir, substring(patient, 1, 7)))

    heatmap <- mutHeatmap(
        maf,
        min.ccf = 0.04,
        use.ccf = FALSE, patient.id = patient,
        use.tumorSampleLabel = TRUE
    )

    ggpubr::ggexport(
        heatmap,
        filename = paste0(patient, "_heatmap_publication.pdf"),
        width = 6,
        height = 6
    )
}

#################
#### PD53330 ####
#################

maf <- make_readable(maf, "PD53330")

make_tree(maf, "PD53330")
make_heatmap(maf, "PD53330")

#################
#### PD53332 ####
#################

maf <- make_readable(maf, "PD53332")

maf[["PD53332"]]@data <- maf[["PD53332"]]@data[!(maf[["PD53332"]]@data$Tumor_Sample_Label %in% c("AM003c")), ]

make_tree(maf, "PD53332")
make_heatmap(maf, "PD53332")

#################
#### PD53333 ####
#################

maf <- make_readable(maf, "PD53333")

make_tree(maf, "PD53333")
make_heatmap(maf, "PD53333")

#################
#### PD53337 ####
#################

maf <- make_readable(maf, "PD53337")

make_tree(maf, "PD53337")
make_heatmap(maf, "PD53337")

#################
#### PD53343 ####
#################

maf <- make_readable(maf, "PD53343")

maf[["PD53343"]]@data$Patient_ID[maf[["PD53343"]]@data$Tumor_Sample_Barcode %in% c("PD53343h", "PD53343h_hum", "PD53343a", "PD53343a_hum", "PD53343d", "PD53343d_hum", "PD53343e", "PD53343e_hum")] <- "PD53343_(A)"
maf[["PD53343"]]@data$Patient_ID[maf[["PD53343"]]@data$Tumor_Sample_Barcode %in% c("PD53343i", "PD53343i_hum", "PD53343c", "PD53343c_hum", "PD53343f", "PD53343f_hum", "PD53343g", "PD53343g_hum")] <- "PD53343_(B)"

maf[["PD53343_(A)"]] <- maf[["PD53343"]]
maf[["PD53343_(A)"]]@data <- maf[["PD53343_(A)"]]@data |> dplyr::filter(Patient_ID == "PD53343_(A)")

maf[["PD53343_(B)"]] <- maf[["PD53343"]]
maf[["PD53343_(B)"]]@data <- maf[["PD53343_(B)"]]@data |> dplyr::filter(Patient_ID == "PD53343_(B)")

make_tree(maf, "PD53343_(A)")
make_heatmap(maf, "PD53343_(A)")
make_tree(maf, "PD53343_(B)")
make_heatmap(maf, "PD53343_(B)")

#################
#### PD53347 ####
#################

maf <- make_readable(maf, "PD53347")

make_tree(maf, "PD53347")
make_heatmap(maf, "PD53347")

#################
#### PD53349 ####
#################

maf <- make_readable(maf, "PD53349")

maf[["PD53349"]]@data$Patient_ID[maf[["PD53349"]]@data$Tumor_Sample_Label %in% c("AM021a", "AM021a-X1_PDX")] <- "PD53349_(A)"
maf[["PD53349"]]@data$Patient_ID[!(maf[["PD53349"]]@data$Tumor_Sample_Label %in% c("AM021a", "AM021a-X1_PDX"))] <- "PD53349_(B)"

maf[["PD53349_(A)"]] <- maf[["PD53349"]]
maf[["PD53349_(A)"]]@data <- maf[["PD53349_(A)"]]@data |> dplyr::filter(Patient_ID == "PD53349_(A)")

maf[["PD53349_(B)"]] <- maf[["PD53349"]]
maf[["PD53349_(B)"]]@data <- maf[["PD53349_(B)"]]@data |> dplyr::filter(Patient_ID == "PD53349_(B)")

make_tree(maf, "PD53349_(A)")
make_heatmap(maf, "PD53349_(A)")
make_tree(maf, "PD53349_(B)")
make_heatmap(maf, "PD53349_(B)")

#################
#### PD53350 ####
#################

maf <- make_readable(maf, "PD53350")

maf[["PD53350"]]@data$Patient_ID[maf[["PD53350"]]@data$Tumor_Sample_Label %in% c("AM022a", "AM022a-X1_PDX")] <- "PD53350_(A)"
maf[["PD53350"]]@data$Patient_ID[!(maf[["PD53350"]]@data$Tumor_Sample_Label %in% c("AM022a", "AM022a-X1_PDX"))] <- "PD53350_(B)"

maf[["PD53350_(A)"]] <- maf[["PD53350"]]
maf[["PD53350_(A)"]]@data <- maf[["PD53350_(A)"]]@data |> dplyr::filter(Patient_ID == "PD53350_(A)")

maf[["PD53350_(B)"]] <- maf[["PD53350"]]
maf[["PD53350_(B)"]]@data <- maf[["PD53350_(B)"]]@data |> dplyr::filter(Patient_ID == "PD53350_(B)")

make_tree(maf, "PD53350_(A)")
make_heatmap(maf, "PD53350_(A)")
make_tree(maf, "PD53350_(B)")
make_heatmap(maf, "PD53350_(B)")

#################
#### PD53352 ####
#################

maf <- make_readable(maf, "PD53352")

maf[["PD53352"]]@data$Patient_ID[maf[["PD53352"]]@data$Tumor_Sample_Label %in% c("AM025a", "AM025a-X1_PDX")] <- "PD53352_(A)"
maf[["PD53352"]]@data$Patient_ID[!(maf[["PD53352"]]@data$Tumor_Sample_Label %in% c("AM025a", "AM025a-X1_PDX"))] <- "PD53352_(B)"

maf[["PD53352_(A)"]] <- maf[["PD53352"]]
maf[["PD53352_(A)"]]@data <- maf[["PD53352_(A)"]]@data |> dplyr::filter(Patient_ID == "PD53352_(A)")

maf[["PD53352_(B)"]] <- maf[["PD53352"]]
maf[["PD53352_(B)"]]@data <- maf[["PD53352_(B)"]]@data |> dplyr::filter(Patient_ID == "PD53352_(B)")

make_tree(maf, "PD53352_(A)")
make_heatmap(maf, "PD53352_(A)")
make_tree(maf, "PD53352_(B)")
make_heatmap(maf, "PD53352_(B)")

#################
#### PD53355 ####
#################

maf <- make_readable(maf, "PD53355")

make_tree(maf, "PD53355")
make_heatmap(maf, "PD53355")

#################
#### PD53357 ####
#################

maf <- make_readable(maf, "PD53357")

maf[["PD53357"]]@data <- maf[["PD53357"]]@data[maf[["PD53357"]]@data$Tumor_Sample_Label %in% c("AM032b", "AM032b-X1_PDX"), ]

make_tree(maf, "PD53357")
make_heatmap(maf, "PD53357")

#################
#### PD53359 ####
#################

maf <- make_readable(maf, "PD53359")

make_tree(maf, "PD53359")
make_heatmap(maf, "PD53359")

#################
#### PD53364 ####
#################

maf <- make_readable(maf, "PD53364")

make_tree(maf, "PD53364")
make_heatmap(maf, "PD53364")
