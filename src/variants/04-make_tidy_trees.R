#!/usr/bin/env Rscript

library(here)
library(MesKit)
library(dplyr)
library(ggpubr)
library(stringr)
library(logger)
library(cowplot)

logger::log_threshold(logger::INFO)

########################
#### Defining paths ####
########################

data_dir <- here::here("data/")
maf_file <- paste0(data_dir, "variants/6633_2729_3248-filtered_mutations_matched_allTum_keep_KEEPONLY_adjusted_for_meskit.maf")
clinical_file <- paste0(data_dir, "variants/6633_2729_3248-filtered_mutations_matched_allTum_keep_KEEPONLY_clinical_for_meskit.tsv")
cosmic_file <- paste0(data_dir, "cancer_gene_census.v97.csv")

mdata_dir <- here::here("metadata/")
mdata_file <- paste0(mdata_dir, "6633_2729_3248_METADATA_PDX_from_Latin_America_WES.txt")

outdir <- here::here("results/variants/")

######################
#### Loading data ####
######################

logger::log_info("Creating MAF object with MesKit...")
maf <- MesKit::readMaf(mafFile = maf_file, clinicalFile = clinical_file, refBuild = "hg38")

logger::log_info("Loading clinical file...")
clinical <- read.csv(clinical_file, sep = "\t", check.names = FALSE)

logger::log_info("Loading metadata...")
metadata <- read.csv(mdata_file, sep = "\t", check.names = FALSE)

logger::log_info("Loading COSMIC Cancer Gene Census data...")
cgc <- read.csv(cosmic_file, header = TRUE, row.names = NULL, check.names = FALSE)

###################
#### Functions ####
###################

make_tree <- function(maf, patient, samples) {
    logger::log_info(paste0("Making tree for patient ", patient, "..."))

    #  Set seed for reproducibility.
    set.seed(42, kind = "L'Ecuyer-CMRG")

    maf[[patient]]@data <- maf[[patient]]@data |>
        dplyr::filter(Tumor_Sample_Barcode %in% samples)

    #  Generate tree.
    tree <- MesKit::getPhyloTree(maf, patient.id = patient, method = "MP", min.vaf = 0.06)

    # Plot tree.
    tree <- MesKit::plotPhyloTree(tree, use.tumorSampleLabel = TRUE)

    # Save plot to output file.
    ggpubr::ggexport(
        tree,
        filename = paste0(outdir, patient, "/", patient, "_tree__tidy.pdf"),
        width = 8,
        height = 6,
        verbose = FALSE
    )
}

####################
#### Make trees ####
####################

make_tree(maf = maf, patient = "AM003a", samples = c("AM003a", "AM003a-X1", "AM003a-X2"))
make_tree(maf = maf, patient = "AM004", samples = c("AM004", "AM004-X1", "AM004-X2", "AM004-X3", "AM004-X4"))
make_tree(maf = maf, patient = "AM015a", samples = c("AM015a", "AM015a-X1", "AM015a-X2", "AM015a-X3", "AM015a-X3-Cell"))
make_tree(maf = maf, patient = "AM015b", samples = c("AM015b", "AM015b-X1", "AM015b-X2", "AM015b-X3", "AM015b-X2-Cell"))
make_tree(maf = maf, patient = "AM021b", samples = c("AM021b", "AM021b-X1", "AM021b-X2-Cell"))
make_tree(maf = maf, patient = "AM022a", samples = c("AM022a", "AM022a-X1", "AM022a-X1-Cell2"))
make_tree(maf = maf, patient = "AM050", samples = c("AM050", "AM050-X1", "AM050-X2-Cell"))
