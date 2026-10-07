#!/usr/bin/env Rscript

library(here)
library(MesKit)
library(dplyr)
library(ggpubr)
library(stringr)
library(logger)

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

add_n_samples_mutation_is_present <- function(df) {
    logger::log_info("Counting in how many samples each mutation is present...")

    #  Identify, for each mutation, in how many samples the mutation is detected.
    df$n_samples_mutation_is_present <- sapply(1:nrow(df), function(x) {
        #  Subset dataframe to keep only a given event.
        subset <- df[df$event == df$event[x], ]
        #  Count number of samples in which that event is detected.
        n_samples <- length(unique(subset$sample[subset$mutation_status_in_sample == 1]))
        return(n_samples)
    })

    return(df)
}

add_sample_ids_mutation_is_present <- function(df) {
    logger::log_info("Getting IDs of samples in which each mutation is present...")

    df$sample_ids_mutation_is_present <- sapply(1:nrow(df), function(x) {
        #  Subset dataframe to keep only a given event.
        subset <- df[df$event == df$event[x], ]

        #  Get sample IDs of samples that bear the mutation.
        sample_ids <- unique(subset$sample[subset$mutation_status_in_sample == 1]) %>%
            paste(collapse = "/")

        return(sample_ids)
    })

    return(df)
}

add_mutation_type <- function(df) {
    logger::log_info("Determining mutation types...")

    df$mutation_type <- sapply(1:nrow(df), function(x) {
        #  If the mutation was not identified in a given sample, just return NA.
        if (df$mutation_status_in_sample[x] == "0") {
            return(NA)

            #  Otherwise...
        } else {
            #  If all samples bear the mutation, then this is a public mutation.
            if (df$n_samples_mutation_is_present[x] == length(unique(df$sample))) {
                return("Public")

                # If a subset of samples bears the mutation, then this is a shared mutation.
            } else if (df$n_samples_mutation_is_present[x] == 1) {
                return("Private")

                # If only one sample bears the mutation, then this is a private mutation.
            } else if (df$n_samples_mutation_is_present[x] > 1 && df$n_samples_mutation_is_present[x] < length(unique(df$sample))) {
                return("Shared")
            }
        }
    })

    return(df)
}

add_cosmic_info <- function(df) {
    logger::log_info("Adding COSMIC info...")

    #  Add gene's role in cancer based on CGC data.
    df$cosmic_cgc_role_in_cancer <- cgc$`Role in Cancer`[match(df$gene, cgc$`Gene Symbol`)]
    df$cosmic_cgc_tumour_types_somatic <- cgc$`Tumour Types(Somatic)`[match(df$gene, cgc$`Gene Symbol`)]
    df$cosmic_cgc_tumour_types_germline <- cgc$`Tumour Types(Germline)`[match(df$gene, cgc$`Gene Symbol`)]

    df$cosmic_cgc_tumour_types_somatic <- sapply(1:nrow(df), function(x) {
        tumour_types <- df$cosmic_cgc_tumour_types_somatic[x]

        if (!is.na(tumour_types)) {
            tumour_types <- gsub(", ", "/", tumour_types)
        }

        return(tumour_types)
    })

    df$cosmic_cgc_tumour_types_germline <- sapply(1:nrow(df), function(x) {
        tumour_types <- df$cosmic_cgc_tumour_types_germline[x]

        if (!is.na(tumour_types)) {
            tumour_types <- gsub(", ", "/", tumour_types)
        }

        return(tumour_types)
    })

    return(df)
}

process_tree_dataframe <- function(tree, maf, patient) {
    logger::log_info("Extracting mutation info...")

    # Extract information from the tree object.
    # This will be a dataframe with samples as columns and mutations as rows.
    tree_df <- as.data.frame(tree@binary.matrix)

    #  Create an empty dataframe, where rearranged information from tree_df will be stored.
    df <- data.frame()

    #  Iterate over samples...
    for (sample in colnames(tree_df)) {
        #  Iterate over "events" (i.e. a concatenation of gene + location + mutation)...
        for (event in rownames(tree_df)) {
            #  Split event into gene, location and mutation.
            event_split <- unlist(str_split(event, ":"))
            gene <- event_split[1]
            location <- paste(event_split[2], event_split[3], sep = ":")
            mutation <- paste(event_split[4], event_split[5], sep = ":")
            #  Add elements as a row to the dataframe created above.
            df <- rbind(df, c(sample, event, gene, location, mutation, tree_df[event, sample]))
        }
    }

    # Rename columns.
    colnames(df) <- c("sample", "event", "gene", "location", "mutation", "mutation_status_in_sample")

    # Remove column called "NORMAL" as this is just a "fake" column.
    df <- df[df$sample != "NORMAL", ]

    #  Add other columns.
    df <- add_n_samples_mutation_is_present(df)
    df <- add_sample_ids_mutation_is_present(df)
    df <- add_mutation_type(df)
    df <- add_cosmic_info(df)

    # Reorder dataframe columns.
    df <- df[, c(
        "sample", "event", "gene", "cosmic_cgc_role_in_cancer",
        "cosmic_cgc_tumour_types_somatic", "cosmic_cgc_tumour_types_germline",
        "location", "mutation", "mutation_status_in_sample",
        "mutation_type", "n_samples_mutation_is_present",
        "sample_ids_mutation_is_present"
    )]

    # Save dataframes to csv files.

    df_path <- paste0(outdir, patient, "/", patient, "_table.csv")
    df |> readr::write_csv(file = df_path)

    df_filt_path <- paste0(outdir, patient, "/", patient, "_table_filt.csv")
    df |>
        dplyr::filter(mutation_status_in_sample != "0") |>
        readr::write_csv(file = df_filt_path)
}

make_tree <- function(maf, patient) {
    logger::log_info(paste0("Making tree for patient ", patient, "..."))

    #  Set seed for reproducibility.
    set.seed(42, kind = "L'Ecuyer-CMRG")

    #  Generate tree.
    tree <- MesKit::getPhyloTree(maf, patient.id = patient, method = "MP", min.vaf = 0.06)

    #  Process the tree data.
    process_tree_dataframe(
        tree = tree,
        maf = maf,
        patient = patient
    )

    # Plot tree.
    tree <- MesKit::plotPhyloTree(tree, use.tumorSampleLabel = TRUE)

    # Save plot to output file.
    ggpubr::ggexport(
        tree,
        filename = paste0(outdir, patient, "/", patient, "_tree.pdf"),
        width = 8,
        height = 6,
        verbose = FALSE
    )

    return(tree)
}

make_heatmap <- function(maf, patient) {
    logger::log_info(paste0("Making heatmap for patient ", patient, "..."))

    heatmap <- MesKit::mutHeatmap(
        maf,
        min.ccf = 0.04,
        use.ccf = FALSE,
        patient.id = patient,
        use.tumorSampleLabel = TRUE
    )

    ggpubr::ggexport(
        heatmap,
        filename = paste0(outdir, patient, "/", patient, "_heatmap.pdf"),
        width = 6,
        height = 6,
        verbose = FALSE
    )

    return(heatmap)
}

plot_tree_and_heatmap <- function(maf, patient) {
    logger::log_info(paste0("Plotting tree and heatmap for patient ", patient, "..."))

    tree <- make_tree(maf, patient)
    heatmap <- make_heatmap(maf, patient)

    #  Combine tree and heatmap into a single plot.
    tree_plus_heatmap <- cowplot::plot_grid(tree, heatmap, nrow = 1, rel_widths = c(1.5, 1))

    ggpubr::ggexport(
        tree_plus_heatmap,
        filename = paste0(outdir, patient, "/", patient, "_tree_plus_heatmap.pdf"),
        width = 12,
        height = 6,
        verbose = FALSE
    )
}

process_patient <- function(maf, patient) {
    logger::log_info(paste0("Processing data of patient ", patient, "..."))

    dir.create(paste0(outdir, patient))
    plot_tree_and_heatmap(maf, patient)

    logger::log_info(paste0("Finished processing data of patient ", patient, "!"))
    cat("\n")
}

####################
#### Make trees ####
####################

for (patient in unique(clinical$Patient_ID)) {
    process_patient(maf = maf, patient = patient)
}


##########################
#### Collate patients ####
##########################

files <- fs::dir_ls(outdir, recurse = TRUE, glob = "*_filt.csv$")
all_samples <- readr::read_csv(files)
all_samples$sample <- gsub("_X", "-X", all_samples$sample)

all_samples |>
    dplyr::mutate(patient = stringr::str_extract(sample, "^AM[0-9]+[ab]?")) |>
    dplyr::relocate(patient, .before = sample) |>
    dplyr::arrange(patient, sample) |>
    readr::write_csv(paste0(outdir, "supplementary_table_6.csv"))
