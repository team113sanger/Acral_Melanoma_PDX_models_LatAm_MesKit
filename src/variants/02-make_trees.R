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
maf_file <- paste0(data_dir, "variants/6633_2729_3248-filtered_mutations_matched_allTum_keep_adjusted_for_meskit.maf")
clinical_file <- paste0(data_dir, "variants/6633_2729_3248-filtered_mutations_matched_allTum_keep_clinical_for_meskit.tsv")
cosmic_file <- paste0(data_dir, "cancer_gene_census.v97.csv")

mdata_dir <- here::here("metadata/")
mdata_file <- paste0(mdata_dir, "6633_METADATA_table.tsv")

outdir <- here::here("results/variants/")

######################
#### Loading data ####
######################

logger::log_info("Creating MAF object with MesKit...")
maf <- MesKit::readMaf(mafFile = maf_file, clinicalFile = clinical_file, refBuild = "hg38")

logger::log_info("Loading metadata...")
metadata <- read.csv(mdata_file, sep = "\t")

logger::log_info("Loading COSMIC Cancer Gene Census data...")
cgc <- read.csv(cosmic_file, header = TRUE, row.names = NULL, check.names = FALSE)

###################
#### Functions ####
###################

make_readable <- function(maf, patient) {
    logger::log_info("Making sample IDs more readable...")

    maf[[patient]]@data$Tumor_Sample_Label <- metadata$sample_ID_in_COSMIC[match(
        maf[[patient]]@data$Tumor_Sample_Barcode, metadata$Tumor_Sample_Barcode
    )]

    maf[[patient]]@data$Tumor_Sample_Label <- gsub("_tumor_skin_melanoma", "", maf[[patient]]@data$Tumor_Sample_Label)
    maf[[patient]]@data$Tumor_Sample_Label <- gsub("_gDNA", "", maf[[patient]]@data$Tumor_Sample_Label)
    maf[[patient]]@data$Tumor_Sample_Label <- gsub("_Lymph_node", "", maf[[patient]]@data$Tumor_Sample_Label)
    maf[[patient]]@data$Tumor_Sample_Label <- gsub("_Tumour - local_recurrence", "", maf[[patient]]@data$Tumor_Sample_Label)
    maf[[patient]]@data$Tumor_Sample_Label <- gsub("_Tumour - primary", "", maf[[patient]]@data$Tumor_Sample_Label)

    return(maf)
}

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
        sample_ids <- unique(subset$sample_interpretable[subset$mutation_status_in_sample == 1]) %>%
            paste(collapse = ", ")

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

    return(df)
}

process_tree_dataframe <- function(tree_df, maf, patient) {
    logger::log_info("Extracting mutation info...")

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
    colnames(df) <- c("sample", "event", "gene", "location", "mutation", "mutation_status_in_sample")

    #  Remove column called "NORMAL" as this is just a "fake" column.
    df <- df[df$sample != "NORMAL", ]

    # Add interpretable sample IDs to the dataframe.
    df$sample_interpretable <- maf[[patient]]@data$Tumor_Sample_Label[match(df$sample, maf[[patient]]@data$Tumor_Sample_Barcode)]

    df <- add_n_samples_mutation_is_present(df)
    df <- add_sample_ids_mutation_is_present(df)
    df <- add_mutation_type(df)
    df <- add_cosmic_info(df)

    # Reorder dataframe columns.
    df <- df[, c(
        "sample", "sample_interpretable",
        "event", "gene", "cosmic_cgc_role_in_cancer",
        "cosmic_cgc_tumour_types_somatic", "cosmic_cgc_tumour_types_germline",
        "location", "mutation", "mutation_status_in_sample",
        "mutation_type", "n_samples_mutation_is_present",
        "sample_ids_mutation_is_present"
    )]

    return(df)
}

make_tree <- function(maf, patient) {
    set.seed(42, kind = "L'Ecuyer-CMRG")
    tree <- MesKit::getPhyloTree(maf, patient.id = patient, method = "MP", min.vaf = 0.06)

    # Extract information from tree object. This will be a dataframe with samples as columns and mutations as rows.
    bin_df <- as.data.frame(tree@binary.matrix)
    bin_df <- process_tree_dataframe(bin_df, maf, patient)

    # Save dataframe to csv file.
    write.csv(bin_df, paste0(outdir, patient, "/", patient, "_table.csv"))
    write.csv(bin_df[bin_df$mutation_status_in_sample != "0", ], paste0(outdir, patient, "/", patient, "_table_filt.csv"))

    logger::log_info(paste0("Making tree for patient ", patient, "..."))

    tree <- MesKit::plotPhyloTree(tree, use.tumorSampleLabel = TRUE)

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

split_a_and_b <- function(maf, patient, samples_in_a) {
    logger::log_info(paste0("Splitting data for patient ", patient, "..."))

    id_A <- paste(patient, "(A)", sep = "_")
    id_B <- paste(patient, "(B)", sep = "_")

    if (TRUE %in% grep("PD", samples_in_a)) {
        maf[[patient]]@data$Patient_ID[maf[[patient]]@data$Tumor_Sample_Barcode %in% samples_in_a] <- id_A
        maf[[patient]]@data$Patient_ID[!(maf[[patient]]@data$Tumor_Sample_Barcode %in% samples_in_a)] <- id_B
    } else if (TRUE %in% grep("AM", samples_in_a)) {
        maf[[patient]]@data$Patient_ID[maf[[patient]]@data$Tumor_Sample_Label %in% samples_in_a] <- id_A
        maf[[patient]]@data$Patient_ID[!(maf[[patient]]@data$Tumor_Sample_Label %in% samples_in_a)] <- id_B
    } else {
        error_message <- paste0("Invalid sample IDs were provided to the function split_a_and_b()")
        logger::log_fatal(error_message)
        stop(error_message)
    }

    maf[[id_A]] <- maf[[patient]]
    maf[[id_A]]@data <- maf[[id_A]]@data |> dplyr::filter(Patient_ID == as.name(id_A))

    maf[[id_B]] <- maf[[patient]]
    maf[[id_B]]@data <- maf[[id_B]]@data |> dplyr::filter(Patient_ID == as.name(id_B))

    return(list(maf, id_A, id_B))
}

process_patient <- function(maf, patient, action, column_for_split, samples_in_a, filter_out, keep) {
    logger::log_info(paste0("Processing data of patient ", patient, "..."))

    maf <- make_readable(maf, patient)

    if (action == "split") {
        post_split <- split_a_and_b(maf, patient, samples_in_a)

        maf <- post_split[[1]]
        id_A <- post_split[[2]]
        id_B <- post_split[[3]]

        for (id in c(id_A, id_B)) {
            dir.create(paste0(outdir, id))
            plot_tree_and_heatmap(maf, id)
        }
    } else {
        dir.create(paste0(outdir, patient))
        if (action == "filter_out") {
            maf[[patient]]@data <- maf[[patient]]@data[!(maf[[patient]]@data$Tumor_Sample_Label %in% filter_out), ]
        } else if (action == "keep") {
            maf[[patient]]@data <- maf[[patient]]@data[maf[[patient]]@data$Tumor_Sample_Label %in% keep, ]
        }
        plot_tree_and_heatmap(maf, patient)
    }

    logger::log_info(paste0("Finished processing data of patient ", patient, "!"))
    cat("\n")
}

#########################
#### Organising data ####
#########################

patients <- list()
patients[["PD53330"]] <- list(action = "none", column_for_split = NA, samples_in_a = NA, filter_out = NA, keep = NA)
patients[["PD53332"]] <- list(action = "filter_out", column_for_split = NA, samples_in_a = NA, filter_out = "AM003c", keep = NA)
patients[["PD53333"]] <- list(action = "none", column_for_split = NA, samples_in_a = NA, filter_out = NA, keep = NA)
patients[["PD53337"]] <- list(action = "none", column_for_split = NA, samples_in_a = NA, filter_out = NA, keep = NA)
patients[["PD53343"]] <- list(
    action = "split",
    column_for_split = "Tumor_Sample_Barcode",
    samples_in_a = c("PD53343h", "PD53343h_hum", "PD53343a", "PD53343a_hum", "PD53343d", "PD53343d_hum", "PD53343e", "PD53343e_hum"),
    filter_out = NA,
    keep = NA
)
patients[["PD53347"]] <- list(action = "none", column_for_split = NA, samples_in_a = NA, filter_out = NA, keep = NA)
patients[["PD53349"]] <- list(
    action = "split",
    column_for_split = "Tumor_Sample_Label",
    samples_in_a = c("AM021a", "AM021a-X1_PDX"),
    filter_out = NA,
    keep = NA
)
patients[["PD53350"]] <- list(
    action = "split",
    column_for_split = "Tumor_Sample_Label",
    samples_in_a = c("AM022a", "AM022a-X1_PDX"),
    filter_out = NA,
    keep = NA
)
patients[["PD53352"]] <- list(
    action = "split",
    column_for_split = "Tumor_Sample_Label",
    samples_in_a = c("AM025a", "AM025a-X1_PDX"),
    filter_out = NA,
    keep = NA
)
patients[["PD53357"]] <- list(
    action = "keep",
    column_for_split = NA,
    samples_in_a = NA,
    filter_out = NA,
    keep = c("AM032b", "AM032b-X1_PDX")
)
patients[["PD53359"]] <- list(action = "none", column_for_split = NA, samples_in_a = NA, filter_out = NA, keep = NA)
patients[["PD53364"]] <- list(action = "none", column_for_split = NA, samples_in_a = NA, filter_out = NA, keep = NA)

for (patient in names(patients)) {
    process_patient(
        maf,
        patient,
        patients[[patient]]$action,
        patients[[patient]]$column_for_split,
        patients[[patient]]$samples_in_a,
        patients[[patient]]$filter_out,
        patients[[patient]]$keep
    )
}
