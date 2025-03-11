library(MesKit)
library(dplyr)
library(ggpubr)
library(stringr)

# Set project directory.
proj_dir <- "/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES/analysis/release_v2/combined_2729_3248/"
setwd(proj_dir)

# Load manifest.
manifest <- read.csv("2729_3248_combined_manifest.txt", sep = "\t")

# Load COSMIC's Cancer Gene Census (CGC) file.
cgc <- read.csv(
    "/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES/resources/COSMIC/cancer_gene_census.v97.csv",
    header = TRUE,
    row.names = NULL,
    check.names = FALSE
)

# Create MAF object.
maf <- readMaf(
    mafFile = "2729_matched_tum_3248_matched_tum_combined_filtered_for_meskit.maf",
    clinicalFile = "2729_matched_tum_3248_matched_tum_metadata_for_meskit.tsv",
    refBuild = "hg38"
)

# Iterate over patients...
for (patient in names(maf)) {

    cat("Running", patient, "\n")

    # Create directory to store patient-level results.
    dir.create(patient)

    # Replace Tumor_Sample_Label with the interpretable sample IDs.
    maf[[patient]]@data$Tumor_Sample_Label <- manifest$tumour_sample_ID_in_COSMIC[match(
        maf[[patient]]@data$Tumor_Sample_Barcode, manifest$Tumour.Sample
    )]

    # Adjust Tumor_Sample_Label IDs to be shorter and fit better in the plots later on.
    maf[[patient]]@data$Tumor_Sample_Label <- gsub("_tumor_skin_melanoma", "", maf[[patient]]@data$Tumor_Sample_Label)
    maf[[patient]]@data$Tumor_Sample_Label <- gsub("_gDNA", "", maf[[patient]]@data$Tumor_Sample_Label)
    maf[[patient]]@data$Tumor_Sample_Label <- gsub("_Lymph_node", "_LN_Met", maf[[patient]]@data$Tumor_Sample_Label)
    maf[[patient]]@data$Tumor_Sample_Label <- gsub("_Tumour - local_recurrence", "_Loc_Rec", maf[[patient]]@data$Tumor_Sample_Label)
    maf[[patient]]@data$Tumor_Sample_Label <- gsub("_Tumour - primary", "_Primary", maf[[patient]]@data$Tumor_Sample_Label)

    # Iterate over tree generation methods...
    for (tree_method in c("NJ", "MP")) {

        cat("\tGenerating tree with method", tree_method, "\n")

        # Generate tree object with a given method for a given patient.
        phyloTree <- getPhyloTree(maf, patient.id = patient, method = tree_method, min.vaf = 0.06)

        # Extract information from tree object. This will be a dataframe with samples as columns and mutations as rows.
        bin_df <- phyloTree@binary.matrix %>% as.data.frame

        # Create an empty dataframe, where rearranged information from bin_df will be stored.
        df <- data.frame()

        # Iterate over samples...
        for (sample in colnames(bin_df)) {

            # Iterate over "events" (i.e. a concatenation of gene + location + mutation)...
            for (event in rownames(bin_df)) {

                # Split event into gene, location and mutation.
                event_split <- unlist(str_split(event, ":"))
                gene <- event_split[1]
                location <- paste(event_split[2], event_split[3], sep = ":")
                mutation <- paste(event_split[4], event_split[5], sep = ":")

                # Add elements as a row to the dataframe created above.
                df <- rbind(df, c(sample, event, gene, location, mutation, bin_df[event, sample]))
            }
        }

        # Adjust column names.
        colnames(df) <- c("sample", "event", "gene", "location", "mutation", "mutation_status_in_sample")

        print("ok1")

        # Remove column called "NORMAL" as this is just a "fake" column.
        df <- df[df$sample != "NORMAL", ]

        # Add interpretable sample IDs to the dataframe.
        df$sample_interpretable <- maf[[patient]]@data$Tumor_Sample_Label[match(df$sample, maf[[patient]]@data$Tumor_Sample_Barcode)] 

        # Identify, for each mutation, in how many samples the mutation is detected.
        df$n_samples_mutation_is_present <- sapply(1:nrow(df), function(x) {

            # Subset dataframe to keep only a given event.
            subset <- df[df$event == df$event[x], ]

            # Count number of samples in which that event is detected.
            n_samples <- subset$sample[subset$mutation_status_in_sample == 1] %>% unique %>% length

            return(n_samples)
        })

        print("ok2")

        # Identify, for each mutation, the samples IDs of samples that bear the mutation.
        df$sample_ids_mutation_is_present <- sapply(1:nrow(df), function(x) {

            # Subset dataframe to keep only a given event.
            subset <- df[df$event == df$event[x], ]

            # Get sample IDs of samples that bear the mutation.           
            sample_ids <- subset$sample_interpretable[subset$mutation_status_in_sample == 1] %>% unique %>% paste(collapse = ", ")

            return(sample_ids)
        })

        print("ok3")

        # For each mutation, define the mutation type (i.e. public, shared or private).
        df$mutation_type <- sapply(1:nrow(df), function(x) {

            # If the mutation was not identified in a given sample, just return NA.
            if (df$mutation_status_in_sample[x] == "0") {
                return(NA)

            # Otherwise...
            } else {

                # If all samples bear the mutation, then this is a public mutation.
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

        # Add gene's role in cancer based on CGC data.
        df$cosmic_cgc_role_in_cancer <- cgc$`Role in Cancer`[match(df$gene, cgc$`Gene Symbol`)]
        df$cosmic_cgc_tumour_types_somatic <- cgc$`Tumour Types(Somatic)`[match(df$gene, cgc$`Gene Symbol`)]
        df$cosmic_cgc_tumour_types_germline <- cgc$`Tumour Types(Germline)`[match(df$gene, cgc$`Gene Symbol`)]

        print("ok4")

        # Reorder dataframe columns.
        df <- df[, c(
            "sample", "sample_interpretable", 
            "event", "gene", "cosmic_cgc_role_in_cancer", 
            "cosmic_cgc_tumour_types_somatic", "cosmic_cgc_tumour_types_germline",
            "location", "mutation", "mutation_status_in_sample", 
            "mutation_type", "n_samples_mutation_is_present", 
            "sample_ids_mutation_is_present"
        )]

        # Save dataframe to csv file.
        write.csv(df, paste0(patient, "/", patient, "_", tree_method, "_table.csv"))
        write.csv(df[df$mutation_status_in_sample != "0", ], paste0(patient, "/", patient, "_", tree_method, "_table_filt.csv"))

        # Plot tree.
        phylotree <- plotPhyloTree(phyloTree, use.tumorSampleLabel = TRUE)
        phylotree <- phylotree + theme(legend.position = "none")

        # Plot heatmap.
        binary_heatmap <- mutHeatmap(maf, min.ccf = 0.04, use.ccf = FALSE, patient.id = patient, use.tumorSampleLabel = TRUE)

        # Combine tree and heatmap into a single plot.
        plot <- cowplot::plot_grid(phylotree, binary_heatmap, nrow = 1, rel_widths = c(1.5, 1))

        # Save combined plot to PDF file.
        ggexport(plot, filename = paste0(patient, "/", patient, "_", tree_method, "_tree.pdf"), width = 12, height = 6)
    }
}