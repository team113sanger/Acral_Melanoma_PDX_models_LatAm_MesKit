library(here)
library(fs)
library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)

# Define the directory with MesKit outputs for variant calling data and load
# filtered results.
meskit_outdir <- here::here("results/variants/")
files <- fs::dir_ls(meskit_outdir, recurse = TRUE, glob = "*_filt.csv$")
all_samples <- readr::read_csv(files)
all_samples <- all_samples |> dplyr::select(!`...1`)

# Get patient IDs.
patient_ids <- all_samples$sample_interpretable[!grepl("-X", all_samples$sample_interpretable)] |> unique()

# Create empty dataframe to store data.
patient_vs_x1_df_detailed <- data.frame(matrix(nrow = 0, ncol = 11))
patient_vs_x1_df_simplified <- data.frame(matrix(nrow = 0, ncol = 5))

# Iterate over patients...
for (patient in patient_ids) {
  
    print(patient)
  
    # Get the sample IDs for PDXs available for the given patient.
    pdx_ids_for_patient <- all_samples$sample_interpretable[grepl(paste0(patient, "-X"), all_samples$sample_interpretable)] |>
        unique()
  
    # Get the number of PDX samples available for the given patient.
    n_pdx_for_patient <- length(pdx_ids_for_patient)
    
    # Set the sample ID for PDX-X1.
    x1_id <- paste0(patient, "-X1_PDX")
    
    # If there isn't a PDX-X1 for the given patient, move on to the next patient.
    if (!(x1_id %in% all_samples$sample_interpretable)) next

    # Subset original dataframe to keep only results related to the given patient and its PDX-X1.
    patient_and_x1 <- all_samples |> dplyr::filter(sample_interpretable %in% c(patient, x1_id))
    
    ############################################
    #### Variants only found in the patient ####
    ############################################

    # Count number of private mutations in the patient (i.e. those that are only present in the patient).
    n_private_patient <- patient_and_x1 |>
        # Subset the dataframe to keep only rows related to the given patient and that contain private
        # mutations (only present in the patient).
        dplyr::filter(sample_interpretable == patient & mutation_type == "Private") |>
        # Count the mutations.
        dplyr::count(mutation_type) |>
        dplyr::select(n) |>
        as.numeric()
    # If there are no private mutations for the patient, save that as a zero.
    if (is.na(n_private_patient)) n_private_patient <- 0
    
    #######################################
    #### Variants only found in PDX-X1 ####
    #######################################

    # Count number of private mutations in the PDX-X1 (i.e. those that are only present in the PDX-X1).
    n_private_x1 <- patient_and_x1 |>
        # Subset the dataframe to keep only rows related to the given PDX-X1 and that contain private
        # mutations (only present in PDX-X1).
        dplyr::filter(sample_interpretable == x1_id & mutation_type == "Private") |>
        # Count the mutations.
        dplyr::count(mutation_type) |>
        dplyr::select(n) |>
        as.numeric()
    # If there are no private mutations for PDX-X1, save that as a zero.
    if (is.na(n_private_x1)) n_private_x1 <- 0
    
    ################################################################
    #### Variants exclusively shared between patient and PDX-X1 ####
    ################################################################

    # Count number of mutations that are exclusively shared between the patient and PDX-X1.
    n_shared_patient_and_x1 <- patient_and_x1 |>
        # Subset the dataframe to keep only rows related to the given patient and that contain shared
        # mutations.
        dplyr::filter(sample_interpretable == patient & mutation_type == "Shared") |>
        # Subset the dataframe again to keep only rows that contain mutations exclusively shared between patient and PDX-X1.
        dplyr::filter(sample_ids_mutation_is_present %in% c(paste0(patient, ", ", x1_id), paste0(x1_id, ", ", patient))) |>
        # Count the mutations.
        dplyr::count(mutation_type) |>
        dplyr::select(n) |>
        as.numeric()
    # If there are none, save that as a zero.
    if (is.na(n_shared_patient_and_x1)) n_shared_patient_and_x1 <- 0
    
    ##################################################################################
    #### Variants shared between patient, PDX-X1 and other PDX(s), but not public ####
    ##################################################################################
    # These mutations contribute towards the total number of variants that are shared between patient and X1!
    
    n_shared_patient_x1_and_others_but_not_public <- 0
    
    sub_df <- patient_and_x1 |>
      # Subset the dataframe to keep only rows related to the given patient and that contain shared
      # mutations.
      dplyr::filter(sample_interpretable == patient & mutation_type == "Shared")
    
    # Each row is a different shared mutation for a set of samples; iterate over each...
    for (sample_list in sub_df$sample_ids_mutation_is_present) {
      
      sample_ids <- strsplit(sample_list, ",\\s*")[[1]]
      
      # If both patient and PDX-X1 are present in the set of samples...
      if (patient %in% sample_ids && x1_id %in% sample_ids) {
        
        # If there are more than 2 samples in the set of samples, considering
        # that both patient and PDX-X1 have been spotted, that means the mutation
        # is also shared with another sample(s) but is not public.
        if (length(sample_ids) > 2) {
          n_shared_patient_x1_and_others_but_not_public <- n_shared_patient_x1_and_others_but_not_public + 1
        }
      }
    }
    
    ###################################################################################
    #### Variants shared between patient and other PDX(s), but not found in PDX-X1 ####
    ###################################################################################
    # These mutations contribute towards the total number of variants that are present in the patient!
    
    n_shared_patient_and_others_but_not_in_x1 <- 0
    
    sub_df <- patient_and_x1 |>
      # Subset the dataframe to keep only rows related to the given patient and that contain shared
      # mutations.
      dplyr::filter(sample_interpretable == patient & mutation_type == "Shared")
    
    # Each row is a different shared mutation for a set of samples; iterate over each...
    for (sample_list in sub_df$sample_ids_mutation_is_present) {
      
      sample_ids <- strsplit(sample_list, ",\\s*")[[1]]
      
      # If the patient is present in the set of samples, but PDX-X1 is not...
      if (patient %in% sample_ids && !(x1_id %in% sample_ids)) {
        n_shared_patient_and_others_but_not_in_x1 <- n_shared_patient_and_others_but_not_in_x1 + 1
      }
    }
    
    #######################################################################################
    #### Variants shared between PDX-X1 and other PDX(s), but not found in the patient ####
    #######################################################################################
    # These mutations contribute towards the total number of variants that are present in PDX-X1!
    
    n_shared_x1_and_others_but_not_in_patient <- 0
    
    sub_df <- patient_and_x1 |>
      # Subset the dataframe to keep only rows related to the given patient and that contain shared
      # mutations.
      dplyr::filter(sample_interpretable == x1_id & mutation_type == "Shared")
    
    # Each row is a different shared mutation for a set of samples; iterate over each...
    for (sample_list in sub_df$sample_ids_mutation_is_present) {
      
      sample_ids <- strsplit(sample_list, ",\\s*")[[1]]
      
      # If PDX-X1 is present in the set of samples, but the patient is not...
      if (!(patient %in% sample_ids) && x1_id %in% sample_ids) {
        n_shared_x1_and_others_but_not_in_patient <- n_shared_x1_and_others_but_not_in_patient + 1
      }
    }
    
    #########################
    #### Public variants ####
    #########################
      
    n_public <- patient_and_x1 |>
        dplyr::filter(sample_interpretable == patient & mutation_type == "Public") |>
        dplyr::count(mutation_type) |>
        dplyr::select(n) |>
        as.numeric()
    if (is.na(n_public)) n_public <- 0
    
    #####################################################
    
    n_total_found_in_patient_but_not_in_x1 <- (n_private_patient + n_shared_patient_and_others_but_not_in_x1)
    n_total_found_in_x1_but_not_in_patient <- (n_private_x1 + n_shared_x1_and_others_but_not_in_patient)
    n_total_shared_between_patient_and_x1 <- (n_shared_patient_and_x1 + n_shared_patient_x1_and_others_but_not_public + n_public)
    
    #####################################################
  
    print("Appending to detailed dataframe")
    
    patient_vs_x1_df_detailed <- rbind(
        patient_vs_x1_df_detailed,
        c(
          patient, 
          n_pdx_for_patient, 
          paste0(pdx_ids_for_patient, collapse = ","),
          x1_id, 
          n_private_patient, 
          n_private_x1, 
          n_shared_patient_and_x1,
          n_shared_x1_and_others_but_not_in_patient,
          n_shared_patient_and_others_but_not_in_x1,
          n_shared_patient_x1_and_others_but_not_public,
          n_public
        )
    )
    
    #####################################################
    
    print("Appending to simplified dataframe")
    
    patient_vs_x1_df_simplified <- rbind(
      patient_vs_x1_df_simplified,
      c(
        patient, 
        x1_id, 
        n_total_found_in_patient_but_not_in_x1, 
        n_total_found_in_x1_but_not_in_patient, 
        n_total_shared_between_patient_and_x1
      )
    )
}

colnames(patient_vs_x1_df_detailed) <- c(
    "patient",
    "n_pdx_for_patient",
    "pdx_sample_ids",
    "x1",
    "n_private_patient",
    "n_private_x1",
    "n_shared_patient_and_x1",
    "n_shared_x1_and_others_but_not_in_patient",
    "n_shared_patient_and_others_but_not_in_x1",
    "n_shared_patient_x1_and_others_but_not_public",
    "n_public"
)

colnames(patient_vs_x1_df_simplified) <- c(
  "patient",
  "x1",
  "n_variants_patient",
  "n_variants_x1",
  "n_variants_shared"
)

write.csv(patient_vs_x1_df_detailed, paste0(meskit_outdir, "patient_vs_x1_detailed.csv"))
write.csv(patient_vs_x1_df_simplified, paste0(meskit_outdir, "patient_vs_x1_simplified.csv"))

p1 <- patient_vs_x1_df_detailed |>
  dplyr::select(!c(n_pdx_for_patient, pdx_sample_ids)) |>
  tidyr::pivot_longer(
    cols = starts_with("n_"),
    names_to = "category",
    values_to = "n_variants"
  ) |>
  ggplot(
    aes(x = patient, y = as.numeric(n_variants), fill = category)
  ) +
    geom_bar(stat = "identity", position = "stack") +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 10, angle = 45, vjust = 1.075, hjust = 1),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      plot.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12)
    ) +
  ggtitle("Patient vs PDX-X1, detailed") +
  xlab("Patient") +
  ylab("# Variants") +
  scale_fill_brewer(
    palette = "Paired",
    labels = c(
      "Only in patient",
      "Only in PDX-X1",
      "Public",
      "Shared by patient and PDX(s) other than PDX-X1",
      "Shared by patient and PDX-X1",
      "Shared by patient, PDX-X1 and other PDX(s) but not public",
      "Shared by PDX-X1 and other PDX(s) but not in patient"
    ),
    name = "Category"
  )
  

p2 <- patient_vs_x1_df_simplified |>
  tidyr::pivot_longer(
    cols = starts_with("n_"),
    names_to = "category",
    values_to = "n_variants"
  ) |>
  ggplot(
    aes(x = patient, y = as.numeric(n_variants), fill = category)
  ) +
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, vjust = 1.075, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.justification = "left"
  ) +
  ggtitle("Patient vs PDX-X1, simplified") +
  xlab("Patient") +
  ylab("# Variants") +
  scale_fill_manual(
    values = c("#B2DF8A", "#1F78B4", "#E31A1C"),
    labels = c("Patient", "Shared", "PDX-X1"),
    name = "Category"
  )

plots_combined <- cowplot::plot_grid(
  plotlist = list(p1, p2),
  ncol = 1,
  align = "hv"
)

ggsave(
  plot = plots_combined,
  filename = paste0(meskit_outdir, "patient_vs_x1_plots.png"),
  dpi = 600,
  width = 10.0,
  height = 8.0
)
  
    




