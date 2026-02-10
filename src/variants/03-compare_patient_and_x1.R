library(here)
library(fs)
library(readr)
library(dplyr)

# Define the directory with MesKit outputs for variant calling data and load
# filtered results.
meskit_outdir <- here::here("results/variants/")
files <- fs::dir_ls(meskit_outdir, recurse = TRUE, glob = "*_filt.csv$")
all_samples <- readr::read_csv(files)
all_samples <- all_samples |> dplyr::select(!`...1`)

# Get patient IDs.
patient_ids <- all_samples$sample_interpretable[!grepl("-X", all_samples$sample_interpretable)] |> unique()

# Create empty dataframe to store data.
patient_vs_x1_df <- data.frame(matrix(nrow = 0, ncol = 9))

# Iterate over patients...
for (patient in patient_ids) {
  
    # Get the sample IDs for PDXs available for the given patient.
    pdx_ids_for_patient <- all_samples$sample_interpretable[grepl(paste0(patient, "-X"), all_samples$sample_interpretable)] |>
        unique()
  
    # Get the number of PDX samples available for the given patient.
    n_pdx_for_patient <- length(pdx_ids_for_patient)
    
    # Set the sample ID for PDX-X1.
    x1_id <- paste0(patient, "-X1_PDX")
    
    # If there isn't a PDX-X1 for the given patient, move on to the next patient.
    if (!(x1_id %in% all_samples$sample_interpretable)) next

    # Subset original dataframe to keep only results related to the given patient and its PDX-X1.
    patient_and_x1 <- all_samples |> dplyr::filter(sample_interpretable %in% c(patient, x1_id))

    # Count number of private mutations in the patient (i.e. those that are only present in the patient).
    n_private_patient <- patient_and_x1 |>
        # Subset the dataframe to keep only rows related to the given patient and that contain private
        # mutations (only present in the patient).
        dplyr::filter(sample_interpretable == patient & mutation_type == "Private") |>
        # Count the mutations.
        dplyr::count(mutation_type) |>
        dplyr::select(n) |>
        as.numeric()
    # If there are no private mutations for the patient, save that as a zero.
    if (is.na(n_private_patient)) n_private_patient <- 0

    # Count number of private mutations in the PDX-X1 (i.e. those that are only present in the PDX-X1).
    n_private_x1 <- patient_and_x1 |>
        # Subset the dataframe to keep only rows related to the given PDX-X1 and that contain private
        # mutations (only present in PDX-X1).
        dplyr::filter(sample_interpretable == x1_id & mutation_type == "Private") |>
        # Count the mutations.
        dplyr::count(mutation_type) |>
        dplyr::select(n) |>
        as.numeric()
    # If there are no private mutations for PDX-X1, save that as a zero.
    if (is.na(n_private_x1)) n_private_x1 <- 0

    # Count number of mutations that are exclusively shared between the patient and PDX-X1.
    n_shared_patient_and_x1 <- patient_and_x1 |>
        # Subset the dataframe to keep only rows related to the given patient and that contain shared
        # mutations.
        dplyr::filter(sample_interpretable == patient & mutation_type == "Shared") |>
        # Subset the dataframe again to keep only rows that contain mutations exclusively shared between patient and PDX-X1.
        dplyr::filter(sample_ids_mutation_is_present %in% c(paste0(patient, ", ", x1_id), paste0(x1_id, ", ", patient))) |>
        # Count the mutations.
        dplyr::count(mutation_type) |>
        dplyr::select(n) |>
        as.numeric()
    # If there are none, save that as a zero.
    if (is.na(n_shared_patient_and_x1)) n_shared_patient_and_x1 <- 0
    
    n_shared_patient_x1_and_others_but_not_public <- 0
    
    sub_df <- patient_and_x1 |>
      # Subset the dataframe to keep only rows related to the given patient and that contain shared
      # mutations.
      dplyr::filter(sample_interpretable == patient & mutation_type == "Shared")
    
    # Each row is a different shared mutation for a set of samples; iterate over each...
    for (sample_list in sub_df$sample_ids_mutation_is_present) {
      
      # If both patient and PDX-X1 are present in the set of samples...
      if (grepl(patient, sample_list) && grepl(x1_id, sample_list)) {
        
        sample_ids <- as.vector(strsplit(sample_list, ",\\s*"))
        
        # If there are more than 2 samples in the set of samples, considering
        # that both patient and PDX-X1 have been spotted, that means the mutation
        # is also shared with another sample(s) but is not public.
        if (length(sample_ids[[1]]) > 2) {
          n_shared_patient_x1_and_others_but_not_public <- n_shared_patient_x1_and_others_but_not_public + 1
        }
      }
    }
      
    n_public <- patient_and_x1 |>
        dplyr::filter(sample_interpretable == patient & mutation_type == "Public") |>
        dplyr::count(mutation_type) |>
        dplyr::select(n) |>
        as.numeric()
    if (is.na(n_public)) n_public <- 0

    n_total <- (n_private_patient + n_private_x1 + n_shared_patient_and_x1 + n_shared_patient_x1_and_others_but_not_public + n_public)

    patient_vs_x1_df <- rbind(
        patient_vs_x1_df,
        c(
          patient, 
          n_pdx_for_patient, 
          paste0(pdx_ids_for_patient, collapse = ","),
          x1_id, 
          n_private_patient, 
          n_private_x1, 
          n_shared_patient_and_x1, 
          n_shared_patient_x1_and_others_but_not_public, 
          n_public, 
          n_total
        )
    )
}
colnames(patient_vs_x1_df) <- c(
    "patient",
    "n_pdx_for_patient",
    "pdx_sample_ids",
    "x1",
    "n_private_patient",
    "n_private_x1",
    "n_shared_patient_x1",
    "n_shared_patient_x1_and_other_but_not_public",
    "n_public",
    "n_total"
)

write.csv(patient_vs_x1_df, paste0(meskit_outdir, "patient_vs_x1.csv"))
