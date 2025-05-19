files <- fs::dir_ls(
    "/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES/analysis/somatic_variants/release_v4/meskit/plots",
    recurse = TRUE,
    glob = "*_filt.csv*"
)

all_samples <- readr::read_csv(files)
all_samples <- all_samples |> dplyr::select(!`...1`)

patient_ids <- all_samples$sample_interpretable[!grepl("-X", all_samples$sample_interpretable)] |> unique()


patient_vs_x1_df <- data.frame(matrix(nrow = 0, ncol = 7))
for (patient in patient_ids) {
    n_pdx_for_patient <- all_samples$sample_interpretable[grepl(paste0(patient, "-X"), all_samples$sample_interpretable)] |>
        unique() |>
        length()

    x1_id <- paste0(patient, "-X1_PDX")
    if (!(x1_id %in% all_samples$sample_interpretable)) next

    patient_and_x1 <- all_samples |> dplyr::filter(sample_interpretable %in% c(patient, x1_id))

    n_private_patient <- patient_and_x1 |>
        dplyr::filter(sample_interpretable == patient & mutation_type == "Private") |>
        dplyr::count(mutation_type) |>
        dplyr::select(n) |>
        as.numeric()
    if (is.na(n_private_patient)) n_private_patient <- 0

    n_private_x1 <- patient_and_x1 |>
        dplyr::filter(sample_interpretable == x1_id & mutation_type == "Private") |>
        dplyr::count(mutation_type) |>
        dplyr::select(n) |>
        as.numeric()
    if (is.na(n_private_x1)) n_private_x1 <- 0

    n_shared_patient_and_x1 <- patient_and_x1 |>
        dplyr::filter(sample_interpretable == patient & mutation_type == "Shared") |>
        dplyr::filter(sample_ids_mutation_is_present %in% c(paste0(patient, ", ", x1_id), paste0(x1_id, ", ", patient))) |>
        dplyr::count(mutation_type) |>
        dplyr::select(n) |>
        as.numeric()
    if (is.na(n_shared_patient_and_x1)) n_shared_patient_and_x1 <- 0

    n_public <- patient_and_x1 |>
        dplyr::filter(sample_interpretable == patient & mutation_type == "Public") |>
        dplyr::count(mutation_type) |>
        dplyr::select(n) |>
        as.numeric()
    if (is.na(n_public)) n_public <- 0

    n_total <- (n_private_patient + n_private_x1 + n_shared_patient_and_x1 + n_public)

    patient_vs_x1_df <- rbind(
        patient_vs_x1_df,
        c(patient, n_pdx_for_patient, x1_id, n_private_patient, n_private_x1, n_shared_patient_and_x1, n_public, n_total)
    )
}
colnames(patient_vs_x1_df) <- c("patient", "n_pdx_for_patient", "x1", "n_private_patient", "n_private_x1", "n_shared_patient_x1", "n_public", "n_total")

write.csv(
    patient_vs_x1_df,
    "/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES/analysis/somatic_variants/release_v4/meskit/plots/patient_vs_x1.csv"
)
