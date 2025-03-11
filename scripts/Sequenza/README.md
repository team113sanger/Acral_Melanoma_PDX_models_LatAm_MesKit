**Generating segment files for 1stModel samples**

```bash
PROJDIR="/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES"
DERMATLAS_RSCRIPT="/software/team113/dermatlas/R/R-4.2.2/bin/Rscript"

# Without normalisation
${DERMATLAS_RSCRIPT} "${PROJDIR}/scripts/MESKIT/acral-melanoma-meskit/scripts/Sequenza/prepare_sequenza_for_meskit.R" \
    --human "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/Braz_Sequenza_1stModelsAll_Human" \
    --pdx "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/Braz_Sequenza_1stModelsAll_PDXs" \
    --outdir "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/MesKit_analysis_without_normalisation" \
    --prefix "no_normalisation_1stModels"

# With normalisation
${DERMATLAS_RSCRIPT} "${PROJDIR}/scripts/MESKIT/acral-melanoma-meskit/scripts/Sequenza/prepare_sequenza_for_meskit.R" \
    --human "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/Braz_Sequenza_1stModelsAll_Human" \
    --pdx "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/Braz_Sequenza_1stModelsAll_PDXs" \
    --outdir "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/MesKit_analysis_with_normalisation" \
    --prefix "normalised_1stModels" \
    --ploidy "${PROJDIR}/scripts/MESKIT/acral-melanoma-meskit/data/Sequenza/ploidy_table_jan2025_data.tsv"
```

**Generating segment files for AltModel samples**

```bash
PROJDIR="/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES"
DERMATLAS_RSCRIPT="/software/team113/dermatlas/R/R-4.2.2/bin/Rscript"

# Without normalisation
${DERMATLAS_RSCRIPT} "${PROJDIR}/scripts/MESKIT/acral-melanoma-meskit/scripts/Sequenza/prepare_sequenza_for_meskit.R" \
    --human "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/Braz_Sequenza_AltModels_Human" \
    --pdx "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/Braz_Sequenza_AltModels_PDXs" \
    --outdir "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/MesKit_analysis_without_normalisation" \
    --prefix "no_normalisation_AltModels"

# With normalisation
${DERMATLAS_RSCRIPT} "${PROJDIR}/scripts/MESKIT/acral-melanoma-meskit/scripts/Sequenza/prepare_sequenza_for_meskit.R" \
    --human "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/Braz_Sequenza_AltModels_Human" \
    --pdx "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/Braz_Sequenza_AltModels_PDXs" \
    --outdir "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/MesKit_analysis_with_normalisation" \
    --prefix "normalised_AltModels"
```

**Concatenating 1stModel and AltModel samples**

```bash
PROJDIR="/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES"

# Results without normalisation
cd "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/MesKit_analysis_without_normalisation"
cat no_normalisation_1stModels_segments_for_meskit.tsv > no_normalisation_segments_for_meskit.tsv
tail --lines=+2 no_normalisation_AltModels_segments_for_meskit.tsv >> no_normalisation_segments_for_meskit.tsv

# Normalised results
cd "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/MesKit_analysis_with_normalisation"
cat normalised_1stModels_segments_for_meskit.tsv > normalised_segments_for_meskit.tsv
tail --lines=+2 normalised_AltModels_segments_for_meskit.tsv >> normalised_segments_for_meskit.tsv
```

**Generating heatmaps with MesKit**
```bash
PROJDIR="/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES"
DERMATLAS_RSCRIPT="/software/team113/dermatlas/R/R-4.2.2/bin/Rscript"

# Results without normalisation
cd "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/MesKit_analysis_without_normalisation"
${DERMATLAS_RSCRIPT} ${PROJDIR}/scripts/MESKIT/acral-melanoma-meskit/scripts/Sequenza/cna_analysis.R \
    --segments "no_normalisation_segments_for_meskit.tsv" \
    --manifest "${PROJDIR}/analysis/release_v2/combined_2729_3248/2729_3248_combined_manifest.txt" \
    --outdir "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/MesKit_analysis_without_normalisation"

# Normalised results
cd "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/MesKit_analysis_with_normalisation"
${DERMATLAS_RSCRIPT} ${PROJDIR}/scripts/MESKIT/acral-melanoma-meskit/scripts/Sequenza/cna_analysis.R \
    --segments "normalised_segments_for_meskit.tsv" \
    --manifest "${PROJDIR}/analysis/release_v2/combined_2729_3248/2729_3248_combined_manifest.txt" \
    --outdir "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/MesKit_analysis_with_normalisation"
```

