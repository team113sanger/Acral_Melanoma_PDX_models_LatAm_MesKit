**Generating segment files for 1stModel samples**

```bash
PROJDIR="/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES"
DERMATLAS_RSCRIPT="/software/team113/dermatlas/R/R-4.2.2/bin/Rscript"

# Copy number results not adjusted
${DERMATLAS_RSCRIPT} "${PROJDIR}/scripts/MESKIT/acral-melanoma-meskit/scripts/copy_number/Sequenza/prepare_sequenza_for_meskit.R" \
    --human "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/Braz_Sequenza_1stModelsAll_Human" \
    --pdx "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/Braz_Sequenza_1stModelsAll_PDXs" \
    --outdir "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/MesKit_analysis_not_adjusted" \
    --prefix "not_adjusted_1stModels"

# Copy number results adjusted
${DERMATLAS_RSCRIPT} "${PROJDIR}/scripts/MESKIT/acral-melanoma-meskit/scripts/copy_number/Sequenza/prepare_sequenza_for_meskit.R" \
    --human "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/Braz_Sequenza_1stModelsAll_Human" \
    --pdx "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/Braz_Sequenza_1stModelsAll_PDXs" \
    --outdir "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/MesKit_analysis_adjusted" \
    --prefix "adjusted_1stModels" \
    --ploidy "${PROJDIR}/scripts/MESKIT/acral-melanoma-meskit/data/Sequenza/ploidy_table_jan2025_data.tsv"
```

**Generating segment files for AltModel samples**

```bash
PROJDIR="/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES"
DERMATLAS_RSCRIPT="/software/team113/dermatlas/R/R-4.2.2/bin/Rscript"

# Copy number results not adjusted
${DERMATLAS_RSCRIPT} "${PROJDIR}/scripts/MESKIT/acral-melanoma-meskit/scripts/copy_number/Sequenza/prepare_sequenza_for_meskit.R" \
    --human "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/Braz_Sequenza_AltModels_Human" \
    --pdx "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/Braz_Sequenza_AltModels_PDXs" \
    --outdir "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/MesKit_analysis_not_adjusted" \
    --prefix "not_adjusted_AltModels"

# Copy number results adjusted
${DERMATLAS_RSCRIPT} "${PROJDIR}/scripts/MESKIT/acral-melanoma-meskit/scripts/copy_number/Sequenza/prepare_sequenza_for_meskit.R" \
    --human "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/Braz_Sequenza_AltModels_Human" \
    --pdx "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/Braz_Sequenza_AltModels_PDXs" \
    --outdir "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/MesKit_analysis_adjusted" \
    --prefix "adjusted_AltModels" \
    --ploidy "${PROJDIR}/scripts/MESKIT/acral-melanoma-meskit/data/Sequenza/ploidy_table_jan2025_data.tsv"
```

**Concatenating 1stModel and AltModel samples**

```bash
PROJDIR="/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES"

# Copy number results not adjusted
cd "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/MesKit_analysis_not_adjusted"
cat not_adjusted_1stModels_segments_for_meskit.tsv > not_adjusted_segments_for_meskit.tsv
tail --lines=+2 not_adjusted_AltModels_segments_for_meskit.tsv >> not_adjusted_segments_for_meskit.tsv

# Copy number results adjusted
cd "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/MesKit_analysis_adjusted"
cat adjusted_1stModels_segments_for_meskit.tsv > adjusted_segments_for_meskit.tsv
tail --lines=+2 adjusted_AltModels_segments_for_meskit.tsv >> adjusted_segments_for_meskit.tsv
```
