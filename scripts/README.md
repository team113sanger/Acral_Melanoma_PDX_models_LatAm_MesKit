# Generating CNA heatmaps with MesKit

## Sequenza

### Adjusted copy number results

```bash
PROJDIR="/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES"
DERMATLAS_RSCRIPT="/software/team113/dermatlas/R/R-4.2.2/bin/Rscript"

cd "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/MesKit_analysis_adjusted"
${DERMATLAS_RSCRIPT} ${PROJDIR}/scripts/MESKIT/acral-melanoma-meskit/scripts/cna_analysis.R \
    --segments "adjusted_segments_for_meskit.tsv" \
    --manifest "${PROJDIR}/analysis/release_v2/combined_2729_3248/2729_3248_combined_manifest.txt" \
    --outdir "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/MesKit_analysis_adjusted"
```

### Not adjusted copy number results

```bash
PROJDIR="/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES"
DERMATLAS_RSCRIPT="/software/team113/dermatlas/R/R-4.2.2/bin/Rscript"

cd "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/MesKit_analysis_not_adjusted"
${DERMATLAS_RSCRIPT} ${PROJDIR}/scripts/MESKIT/acral-melanoma-meskit/scripts/cna_analysis.R \
    --segments "not_adjusted_segments_for_meskit.tsv" \
    --manifest "${PROJDIR}/analysis/release_v2/combined_2729_3248/2729_3248_combined_manifest.txt" \
    --outdir "${PROJDIR}/analysis/Sequenza/BrazilianPop_PairedSamples_Hum-PDX_SequenzaFiles_Jan2025/MesKit_analysis_not_adjusted"
```

## ASCAT

### Adjusted copy number results

```bash
PROJDIR="/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES"
DERMATLAS_RSCRIPT="/software/team113/dermatlas/R/R-4.2.2/bin/Rscript"

cd "${PROJDIR}/analysis/ASCAT/MesKit_analysis_adjusted"
${DERMATLAS_RSCRIPT} ${PROJDIR}/scripts/MESKIT/acral-melanoma-meskit/scripts/cna_analysis.R \
    --segments "adjusted_segments_for_meskit.tsv" \
    --manifest "${PROJDIR}/analysis/release_v2/combined_2729_3248/2729_3248_combined_manifest.txt" \
    --outdir "${PROJDIR}/analysis/ASCAT/MesKit_analysis_adjusted"
```

### Not adjusted copy number results

```bash
PROJDIR="/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES"
DERMATLAS_RSCRIPT="/software/team113/dermatlas/R/R-4.2.2/bin/Rscript"

cd "${PROJDIR}/analysis/ASCAT/MesKit_analysis_not_adjusted"
${DERMATLAS_RSCRIPT} ${PROJDIR}/scripts/MESKIT/acral-melanoma-meskit/scripts/cna_analysis.R \
    --segments "not_adjusted_segments_for_meskit.tsv" \
    --manifest "${PROJDIR}/analysis/release_v2/combined_2729_3248/2729_3248_combined_manifest.txt" \
    --outdir "${PROJDIR}/analysis/ASCAT/MesKit_analysis_not_adjusted"
```