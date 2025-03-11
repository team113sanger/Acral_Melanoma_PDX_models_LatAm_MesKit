 **Generating segment files for MesKit with ASCAT outputs**

```bash
PROJDIR="/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES"
DERMATLAS_RSCRIPT="/software/team113/dermatlas/R/R-4.2.2/bin/Rscript"

# Copy number results not adjusted
${DERMATLAS_RSCRIPT} "${PROJDIR}/scripts/MESKIT/acral-melanoma-meskit/scripts/ASCAT/prepare_ascat_for_meskit.R" \
    --segments "${PROJDIR}/analysis/ASCAT/release_v1/all_samples/plots/all_samp_segments.tsv" \
    --outdir "${PROJDIR}/analysis/ASCAT/MesKit_analysis_not_adjusted" \
    --prefix "not_adjusted"

# Copy number results adjusted
${DERMATLAS_RSCRIPT} "${PROJDIR}/scripts/MESKIT/acral-melanoma-meskit/scripts/ASCAT/prepare_ascat_for_meskit.R" \
    --segments "${PROJDIR}/analysis/ASCAT/release_v1/all_samples/plots/all_samp_segments.tsv" \
    --outdir "${PROJDIR}/analysis/ASCAT/MesKit_analysis_adjusted" \
    --prefix "adjusted" \
    --ploidy
```
