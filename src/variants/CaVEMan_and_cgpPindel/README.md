# Generating phylogenetic trees with MesKit

```bash
DATADIR="/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES/analysis/somatic_variants/release_v4"
MAF="${DATADIR}/Combined_2729_3248/matched_samples/6633_2729_3248-filtered_mutations_matched_allTum_keep.maf"
PREFIX="6633_2729_3248-filtered_mutations_matched_allTum_keep"
OUTDIR="${DATADIR}/meskit"

bash adjust_maf_for_meskit.sh -m "${MAF}" -p "${PREFIX}" -o "${OUTDIR}"
```

```bash
PROJDIR="/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES"
RSCRIPT="/software/team113/dermatlas/R/R-4.2.2/bin/Rscript"
SCRIPTS_DIR="${PROJDIR}/scripts/MESKIT/acral-melanoma-meskit/scripts"

DATADIR="/lustre/scratch124/casm/team113/projects/6633_PDX_models_Latin_America_WES/analysis/somatic_variants/release_v4"
MAF="${DATADIR}/meskit/6633_2729_3248-filtered_mutations_matched_allTum_keep_adjusted_for_meskit.maf"
CLINICAL="${DATADIR}/meskit/6633_2729_3248-filtered_mutations_matched_allTum_keep_clinical_for_meskit.tsv"
MDATA="${DATADIR}/Combined_2729_3248/6633_METADATA_table.tsv"
COSMIC="${DATADIR}/Combined_2729_3248/COSMIC/cancer_gene_census.v97.csv"
OUTDIR="${DATADIR}/meskit"

"${RSCRIPT}" "${SCRIPTS_DIR}/SNVs_and_indels/CaVEMan_and_cgpPindel/make_trees.R" \
    --maf "${MAF}" \
    --clinical "${CLINICAL}" \
    --metadata "${MDATA}" \
    --cosmic "${COSMIC}" \
    --outdir "${OUTDIR}" > ${OUTDIR}/meskit_log.txt 2>&1
```
