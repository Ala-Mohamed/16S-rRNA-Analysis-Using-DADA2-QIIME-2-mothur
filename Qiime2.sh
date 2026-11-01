cat << 'EOF' > run_mocks.sh
#!/usr/bin/env bash
set -e

WORKDIR="/home/alaa/Downloads"
THREADS=4

cd "$WORKDIR"

# Process Mock 18, Mock 16, and Mock Default
MOCKS=(
  " 18|mock_18"
  " 16|mock_16"
  "|mock_default"
)

for ENTRY in "${MOCKS[@]}"; do
  IFS="|" read -r SUFFIX LABEL <<< "$ENTRY"

  if [ -z "$SUFFIX" ]; then
    EXPECTED_FASTA="${WORKDIR}/expected-sequences.fasta"
    FWD_FASTQ="${WORKDIR}/mock-forward-read.fastq.gz"
    REV_FASTQ="${WORKDIR}/mock-reverse-read.fastq.gz"
  else
    EXPECTED_FASTA="${WORKDIR}/expected-sequences${SUFFIX}.fasta"
    FWD_FASTQ="${WORKDIR}/mock-forward-read.fastq${SUFFIX}.gz"
    REV_FASTQ="${WORKDIR}/mock-reverse-read.fastq${SUFFIX}.gz"
  fi

  OUT_DIR="${WORKDIR}/results_${LABEL}"
  mkdir -p "${OUT_DIR}"

  echo "=================================================="
  echo " Processing Community: ${LABEL}"
  echo "=================================================="

  # Step 1: Manifest & Import
  MANIFEST="${OUT_DIR}/manifest.tsv"
  if [ ! -f "${OUT_DIR}/demux.qza" ]; then
    echo "--> Writing manifest and importing raw reads..."
    printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > "${MANIFEST}"
    printf "%s\t%s\t%s\n" "${LABEL}" "${FWD_FASTQ}" "${REV_FASTQ}" >> "${MANIFEST}"

    qiime tools import \
      --type 'SampleData[PairedEndSequencesWithQuality]' \
      --input-path "${MANIFEST}" \
      --output-path "${OUT_DIR}/demux.qza" \
      --input-format PairedEndFastqManifestPhred33V2
  else
    echo "--> demux.qza already exists. Skipping import..."
  fi

  # Step 2: Merge Paired-End Reads
  if [ ! -f "${OUT_DIR}/joined.qza" ]; then
    echo "--> Merging paired-end reads..."
    qiime vsearch merge-pairs \
      --i-demultiplexed-seqs "${OUT_DIR}/demux.qza" \
      --o-merged-sequences "${OUT_DIR}/joined.qza" \
      --o-unmerged-sequences "${OUT_DIR}/unmerged.qza"
  else
    echo "--> joined.qza already exists. Skipping merging..."
  fi

  # Step 3: Dereplicate Reads
  if [ ! -f "${OUT_DIR}/derep-table.qza" ]; then
    echo "--> Dereplicating merged reads..."
    qiime vsearch dereplicate-sequences \
      --i-sequences "${OUT_DIR}/joined.qza" \
      --o-dereplicated-table "${OUT_DIR}/derep-table.qza" \
      --o-dereplicated-sequences "${OUT_DIR}/derep-seqs.qza"
  else
    echo "--> derep-table.qza already exists. Skipping dereplication..."
  fi

  # Step 4: De Novo Chimera Filtering
  if [ ! -f "${OUT_DIR}/nonchimeras-table.qza" ]; then
    echo "--> Identifying chimeras with uchime-denovo..."
    qiime vsearch uchime-denovo \
      --i-table "${OUT_DIR}/derep-table.qza" \
      --i-sequences "${OUT_DIR}/derep-seqs.qza" \
      --o-chimeras "${OUT_DIR}/chimeras.qza" \
      --o-nonchimeras "${OUT_DIR}/nonchimeras-seqs.qza" \
      --o-stats "${OUT_DIR}/uchime-stats.qza"

    echo "--> Filtering chimeras out of the feature table..."
    qiime feature-table filter-features \
      --i-table "${OUT_DIR}/derep-table.qza" \
      --m-metadata-file "${OUT_DIR}/nonchimeras-seqs.qza" \
      --o-filtered-table "${OUT_DIR}/nonchimeras-table.qza"
  else
    echo "--> nonchimeras-table.qza already exists. Skipping chimera filtering..."
  fi

  # Step 5: Import Reference Sequences
  if [ ! -f "${OUT_DIR}/expected-mock-ref.qza" ]; then
    echo "--> Importing reference sequences..."
    qiime tools import \
      --type 'FeatureData[Sequence]' \
      --input-path "${EXPECTED_FASTA}" \
      --output-path "${OUT_DIR}/expected-mock-ref.qza"
  else
    echo "--> expected-mock-ref.qza already exists. Skipping..."
  fi

  # Step 6: VSEARCH Open-Reference OTU Clustering
  if [ ! -f "${OUT_DIR}/vsearch-otu-table.qza" ]; then
    echo "--> Running VSEARCH open-reference OTU clustering on non-chimeric reads..."
    qiime vsearch cluster-features-open-reference \
      --i-table "${OUT_DIR}/nonchimeras-table.qza" \
      --i-sequences "${OUT_DIR}/nonchimeras-seqs.qza" \
      --i-reference-sequences "${OUT_DIR}/expected-mock-ref.qza" \
      --p-perc-identity 0.97 \
      --p-threads ${THREADS} \
      --o-clustered-table "${OUT_DIR}/vsearch-otu-table.qza" \
      --o-clustered-sequences "${OUT_DIR}/vsearch-otu-seqs.qza" \
      --o-new-reference-sequences "${OUT_DIR}/vsearch-unmatched-seqs.qza"
  else
    echo "--> vsearch-otu-table.qza already exists. Skipping VSEARCH clustering..."
  fi

  # Step 7: dbOTU Clustering
  if [ ! -f "${OUT_DIR}/dbotu-table.qza" ]; then
    echo "--> Running dbOTU clustering on non-chimeric reads..."
    qiime dbotu call-otus \
      --i-table "${OUT_DIR}/nonchimeras-table.qza" \
      --i-sequences "${OUT_DIR}/nonchimeras-seqs.qza" \
      --p-gen-crit 0.10 \
      --p-abund-crit 10.0 \
      --p-pval-crit 0.0001 \
      --o-dbotu-table "${OUT_DIR}/dbotu-table.qza" \
      --o-representative-sequences "${OUT_DIR}/dbotu-seqs.qza"
  else
    echo "--> dbotu-table.qza already exists. Skipping dbOTU..."
  fi

  # Step 8: Visualizations
  echo "--> Generating summary visualizations..."
  qiime feature-table summarize \
    --i-table "${OUT_DIR}/vsearch-otu-table.qza" \
    --o-visualization "${OUT_DIR}/vsearch-otu-table.qzv"

  qiime feature-table summarize \
    --i-table "${OUT_DIR}/dbotu-table.qza" \
    --o-visualization "${OUT_DIR}/dbotu-table.qzv"

  echo "Finished processing ${LABEL}."
  echo ""
done

echo "=== Processing complete for all mock samples ==="
EOF