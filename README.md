# 16S rRNA Analysis Using DADA2, QIIME 2 & mothur

A benchmarking project comparing 16S rRNA amplicon-sequencing bioinformatics pipelines — **DADA2** (ASV-based denoising), **QIIME 2**, and **mothur** (OTU-based clustering) — to evaluate differences in their algorithms, outputs, and practical usability. Mock microbial communities with known composition are used as ground truth, with results cross-converted into a common **mothur-style count-table format** so outputs from all three pipelines can be compared directly, like for like.

## Data

Mock communities **Mock-16**, **Mock-18**, and **Mock-19**, obtained from the public [mockrobiota](https://github.com/caporaso-lab/mockrobiota/tree/master/data) repository.

> ⚠️ **Primers:** Pay attention to primer sequences before processing — raw reads must be checked for residual PCR primers (e.g., via `trimmid_primers.sh`) before filtering/clustering, as untrimmed primers can inflate apparent sequence diversity and cause reference-matching failures.

## Repository Structure

| File / Folder | Purpose |
|---|---|
| `DADA2/` | DADA2 (R) pipeline scripts for ASV-based processing of the mock datasets — filtering, error learning, denoising, merging, and chimera removal. |
| `Demultiplex1.sh` | Prepares/pairs raw FASTQ files per sample ahead of processing. |
| `trimmid_primers.sh` | Removes 16S PCR primer sequences from raw reads prior to downstream filtering and clustering. |
| `script uparse.sh` | Runs UPARSE-based OTU clustering (usearch) on quality-filtered, primer-trimmed reads, as an additional OTU-based comparison point alongside mothur's own clustering. |
| `Qiime2.sh` | QIIME 2 command-line pipeline — import, denoising, and feature-table generation. |
| `run_mocks.sh` | Orchestrator script that runs the full pipeline across Mock-16, Mock-18, and Mock-19. |
| `generate_mothur.py` | Converts pipeline output (OTU/ASV table + representative sequences) into a mothur-style `make.shared`-equivalent count table for cross-pipeline comparison. |
| `fix_mothur.py` | Cleans/corrects formatting issues in generated mothur-format files (e.g., label/group header mismatches). |
| `convert_to_count_table_format.sh` | Wraps the conversion of final abundance tables into a standardized count-table format shared across pipelines. |

> **Note on mothur-format output:** mothur's native pipeline builds this table via `make.contigs → screen.seqs → unique.seqs → align.seqs → filter.seqs → pre.cluster → chimera.vsearch → classify.seqs → remove.lineage → cluster.split → make.shared → summary.single/dist.shared`. The scripts above replicate the *output format* of this workflow (the `.shared` count table) so that DADA2 and QIIME 2 results can be compared against mothur-native runs without re-implementing mothur's full clustering pipeline.

## Workflow Overview

1. **Primer trimming** (`trimmid_primers.sh`) and **read preparation** (`Demultiplex1.sh`)
2. **Pipeline-specific processing** — DADA2 (`DADA2/`), QIIME 2 (`Qiime2.sh`), and UPARSE/mothur-style OTU clustering (`script uparse.sh`)
3. **Batch execution across all three mock datasets** (`run_mocks.sh`)
4. **Cross-pipeline standardization** — converting DADA2 and QIIME 2 output into a common mothur-style count table (`generate_mothur.py`, `fix_mothur.py`, `convert_to_count_table_format.sh`) for direct, like-for-like comparison against mothur's native output

## Author

Alaa M. Elmansy — 2026 Research Summer Internship in Bioinformatics, *Benchmarking Analysis of 16S rRNA Amplicon Sequencing Pipelines*.
