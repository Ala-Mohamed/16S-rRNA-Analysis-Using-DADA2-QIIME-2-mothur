library(dada2); packageVersion("dada2")
library(Biostrings)
library(ggplot2)
theme_set(theme_bw())

# ---- 0. Setup Paths & Metadata ----
library(dada2)

path <- "C:/Users/alsaad/Downloads/data"
sample.name <- "mock-18"

fnF <- file.path(path, sample.name, "mock-18_R1.fastq.gz")
fnR <- file.path(path, sample.name, "mock-18_R2.fastq.gz")

# -------------------------------------------------------------------------
# Pass 1: Remove 5' Primers Only (21 bp FWD, 20 bp REV)
# -------------------------------------------------------------------------
noprimers_dir <- file.path(path, sample.name, "noprimers")
dir.create(noprimers_dir, showWarnings = FALSE)

nopF <- file.path(noprimers_dir, "mock18_F_noprimers.fastq.gz")
nopR <- file.path(noprimers_dir, "mock18_R_noprimers.fastq.gz")

# Only use trimLeft here—do NOT set truncLen yet
out_pass1 <- filterAndTrim(fnF, nopF, fnR, nopR,
                           trimLeft = c(21, 20),
                           maxN = 0, compress = TRUE, multithread = FALSE)
head(out_pass1)


# -------------------------------------------------------------------------
# Step 2: Plot Quality Profiles of the Primer-Free Reads
# -------------------------------------------------------------------------
# Position 1 on these plots now represents your actual biological sequence
plotQualityProfile(nopF)
plotQualityProfile(nopR)

# -------------------------------------------------------------------------
# Pass 2: Apply 3' Quality Truncation (truncLen) on Primer-Free Reads
# -------------------------------------------------------------------------
filtered_dir <- file.path(path, sample.name, "filtered")
dir.create(filtered_dir, showWarnings = FALSE)

filtF <- file.path(filtered_dir, "mock18_F_filt.fastq.gz")
filtR <- file.path(filtered_dir, "mock18_R_filt.fastq.gz")

# Apply truncLen chosen from the plots above (e.g., c(220, 140))
# Note: Do NOT set trimLeft again here, as primers were removed in Pass 1
out_pass2 <- filterAndTrim(nopF, filtF, nopR, filtR,
                           truncLen = c(220, 190), # Adjust based on your plot inspection
                           maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE,
                           compress = TRUE, multithread = FALSE)

head(out_pass2)

# ---- 3. Learn Error Rates & Plot ----
errF <- learnErrors(filtF, multithread = TRUE)
errR <- learnErrors(filtR, multithread = TRUE)

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

# ---- 4. Denoising ----
dadaF <- dada(filtF, err = errF, multithread = TRUE)
dadaR <- dada(filtR, err = errR, multithread = TRUE)

# ---- 5. Merge Pairs ----
merger <- mergePairs(dadaF, filtF, dadaR, filtR, verbose = TRUE)

# ---- 6. Sequence Table & Chimera Removal ----
seqtab <- makeSequenceTable(list(merger))
rownames(seqtab) <- sample.name

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)

# ---- 7. Pipeline Tracking ----
getN <- function(x) sum(getUniques(x))
track <- cbind(
  input = out_pass1[1, "reads.in"],
  filtered = out_pass2[1, "reads.out"],
  denoisedF = getN(dadaF),
  denoisedR = getN(dadaR),
  merged = getN(merger),
  nonchim = sum(seqtab.nochim)
)
rownames(track) <- sample.name
print(track)

# =========================================================================
# 8. EXPORT ASVs: FASTA & mothur Shared Format
# =========================================================================
output_dir <- file.path(path, sample.name, "results")
dir.create(output_dir, showWarnings = FALSE)

# Header naming for ASVs (ASV0001, ASV0002, ...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- sprintf("ASV%04d", seq_len(length(asv_seqs)))

# 8A. Save ASVs as FASTA file
asv_fasta <- DNAStringSet(asv_seqs)
names(asv_fasta) <- asv_headers
writeXStringSet(asv_fasta, file.path(output_dir, "ASVs.fasta"))

# 8B. Save as mothur .count_table format
mothur_count_table <- data.frame(
  Representative_Sequence = asv_headers,
  total = colSums(seqtab.nochim),
  t(seqtab.nochim),
  check.names = FALSE
)

# Write to file with tab separation
write.table(mothur_count_table, 
            file = file.path(output_dir, "mock18.count_table"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# 8B. Save Count Table (mothur .shared format)
mothur_shared <- data.frame(
  label = "0.03",
  Group = rownames(seqtab.nochim),
  numOtus = ncol(seqtab.nochim),
  seqtab.nochim,
  check.names = FALSE
)
colnames(mothur_shared)[4:ncol(mothur_shared)] <- asv_headers

write.table(mothur_shared, 
            file = file.path(output_dir, "mock19.shared"), 
            sep = "\t", quote = FALSE, row.names = FALSE)