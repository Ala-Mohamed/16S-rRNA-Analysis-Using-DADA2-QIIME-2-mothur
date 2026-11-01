cd ~/Desktop/bio/lotus-mock16
mkdir -p 01-input-trimmed
cutadapt \
    -g GTGCCAGCMGCCGCGGTAA \
    -G GGACTACHVGGGTWTCTAAT \
    --discard-untrimmed \
    -o 01-input-trimmed/mock-forward-read.fastq.gz \
    -p 01-input-trimmed/mock-reverse-read.fastq.gz \
    01-input/mock-forward-read.fastq.gz \
    01-input/mock-reverse-read.fastq.gz \
    > cutadapt_report.txt
cat cutadapt_report.txt