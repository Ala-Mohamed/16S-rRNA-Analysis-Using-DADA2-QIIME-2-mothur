install -c bioconda -c conda-forge usearsh -y

# Check that USEARCH is available
which usearch
# Link USEARCH to LotuS3
../LotuS3/lotus3 -link_usearch "$(which usearch)"

../LotuS3/lotus3 \
    -i 01-input-trimmed \
    -m mapping.txt \
    -o 03-otus-uparse \
    -p miSeq \
    -amplicon_type SSU \
    -CL uparse \
    -id 0.97 \
    -refDB SLV \
    -taxAligner vsearch \
    -useVsearch 1 \
    -buildPhylo 0 \
    -t 4