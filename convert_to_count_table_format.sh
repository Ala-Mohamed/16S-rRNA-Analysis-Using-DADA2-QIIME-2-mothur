cat << 'EOF' > convert_tables.sh
#!/bin/bash

# Navigate to working directory
cd /home/alaa/Downloads || exit 1

# Loop through the target files
for file in Dada2_mock16.count_table Dada2_mock18.count_table Dada2_mock19.count_table; do
    # Verify file exists before processing
    if [ ! -f "$file" ]; then
        echo "Warning: File $file not found, skipping."
        continue
    fi

    # Extract mock identifier (e.g., mock16, mock18, mock19)
    sample=$(echo "$file" | grep -oP 'mock\d+')
    outfile="${sample}.dada2.count_table"

    # Write formatted mothur headers and process table rows
    {
        echo "#Compressed Format: groupIndex,abundance. For example 1,6 would mean the read has an abundance of 6 for group ${sample}."
        echo "#1,${sample}"
        printf 'Representative_Sequence\ttotal\t%s\n' "${sample}"
        awk 'BEGIN {FS=OFS="\t"} NR>1 && NF>=2 {gsub(/\r/,"",$2); printf "%s\t%d\t1,%d\n",$1,$2,$2}' "$file"
    } > "$outfile"

    echo "Successfully generated $outfile"
done
EOF
