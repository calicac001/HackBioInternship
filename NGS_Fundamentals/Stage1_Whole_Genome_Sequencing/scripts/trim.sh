#!/bin/bash

# Create output dir if it doesn’t exist
mkdir -p trimmed

SAMPLES="samples.txt"

while IFS= read -r sample; do
    # Define input files (with wildcards if needed)
    R1=$(ls dataset/${sample}*_1.fastq.gz 2>/dev/null)
    R2=$(ls dataset/${sample}*_2.fastq.gz 2>/dev/null)

    # Skip if no files found
    if [[ -z "$R1" || -z "$R2" ]]; then
        echo "Skipping $sample — input files not found"
        continue
    fi

    # Define outputs (no wildcards!)
    OUT1="trimmed/${sample}_1.fastq.gz"
    OUT2="trimmed/${sample}_2.fastq.gz"
    HTML="trimmed/${sample}_fastp.html"

    # Run fastp
    echo "Running fastp for $sample..."
    fastp \
        -i "$R1" \
        -I "$R2" \
        -o "$OUT1" \
        -O "$OUT2" \
        --html "$HTML"

done < "$SAMPLES"

