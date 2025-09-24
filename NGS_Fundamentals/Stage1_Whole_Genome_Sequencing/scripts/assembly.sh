#!/bin/bash

# Make output dir
mkdir -p spades_assembly

SAMPLES="samples.txt"

while IFS= read -r sample; do
    spades.py \
        -1 "trimmed_front_back/${sample}_1.trim.fastq.gz" \
        -2 "trimmed_front_back/${sample}_2.trim.fastq.gz" \
        -o "spades_assembly/${sample}" \
        --phred-offset 33
done < "$SAMPLES"
