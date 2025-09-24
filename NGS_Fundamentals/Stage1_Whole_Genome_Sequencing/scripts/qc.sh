#!/bin/bash

SAMPLES="samples.txt"

while IFS= read -r sample; do
    # Find matching FASTQ files for this sample
    files=(dataset/${sample}*_1.fastq.gz dataset/${sample}*_2.fastq.gz)

    # Check if files exist
    if [[ ! -e "${files[0]}" ]]; then
        echo "No files found for $sample, skipping"
        continue
    fi

    # Optionally check if FastQC has already been run
    html_files=(qc_reports/*${sample}*_fastqc.html)
    zip_files=(qc_reports/*${sample}*_fastqc.zip)
    if [[ -e "${html_files[0]}" && -e "${zip_files[0]}" 
&& -e "${html_files[1]}" && -e "${zip_files[1]}" ]]; then
        echo "Skipping $sample — already processed"
        continue
    fi

    # Run FastQC on both paired files
    echo "Running FastQC for $sample..."
    fastqc "${files[@]}" -o qc_reports/

done < "$SAMPLES"

