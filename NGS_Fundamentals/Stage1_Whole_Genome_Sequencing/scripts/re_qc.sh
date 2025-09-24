#!/bin/bash
mkdir -p qc_reports_trimmed

SAMPLES="samples.txt"

while IFS= read -r sample; do
    # Find matching FASTQ files for this sample
    files=(trimmed/${sample}*_1.fastq.gz trimmed/${sample}*_2.fastq.gz)

    # Check if files exist
    if [[ ! -e "${files[0]}" ]]; then
        echo "No files found for $sample, skipping"
        continue
    fi

    # Check if FastQC has already been run
    html_files=(qc_reports_trimmed/*${sample}*_fastqc.html)
    zip_files=(qc_reports_trimmed/*${sample}*_fastqc.zip)
    if [[ -e "${html_files[0]}" && -e "${zip_files[0]}" 
&& -e "${html_files[1]}" && -e "${zip_files[1]}" ]]; then
        echo "Skipping $sample — already processed"
        continue
    fi

    # Run FastQC on both paired files
    echo "Running FastQC for $sample..."
    fastqc "${files[@]}" -o qc_reports_trimmed/

done < "$SAMPLES"

# Run MultiQC
multiqc qc_reports_trimmed/
