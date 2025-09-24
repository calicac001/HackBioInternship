#!/bin/bash

OUT_DIR="blast_results"
SUMMARY="blast_summary_species.tsv"

echo -e "Sample\tSpecies\tE-value" > "$SUMMARY"

# Loop through each sample folder
for sample_dir in "$OUT_DIR"/*/; do
    sample=$(basename "$sample_dir")  # get folder name
    file="$sample_dir/blast.tsv"

    if [[ -s "$file" ]]; then
        # Extract the first hit (collapse isolates)
        evalue=$(awk 'NR==1{print $3}' "$file")
        species=$(awk 'NR==1{for(i=5;i<=NF;i++) printf "%s ", $i; print ""}' "$file" | awk '{print $1, $2}')
        
        # Default if empty
        [[ -z "$species" ]] && species="No_hit"
        [[ -z "$evalue" ]] && evalue="NA"

        echo -e "${sample}\t${species}\t${evalue}" >> "$SUMMARY"
    else
        echo -e "${sample}\tNo_hit\tNA" >> "$SUMMARY"
    fi
done

echo "[DONE] Summary written to $SUMMARY"
