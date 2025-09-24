#!/bin/bash

SUMMARY="blast_summary_species.tsv"
FREQ_FILE="blast_species_frequency.tsv"

# Count total number of samples
TOTAL=$(awk 'NR>1{count++} END{print count}' "$SUMMARY")

echo -e "Species\tCount\tPercentage" > "$FREQ_FILE"

# Tally species counts using all fields from column 2 onward
awk -v total="$TOTAL" 'NR>1 {
	if ($2 == "No_hit") {
		species = $2
	} else {
		species = $2 " " $3
	}
    
    species_count[species]++
}
END {
    for (sp in species_count) {
        pct = (species_count[sp]/total)*100
        printf "%s\t%d\t%.2f%%\n", sp, species_count[sp], pct
    }
}' "$SUMMARY" >> "$FREQ_FILE"

echo "[DONE] Frequency summary written to $FREQ_FILE"
