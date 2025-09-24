#!/bin/bash

ASSEMBLY_DIR="spades_assembly"
OUT_DIR="abricate_results"
mkdir -p "$OUT_DIR"

DBS=("ncbi" "vfdb")   # ncbi = AMR genes, vfdb = virulence/toxins

for sample in $(ls "$ASSEMBLY_DIR"); do
    contigs="$ASSEMBLY_DIR/$sample/contigs.fasta"

    # Skip if contigs missing
    if [[ ! -f "$contigs" ]]; then
        echo "Skipping $sample: contigs.fasta not found"
        continue
    fi

    for db in "${DBS[@]}"; do
        outfile="$OUT_DIR/${sample}_${db}.tab"

        # Skip if abricate already run
        if [[ -f "$outfile" ]]; then
            echo "Skipping $sample ($db): already processed"
            continue
        fi

        # Run abricate
        abricate --db "$db" "$contigs" > "$outfile"
        echo "Finished abricate for $sample ($db)"
    done
done
