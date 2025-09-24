#!/usr/bin/env bash
set -euo pipefail

ASSEMBLY_DIR="spades_assembly"       
OUT_DIR="blast_results"
mkdir -p "$OUT_DIR"

# iterate over assembly directories
for dir in "$ASSEMBLY_DIR"/*; do
    [[ -d "$dir" ]] || continue
    sample=$(basename "$dir")
    contigs="$dir/contigs.fasta"
    SAMPLE_DIR="$OUT_DIR/${sample}"
    mkdir -p "$SAMPLE_DIR"
    query="$SAMPLE_DIR/query_first1k.fa"
    blastout="$SAMPLE_DIR/blast.tsv"

    # skip if contigs missing
    if [[ ! -f "$contigs" ]]; then
        echo "[SKIP] $sample: contigs.fasta not found"
        continue
    fi

    # skip if blast already done
    if [[ -f "$blastout" ]]; then
        echo "[SKIP] $sample: BLAST output already exists"
        continue
    fi

    echo "[INFO] Preparing query for $sample"

    
    # Take the first contig (everything from first '>' to next '>') and limit to 1000 bp
    awk 'BEGIN{RS=">"; ORS=""}
     NR>1 {
        split($0, lines, "\n");
        header = lines[1];
        seq = "";
        for(i=2; i<=length(lines); i++) seq = seq lines[i];
        if(length(seq)>1000) seq = substr(seq,1,1000);
        print ">" header "\n" seq "\n";
        exit
     }' "$contigs" > "$query"

    # quick sanity check: ensure query was written and non-empty
    if [[ ! -s "$query" ]]; then
        echo "[WARN] $sample: extracted query is empty - skipping BLAST"
        rm -f "$query"
        continue
    fi

    echo "[INFO] Running remote BLAST for $sample (this may take a while)"

    blastn -query "$query" -db nt -remote \
        -outfmt "6 qseqid ssciname evalue bitscore stitle" \
        -max_target_seqs 1 -out "$blastout" 2> "$SAMPLE_DIR/blast.log" || {
            echo "[ERROR] BLAST failed for $sample (see $SAMPLE_DIR/blast.log)"
            continue
            }

    echo "[DONE] $sample -> BLAST result: $blastout"

    # sleep 5 seconds to avoid hitting NCBI request limits
    sleep 5
done

echo "[ALL DONE] BLASTs placed in $OUT_DIR"
