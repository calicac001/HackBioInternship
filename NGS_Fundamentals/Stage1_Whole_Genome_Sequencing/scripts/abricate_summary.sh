#!/usr/bin/env bash

# Directory containing abricate results
RESULTS_DIR="abricate_results"

# Output files
AMR_OUT="amr_summary.tsv"
VIR_OUT="vir_summary.tsv"

# Temporary files
AMR_TMP="amr_gene_sample.tsv"
VIR_TMP="vir_gene_sample.tsv"

# --- AMR genes ---
echo "Summarizing AMR genes..."
for f in "$RESULTS_DIR"/*_ncbi.tab; do
    sample=$(basename "$f" _ncbi.tab)
    # get gene names (column 6), skip header
    awk 'NR>1 {print $6}' "$f" | sort -u | \
    awk -v s="$sample" '{print $1"\t"s}'
done | sort -u > "$AMR_TMP"

awk '{genes[$1]++; samples[$2]++}
     END{
         total=length(samples);
         for (g in genes)
             printf "%s\t%d\t%.2f%%\n", g, genes[g], (genes[g]/total*100)
     }' "$AMR_TMP" | sort -k2,2nr > "$AMR_OUT"


# --- Virulence genes ---
echo "Summarizing virulence genes..."
for f in "$RESULTS_DIR"/*_vfdb.tab; do
    sample=$(basename "$f" _vfdb.tab)
    awk 'NR>1 {print $6}' "$f" | sort -u | \
    awk -v s="$sample" '{print $1"\t"s}'
done | sort -u > "$VIR_TMP"

awk '{genes[$1]++; samples[$2]++}
     END{
         total=length(samples);
         for (g in genes)
             printf "%s\t%d\t%.2f%%\n", g, genes[g], (genes[g]/total*100)
     }' "$VIR_TMP" | sort -k2,2nr > "$VIR_OUT"


# Cleanup
rm "$AMR_TMP" "$VIR_TMP"

echo "Done! Results saved to $AMR_OUT and $VIR_OUT"
