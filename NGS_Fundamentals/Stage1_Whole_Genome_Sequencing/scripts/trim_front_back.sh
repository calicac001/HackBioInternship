#!/usr/bin/env bash
set -euo pipefail

SAMPLES="samples.txt"
RAWDIR="dataset"          # directory with raw fastq.gz files
OUTDIR="trimmed_front_back"   # output directory
FASTQC_OUT="qc_trimmed_front_back" # FastQC output

mkdir -p "$OUTDIR" "$FASTQC_OUT"

while IFS= read -r sample || [[ -n "$sample" ]]; do
  # Find paired-end reads (adjust pattern if needed)
  R1=( "$RAWDIR"/${sample}*1.fastq.gz "$RAWDIR"/${sample}*_R1.fastq.gz )
  R2=( "$RAWDIR"/${sample}*2.fastq.gz "$RAWDIR"/${sample}*_R2.fastq.gz )

  if [[ ! -e "${R1[0]}" || ! -e "${R2[0]}" ]]; then
    echo "[WARN] Missing mates for $sample; skipping"
    continue
  fi

  IN1="${R1[0]}"
  IN2="${R2[0]}"
  OUT1="$OUTDIR/${sample}_1.trim.fastq.gz"
  OUT2="$OUTDIR/${sample}_2.trim.fastq.gz"
  HTML="$OUTDIR/${sample}_fastp.html"
  JSON="$OUTDIR/${sample}_fastp.json"

  echo "[INFO] Trimming $sample..."
  fastp \
    --detect_adapter_for_pe \
    --trim_front1 10 --trim_front2 10 \
    --trim_tail1 5  --trim_tail2 5 \
    -i "$IN1" -I "$IN2" \
    -o "$OUT1" -O "$OUT2" \
    --html "$HTML"

done < "$SAMPLES"

# Run FastQC to check results
fastqc "$OUTDIR"/*_1.trim.fastq.gz "$OUTDIR"/*_2.trim.fastq.gz -o "$FASTQC_OUT"

echo "[DONE] Check FastQC outputs in $FASTQC_OUT; optionally run MultiQC to summarize."
