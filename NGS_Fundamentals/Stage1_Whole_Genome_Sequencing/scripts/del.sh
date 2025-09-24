#!/bin/bash

SAMPLES="samples.txt"
TARGET_DIR="qc_reports"  # <-- specify your directory here

# Loop over all files in the target directory
for file in "$TARGET_DIR"/*; do
    # Skip if it's a directory
    [ -d "$file" ] && continue

    keep=false

    # Check each sample in samples.txt
    while IFS= read -r sample; do
        if [[ "$(basename "$file")" == *"$sample"* ]]; then
            keep=true
            break
        fi
    done < "$SAMPLES"

    # Delete file if it doesn't match any sample
    if [ "$keep" = false ]; then
        echo "Deleting $file"
        rm -v "$file"
    fi
done
