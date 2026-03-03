#!/usr/bin/env bash
set -euo pipefail

# Define input directory and output file
input_dir="/wynton/group/sirota/clairedubin/non_pregnancy_datasets/Norenhag2024_PRJEB72778/raw_reads"
output_csv="manifest_pre_qc.csv"

# Write header
echo "specimen,R1,R2" > "$output_csv"

# Loop through all *_1.fq.gz files
for r1_file in "$input_dir"/*_1.fq.gz; do
    # Get base name without _1.fq.gz
    base_name=$(basename "$r1_file" _1.fq.gz)
    r2_file="$input_dir/${base_name}_2.fq.gz"

    # Check if R2 file exists
    if [[ -f "$r2_file" ]]; then
        echo "$base_name,$r1_file,$r2_file" >> "$output_csv"
    else
        echo "Warning: Missing R2 file for specimen $base_name" >&2
    fi
done

