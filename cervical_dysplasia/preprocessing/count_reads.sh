#!/usr/bin/env bash
#$ -S /bin/bash 
#$ -cwd 
#$ -l mem_free=500M     
#$ -l h_rt=00:30:00  
#$ -o FASTQ_read_count.$JOB_ID.out
#$ -e FASTQ_read_count.$JOB_ID.err

set -euo pipefail

FASTQ_list=$1
OUTDIR=$2

SAMPLE_ID=$(awk "NR==$SGE_TASK_ID" "$FASTQ_list" | awk -F',' '{print $1}')
FASTQ_1=$(awk "NR==$SGE_TASK_ID" "$FASTQ_list" | awk -F',' '{print $2}')
FASTQ_2=$(awk "NR==$SGE_TASK_ID" "$FASTQ_list" | awk -F',' '{print $3}')

OUTFILE="${OUTDIR}/${SAMPLE_ID}.csv"

READ_COUNT_1=$(zcat "$FASTQ_1" | wc -l)
READ_COUNT_2=$(zcat "$FASTQ_2" | wc -l)

READ_COUNT_1=$((READ_COUNT_1 / 4))
READ_COUNT_2=$((READ_COUNT_2 / 4))

if [[ "$READ_COUNT_1" -ne "$READ_COUNT_2" ]]; then
    echo "ERROR: Read counts for $SAMPLE_ID are not equal (R1: $READ_COUNT_1, R2: $READ_COUNT_2)" >&2
    exit 1
fi

READ_COUNT=$((READ_COUNT_1 + READ_COUNT_2))

echo "$SAMPLE_ID,$READ_COUNT" > "$OUTFILE"

[[ -n "${JOB_ID:-}" ]] && qstat -j "$JOB_ID"

