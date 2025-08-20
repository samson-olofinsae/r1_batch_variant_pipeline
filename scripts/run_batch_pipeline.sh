#!/bin/bash

# Wrapper script to batch-run split_snv_indel_variants.py on all sample pairs

set -e  # Exit on error

INPUT_DIR="inputs"
SCRIPT_DIR="scripts"
REF_GENOME="ref/hg19_chr8.fa"
LOG_DIR="results/logs"

mkdir -p "$LOG_DIR"

echo " Batch processing started: $(date)" | tee "$LOG_DIR/run_batch.log"

# Loop over all *_R1.fastq.gz files
for R1 in "$INPUT_DIR"/*_R1.fastq.gz; do
    sample_base=$(basename "$R1" _R1.fastq.gz)
    R2="$INPUT_DIR/${sample_base}_R2.fastq.gz"

    # Check if R2 exists
    if [[ -f "$R2" ]]; then
        echo " Processing $sample_base..." | tee -a "$LOG_DIR/run_batch.log"

        # Call the Python script with sample-specific arguments
        python3 "$SCRIPT_DIR/split_snv_indel_variants.py" \
            "$R1" "$R2" "$REF_GENOME" > "$LOG_DIR/${sample_base}.log" 2> "$LOG_DIR/${sample_base}.err"

        # Check exit status
        if [[ $? -eq 0 ]]; then
            echo " Completed $sample_base" | tee -a "$LOG_DIR/run_batch.log"
        else
            echo " Error in $sample_base — check ${sample_base}.err" | tee -a "$LOG_DIR/run_batch.log"
        fi
    else
        echo " Missing R2 pair for $sample_base, skipping." | tee -a "$LOG_DIR/run_batch.log"
    fi
done

echo " Batch run finished: $(date)" | tee -a "$LOG_DIR/run_batch.log"
