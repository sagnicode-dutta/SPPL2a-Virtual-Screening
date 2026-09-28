#!/bin/bash

# Ensure the scratch directory environment variable is explicitly set for this session
export GAUSS_SCRDIR=/home/dell/g16/g16_scratch

# List of your input files to process sequentially
INPUT_FILES=("lig_control.com" "lig_enaminetop2.com" "lig_zinctop1.com")

echo "========================================="
echo " Starting Gaussian 16 Sequential Batch Queue"
echo "========================================="

for file in "${INPUT_FILES[@]}"; do
    # Extract the base name (e.g., "lig_control" from "lig_control.com")
    base_name="${file%.com}"
    
    echo "Current time: $(date)"
    echo "Running calculation for: ${file}..."
    
    # Run Gaussian natively and wait until it completes entirely before moving on
    g16 < "$file" > "${base_name}.log"
    
    echo "Finished ${file}. Log saved to ${base_name}.log"
    echo "-----------------------------------------"
done

echo "All calculations in the batch queue have completed successfully!"