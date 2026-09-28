#!/bin/bash
set -e

BASE_DIR=$(pwd)
RESULTS_DIR="$BASE_DIR/Analysis"

systems=(
    "control" 
    "aurigene_top1" "aurigene_top2" "aurigene_top3"
    "enamine_top1" "enamine_top2" "enamine_top3"
    "zinc_top1" "zinc_top2" "zinc_top3"
)

for sys in "${systems[@]}"; do
    echo "Calculating C-alpha RMSF for $sys..."
    SYS_DIR="$BASE_DIR/$sys"
    OUT_DIR="$RESULTS_DIR/$sys"
    
    # Use C-alpha group (3) for calculation
    # We use the PBC-corrected trajectory and fitting on C-alpha
    echo "3" | gmx rmsf -f "$OUT_DIR/no_pbc.xtc" -s "$SYS_DIR/step7_1.tpr" -n "$OUT_DIR/index.ndx" \
        -o "$OUT_DIR/rmsf.xvg" -res -quiet
    
    echo "Finished $sys RMSF"
done
