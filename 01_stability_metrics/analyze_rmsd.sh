#!/bin/bash
set -e

BASE_DIR="/mnt/d/Project_data/SPPL2A/Virtual_Screening_SPPL2a"
RESULTS_DIR="$BASE_DIR/analysis_results"

systems=(
    "control" 
    "aurigene_top1" "aurigene_top2" "aurigene_top3"
    "enamine_top1" "enamine_top2" "enamine_top3"
    "zinc_top1" "zinc_top2" "zinc_top3"
)

ligands=(
    "FTO" 
    "UNK" "UNK" "UNK"
    "UNK" "UNK" "UNK"
    "UNK" "UNK" "UNK"
)

for i in "${!systems[@]}"; do
    sys=${systems[$i]}
    lig=${ligands[$i]}
    
    echo "Processing $sys with PBC removal..."
    SYS_DIR="$BASE_DIR/$sys"
    OUT_DIR="$RESULTS_DIR/$sys"
    mkdir -p "$OUT_DIR"

    # 1. Create Index File
    if [ ! -f "$OUT_DIR/index.ndx" ]; then
        printf "q\n" | gmx make_ndx -f "$SYS_DIR/step7_1.tpr" -o "$OUT_DIR/index.ndx" -quiet
    fi

    # 2. PBC Removal: Center Protein and make molecules whole
    # We use -pbc mol -center and select Protein for centering
    if [ ! -f "$OUT_DIR/no_pbc.xtc" ]; then
        echo "1 0" | gmx trjconv -f "$SYS_DIR/step7_1.xtc" -s "$SYS_DIR/step7_1.tpr" -n "$OUT_DIR/index.ndx" \
            -pbc mol -center -o "$OUT_DIR/no_pbc.xtc" -quiet
    fi

    # 3. Find ligand group
    LIG_GROUP=$(grep "\[" "$OUT_DIR/index.ndx" | grep -n " $lig " | cut -d: -f1)
    if [ ! -z "$LIG_GROUP" ]; then
        LIG_GROUP=$((LIG_GROUP-1))
    fi

    # 4. Protein Backbone RMSD (Using PBC-corrected trajectory)
    echo "4 4" | gmx rms -f "$OUT_DIR/no_pbc.xtc" -s "$SYS_DIR/step7_1.tpr" -n "$OUT_DIR/index.ndx" -o "$OUT_DIR/rmsd_backbone.xvg" -tu ns -quiet
    
    # 5. Ligand RMSD
    if [ ! -z "$LIG_GROUP" ]; then
        echo "4 $LIG_GROUP" | gmx rms -f "$OUT_DIR/no_pbc.xtc" -s "$SYS_DIR/step7_1.tpr" -n "$OUT_DIR/index.ndx" -o "$OUT_DIR/rmsd_ligand.xvg" -tu ns -quiet
    fi
    
    echo "Finished $sys RMSD"
done
