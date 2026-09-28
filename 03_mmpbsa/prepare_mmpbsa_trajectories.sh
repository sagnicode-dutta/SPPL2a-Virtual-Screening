#!/bin/bash
# Script to prepare systems for MMPBSA analysis
# Extracts last 200ns (300-500ns) and creates reduced TPR/TOP

BASE_DIR="/mnt/d/Project_data/SPPL2A/Virtual_Screening_SPPL2a"
OUT_BASE="$BASE_DIR/Analysis/mmpbsa"
MINI_MDP="$OUT_BASE/mini.mdp"

# System names
systems=("control" "aurigene_top1" "aurigene_top2" "aurigene_top3" "enamine_top1" "enamine_top2" "enamine_top3" "zinc_top1" "zinc_top2" "zinc_top3")
# Corresponding ligand residues
ligands=("FTO" "UNK" "UNK" "UNK" "UNK" "UNK" "UNK" "UNK" "UNK" "UNK")

mkdir -p "$OUT_BASE"

for i in "${!systems[@]}"; do
    sys=${systems[$i]}
    lig=${ligands[$i]}
    
    echo "================================================================"
    echo "Processing $sys (Ligand: $lig)..."
    echo "================================================================"
    
    SYS_DIR="$BASE_DIR/$sys"
    OUT_DIR="$OUT_BASE/$sys"
    mkdir -p "$OUT_DIR"
    cp -r "$SYS_DIR/toppar" "$OUT_DIR/"

    # 1. Create index file
    # We create a group for the ligand and then the complex (Protein + Ligand)
    # Using 'q' at the end to save and exit
    printf "r $lig\n1 | r $lig\nname 19 Complex\nq\n" | gmx make_ndx -f "$SYS_DIR/step7_1.tpr" -o "$OUT_DIR/mmpbsa.ndx" -quiet

    # Identify group numbers by searching for names in the generated index file
    # We need the group index (0-based)
    PROT_GRP=$(grep "\[" "$OUT_DIR/mmpbsa.ndx" | grep -n "\[ Protein \]" | head -n 1 | cut -d: -f1)
    PROT_GRP=$((PROT_GRP-1))
    
    LIG_GRP=$(grep "\[" "$OUT_DIR/mmpbsa.ndx" | grep -n "\[ $lig \]" | head -n 1 | cut -d: -f1)
    LIG_GRP=$((LIG_GRP-1))
    
    COMP_GRP=$(grep "\[" "$OUT_DIR/mmpbsa.ndx" | grep -n "\[ Complex \]" | head -n 1 | cut -d: -f1)
    COMP_GRP=$((COMP_GRP-1))

    echo "Groups: Protein=$PROT_GRP, Ligand=$LIG_GRP, Complex=$COMP_GRP"

    # 2. Extract last 200ns trajectory for the Complex group ONLY
    echo "0" | gmx trjconv -f "$SYS_DIR/step7_1.xtc" -s "$SYS_DIR/step7_1.tpr" -n "$OUT_DIR/mmpbsa.ndx" \
        -o "$OUT_DIR/tmp_nojump.xtc" -b 300000 -e 500000 -dt 1000 -pbc nojump -quiet

    # Second pass: Center on Protein, output Complex
    echo "$PROT_GRP $COMP_GRP" | gmx trjconv -f "$OUT_DIR/tmp_nojump.xtc" -s "$SYS_DIR/step7_1.tpr" -n "$OUT_DIR/mmpbsa.ndx" \
        -o "$OUT_DIR/${sys}_mmpbsa_200ns.xtc" -center -pbc mol -ur compact -quiet

    rm "$OUT_DIR/tmp_nojump.xtc"

    # 3. Create reduced complex.gro
    echo "$COMP_GRP" | gmx trjconv -f "$SYS_DIR/step7_1.gro" -s "$SYS_DIR/step7_1.tpr" -n "$OUT_DIR/mmpbsa.ndx" \
        -o "$OUT_DIR/complex.gro" -quiet

    # 4. Create reduced complex.top
    # We take the original topol.top and keep only PROA and the ligand in [ molecules ]
    sed '/\[ molecules \]/,$d' "$SYS_DIR/topol.top" > "$OUT_DIR/complex.top"
    echo "[ molecules ]" >> "$OUT_DIR/complex.top"
    echo "; Compound      #mols" >> "$OUT_DIR/complex.top"
    echo "PROA               1" >> "$OUT_DIR/complex.top"
    echo "$lig                1" >> "$OUT_DIR/complex.top"

    # 5. Generate reduced complex.tpr
    gmx grompp -f "$MINI_MDP" -c "$OUT_DIR/complex.gro" -p "$OUT_DIR/complex.top" -o "$OUT_DIR/complex.tpr" -maxwarn 5 -quiet

    # 6. Generate a fresh index file for the reduced complex
    # This ensures the group indices match the reduced system
    echo "q" | gmx make_ndx -f "$OUT_DIR/complex.tpr" -o "$OUT_DIR/complex.ndx" -quiet

    echo "Successfully prepared $sys"
done
