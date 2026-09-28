#!/bin/bash
systems=("control" "aurigene_top1" "aurigene_top2" "aurigene_top3" "enamine_top1" "enamine_top2" "enamine_top3" "zinc_top1" "zinc_top2" "zinc_top3")
REF="Analysis/pca/universal_ref.pdb"
BACKBONE_TPR="Analysis/pca/combined_ref.tpr"
mkdir -p Analysis/pca/combined

echo "=== Step 1: Fitting all trajectories to Universal Reference ==="
for sys in "${systems[@]}"; do
    mkdir -p "Analysis/pca/$sys"
    if [ ! -f "Analysis/pca/$sys/uni_fitted.xtc" ]; then
        echo "PBC correction for $sys..."
        echo "0 1" | gmx trjconv -s "$sys/step7_1.tpr" -f "$sys/step7_1.xtc" -o "Analysis/pca/$sys/pbc_corrected.xtc" -pbc nojump -center
        echo "Fitting $sys to universal reference..."
        # Use Backbone-only TPR as structure template (-s)
        echo "0 0" | gmx trjconv -s "$BACKBONE_TPR" -f "Analysis/pca/$sys/pbc_corrected.xtc" -o "Analysis/pca/$sys/uni_fitted.xtc" -fit rot+trans
        rm "Analysis/pca/$sys/pbc_corrected.xtc"
    fi
done

echo "=== Step 2: Concatenating trajectories ==="
if [ ! -f "Analysis/pca/combined/mega_traj.xtc" ]; then
    traj_list=""
    for sys in "${systems[@]}"; do
        traj_list="$traj_list Analysis/pca/$sys/uni_fitted.xtc"
    done
    gmx trjcat -f $traj_list -o Analysis/pca/combined/mega_traj.xtc -cat
fi

echo "=== Step 3: Calculating Shared Covariance Matrix ==="
if [ ! -f "Analysis/pca/combined/shared_eigenvec.trr" ]; then
    echo "Calculating shared covariance..."
    # Select 0 0 (System is Backbone)
    echo "0 0" | gmx covar -s "$BACKBONE_TPR" -f "Analysis/pca/combined/mega_traj.xtc" \
        -o Analysis/pca/combined/shared_eigenval.xvg \
        -v Analysis/pca/combined/shared_eigenvec.trr \
        -av Analysis/pca/combined/shared_average.pdb \
        -l Analysis/pca/combined/shared_covar.log
fi

echo "=== Step 4: Individual Projections onto Shared Axes ==="
for sys in "${systems[@]}"; do
    if [ ! -f "Analysis/pca/$sys/shared_2dproj.xvg" ]; then
        echo "Projecting $sys..."
        echo "0 0" | gmx anaeig -s "$BACKBONE_TPR" \
            -f "Analysis/pca/$sys/uni_fitted.xtc" \
            -v Analysis/pca/combined/shared_eigenvec.trr \
            -2d "Analysis/pca/$sys/shared_2dproj.xvg" -first 1 -last 2
    fi
done
