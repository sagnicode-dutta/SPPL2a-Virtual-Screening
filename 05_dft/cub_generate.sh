#!/bin/bash

# Ensure the script stops if any individual command fails
set -e

echo "========================================================="
echo "Starting Master Cubegen Pipeline for Structural Analysis"
echo "========================================================="

# --- 1. CONTROL LIGAND ---
echo "--> Processing: lig_control..."
cubegen 0 MO=181 lig_control.fchk control_homo.cub
cubegen 0 MO=182 lig_control.fchk control_lumo.cub
cubegen 0 Density=SCF lig_control.fchk control_density.cub
cubegen 0 Potential=SCF lig_control.fchk control_mep.cub

# --- 2. ENAMINE TOP 2 ---
echo "--> Processing: lig_enaminetop2..."
cubegen 0 MO=116 lig_enaminetop2.fchk enamine_homo.cub
cubegen 0 MO=117 lig_enaminetop2.fchk enamine_lumo.cub
cubegen 0 Density=SCF lig_enaminetop2.fchk enamine_density.cub
cubegen 0 Potential=SCF lig_enaminetop2.fchk enamine_mep.cub

# --- 3. ZINC TOP 1 ---
echo "--> Processing: lig_zinctop1..."
cubegen 0 MO=128 lig_zinctop1.fchk zinc_homo.cub
cubegen 0 MO=129 lig_zinctop1.fchk zinc_lumo.cub
cubegen 0 Density=SCF lig_zinctop1.fchk zinc_density.cub
cubegen 0 Potential=SCF lig_zinctop1.fchk zinc_mep.cub

echo "========================================================="
echo "SUCCESS! All 12 grid cubes generated cleanly."
echo "You can now safely transfer these to your laptop."
echo "========================================================="