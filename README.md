# SPPL2a Virtual Screening & Lead Analysis Pipeline

This repository contains the full molecular dynamics (MD), binding energetics, interaction mapping, and quantum chemical (DFT) analysis scripts for the virtual screening and lead characterization study targeting human **SPPL2a** (Signal Peptide Peptidase-Like 2a).

The analysis pipeline characterizes top hit compounds from **ZINC15**, **Enamine-REAL**, and **Aurigene** libraries in comparison with a reference control compound:
* **Control Compound**
* **Aurigene Hits:** `aurigene_top1`, `aurigene_top2`, `aurigene_top3`
* **Enamine Hits:** `enamine_top1`, `enamine_top2`, `enamine_top3`
* **ZINC Hits:** `zinc_top1`, `zinc_top2`, `zinc_top3`

---

## Repository Structure

The scripts are organized into modular, workflow-based directories:

### `01_stability_metrics/`
Calculations and multi-panel plotting for global structural stability and equilibration of protein–ligand complexes.
* `analyze_rmsd.sh`: GROMACS batch script computing protein backbone and ligand root-mean-square deviation (RMSD).
* `analyze_rmsf.sh`: Calculates per-residue root-mean-square fluctuation (RMSF) of SPPL2a.
* `analyze_sasa_rg.sh`: Computes Radius of Gyration ($R_g$) and Solvent Accessible Surface Area (SASA) time-series.
* `plot_rmsd.py`: Generates publication-ready comparative RMSD plots across systems.
* `plot_metrics.py`: Multi-panel plotting for $R_g$, SASA, and RMSF comparisons.

---

### `02_pca_fel/`
Essential dynamics and 2D Free Energy Landscape (FEL) contour mapping.
* `analyze_pca_combined.sh`: Fits all trajectories to a universal reference, constructs joint covariance matrices, and projects conformations onto principal components (PC1/PC2).
* `plot_pca_combined.py`: Constructs 2D Gibbs free energy landscapes ($\Delta G = -k_B T \ln P$) and maps conformational energy minima.

---

### `03_mmpbsa/`
End-point binding free energy and per-residue hotspot energy decomposition using MM-PBSA.
* `prepare_mmpbsa_trajectories.sh`: Extracts equilibrated trajectory snapshots and prepares complex, receptor, and ligand topologies for `gmx_MMPBSA`.
* `get_mmpbsa_table.py`: Parses and compiles binding free energy ($\Delta G_{\text{bind}}$) with individual energetic components (vdW, electrostatics, polar, and apolar solvation).
* `get_decomp_top.py`: Extracts per-residue energy decomposition to profile key catalytic dyad and pocket hotspot interactions.

---

### `04_interactions/`
Dynamic protein–ligand interaction profiling, hydrogen bond occupancies, and bipartite chord network visualization.
* `analyze_hbonds.py`: Time-resolved donor–acceptor hydrogen bond detection.
* `hbond_occupancy_analysis.py`: Calculates simulation-wide percentage occupancy and persistence of specific H-bonds.
* `extract_all_interactions.py`: Comprehensive extraction of non-covalent interactions (hydrogen bonds, hydrophobic contacts, $\pi$-stacking, and ionic interactions).
* `generate_bipartite_circos.py`: Generates bipartite Circos chord diagrams linking compound pharmacophore motifs to SPPL2a transmembrane helices.

---

### `05_dft/`
Quantum chemical calculations and electronic property profiling using Gaussian 16.
* `run_all.sh`: Batch execution script for DFT geometry optimization and vibrational frequency calculations.
* `cub_generate.sh`: Automated `cubegen` pipeline generating 3D volumetric surfaces (HOMO, LUMO, total electron density, and Molecular Electrostatic Potential / MEP).
* `extract_data.py`: Parses Gaussian output logs to extract HOMO/LUMO energies, band gap ($\Delta E$), dipole moments, and polarizabilities.

---

## Prerequisites & Dependencies

* **MD Simulation & Analysis:**
  * [GROMACS](http://www.gromacs.org/) ($\ge$ 2021)
  * [gmx_MMPBSA](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/) ($\ge$ v1.5)
* **Quantum Chemistry:**
  * [Gaussian 16](https://gaussian.com/g16/) / GaussView
* **Python Environment:**
  * Python $\ge$ 3.8
  * `numpy`, `scipy`, `pandas`, `matplotlib`, `seaborn`
  * `MDAnalysis`, `mdtraj`, `rdkit`
  * `pyCircos` (for bipartite chord network plotting)

---
## Data Availability

The reduced trajectories and structure files required to run these scripts are available on Zenodo: '10.5281/zenodo.23011361'
