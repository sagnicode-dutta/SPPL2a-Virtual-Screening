import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbondanalysis import HydrogenBondAnalysis
import pandas as pd
import os

# Analysis configuration
systems = [
    "control",
    "aurigene_top1", "aurigene_top2", "aurigene_top3",
    "enamine_top1", "enamine_top2", "enamine_top3",
    "zinc_top1", "zinc_top2", "zinc_top3"
]

ligands = {
    "control": "resname FTO",
    "aurigene_top1": "resname UNK",
    "aurigene_top2": "resname UNK",
    "aurigene_top3": "resname UNK",
    "enamine_top1": "resname UNK",
    "enamine_top2": "resname UNK",
    "enamine_top3": "resname UNK",
    "zinc_top1": "resname UNK",
    "zinc_top2": "resname UNK",
    "zinc_top3": "resname UNK"
}

base_dir = "."
results_dir = "analysis_results"

def analyze_system(sys_name):
    print(f"Analyzing H-bonds for {sys_name}...")
    tpr = os.path.join(base_dir, sys_name, "step7_1.tpr")
    # Using the PBC-corrected trajectory generated in previous steps
    xtc = os.path.join(results_dir, sys_name, "no_pbc.xtc")
    
    if not (os.path.exists(tpr) and os.path.exists(xtc)):
        print(f"Error: Files missing for {sys_name}")
        return

    u = mda.Universe(tpr, xtc)
    
    lig_selection = ligands[sys_name]
    prot_selection = "protein"
    
    # Initialize H-bond analysis
    # Distance cutoff: 3.5 A, Angle cutoff: 150 degrees
    hb = HydrogenBondAnalysis(
        universe=u,
        donors_sel=None, # Automatically detected
        acceptors_sel=None, # Automatically detected
        d_a_cutoff=3.5,
        d_h_a_angle_cutoff=150,
        update_selections=True
    )
    
    # Run analysis for Protein -> Ligand and Ligand -> Protein
    # We define the selections for donors and acceptors specifically
    hb.run(
        donors_sel=f"({prot_selection}) or ({lig_selection})",
        acceptors_sel=f"({prot_selection}) or ({lig_selection})"
    )
    
    # Filter only Protein-Ligand interactions
    # hb.results.hbonds contains [frame, donor_index, hydrogen_index, acceptor_index, distance, angle]
    data = []
    n_frames = len(u.trajectory)
    
    # Convert results to a temporary dataframe for filtering
    res_df = pd.DataFrame(hb.results.hbonds, columns=['frame', 'donor_idx', 'h_idx', 'acc_idx', 'dist', 'angle'])
    
    # Get atom information
    atoms = u.atoms
    
    # Track bond counts: {(donor_idx, acc_idx): count}
    bond_counts = {}
    
    for _, row in res_df.iterrows():
        d_idx, a_idx = int(row['donor_idx']), int(row['acc_idx'])
        d_atom = atoms[d_idx]
        a_atom = atoms[a_idx]
        
        # Check if one is protein and the other is ligand
        is_prot_lig = (
            (d_atom.resname != "UNK" and d_atom.resname != "FTO" and (a_atom.resname == "UNK" or a_atom.resname == "FTO")) or
            ((d_atom.resname == "UNK" or d_atom.resname == "FTO") and a_atom.resname != "UNK" and a_atom.resname != "FTO")
        )
        
        if is_prot_lig:
            bond_key = (d_idx, a_idx)
            bond_counts[bond_key] = bond_counts.get(bond_key, 0) + 1
            
    # Calculate occupancy and format results
    output = []
    for (d_idx, a_idx), count in bond_counts.items():
        d_atom = atoms[d_idx]
        a_atom = atoms[a_idx]
        occupancy = (count / n_frames) * 100
        
        # Identify which is the protein residue and which is the ligand atom
        if d_atom.resname in ["UNK", "FTO"]:
            lig_atom = f"{d_atom.name}"
            prot_res = f"{d_atom.universe.atoms[a_idx].resname}{d_atom.universe.atoms[a_idx].resid}"
            prot_atom = f"{d_atom.universe.atoms[a_idx].name}"
            bond_type = "Ligand-Donor"
        else:
            prot_res = f"{d_atom.resname}{d_atom.resid}"
            prot_atom = f"{d_atom.name}"
            lig_atom = f"{d_atom.universe.atoms[a_idx].name}"
            bond_type = "Protein-Donor"
            
        output.append({
            "Residue": prot_res,
            "Prot_Atom": prot_atom,
            "Lig_Atom": lig_atom,
            "Type": bond_type,
            "Occupancy_Pct": round(occupancy, 2)
        })
    
    if output:
        out_df = pd.DataFrame(output).sort_values("Occupancy_Pct", ascending=False)
        out_path = os.path.join(results_dir, sys_name, "hbonds_occupancy.csv")
        out_df.to_csv(out_path, index=False)
        print(f"Saved {len(output)} bonds to {out_path}")
    else:
        print(f"No Protein-Ligand H-bonds found for {sys_name}")

for sys in systems:
    analyze_system(sys)
