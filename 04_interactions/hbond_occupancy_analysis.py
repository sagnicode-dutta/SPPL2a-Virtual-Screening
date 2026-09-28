import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
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

def get_pretty_name(sys_name):
    if sys_name == "control":
        return "control"
    parts = sys_name.split("_")
    if len(parts) == 2:
        # e.g., aurigene_top1 -> top1aurigene
        return f"{parts[1]}{parts[0]}"
    return sys_name

base_dir = "."
results_dir = "analysis_results"

def analyze_system(sys_name):
    pretty_name = get_pretty_name(sys_name)
    print(f"Analyzing H-bonds for {sys_name} ({pretty_name})...")
    tpr = os.path.join(base_dir, sys_name, "step7_1.tpr")
    xtc = os.path.join(results_dir, sys_name, "no_pbc.xtc")
    
    if not (os.path.exists(tpr) and os.path.exists(xtc)):
        print(f"Error: Files missing for {sys_name}")
        return None

    try:
        u = mda.Universe(tpr, xtc)
    except Exception as e:
        print(f"Error loading universe for {sys_name}: {e}")
        return None
    
    lig_selection = ligands[sys_name]
    prot_selection = "protein"
    
    # Initialize H-bond analysis
    hb = HydrogenBondAnalysis(
        universe=u,
        donors_sel=f"({prot_selection}) or ({lig_selection})",
        acceptors_sel=f"({prot_selection}) or ({lig_selection})",
        d_a_cutoff=3.5,
        d_h_a_angle_cutoff=150,
        update_selections=True
    )
    
    # Run analysis with step=10 to speed up
    print(f"Running analysis for {len(u.trajectory)} frames (step=10)...")
    hb.run(step=10)
    
    n_analyzed_frames = len(u.trajectory[::10])
    
    # Filter only Protein-Ligand interactions
    res_df = pd.DataFrame(hb.results.hbonds, columns=['frame', 'donor_idx', 'h_idx', 'acc_idx', 'dist', 'angle'])
    atoms = u.atoms
    bond_counts = {}
    
    for _, row in res_df.iterrows():
        d_idx, a_idx = int(row['donor_idx']), int(row['acc_idx'])
        d_atom = atoms[d_idx]
        a_atom = atoms[a_idx]
        
        # Check if one is protein and the other is ligand
        is_lig = lambda a: a.resname in ["UNK", "FTO"]
        # In MDA, protein atoms usually have 'protein' in their segid or belong to a protein residue
        # Simple check: if it's not water/ion/ligand, it's protein
        is_prot = lambda a: a.resname not in ["UNK", "FTO", "SOL", "WAT", "HOH", "TIP3", "NA", "CL", "K"]
        
        if (is_prot(d_atom) and is_lig(a_atom)) or (is_lig(d_atom) and is_prot(a_atom)):
            bond_key = (d_idx, a_idx)
            bond_counts[bond_key] = bond_counts.get(bond_key, 0) + 1
            
    output = []
    for (d_idx, a_idx), count in bond_counts.items():
        d_atom = atoms[d_idx]
        a_atom = atoms[a_idx]
        occupancy = (count / n_analyzed_frames) * 100
        
        if d_atom.resname in ["UNK", "FTO"]:
            lig_atom = f"{d_atom.name}"
            prot_res = f"{a_atom.resname}{a_atom.resid}"
            prot_atom = f"{a_atom.name}"
            bond_type = "Ligand-Donor"
        else:
            prot_res = f"{d_atom.resname}{d_atom.resid}"
            prot_atom = f"{d_atom.name}"
            lig_atom = f"{a_atom.name}"
            bond_type = "Protein-Donor"
            
        output.append({
            "System": pretty_name,
            "Residue": prot_res,
            "Prot_Atom": prot_atom,
            "Lig_Atom": lig_atom,
            "Type": bond_type,
            "Occupancy_Pct": round(occupancy, 2)
        })
    
    if output:
        out_df = pd.DataFrame(output).sort_values("Occupancy_Pct", ascending=False)
        out_path = os.path.join(results_dir, sys_name, "hbonds_occupancy.csv")
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        out_df.to_csv(out_path, index=False)
        print(f"Saved {len(output)} bonds to {out_path}")
        return out_df
    else:
        print(f"No Protein-Ligand H-bonds found for {sys_name}")
        return None

all_results = []
for sys in systems:
    df = analyze_system(sys)
    if df is not None:
        all_results.append(df)

if all_results:
    summary_df = pd.concat(all_results)
    summary_df.to_csv(os.path.join(results_dir, "all_hbonds_occupancy.csv"), index=False)
    print(f"\nCombined results saved to {results_dir}/all_hbonds_occupancy.csv")
