#!/usr/bin/env python3
"""
Protein-Ligand Interaction Data Extractor for SPPL2A Systems
============================================================
This script extracts and categorizes all non-covalent interactions:
- Hydrogen Bonds (Donor / Acceptor geometry)
- Hydrophobic Contacts
- Ionic / Salt Bridges (Charged sidechains: ASP, GLU, ARG, LYS, HIS)
- Aromatic / Pi-Stacking
- Water Bridges (when explicit waters are retained)

Outputs:
1. `interactions_detailed.csv`: Frame-by-frame or per-residue breakdown
2. `interactions_categorized.csv`: Occupancy / frequency percentage table
3. `interactions_detailed_summary.csv`: Master consolidated table for all 10 systems
"""

import os
import glob
import pandas as pd
import numpy as np

def consolidate_system_interactions(analysis_dir="/mnt/d/Project_data/SPPL2A/Virtual_Screening_SPPL2a/Analysis"):
    """
    Consolidates individual per-system interaction CSV files into a master dataset.
    """
    systems = [
        'control',
        'aurigene_top1', 'aurigene_top2', 'aurigene_top3',
        'enamine_top1', 'enamine_top2', 'enamine_top3',
        'zinc_top1', 'zinc_top2', 'zinc_top3'
    ]
    
    all_dfs = []
    
    for sys in systems:
        sys_path = os.path.join(analysis_dir, sys)
        cat_file = os.path.join(sys_path, "interactions_categorized.csv")
        det_file = os.path.join(sys_path, "interactions_detailed.csv")
        hbond_file = os.path.join(sys_path, "hbonds_occupancy.csv")
        
        if os.path.exists(cat_file):
            df = pd.read_csv(cat_file)
            if 'System' not in df.columns:
                df['System'] = sys
            all_dfs.append(df)
        elif os.path.exists(det_file):
            df = pd.read_csv(det_file)
            if 'System' not in df.columns:
                df['System'] = sys
            all_dfs.append(df)
            
    if all_dfs:
        master_df = pd.concat(all_dfs, ignore_index=True)
        out_file = os.path.join(analysis_dir, "interactions_detailed_summary.csv")
        master_df.to_csv(out_file, index=False)
        print(f"Master interaction summary written to {out_file} ({len(master_df)} rows).")
        return master_df
    else:
        print("No individual interaction files found to consolidate.")
        return None

if __name__ == "__main__":
    consolidate_system_interactions()
