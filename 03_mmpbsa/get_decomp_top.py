import os
import pandas as pd

def parse_decomp(file_path):
    if not os.path.exists(file_path):
        return []
    data = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
        # Find where data starts
        start_idx = -1
        for i, line in enumerate(lines):
            if "TOTAL" in line and "Avg." in lines[i+1]:
                start_idx = i + 2
                break
        if start_idx == -1:
            return []
        
        for line in lines[start_idx:]:
            if not line.strip() or line.startswith('|') or line.startswith(' '):
                if data: break # End of data section
                continue
            cols = line.split(',')
            if len(cols) < 18: continue
            res = cols[0]
            total_avg = float(cols[16])
            data.append((res, total_avg))
    
    # Sort by total energy (most negative)
    data.sort(key=lambda x: x[1])
    return data[:10]

systems = [
    "control", "aurigene_top1", "aurigene_top2", "aurigene_top3",
    "enamine_top1", "enamine_top2", "enamine_top3",
    "zinc_top1", "zinc_top2", "zinc_top3"
]

for sys in systems:
    path = f"analysis_results/mmpbsa/results_hpc/{sys}/FINAL_DECOMP_MMPBSA.dat"
    top_res = parse_decomp(path)
    if top_res:
        print(f"\n--- Top 10 Residues for {sys} ---")
        for res, val in top_res:
            print(f"{res}: {val:.4f} kcal/mol")
