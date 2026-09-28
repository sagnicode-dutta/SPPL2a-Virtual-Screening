import os

def get_mmpbsa_breakdown(sys):
    path = f"analysis_results/mmpbsa/results_hpc/{sys}/FINAL_RESULTS_MMPBSA.dat"
    if not os.path.exists(path): return None
    
    results = {}
    with open(path, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if "GENERALIZED BORN" in line or "POISSON BOLTZMANN" in line:
                # Find DELTA section
                for j in range(i, len(lines)):
                    if "Delta:" in lines[j]:
                        # Read components
                        for k in range(j+1, len(lines)):
                            l = lines[k]
                            if "VDWAALS" in l: results['vdW'] = float(l.split()[1])
                            if "EEL" in l and "1-4" not in l: results['Elec'] = float(l.split()[1])
                            if "EPB" in l or "EGB" in l: results['Pol'] = float(l.split()[1])
                            if "ENPOLAR" in l or "ESURF" in l: results['NonPol'] = float(l.split()[1])
                            if "TOTAL" in l: 
                                results['Total'] = float(l.split()[1])
                                return results
    return None

systems = [
    "control", "aurigene_top1", "aurigene_top2", "aurigene_top3",
    "enamine_top1", "enamine_top2", "enamine_top3",
    "zinc_top1", "zinc_top2", "zinc_top3"
]

print("System,vdW,Electrostatic,Polar_Solv,NonPolar_Solv,TOTAL_DeltaH")
for sys in systems:
    res = get_mmpbsa_breakdown(sys)
    if res:
        print(f"{sys},{res['vdW']:.2f},{res['Elec']:.2f},{res['Pol']:.2f},{res['NonPol']:.2f},{res['Total']:.2f}")
