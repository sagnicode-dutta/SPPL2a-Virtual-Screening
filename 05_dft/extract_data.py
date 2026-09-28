import re
import csv

def extract_properties(file_path):
    homo = None
    lumo = None
    dipole = None
    polar = None
    
    occ_pattern = re.compile(r"Alpha\s+occ\.\s+eigenvalues\s+--\s+(.*)")
    virt_pattern = re.compile(r"Alpha\s+virt\.\s+eigenvalues\s+--\s+(.*)")
    dipole_pattern = re.compile(r"Dipole moment \(field-independent basis, Debye\):")
    polar_pattern = re.compile(r"Isotropic polarizability for W=.*?\s+([\d\.]+)\s+Bohr\*\*3")
    exact_polar_pattern = re.compile(r"Exact polarizability:\s+([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\.\-]+)")

    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
            
            # Find HOMO/LUMO (last set)
            last_occ_line_idx = -1
            for i, line in enumerate(lines):
                if "Alpha  occ. eigenvalues" in line:
                    last_occ_line_idx = i
            
            if last_occ_line_idx != -1:
                all_occ = []
                i = last_occ_line_idx
                while i >= 0 and "Alpha  occ. eigenvalues" in lines[i]:
                    all_occ.append(lines[i])
                    i -= 1
                all_occ.reverse()
                homo = float(all_occ[-1].split("--")[1].split()[-1])
                
                for j in range(last_occ_line_idx + 1, len(lines)):
                    if "Alpha virt. eigenvalues" in lines[j]:
                        lumo = float(lines[j].split("--")[1].split()[0])
                        break

            # Find Dipole Moment (last set)
            for i in range(len(lines) - 1, 0, -1):
                if dipole_pattern.search(lines[i]):
                    # The next line contains the values
                    next_line = lines[i+1]
                    if "Tot=" in next_line:
                        dipole = float(next_line.split("Tot=")[1].strip())
                        break

            # Find Isotropic Polarizability (last set)
            for i in range(len(lines) - 1, 0, -1):
                m = polar_pattern.search(lines[i])
                if m:
                    polar = float(m.group(1))
                    break
            
            # If not found, try Exact Polarizability
            if polar is None:
                for i in range(len(lines) - 1, 0, -1):
                    m = exact_polar_pattern.search(lines[i])
                    if m:
                        xx = float(m.group(1))
                        yy = float(m.group(3))
                        zz = float(m.group(6))
                        polar = (xx + yy + zz) / 3.0
                        break
                    
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        
    return homo, lumo, dipole, polar

molecules = [
    ("Control", "lig_control.log"),
    ("Zinctop1", "lig_zinctop1.log"),
    ("Enaminetop2", "lig_enaminetop2.log")
]

results = []
conversion_factor = 27.2114

for name, file_path in molecules:
    homo_hartree, lumo_hartree, dipole, polar = extract_properties(file_path)
    if homo_hartree is not None and lumo_hartree is not None:
        homo_ev = homo_hartree * conversion_factor
        lumo_ev = lumo_hartree * conversion_factor
        
        I = -homo_ev
        A = -lumo_ev
        gap = lumo_ev - homo_ev
        chi = (I + A) / 2
        eta = (I - A) / 2
        S = 1 / (2 * eta) if eta != 0 else 0
        omega = (chi**2) / (2 * eta) if eta != 0 else 0
        
        results.append({
            "Molecule Name": name,
            "EHOMO (eV)": round(homo_ev, 4),
            "ELUMO (eV)": round(lumo_ev, 4),
            "Energy Gap ΔE (eV)": round(gap, 4),
            "Electronegativity χ (eV)": round(chi, 4),
            "Hardness η (eV)": round(eta, 4),
            "Global Softness S (eV−1)": round(S, 4),
            "Electrophilicity ω (eV)": round(omega, 4),
            "Dipole Moment (Debye)": round(dipole, 4) if dipole is not None else "N/A",
            "Polarizability (Bohr^3)": round(polar, 4) if polar is not None else "N/A"
        })


if results:
    keys = results[0].keys()
    with open('molecular_properties.csv', 'w', newline='') as output_file:
        dict_writer = csv.DictWriter(output_file, fieldnames=keys)
        dict_writer.writeheader()
        dict_writer.writerows(results)
    print("Results saved to molecular_properties.csv")
    for res in results:
        print(res)
else:
    print("No results found.")
