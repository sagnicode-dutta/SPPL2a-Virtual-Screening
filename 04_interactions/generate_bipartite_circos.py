#!/usr/bin/env python3
"""
Master Bipartite Circos Plot Generator for SPPL2A Protein-Ligand Interactions
=============================================================================
Directly powered by the 1000-frame MD trajectory interaction extraction table:
`Analysis/interaction/master_interactions_all_systems_1000frames.csv`

Features:
- Layering & Z-Order: H-Bonds rendered prominently OVER Aromatic and Hydrophobic chords.
- Residue Orientation: Radially/vertically oriented along the outer circle to prevent any text overlapping.
- Residue Translation: 1-letter canonical amino acid code (e.g., D351, D412, Y350, I372).
- Catalytic Dyad: Exclusively highlights D351 and D412 with distinct red arc sectors.
- Color Scheme:
    * Hydrophobic: Crisp emerald green ribbons (#388e3c, base layer)
    * Aromatic: Amber gold ribbons (#f57c00, zorder=2)
    * H-Bond: Vivid cobalt blue ribbons (#0277bd, zorder=4, layered over Aromatic)
    * Ionic: Intense crimson red ribbons (#d32f2f, zorder=5)
- Noise Reduction: Filtered at >= 20% simulation occupancy.
- Canvas Dimensions: Formatted for vertical A4 page (width <= 8.27 in).
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path

# Strict catalytic residues (Canonical SPPL2A sequence numbering)
CATALYTIC_RESIDUES = {'D351', 'D412'}

# Refined color hierarchy and balanced line widths (H-Bond strictly prioritized over Aromatic)
STYLE_CONFIG = {
    'Hydrophobic': {
        'color': '#388e3c',      # Crisp Emerald Green
        'zorder': 1,
        'alpha_base': 0.35,
        'alpha_scale': 0.25,
        'lw_base': 0.45,
        'lw_scale': 0.65,
    },
    'Aromatic': {
        'color': '#f57c00',      # Amber Gold
        'zorder': 2,
        'alpha_base': 0.65,
        'alpha_scale': 0.20,
        'lw_base': 0.60,
        'lw_scale': 0.80,
    },
    'H-Bond': {
        'color': '#0277bd',      # Vivid Cobalt / Royal Blue (Layered over Aromatic)
        'zorder': 4,
        'alpha_base': 0.88,
        'alpha_scale': 0.12,
        'lw_base': 0.85,
        'lw_scale': 1.00,
    },
    'Ionic': {
        'color': '#d32f2f',      # Intense Crimson Red
        'zorder': 5,
        'alpha_base': 0.90,
        'alpha_scale': 0.10,
        'lw_base': 0.90,
        'lw_scale': 1.10,
    }
}

def parse_residue_number(res_str):
    """Extract numeric integer for sequential sorting."""
    match = re.search(r'\d+', str(res_str))
    return int(match.group()) if match else 9999

def draw_open_bezier_chord(ax, p1, p2, color, alpha=0.5, linewidth=1.0, zorder=1):
    """Draws an aesthetically routed Bézier ribbon that keeps the center hollow."""
    ctrl_scale = 0.44
    c1 = p1 * ctrl_scale
    c2 = p2 * ctrl_scale
    
    verts = [
        (p1[0], p1[1]),
        (c1[0], c1[1]),
        (c2[0], c2[1]),
        (p2[0], p2[1]),
    ]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor='none', edgecolor=color,
                              linewidth=linewidth, alpha=alpha, capstyle='round',
                              zorder=zorder)
    ax.add_patch(patch)

def generate_publication_circos(base_dir="/mnt/d/Project_data/SPPL2A/Virtual_Screening_SPPL2a",
                                 min_occupancy=20.0):
    input_csv = os.path.join(base_dir, "Analysis", "interaction", "master_interactions_all_systems_1000frames.csv")
    plots_dir = os.path.join(base_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)
    
    if not os.path.exists(input_csv):
        raise FileNotFoundError(f"Missing master interaction dataset: {input_csv}")
        
    master_df = pd.read_csv(input_csv)
    
    # Standardize column for residue 1-letter code
    if 'Residue_1Letter' in master_df.columns:
        master_df['Residue'] = master_df['Residue_1Letter']
        
    # Standard systems ordering
    systems_order = [
        'control',
        'aurigene_top1', 'aurigene_top2', 'aurigene_top3',
        'enamine_top1', 'enamine_top2', 'enamine_top3',
        'zinc_top1', 'zinc_top2', 'zinc_top3'
    ]
    
    # Filter by Occupancy (>= min_occupancy %)
    filtered_df = master_df[master_df['Frequency'] >= min_occupancy].copy()
    
    available_systems = [s for s in systems_order if s in filtered_df['System'].unique()]
    unique_residues = sorted(filtered_df['Residue'].unique(), key=parse_residue_number)
    
    n_sys = len(available_systems)
    n_res = len(unique_residues)
    
    print(f" Plotting Circos: {len(filtered_df)} links | {n_sys} systems | {n_res} residues (Filtered >= {min_occupancy}%)")
    
    # Canvas Setup - Strict Portrait A4 Dimensions
    fig, ax = plt.subplots(figsize=(8.2, 9.6), subplot_kw={'aspect': 'equal'})
    ax.set_xlim(-1.44, 1.44)
    ax.set_ylim(-1.38, 1.38)
    ax.axis('off')
    
    R_inner = 0.98
    R_outer = 1.05
    R_text = 1.075
    
    # 4. Map Systems onto Left Arc (98° to 262°)
    sys_angles_deg = np.linspace(98, 262, n_sys)
    sys_positions = {}
    for i, sys_name in enumerate(available_systems):
        angle_rad = np.deg2rad(sys_angles_deg[i])
        sys_positions[sys_name] = angle_rad
        
        # Sector wedge
        arc_width_deg = (164 / n_sys) * 0.72
        theta1 = sys_angles_deg[i] - arc_width_deg / 2
        theta2 = sys_angles_deg[i] + arc_width_deg / 2
        
        wedge = patches.Wedge((0, 0), R_outer, theta1, theta2, width=R_outer-R_inner,
                              facecolor='#1b3b5f', edgecolor='#ffffff', linewidth=1.2, zorder=2)
        ax.add_patch(wedge)
        
        # System Label (Horizontal on left side)
        tx = (R_text + 0.03) * np.cos(angle_rad)
        ty = (R_text + 0.03) * np.sin(angle_rad)
        display_name = sys_name.replace('_', ' ').title()
        ax.text(tx, ty, display_name, ha='right', va='center',
                fontsize=8.5, fontweight='bold', color='#0a192f')

    # 5. Map Residues onto Right Arc (82° down to -82°) - Radially/Vertically Oriented
    res_angles_deg = np.linspace(82, -82, n_res)
    res_positions = {}
    for j, res_name in enumerate(unique_residues):
        angle_rad = np.deg2rad(res_angles_deg[j])
        res_positions[res_name] = angle_rad
        
        is_catalytic = (res_name in CATALYTIC_RESIDUES)
        
        # Sector wedge
        arc_width_deg = (164 / n_res) * 0.82
        theta1 = res_angles_deg[j] - arc_width_deg / 2
        theta2 = res_angles_deg[j] + arc_width_deg / 2
        
        face_col = '#d32f2f' if is_catalytic else '#2e7d32'
        edge_col = '#7f0000' if is_catalytic else '#1b5e20'
        
        wedge = patches.Wedge((0, 0), R_outer + (0.025 if is_catalytic else 0),
                              theta1, theta2, width=R_outer-R_inner,
                              facecolor=face_col, edgecolor=edge_col,
                              linewidth=1.8 if is_catalytic else 0.5, zorder=2)
        ax.add_patch(wedge)
        
        # Radial/Vertical Text Placement
        tx = (R_text + (0.02 if is_catalytic else 0)) * np.cos(angle_rad)
        ty = (R_text + (0.02 if is_catalytic else 0)) * np.sin(angle_rad)
        
        rot_angle = res_angles_deg[j]
        
        ax.text(tx, ty, res_name, ha='left', va='center',
                rotation=rot_angle, rotation_mode='anchor',
                fontsize=7.8 if not is_catalytic else 9.2,
                fontweight='bold' if is_catalytic else 'normal',
                color='#b71c1c' if is_catalytic else '#212121')

    # 6. Draw Chords strictly in ascending Z-Order: Hydrophobic (1) -> Aromatic (2) -> H-Bond (4) -> Ionic (5)
    filtered_df['zorder'] = filtered_df['Type'].apply(lambda t: STYLE_CONFIG.get(t, {}).get('zorder', 2))
    sorted_chords = filtered_df.sort_values(by=['zorder', 'Frequency'], ascending=[True, True])
    
    for _, row in sorted_chords.iterrows():
        sys = row['System']
        res = row['Residue']
        itype = row['Type']
        freq = float(row['Frequency'])
        
        if sys not in sys_positions or res not in res_positions:
            continue
            
        th_sys = sys_positions[sys]
        th_res = res_positions[res]
        
        p_sys = np.array([R_inner * np.cos(th_sys), R_inner * np.sin(th_sys)])
        p_res = np.array([R_inner * np.cos(th_res), R_inner * np.sin(th_res)])
        
        st = STYLE_CONFIG.get(itype, {
            'color': '#757575', 'zorder': 2,
            'alpha_base': 0.4, 'alpha_scale': 0.2,
            'lw_base': 0.5, 'lw_scale': 0.5
        })
        
        lw = st['lw_base'] + (freq / 100.0) * st['lw_scale']
        alpha = st['alpha_base'] + (freq / 100.0) * st['alpha_scale']
        
        draw_open_bezier_chord(ax, p_sys, p_res,
                                color=st['color'],
                                alpha=min(1.0, alpha),
                                linewidth=lw,
                                zorder=st['zorder'])

    # 7. Legend & Title
    plt.title("SPPL2A Protein-Ligand Interaction Network (1,000 Frames MD)",
              fontsize=12, fontweight='bold', pad=10, color='#0a192f')
    
    legend_patches = []
    # Ordered display in legend: Hydrophobic, Aromatic, H-Bond, Ionic
    ordered_legend_types = ['Hydrophobic', 'Aromatic', 'H-Bond', 'Ionic']
    for t in ordered_legend_types:
        if t in filtered_df['Type'].unique() and t in STYLE_CONFIG:
            legend_patches.append(
                patches.Patch(facecolor=STYLE_CONFIG[t]['color'], edgecolor='none', label=t)
            )
    
    ax.legend(handles=legend_patches, loc='lower center', bbox_to_anchor=(0.5, -0.06),
              ncol=len(legend_patches), frameon=True, fontsize=9.0,
              title="Interaction Types",
              title_fontsize=10.0, edgecolor='#b0bec5')

    # 8. Save PNG and PDF to plots folder
    suffix = f"_cutoff{int(min_occupancy)}" if min_occupancy != 20.0 else ""
    png_path = os.path.join(plots_dir, f"SPPL2A_bipartite_circos{suffix}.png")
    pdf_path = os.path.join(plots_dir, f"SPPL2A_bipartite_circos{suffix}.pdf")
    
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, format='pdf', bbox_inches='tight')
    plt.close()
    
    print(f" Saved updated plots to:\n  - {png_path}\n  - {pdf_path}")

def update_excel_summary(base_dir="/mnt/d/Project_data/SPPL2A/Virtual_Screening_SPPL2a"):
    """Updates the publication Excel summary workbook."""
    input_csv = os.path.join(base_dir, "Analysis", "interaction", "master_interactions_all_systems_1000frames.csv")
    out_xlsx = os.path.join(base_dir, "plots", "circos-interaction-summary_manuscript.xlsx")
    
    df = pd.read_csv(input_csv)
    df_20 = df[df['Frequency'] >= 20.0]
    
    systems = [
        'control',
        'aurigene_top1', 'aurigene_top2', 'aurigene_top3',
        'enamine_top1', 'enamine_top2', 'enamine_top3',
        'zinc_top1', 'zinc_top2', 'zinc_top3'
    ]
    
    rows = []
    for s in systems:
        sdf = df_20[df_20['System'] == s]
        counts = sdf['Type'].value_counts().to_dict()
        
        # Check catalytic dyad contacts
        d351_types = sdf[sdf['Residue_1Letter'] == 'D351']['Type'].tolist()
        d412_types = sdf[sdf['Residue_1Letter'] == 'D412']['Type'].tolist()
        
        rows.append({
            'System': s,
            'Total Contacts (>=20%)': len(sdf),
            'Hydrogen Bonds': counts.get('H-Bond', 0),
            'Hydrophobic Contacts': counts.get('Hydrophobic', 0),
            'Aromatic Contacts': counts.get('Aromatic', 0),
            'Ionic / Salt Bridges': counts.get('Ionic', 0),
            'Asp351 (D351) Contact': ', '.join(d351_types) if d351_types else 'None',
            'Asp412 (D412) Contact': ', '.join(d412_types) if d412_types else 'None',
            'Interacting Residues List': ', '.join(sorted(sdf['Residue_1Letter'].unique(), key=parse_residue_number))
        })
        
    summary_df = pd.DataFrame(rows)
    with pd.ExcelWriter(out_xlsx, engine='openpyxl') as writer:
        summary_df.to_excel(writer, sheet_name='Interaction_Summary_20pct', index=False)
        df.to_excel(writer, sheet_name='All_Interactions_Raw', index=False)
        
    print(f" Manuscript Excel updated: {out_xlsx}")

def main():
    base_dir = "/mnt/d/Project_data/SPPL2A/Virtual_Screening_SPPL2a"
    generate_publication_circos(base_dir=base_dir, min_occupancy=20.0)
    generate_publication_circos(base_dir=base_dir, min_occupancy=30.0)
    update_excel_summary(base_dir=base_dir)

if __name__ == "__main__":
    main()
