import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib

# High-Visibility Settings (Adjusted for A4 Width)
matplotlib.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['xtick.major.size'] = 6
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['legend.fontsize'] = 10

systems = ['control', 'aurigene_top1', 'aurigene_top2', 'aurigene_top3', 
           'enamine_top1', 'enamine_top2', 'enamine_top3', 
           'zinc_top1', 'zinc_top2', 'zinc_top3']

display_names = ['Control', 'Aurigene T1', 'Aurigene T2', 'Aurigene T3',
                 'Enamine T1', 'Enamine T2', 'Enamine T3',
                 'Zinc T1', 'Zinc T2', 'Zinc T3']

colors = {
    'control': '#000000',
    'aurigene_top1': '#068438', 'aurigene_top2': '#e90000', 'aurigene_top3': '#0000f9',
    'enamine_top1': '#8a45e2', 'enamine_top2': '#00d500', 'enamine_top3': '#ffa800',
    'zinc_top1': '#ff3cda', 'zinc_top2': '#9abeff', 'zinc_top3': '#00d1a2'
}

base_dir = "Analysis/pca"
output_dir = "plots"
os.makedirs(output_dir, exist_ok=True)

def read_proj(file_path):
    # gmx anaeig -2d outputs 2 columns: PC1 and PC2
    if not os.path.exists(file_path): return np.array([]), np.array([])
    try:
        data = np.loadtxt(file_path, comments=['#', '@'])
        if data.shape[1] >= 2:
            return data[:, 0], data[:, 1]
    except:
        pass
    return np.array([]), np.array([])

def smooth_2d(data, sigma=1.5):
    """Pure numpy 2D Gaussian smoothing."""
    size = int(2 * sigma) * 2 + 1
    x = np.arange(-size // 2 + 1, size // 2 + 1)
    y = np.arange(-size // 2 + 1, size // 2 + 1)
    x, y = np.meshgrid(x, y)
    kernel = np.exp(-(x**2 + y**2) / (2 * sigma**2))
    kernel /= kernel.sum()
    
    # Pad and convolve
    padded = np.pad(data, size // 2, mode='edge')
    output = np.zeros_like(data)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            output[i, j] = np.sum(padded[i:i+size, j:j+size] * kernel)
    return output

def plot_combined_FEL():
    kT = 0.596 # kcal/mol at 300K
    
    # New specific order requested by user
    custom_order = [
        'zinc_top1', 'zinc_top2', 'zinc_top3',
        'enamine_top1', 'enamine_top2', 'enamine_top3',
        'aurigene_top1', 'aurigene_top2', 'aurigene_top3',
        'control'
    ]
    
    names_dict = {
        'zinc_top1': 'Zinc T1', 'zinc_top2': 'Zinc T2', 'zinc_top3': 'Zinc T3',
        'enamine_top1': 'Enamine T1', 'enamine_top2': 'Enamine T2', 'enamine_top3': 'Enamine T3',
        'aurigene_top1': 'Aurigene T1', 'aurigene_top2': 'Aurigene T2', 'aurigene_top3': 'Aurigene T3',
        'control': 'Control'
    }

    # 4 rows x 3 columns grid
    fig, axes = plt.subplots(4, 3, figsize=(8.27, 11.5), facecolor='white')
    
    all_pc1 = []
    all_pc2 = []
    for sys in custom_order:
        proj_file = os.path.join(base_dir, sys, "shared_2dproj.xvg")
        pc1, pc2 = read_proj(proj_file)
        if len(pc1) > 0:
            # Convert nm to Angstrom
            all_pc1.extend(pc1 * 10)
            all_pc2.extend(pc2 * 10)
    
    if not all_pc1: return
    pc1_min, pc1_max = min(all_pc1), max(all_pc1)
    pc2_min, pc2_max = min(all_pc2), max(all_pc2)

    # Mapping custom order to 4x3 grid indices
    grid_indices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 10] 

    for i, sys in zip(grid_indices, custom_order):
        ax = axes.flat[i]
        name = names_dict[sys]
        proj_file = os.path.join(base_dir, sys, "shared_2dproj.xvg")
        pc1, pc2 = read_proj(proj_file)
        
        if len(pc1) == 0: continue
        
        # Convert nm to Angstrom
        pc1_a = pc1 * 10
        pc2_a = pc2 * 10
        
        bins = 60
        H, xedges, yedges = np.histogram2d(pc1_a, pc2_a, bins=bins, 
                                          range=[[pc1_min, pc1_max], [pc2_min, pc2_max]])
        
        H = smooth_2d(H, sigma=3.0)
        H = H / H.sum()
        G = -kT * np.log(H + 1e-10)
        G = G - np.nanmin(G)
        G_masked = np.ma.masked_where(G > 4.0, G)
        
        im = ax.contourf(0.5*(xedges[:-1]+xedges[1:]),
                         0.5*(yedges[:-1]+yedges[1:]),
                         G_masked.T, levels=20, cmap='jet')
        
        ax.contour(0.5*(xedges[:-1]+xedges[1:]),
                   0.5*(yedges[:-1]+yedges[1:]),
                   G_masked.T, levels=6, colors='black', linewidths=0.3, alpha=0.3)

        ax.set_title(name, fontweight='bold', fontsize=11)
        ax.set_aspect('equal')
        ax.set_xlim(pc1_min, pc1_max)
        ax.set_ylim(pc2_min, pc2_max)
        
        for spine in ax.spines.values():
            spine.set_linewidth(1.0)
            
        ax.tick_params(axis='both', which='major', labelsize=8, width=1.0, size=4)
        
        # Add X-labels to the bottom of each column
        # Col 0 bottom: 6 (Aurigene T1), Col 1 bottom: 10 (Control), Col 2 bottom: 8 (Aurigene T3)
        if i in [6, 8, 10]:
            ax.set_xlabel('PC1 (Å)', fontweight='bold', fontsize=9)
        if i % 3 == 0:
            ax.set_ylabel('PC2 (Å)', fontweight='bold', fontsize=9)

    # Hide the empty bottom-left (9) and bottom-right (11) slots
    axes.flat[9].axis('off')
    axes.flat[11].axis('off')

    # Reduced wspace to 0.3 for compactness while attempting to avoid label overlap
    fig.subplots_adjust(bottom=0.14, hspace=0.35, wspace=0.3)
    # Move colorbar up and adjust spacing
    cbar_ax = fig.add_axes([0.2, 0.10, 0.6, 0.012])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar.outline.set_linewidth(1.0)
    cbar.set_label('ΔG (kcal/mol)', fontweight='bold', fontsize=10)
    cbar.ax.tick_params(labelsize=9)
    
    plt.savefig(os.path.join(output_dir, 'shared_PCA_FEL_grid_Angstrom.pdf'), dpi=300, bbox_inches='tight')
    print("Generated shared_PCA_FEL_grid_Angstrom.pdf")
    plt.close()

def plot_combined_scatter():
    fig, ax = plt.subplots(figsize=(8.27, 8), facecolor='white')
    
    # Use the same systems/order as FEL
    ordered_systems = [
        'control', 'zinc_top1', 'zinc_top2', 'zinc_top3',
        'enamine_top1', 'enamine_top2', 'enamine_top3',
        'aurigene_top1', 'aurigene_top2', 'aurigene_top3'
    ]
    
    # Plot individual clouds (subsampled for clarity)
    for sys in ordered_systems:
        proj_file = os.path.join(base_dir, sys, "shared_2dproj.xvg")
        pc1, pc2 = read_proj(proj_file)
        if len(pc1) > 0:
            # Convert to Angstrom
            pc1_a = pc1 * 10
            pc2_a = pc2 * 10
            ax.scatter(pc1_a[::5], pc2_a[::5],
                       alpha=0.4, s=3, color=colors[sys], label='_nolegend_')

    # Plot centroids with large circles (using Peak Density/Mode)
    for sys in ordered_systems:
        name = [d for s, d in zip(systems, display_names) if s == sys][0]
        proj_file = os.path.join(base_dir, sys, "shared_2dproj.xvg")
        pc1, pc2 = read_proj(proj_file)
        if len(pc1) > 0:
            pc1_a = pc1 * 10
            pc2_a = pc2 * 10
            
            # Find the peak density (Mode) using a 2D histogram
            H, xedges, yedges = np.histogram2d(pc1_a, pc2_a, bins=50)
            H_smooth = smooth_2d(H, sigma=2.0) # Smooth to find the true basin center
            
            # Get index of the max density bin
            idx = np.unravel_index(np.argmax(H_smooth), H_smooth.shape)
            
            # Convert bin index back to coordinates (using midpoints)
            centroid_x = (xedges[idx[0]] + xedges[idx[0]+1]) / 2
            centroid_y = (yedges[idx[1]] + yedges[idx[1]+1]) / 2
            
            ax.scatter(centroid_x, centroid_y,
                       s=150, color=colors[sys], marker='o', 
                       edgecolor='none', label=name, zorder=10)

    ax.set_xlabel('Principal Component 1 (Å)', fontweight='bold', fontsize=12)
    ax.set_ylabel('Principal Component 2 (Å)', fontweight='bold', fontsize=12)
    ax.set_title('PCA Conformational Space', pad=20, fontweight='bold', fontsize=14)
    ax.legend(markerscale=0.8, fontsize=9, loc='upper right', frameon=True, edgecolor='black')
    ax.grid(True, linestyle='--', alpha=0.3)
    
    # Enforce square aspect for accuracy
    ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'shared_PCA_centroids_Angstrom.pdf'), dpi=300)
    print("Generated shared_PCA_centroids_Angstrom.pdf")
    plt.close()

if __name__ == "__main__":
    plot_combined_FEL()
    plot_combined_scatter()
