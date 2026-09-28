import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib

# High-Visibility Settings
matplotlib.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.major.width'] = 3
plt.rcParams['ytick.major.width'] = 3
plt.rcParams['xtick.major.size'] = 10
plt.rcParams['ytick.major.size'] = 10
plt.rcParams['xtick.labelsize'] = 24
plt.rcParams['ytick.labelsize'] = 24
plt.rcParams['axes.labelsize'] = 28
plt.rcParams['axes.titlesize'] = 32
plt.rcParams['legend.fontsize'] = 20

# Database Configurations
configs = {
    'aurigene': {
        'systems': ['control', 'aurigene_top1', 'aurigene_top2', 'aurigene_top3'],
        'colors': {
            'control': '#000000',
            'aurigene_top1': '#068438', # Keeping previous Aurigene colors as they were approved
            'aurigene_top2': '#e90000',
            'aurigene_top3': '#0000f9'
        },
        'labels': {
            'control': 'Control',
            'aurigene_top1': 'Aurigene Top 1',
            'aurigene_top2': 'Aurigene Top 2',
            'aurigene_top3': 'Aurigene Top 3'
        }
    },
    'enamine': {
        'systems': ['control', 'enamine_top1', 'enamine_top2', 'enamine_top3'],
        'colors': {
            'control': '#000000',
            'enamine_top1': '#8a45e2',
            'enamine_top2': '#00d500',
            'enamine_top3': '#ffa800'
        },
        'labels': {
            'control': 'Control',
            'enamine_top1': 'Enamine Top 1',
            'enamine_top2': 'Enamine Top 2',
            'enamine_top3': 'Enamine Top 3'
        }
    },
    'zinc': {
        'systems': ['control', 'zinc_top1', 'zinc_top2', 'zinc_top3'],
        'colors': {
            'control': '#000000',
            'zinc_top1': '#ff3cda',
            'zinc_top2': '#9abeff',
            'zinc_top3': '#00d1a2'
        },
        'labels': {
            'control': 'Control',
            'zinc_top1': 'ZINC15 Top 1',
            'zinc_top2': 'ZINC15 Top 2',
            'zinc_top3': 'ZINC15 Top 3'
        }
    }
}

base_dir = "analysis_results"
output_dir = "plots"
os.makedirs(output_dir, exist_ok=True)

def moving_average(data, window=50):
    if len(data) < window:
        return data
    return np.convolve(data, np.ones(window)/window, mode='same')

def read_xvg(file_path):
    time, data = [], []
    if not os.path.exists(file_path): return np.array([]), np.array([])
    with open(file_path, 'r') as f:
        for line in f:
            if not line.startswith(('@', '#')):
                parts = line.split()
                if len(parts) >= 2:
                    time.append(float(parts[0]))
                    data.append(float(parts[1]))
    return np.array(time), np.array(data)

def plot_database(db_name):
    config = configs[db_name]
    metrics = [
        ('rmsd_backbone.xvg', f'rmsd_{db_name}_backbone.pdf', 'Protein Backbone RMSD'),
        ('rmsd_ligand.xvg', f'rmsd_{db_name}_ligand.pdf', 'Ligand RMSD')
    ]

    for metric_file, output_pdf, title in metrics:
        plt.figure(figsize=(16, 12))
        valid_plot = False
        
        for sys in config['systems']:
            path = os.path.join(base_dir, sys, metric_file)
            t, d = read_xvg(path)
            
            if len(d) > 0:
                color = config['colors'][sys]
                label = config['labels'][sys]
                
                # Calculate moving average
                d_smooth = moving_average(d, window=50)
                
                # Plot faded raw data
                plt.plot(t, d, color=color, alpha=0.3, linewidth=1.2, label='_nolegend_')
                
                # Plot bold smoothed data
                plt.plot(t, d_smooth, color=color, label=label, linewidth=3, alpha=1.0)
                valid_plot = True
        
        if valid_plot:
            plt.xlabel('Time (ns)', fontweight='bold')
            plt.ylabel('RMSD (nm)', fontweight='bold')
            plt.title(title, pad=20, fontweight='bold')
            plt.legend(frameon=True, edgecolor='black', loc='upper left', fancybox=False)
            plt.grid(True, linestyle='--', alpha=0.3, linewidth=2)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, output_pdf), format='pdf', dpi=300)
            print(f"Generated {output_pdf}")
        
        plt.close()

# Run for all databases
for db in ['aurigene', 'enamine', 'zinc']:
    plot_database(db)
