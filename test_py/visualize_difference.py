import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import yaml

def load_config():
    yaml_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "ring_test.yaml"))
    if not os.path.exists(yaml_path):
        yaml_path = "ring_test.yaml"
    with open(yaml_path, "r") as f:
        config = yaml.safe_load(f)
    return config

if __name__ == "__main__":
    config = load_config()
    output_dir = config["output_parameters"]["output_directory_name"]
    r0 = config["test"]["ring_center_au"]
    nu = config["test"]["fixed_ring_viscosity"]
    
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    logs_dir = os.path.join(base_dir, output_dir, "LOGS")
    an_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "analytical_outputs"))
    
    # 3 paneles ábra létrehozása egymás mellett
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), dpi=120, sharex=True, sharey=True)
    ax_an, ax_sim, ax_comp = axes
    
    # --- 1. PANEL: Analitikus referenciák ---
    ax_an.set_title("1. Analytical References", fontsize=12)
    an_files = sorted(glob.glob(os.path.join(an_dir, "analytical_T_*.dat")))
    colors_an = plt.cm.viridis(np.linspace(0, 1, max(len(an_files), 1)))
    
    for idx, an_path in enumerate(an_files):
        data = np.loadtxt(an_path, comments="#")
        filename = os.path.basename(an_path)
        ax_an.plot(data[:, 0], data[:, 1], '-', color=colors_an[idx], lw=2, label=filename.replace(".dat", ""))
    
    ax_an.set_xlabel("Radius [AU]")
    ax_an.set_ylabel(r"Surface Density [$M_\odot / \text{AU}^2$]")
    ax_an.legend(fontsize=8, loc="upper right")
    ax_an.grid(True, linestyle="--", alpha=0.6)
    
    # --- 2. PANEL: C kód simulációk (Javítva: r = 2. oszlop [1], sigma = 3. oszlop [2]) ---
    ax_sim.set_title("2. C Simulation Snapshots", fontsize=12)
    sim_files = sorted(glob.glob(os.path.join(logs_dir, "density_profile_*.dat")))
    colors_sim = plt.cm.plasma(np.linspace(0, 1, max(len(sim_files), 1)))
    
    for idx, sim_path in enumerate(sim_files):
        data = np.loadtxt(sim_path, comments="#")
        filename = os.path.basename(sim_path)
        
        r_sim = data[:, 1]     # 2. oszlop
        sigma_sim = data[:, 2] # 3. oszlop
        
        show_label = (idx % max(1, len(sim_files)//5) == 0)
        ax_sim.plot(r_sim, sigma_sim, '--', color=colors_sim[idx], lw=1.5, 
                    label=filename if show_label else "")
                    
    ax_sim.set_xlabel("Radius [AU]")
    ax_sim.legend(fontsize=8, loc="upper right")
    ax_sim.grid(True, linestyle="--", alpha=0.6)
    
    # --- 3. PANEL: Összehasonlítás ---
    ax_comp.set_title("3. Direct Comparison", fontsize=12)
    
    for idx, an_path in enumerate(an_files):
        data = np.loadtxt(an_path, comments="#")
        ax_comp.plot(data[:, 0], data[:, 1], '-', color="black", alpha=0.4, lw=1.5)
        
    for idx, sim_path in enumerate(sim_files):
        data = np.loadtxt(sim_path, comments="#")
        r_sim = data[:, 1]
        sigma_sim = data[:, 2]
        ax_comp.plot(r_sim, sigma_sim, 'o--', color="red", alpha=0.6, markersize=3)
        
    ax_comp.plot([], [], '-', color="black", lw=2, label="Analytical")
    ax_comp.plot([], [], 'o--', color="red", label="C Simulation")
    
    ax_comp.set_xlabel("Radius [AU]")
    ax_comp.legend(fontsize=10, loc="upper right")
    ax_comp.grid(True, linestyle="--", alpha=0.6)
    
    for ax in axes:
        ax.set_xlim(config["disk_parameters"]["inner_radius_au"], config["disk_parameters"]["outer_radius_au"])
        ax.set_ylim(bottom=0)
        
    plt.tight_layout()
    plt.show()