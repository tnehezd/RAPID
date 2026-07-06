import os
import glob
import time
import h5py
import matplotlib.pyplot as plt
import numpy as np

# --- CONFIGURATION ---
# Physical constants and conversion factors
MSUN = 1.98847e33
AU = 1.495978707e13
MSUN_PER_AU2_TO_CGS = MSUN / (AU**2)

# Path settings
logs_dir = "output/LOGS"  # Ensure this matches your actual output directory
file_pattern = os.path.join(logs_dir, "snapshot_*.h5")

# --- REFERENCE SNAPSHOT (t=0) ---
# Load initial conditions to calculate relative ratios
ref_file = os.path.join(logs_dir, "snapshot_00000000.h5")
with h5py.File(ref_file, "r") as f0:
    r0_gas = f0["/gas_grid/radial_grid"][:]
    sigma0_gas = f0["/gas_grid/surface_density"][:] * MSUN_PER_AU2_TO_CGS
    
    # Calculate reference dust densities based on dust-to-gas ratio (eps=0.01)
    # and population ratio (e.g., 0.85 for primary, 0.15 for micron)
    sigma0_dust = 0.01 * 0.85 * sigma0_gas 
    sigma0_mic = 0.01 * 0.15 * sigma0_gas  

# Initialize interactive plot
plt.ion()
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

while True:
    # Get all snapshot files, sort them, and select the latest one
    all_snapshots = sorted(glob.glob(file_pattern))
    if not all_snapshots:
        time.sleep(5)
        continue

    latest = all_snapshots[-1]

    # Clear axes for next update
    ax1.clear()
    ax2.clear()

    with h5py.File(latest, "r") as f:
        # --- GAS GRID ---
        # Radial grid is shared across all components
        r_gas = f["/gas_grid/radial_grid"][:]
        sigma_gas = f["/gas_grid/surface_density"][:] * MSUN_PER_AU2_TO_CGS
        ax1.loglog(r_gas, sigma_gas, color="teal", linewidth=2, label="Gas")
        ax2.semilogx(r_gas, sigma0_gas / sigma_gas, color="teal", linewidth=2, label="Gas (Ratio)")

        # --- PRIMARY DUST (CM) ---
        if "/dust_grid" in f:
            sigma_dust = f["/dust_grid/surface_density"][:] * MSUN_PER_AU2_TO_CGS
            ax1.loglog(r_gas, sigma_dust, color="crimson", linestyle="--", linewidth=2, label="Dust (cm)")
            
            # Use shared gas grid for ratio calculation
            ax2.semilogx(r_gas, sigma_dust / sigma0_dust, color="crimson", linestyle="--", linewidth=2, label="Dust (cm) / Init")

        # --- MICRON DUST ---
        if "/micron_grid" in f:
            sigma_mic = f["/micron_grid/surface_density"][:] * MSUN_PER_AU2_TO_CGS
            ax1.loglog(r_gas, sigma_mic, color="purple", linestyle=":", linewidth=2, label="Dust (micron)")
            
            # Use shared gas grid for ratio calculation
            ax2.semilogx(r_gas, sigma_mic / sigma0_mic, color="purple", linestyle=":", linewidth=2, label="Dust (micron) / Init")

    # --- FORMATTING AND LABELS ---
    ax1.set_xlabel("Radius [AU]")
    ax1.set_ylabel(r"Surface Density [g cm$^{-2}$]")
    ax1.set_title("Absolute Surface Density")
    ax1.grid(True, which="both", ls="--", alpha=0.4)
    ax1.legend()

    ax2.set_xlabel("Radius [AU]")
    ax2.set_ylabel("Relative Ratio")
    ax2.set_title("Evolution Relative to Initial Snapshot")
    ax2.set_yscale('log')
    ax2.set_ylim(1e-2, 20)
    ax2.grid(True, which="both", ls="--", alpha=0.4)
    ax2.legend()

    plt.tight_layout()
    plt.pause(5)