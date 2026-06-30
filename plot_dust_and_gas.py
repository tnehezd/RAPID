import os
import glob
import time
import h5py
import matplotlib.pyplot as plt
import numpy as np

# --- BEÁLLÍTÁSOK ---
MSUN = 1.98847e33
AU = 1.495978707e13
MSUN_PER_AU2_TO_CGS = MSUN / (AU**2)
logs_dir = "output_0043/LOGS"
file_pattern = os.path.join(logs_dir, "snapshot_*.h5")

# --- REFERENCIA SNAPSHOT (t=0) ---
ref_file = os.path.join(logs_dir, "snapshot_00000000.h5")
with h5py.File(ref_file, "r") as f0:
    r0_gas = f0["/gas_grid/radial_grid"][:]
    sigma0_gas = f0["/gas_grid/surface_density"][:] * MSUN_PER_AU2_TO_CGS

    # Por kezdeti értéke = 0.01 * gáz kezdeti
    sigma0_dust = 0.01 * sigma0_gas

plt.ion()
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

while True:
    all_snapshots = sorted(glob.glob(file_pattern))
    if not all_snapshots:
        time.sleep(5)
        continue

    latest = all_snapshots[-1]   # <<< CSAK A LEGUTOLSÓ

    ax1.clear()
    ax2.clear()

    with h5py.File(latest, "r") as f:
        # --- GÁZ ---
        r_gas = f["/gas_grid/radial_grid"][:]
        sigma_gas = f["/gas_grid/surface_density"][:] * MSUN_PER_AU2_TO_CGS

        # BAL PANEL: aktuális gáz
        ax1.loglog(r_gas, sigma_gas, color="teal", linewidth=2, label="Gas")

        # JOBB PANEL: gáz0 / gáz(t)
        gas_ratio = sigma0_gas / sigma_gas
        ax2.semilogx(r_gas, gas_ratio, color="teal", linewidth=2, label="Gas")

        # --- POR ---
        if "/dust_grid" in f:
            r_dust = f["/dust_grid/radial_grid"][:]
            sigma_dust = f["/dust_grid/surface_density"][:] * MSUN_PER_AU2_TO_CGS

            # BAL PANEL: aktuális por
            ax1.loglog(r_dust, sigma_dust, color="crimson", linestyle="--", linewidth=2, label="Dust")

            # Interpoláljuk a por rácsot a gáz rácsára
            sigma_dust_interp = np.interp(r_gas, r_dust, sigma_dust)

            # JOBB PANEL: por(t) / (gáz0 * 0.01)
            dust_ratio = sigma_dust_interp / sigma0_dust
            ax2.semilogx(r_gas, dust_ratio, color="crimson", linestyle="--", linewidth=2, label="Dust")

    # --- BAL PANEL FORMÁZÁS ---
    ax1.set_xlabel("Radius [AU]")
    ax1.set_ylabel(r"Surface Density [g cm$^{-2}$]")
    ax1.set_title("Absolute Surface Density")
    ax1.grid(True, which="both", ls="--", alpha=0.4)
    ax1.legend()

    # --- JOBB PANEL FORMÁZÁS ---
    ax2.set_xlabel("Radius [AU]")
    ax2.set_ylabel(r"$\Sigma(0) / \Sigma(t)$  and  $\Sigma_{\rm dust}(t)/(0.01\Sigma_{\rm gas}(0))$")
    ax2.set_title("Evolution Relative to Initial Snapshot")
    ax2.set_yscale('log')
    ax2.grid(True, which="both", ls="--", alpha=0.4)
    ax2.legend()

    plt.tight_layout()
    plt.pause(5)
