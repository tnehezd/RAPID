import os
import glob
import time
import h5py
import matplotlib.pyplot as plt
import matplotlib as mpl

# --- CONSTANTS ---
MSUN = 1.98847e33
AU = 1.495978707e13
MSUN_PER_AU2_TO_CGS = MSUN / (AU**2)

logs_dir = "test_output/LOGS"
snapshot_step = 15000

file_pattern = os.path.join(logs_dir, "snapshot_*.h5")

plt.ion()  # interactive mode

fig, ax = plt.subplots(figsize=(9,6))

while True:

    all_snapshots = sorted(glob.glob(file_pattern))
    profiles = []

    for file_path in all_snapshots:
        with h5py.File(file_path, "r") as f:
            t = f["/frame"]["time"][0]

            if t == 0 or abs(t % snapshot_step) < 1e-3:
                gas = f["/gas_grid"]

                r = gas["radial_grid"][:]
                sigma = gas["surface_density"][:] * MSUN_PER_AU2_TO_CGS

                profiles.append({"time": t, "r": r, "sigma": sigma})

    if not profiles:
        time.sleep(2)
        continue

    ax.clear()

    times = [p["time"] for p in profiles]
    norm = mpl.colors.Normalize(min(times), max(times))
    cmap = plt.get_cmap("viridis")

    for p in profiles:
        ax.semilogy(
            p["r"],
            p["sigma"],
            color=cmap(norm(p["time"])),
            lw=2
        )


    ax.set_xlabel("Radius [AU]")
    ax.set_ylabel(r"Gas Surface Density [g cm$^{-2}$]")
    
    ax.set_title("Live Gas Surface Density Evolution")
    ax.grid(True, which="both", ls="--", alpha=0.6)

    plt.pause(10)  # frissítés 2 másodpercenként