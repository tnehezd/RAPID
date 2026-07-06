import numpy as np
import matplotlib.pyplot as plt

# --- Timesteps to plot ---
steps = [0, 200, 400, 600, 800]

# --- Colors for each timestep ---
colors = ["#1f77b4", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b"]

plt.figure(figsize=(10, 6))

for step, col in zip(steps, colors):
    dust_file   = f"output_0052/LOGS/dust_density_profile_{step:08d}.dat"
    micron_file = f"output_0052/LOGS/dust_micron_density_profile_{step:08d}.dat"

    # Load data (col 2 = radius, col 3 = Sigma)
    r_dust, sigma_dust = np.loadtxt(dust_file, usecols=(1,2), unpack=True)
    r_mic, sigma_mic   = np.loadtxt(micron_file, usecols=(1,2), unpack=True)

    # Plot dust (solid)
    plt.plot(r_dust, sigma_dust, color=col, lw=2, label=f"dust {step}")

    # Plot micron (dashed)
    plt.plot(r_mic, sigma_mic, color=col, lw=2, ls="--", label=f"micron {step}")

plt.xlabel("Radius [AU]")
plt.ylabel("Surface density [g/cm²]")
plt.title("Dust vs Micron Dust – timesteps 0..800")
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
