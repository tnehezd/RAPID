import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import os
import time

LOGDIR = "output_0069/LOGS"

def extract_step(fname):
    m = re.search(r"(\d{8})", fname)
    return int(m.group(1)) if m else None

def load_pairs():
    dust_files   = sorted(glob.glob(os.path.join(LOGDIR, "dust_density_profile_*.dat")))
    micron_files = sorted(glob.glob(os.path.join(LOGDIR, "dust_micron_density_profile_*.dat")))

    pairs = {}
    for f in dust_files:
        step = extract_step(f)
        if step is not None:
            pairs.setdefault(step, {})["dust"] = f

    for f in micron_files:
        step = extract_step(f)
        if step is not None:
            pairs.setdefault(step, {})["micron"] = f

    return {s: p for s, p in pairs.items() if "dust" in p and "micron" in p}

plt.ion()
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6), sharey=False)

last_steps = set()

while True:
    pairs = load_pairs()
    steps = sorted(pairs.keys())

    if set(steps) != last_steps:
        ax1.clear()
        ax2.clear()

        colors = plt.cm.plasma(np.linspace(0, 1, len(steps)))

        for color, step in zip(colors, steps):
            dust_file   = pairs[step]["dust"]
            micron_file = pairs[step]["micron"]

            r_dust, sigma_dust = np.loadtxt(dust_file, usecols=(1,2), unpack=True)
            r_mic, sigma_mic   = np.loadtxt(micron_file, usecols=(1,2), unpack=True)

            # --- Panel 1: dust + micron ---
            ax1.plot(r_dust, sigma_dust, color=color, lw=2)
            ax1.plot(r_mic, sigma_mic, color=color, lw=2, ls="--")

            # --- Panel 2: Σ_total ---
            sigma_mic_interp = np.interp(r_dust, r_mic, sigma_mic)
            sigma_total = sigma_dust + sigma_mic_interp

            ax2.plot(r_dust, sigma_total, color=color, lw=2)

        # --- Formatting ---
        ax1.set_xlabel("Radius [AU]")
        ax1.set_ylabel("Surface density [g/cm²]")
        ax1.set_title("Dust + Micron (log scale)")
        ax1.set_yscale("log")
        ax1.set_ylim(1e-10,1e-6)
        ax1.grid(True)

        ax2.set_xlabel("Radius [AU]")
        ax2.set_ylabel("Σ_total [g/cm²]")
        ax2.set_title("Total Σ (log scale)")
        ax2.set_yscale("log")
        ax2.set_ylim(1e-10,1e-6)

        ax2.grid(True)


        fig.tight_layout()
        fig.canvas.draw()
        fig.canvas.flush_events()

        last_steps = set(steps)

    time.sleep(1.0)
