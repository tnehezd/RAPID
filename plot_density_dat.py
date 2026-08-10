import numpy as np
import matplotlib.pyplot as plt

steps = [0, 500, 1000, 5000, 10000, 50000, 99000]

plt.figure(figsize=(10, 6))

for step in steps:

    filename = f"test_output_0004/LOGS/density_profile_{step:08d}.dat"

    # nézd meg először: col 1 = radius, col 2 = Sigma?
    r, sigma = np.loadtxt(filename, usecols=(1, 2), unpack=True)

    plt.plot(
        r,
        sigma,
        lw=2,
        label=f"step {step}"
    )

plt.xlabel("Radius [AU]")
plt.ylabel(r"$\Sigma_{\rm gas}$ [g/cm$^2$]")
plt.title("Viscous ring test - gas surface density")
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()