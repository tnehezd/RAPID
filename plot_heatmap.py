import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

filename = "all_particles_vertical_profiles.dat"

# Beolvasás
r_list = []
z_list = []
rho_list = []

with open(filename) as f:
    z_vals, rho_vals, r_val = [], [], None
    for line in f:
        line = line.strip()
        if line.startswith("# Particle"):
            if z_vals:
                z_list.append(np.array(z_vals))
                rho_list.append(np.array(rho_vals))
                r_list.append(r_val)
                z_vals, rho_vals = [], []
            # r érték kinyerése
            parts = line.split(",")
            r_part = [p for p in parts if "r =" in p][0]
            r_val = float(r_part.split("=")[1].split()[0])
        elif line.startswith("#") or not line:
            continue
        else:
            z, rho = map(float, line.split())
            z_vals.append(z)
            rho_vals.append(rho)
    if z_vals:
        z_list.append(np.array(z_vals))
        rho_list.append(np.array(rho_vals))
        r_list.append(r_val)

# Pontok létrehozása scatterhez
r_points = []
z_points = []
rho_points = []

for i in range(len(r_list)):
    r_array = np.full_like(z_list[i], r_list[i])
    r_points.extend(r_array)
    z_points.extend(z_list[i])
    rho_points.extend(rho_list[i])

r_points = np.array(r_points)
z_points = np.array(z_points)
rho_points = np.array(rho_points)

# Scatter plot log színskálával
plt.figure(figsize=(10,6))
plt.scatter(r_points, z_points, c=rho_points, norm=LogNorm(), cmap='viridis', s=5)
plt.colorbar(label=r'$\rho(z)$ [g/cm$^2$]')
plt.xlabel("r [AU]")
#plt.yscale("log")   # log skála a z tengelyre

plt.ylabel("z [AU]")
plt.title("Vertical dust density profiles (cone shape)")
plt.show()
