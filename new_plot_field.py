import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os

# 1. Elérési utak
path = "output_0078/config"
output_dir = "plots"
os.makedirs(output_dir, exist_ok=True)

gas_file = f"{path}/initial_gas_density_2d.dat"
gas_grid_file = f"{path}/initial_gas_grid_2d.dat"
dust_mass_file = f"{path}/initial_dust2d_mass.dat"
dust_grid_file = f"{path}/initial_dust2d_grid.dat"

# Adatok beolvasása
gas_data = np.loadtxt(gas_file)
gas_grid = np.loadtxt(gas_grid_file)
dust_data = np.loadtxt(dust_mass_file)
dust_grid = np.loadtxt(dust_grid_file)

# Mátrixok szétválasztása
gas_matrix = gas_data[:, 1:]
R_gas = gas_grid[:, 0::2]
Z_gas = gas_grid[:, 1::2]

# Por felismerése (r oszlop kezelése)
dust_matrix = dust_data[:, 1:] if dust_data.shape[1] > (dust_grid.shape[1]//2) else dust_data
R_dust = dust_grid[:, 0::2]
Z_dust = dust_grid[:, 1::2]

# --- VIZUÁLIS ELLENŐRZÉS (A valóság felfedése) ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 12), sharex=True)

# Gáz plot gridvonalakkal
im1 = ax1.pcolormesh(R_gas, Z_gas, gas_matrix, 
                     shading='auto', 
                     norm=LogNorm(vmin=gas_matrix[gas_matrix>0].min(), vmax=gas_matrix.max()), 
                     cmap='viridis',
                     edgecolors='white', linewidths=0.05) # Itt látszanak a cellák!
fig.colorbar(im1, ax=ax1, label='Gas Density [g/cm³]')
ax1.set_title('Gas Density - Grid Structure Verification')
ax1.set_ylabel('z [AU]')

# Por plot
im2 = ax2.pcolormesh(R_dust, Z_dust, dust_matrix, 
                     shading='auto', 
                     norm=LogNorm(vmin=dust_matrix[dust_matrix>0].min(), vmax=dust_matrix.max()), 
                     cmap='magma')
fig.colorbar(im2, ax=ax2, label='Dust Mass [g]')
ax2.set_title('Dust Distribution (Should be thin, but resolved)')
ax2.set_xlabel('r [AU]')
ax2.set_ylabel('z [AU]')

# Azonos skála a porhoz is, hogy lásd a gázhoz képest hol van
ax2.set_ylim(ax1.get_ylim())

plt.tight_layout()
plt.savefig(f"{output_dir}/grid_verification_map.png", dpi=300)

# --- A "VALÓSÁG" DIAGNOSZTIKA (Nincs csalás) ---
r_idx = len(R_gas) // 2
dz = np.diff(Z_gas[r_idx, :])

fig2, ax_dz = plt.subplots(figsize=(8, 5))
ax_dz.plot(dz, 'r-o', markersize=3)
ax_dz.set_title(f"Cell thickness (dz) at r={R_gas[r_idx,0]:.2f} AU")
ax_dz.set_xlabel("Cell index (midplane is in the middle)")
ax_dz.set_ylabel("Cella vastagság [AU]")
ax_dz.grid(True, which='both', linestyle='--', alpha=0.5)

# Számadatok kiírása az ábrára a bizonyítékhoz
mid_idx = len(dz) // 2
text = (f"Cella méret a szélen: {dz[0]:.6f} AU\n"
        f"Cella méret középen: {dz[mid_idx]:.6f} AU\n"
        f"Sűrítési arány: {dz[0]/dz[mid_idx]:.1f}x")
ax_dz.text(0.05, 0.95, text, transform=ax_dz.transAxes, verticalalignment='top', 
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.savefig(f"{output_dir}/grid_cell_size_audit.png")
print(f"Ábrák mentve a '{output_dir}' mappába.")
plt.show()




import numpy as np
import matplotlib.pyplot as plt

xi = np.linspace(-1, 1, 100)
s = 3.0

# Különböző leképezések (z koordináták)
z_sinh = np.sinh(s * xi) / np.sinh(s)
z_asinh = np.arcsinh(s * xi) / np.arcsinh(s)
z_para = np.sign(xi) * xi**2

# Cellavastagságok (dz)
dz_sinh = np.diff(z_sinh)
dz_asinh = np.diff(z_asinh)
dz_para = np.diff(z_para)

plt.figure(figsize=(10, 6))
plt.plot(dz_sinh, label='Sinh (Rossz: középen ritka)', color='red')
plt.plot(dz_asinh, label='Asinh (Jó: középen sűrű völgy)', color='green', lw=3)
plt.plot(dz_para, label='Parabola (V-alakú)', color='blue', linestyle='--')

plt.title("Melyik rács sűrít középen? (dz eloszlás)")
plt.xlabel("Cella index")
plt.ylabel("Cella vastagság (dz)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()