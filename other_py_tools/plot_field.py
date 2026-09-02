import numpy as np
import matplotlib.pyplot as plt

# beolvassuk a tömeg mátrixot
mass_matrix = np.loadtxt("output_0031/config/mass_field.dat")

# beolvassuk a z mátrixot (első oszlop = r, többi = z pozíciók)
z_data = np.loadtxt("output_0031/config/grid_field.dat")
r = z_data[:, 0]                   # r pozíciók AU-ban
z_matrix = z_data[:, 1:]            # z értékek AU-ban vagy H-ben

# készítünk meshgrid-et plotoláshoz
R, Z = np.meshgrid(np.arange(mass_matrix.shape[1]), r)

# ábrázolás
plt.figure(figsize=(8,6))
plt.pcolormesh(Z, z_matrix, mass_matrix, shading='auto', cmap='viridis')
plt.colorbar(label='Dust mass [g]')
plt.xlabel('Z position [AU or H]')
plt.ylabel('Radial distance r [AU]')
plt.title('Dust mass distribution in r-z plane')
plt.tight_layout()
plt.show()
