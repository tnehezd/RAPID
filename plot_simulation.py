import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# --- Adatfájl útvonala ---
# Győződj meg róla, hogy ez az útvonal helyes a szkript futtatási helyéhez képest!
# Feltételezi, hogy a szkript a projekt gyökérkönyvtárában van.
file_path = 'output_0015/LOGS/dust_particle_evolution.dat'

# --- Adatok beolvasása ---
# A Gnuplot példához hasonlóan, ignoráljuk a komment sorokat (# jellel kezdődőket).
try:
    df = pd.read_csv(file_path, sep=r'\s+', comment='#', header=None,
                     names=['Time_Step', 'Particle_ID', 'Radius_AU'])
except FileNotFoundError:
    print(f"Hiba: A fájl nem található a következő útvonalon: {file_path}")
    print("Kérlek, ellenőrizd az útvonalat és győződj meg róla, hogy a szkript a projekt gyökérkönyvtárában van.")
    exit()
except Exception as e:
    print(f"Hiba az adatok beolvasása közben: {e}")
    exit()

# --- Y-tengely tartományának meghatározása a T=0 adatok alapján ---
# Szűrjük az adatokat a Time_Step == 0 értékre
df_t0 = df[df['Time_Step'] == 0]

# Ellenőrizzük, hogy van-e adat T=0-ban
if df_t0.empty:
    print("Figyelmeztetés: Nincs adat a T=0 időpillanatban. Az Y-tengely tartománya automatikus lesz.")
    y_min = df['Radius_AU'].min() * 0.95
    y_max = df['Radius_AU'].max() * 1.05
else:
    # Kiszámoljuk a sugár minimumát és maximumát T=0-ban
    y_min_t0 = df_t0['Radius_AU'].min()
    y_max_t0 = df_t0['Radius_AU'].max()

    # Kicsit kibővítjük a tartományt, hogy legyen "levegő" a pontok körül
    y_min = y_min_t0 * 0.95
    y_max = y_max_t0 * 1.05

# Ha az egész adatállományban van a T=0 tartományon kívüli adat, és azt is látni akarjuk,
# akkor az y_min és y_max értékeket módosíthatjuk úgy, hogy ne vágjuk le az ábrát.
# De a feladat szerint fix T=0 tartomány kell.

# --- Ábrázolás Matplotlib és Seaborn segítségével ---

# Ábra és tengelyek létrehozása
plt.figure(figsize=(12, 8)) # Ábra mérete (szélesség, magasság inch-ben)

# Seaborn scatterplot használata a színezéshez
# A hue='Particle_ID' paraméter automatikusan a Particle_ID szerint színez
# A palette='viridis' a színséma (másokat is választhatsz, pl. 'magma', 'plasma', 'jet', 'cividis')
# A s=20 a pontok mérete
sns.scatterplot(data=df, x='Time_Step', y='Radius_AU', hue='Particle_ID',
                palette='viridis', s=20, alpha=0.7,  # alpha a pontok átlátszósága
                legend='full') # Teljes színskála legenda megjelenítése

# Címek és feliratok
plt.title("Porrészecskék sugarának evolúciója (Részecske ID szerint színezve)", fontsize=16)
plt.xlabel("Idő [év]", fontsize=14)
plt.ylabel("Sugár [AU]", fontsize=14)

# Y-tengely tartományának beállítása
plt.ylim(y_min, y_max)

# Színskála címe
# Megtaláljuk a színskála objektumot és beállítjuk a címét
handles, labels = plt.gca().get_legend_handles_labels()
# A scatterplot által generált legenda az első, amit módosítani kell.
# Általában a 'hue' legenda a legutolsó a listában.
if plt.gca().get_legend() is not None:
    plt.gca().get_legend().set_title("Részecske ID")

# Rács hozzáadása a jobb olvashatóság érdekében
plt.grid(True, linestyle='--', alpha=0.6)

# Elrendezés optimalizálása (hogy ne lógjanak le a címkék)
plt.tight_layout()

# Ábra mentése fájlba
output_image_path = 'dust_radius_evolution_python.png' # PNG kimenet
# output_image_path = 'dust_radius_evolution_python.pdf' # Vagy PDF-be
plt.savefig(output_image_path, dpi=300) # dpi a felbontás

print(f"Az ábra elmentve ide: {output_image_path}")

# Ábra megjelenítése (opcionális, ha lokálisan futtatod)
plt.show()