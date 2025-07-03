import pandas as pd
import numpy as np

# --- Beállítások ---
input_file_path = 'output_0001/LOGS/dust_particle_evolution.dat'
output_file_path = 'output_0001/LOGS/dust_particle_evolution_prepared.dat' # Ideiglenes fájl a Gnuplotnak

# --- Adatok beolvasása ---
try:
    # A 'comment'='#': figyelmen kívül hagyja a '#' jellel kezdődő sorokat
    # 'header=None': nincs fejlécsor
    # 'names': oszlopnevek definiálása
    # 'sep=r'\s+'': a szóközök (vagy több szóköz) a szeparátorok
    df = pd.read_csv(input_file_path, sep=r'\s+', comment='#', header=None,
                     names=['Time_Step', 'Particle_ID', 'Radius_AU'])
    print(f"Sikeresen beolvasva: {input_file_path}")

except FileNotFoundError:
    print(f"HIBA: A bemeneti fájl nem található: {input_file_path}")
    print("Kérlek, ellenőrizd az útvonalat és a fájlnevet.")
    exit()
except pd.errors.EmptyDataError:
    print(f"HIBA: A bemeneti fájl üres: {input_file_path}")
    print("Nincs adat a feldolgozáshoz.")
    exit()
except Exception as e:
    print(f"HIBA az adatbeolvasás közben: {e}")
    exit()

# --- Adatok átrendezése Gnuplot-barát formátumba ---
# Cél: minden Particle_ID adatai egy blokkban, üres sorokkal elválasztva.
# Rendezés Particle_ID, majd Time_Step szerint.
df_sorted = df.sort_values(by=['Particle_ID', 'Time_Step'])

# Ideiglenes fájl írása
with open(output_file_path, 'w') as f:
    current_particle_id = None
    for index, row in df_sorted.iterrows():
        # Ha új Particle_ID-t látunk, és ez nem az első blokk, írjunk egy üres sort
        if current_particle_id is not None and row['Particle_ID'] != current_particle_id:
            f.write("\n\n") # Két új sor = egy üres sor a Gnuplot számára
        
        # Írjuk ki az adatokat a fájlba
        f.write(f"{row['Time_Step']} {row['Particle_ID']} {row['Radius_AU']}\n")
        
        # Frissítjük az aktuális Particle_ID-t
        current_particle_id = row['Particle_ID']

print(f"Adatok sikeresen előkészítve és elmentve ide: {output_file_path}")

# Meghatározzuk a maximális ID-t, ami hasznos a Gnuplot szkriptnek
max_particle_id = df['Particle_ID'].max()
print(f"Maximális részecske ID: {max_particle_id}")

# Megjegyezzük a maximális ID-t egy kis segédfájlba a Gnuplot számára (opcionális, de hasznos)
with open('max_particle_id.tmp', 'w') as f_id:
    f_id.write(str(int(max_particle_id)))