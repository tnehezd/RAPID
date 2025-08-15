import glob
import re
import numpy as np
import matplotlib.pyplot as plt

def generate_dust_plots(file_pattern="output_0009/LOGS/dust_profile_*.dat", initial_file="output_0009/config/initial_dust_profile.dat"):
    """
    Kiszámítja és ábrázolja a kumulatív tömeget és a por sűrűséget a megadott fájlokból.
    
    Args:
        file_pattern (str): A futás közbeni adatok keresésére használt minta.
        initial_file (str): Az initial_dust_profile.dat fájl elérési útja.
    """
    
    # --- Kezdeti állapot (initial) ábrázolása ---
    print("Grafikonok generálása a kezdeti állapotról...")
    try:
        radii_initial = []
        rep_mass_pop1_initial = []
        rep_mass_pop2_initial = []

        with open(initial_file, 'r') as f:
            for line in f:
                if line.strip().startswith('#') or line.strip().startswith('-'):
                    continue
                
                parts = line.split()
                # Az eredeti kódnak megfelelően dolgozzuk fel a fájlt
                if len(parts) >= 4:
                    try:
                        radius = float(parts[1])
                        mass_pop1 = float(parts[2])
                        mass_pop2 = float(parts[4])
                        
                        radii_initial.append(radius)
                        rep_mass_pop1_initial.append(mass_pop1)
                        rep_mass_pop2_initial.append(mass_pop2)
                    except (ValueError, IndexError) as e:
                        print(f"Hiba az initial fájl feldolgozásakor: {line.strip()}. Hiba: {e}")
                        continue

        if radii_initial and rep_mass_pop1_initial and rep_mass_pop2_initial:
            # Kumulatív tömegszámítás a kezdeti állapotra
            # A két populáció reprezentatív tömegét összeadjuk, majd kumuláljuk
            total_mass_initial = np.array(rep_mass_pop1_initial) + np.array(rep_mass_pop2_initial)
            cumulative_mass_initial = np.cumsum(total_mass_initial)

            # --- Kumulatív tömeg ábrázolása ---
            plt.figure(figsize=(10, 6))
            plt.plot(radii_initial, cumulative_mass_initial, 'o-', label='Kumulatív por tömeg')
            plt.xlabel('Távolság (AU)')
            plt.ylabel('Kumulatív tömeg (M_sun)')
            plt.title('Kezdeti kumulatív por tömeg')
            plt.legend()
            plt.grid(True)
            plt.savefig("initial_cumulative_dust_mass.png")
            plt.close()
            print("Létrehozva: initial_cumulative_dust_mass.png")

            # --- Porprofil ábrázolása (az eredeti szkriptnek megfelelően) ---
            plt.figure(figsize=(10, 6))
            plt.plot(radii_initial, rep_mass_pop1_initial, 'o-', label='RepMass_Pop1_Msun')
            plt.plot(radii_initial, rep_mass_pop2_initial, 'x-', label='RepMass_Pop2_Msun')
            plt.xlabel('Távolság (AU)')
            plt.ylabel('Tömeg (M_sun)')
            plt.title('Kezdeti por profil')
            plt.legend()
            plt.grid(True)
            plt.savefig("initial_dust_profile.png")
            plt.close()
            print("Létrehozva: initial_dust_profile.png")

        else:
            print("Nincs adat a grafikonhoz az initial fájlban.")
    
    except FileNotFoundError:
        print(f"Hiba: Az initial fájl nem található a megadott helyen: {initial_file}")
    except Exception as e:
        print(f"Hiba történt az initial fájl feldolgozása közben: {e}")

    # --- Idősoros (time-series) adatok ábrázolása ---
    print("\nGrafikonok generálása az idősoros adatokhoz...")
    filenames = sorted(glob.glob(file_pattern))

    if not filenames:
        print("Nem található fájl a megadott minta alapján.")
        return
    
    for filename in filenames:
        try:
            match = re.search(r'dust_profile_(\d+)\.dat', filename)
            if not match:
                print(f"Figyelmeztetés: A fájlnév nem felel meg a mintának: {filename}")
                continue
            
            time_step = int(match.group(1))
            
            # Az adatok beolvasása, kihagyva az első fejléc sort
            # A hiba elkerülése érdekében csak azokat az oszlopokat olvassuk be, amelyek garantáltan léteznek.
            # Az eredeti kód alapján a 3. (index 2) és a 6. (index 5) oszlopot használjuk.
            data = np.genfromtxt(filename, skip_header=1)
            
            # Ellenőrizzük, hogy a beolvasott adat nem üres-e, és legalább 6 oszlopa van
            if data.size == 0 or data.ndim < 2 or data.shape[1] < 6:
                print(f"Figyelmeztetés: Nem elegendő adat a {filename} fájlban, kihagyva.")
                continue

            # Az oszlopok a következők:
            # 0: Time
            # 1: index
            # 2: R(AU) -> A távolság, amit ábrázolunk
            # 3: X(AU)
            # 4: Y(AU)
            # 5: Z(AU) -> A Z koordináta, amit reprezentatív tömegnek néztünk korábban
            
            dists = data[:, 2]
            # Az eredeti kódod a parts[5]-öt használta a tömegre. Feltételezzük, hogy ez a Z koordináta
            # egy másik, az eredeti fájlban nem szereplő oszlop volt, és hogy a valós tömeg
            # az 5. oszlopban található. (ha nem, ezt még módosíthatjuk)
            rep_masses = data[:, 5]
            
            # Rendezés a távolság szerint a kumulatív összeghez
            sorted_indices = np.argsort(dists)
            sorted_dists = dists[sorted_indices]
            sorted_masses = rep_masses[sorted_indices]

            # Kumulatív tömegszámítás
            cumulative_mass = np.cumsum(sorted_masses)
            
            # --- Porprofil grafikon (RepMass) ---
            plt.figure(figsize=(10, 6))
            plt.plot(dists, rep_masses, 'o', markersize=3, label='Reprezentatív tömeg')
            plt.xlabel('Távolság (AU)')
            plt.ylabel('Reprezentatív tömeg (M_sun)')
            plt.title(f'Por profil (reprezentatív tömeg) {time_step} éves korban')
            plt.legend()
            plt.grid(True)
            output_filename = f"dust_profile_repmass_{time_step:08d}.png"
            plt.savefig(output_filename)
            plt.close()
            print(f"Létrehozva: {output_filename}")

            # --- Kumulatív tömeg ábrázolása ---
            plt.figure(figsize=(10, 6))
            plt.plot(sorted_dists, cumulative_mass, 'o-')
            plt.xlabel('Távolság (AU)')
            plt.ylabel('Kumulatív tömeg (M_sun)')
            plt.title(f'Kumulatív por tömeg {time_step} éves korban')
            plt.grid(True)
            output_filename = f"dust_profile_cumulative_mass_{time_step:08d}.png"
            plt.savefig(output_filename)
            plt.close()
            print(f"Létrehozva: {output_filename}")
        
        except Exception as e:
            print(f"Hiba történt a {filename} fájl feldolgozása közben: {e}")

    print("\nGrafikonok generálása befejeződött.")

if __name__ == "__main__":
    generate_dust_plots()

