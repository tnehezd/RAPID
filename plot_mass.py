import glob
import re
import numpy as np
import matplotlib.pyplot as plt

def generate_dust_plots(file_pattern="output_0008/LOGS/dust_profile_*.dat", initial_file="output_0008/config/initial_dust_profile.dat"):
    """
    Kiszámítja a kumulatív tömeget és generál grafikonokat minden egyes fájlból.

    Args:
        file_pattern (str): A futás közbeni adatok keresésére használt minta.
        initial_file (str): Az initial_dust_profile.dat fájl elérési útja.
    """
    
    # --- Kezdeti állapot (initial) ábrázolása ---
    
    try:
        radii = []
        rep_mass_pop1 = []
        rep_mass_pop2 = []

        with open(initial_file, 'r') as f:
            for line in f:
                if line.strip().startswith('#') or line.strip().startswith('-'):
                    continue
                
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        radius = float(parts[1])
                        mass_pop1 = float(parts[2])
                        mass_pop2 = float(parts[4])
                        radii.append(radius)
                        rep_mass_pop1.append(mass_pop1)
                        rep_mass_pop2.append(mass_pop2)
                    except (ValueError, IndexError) as e:
                        print(f"Hiba az initial fájl feldolgozásakor: {line.strip()}. Hiba: {e}")
                        continue

        if radii and rep_mass_pop1 and rep_mass_pop2:
            plt.figure(figsize=(10, 6))
            plt.plot(radii, rep_mass_pop1, 'o-', label='RepMass_Pop1_Msun')
            plt.plot(radii, rep_mass_pop2, 'x-', label='RepMass_Pop2_Msun')
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
            distances = []
            rep_masses = []

            with open(filename, 'r') as f:
                for line in f:
                    if line.strip().startswith('#') or line.strip().startswith('-'):
                        continue
                    
                    parts = line.split()
                    if len(parts) >= 6:
                        try:
                            distance = float(parts[2])
                            rep_mass = float(parts[5])
                            distances.append(distance)
                            rep_masses.append(rep_mass)
                        except (ValueError, IndexError) as e:
                            print(f"Hiba a sor feldolgozásakor a {filename} fájlban: {line.strip()}. Hiba: {e}")
                            continue

            if distances and rep_masses:
                plt.figure(figsize=(10, 6))
                plt.plot(distances, rep_masses, 'o-')
                plt.xlabel('Távolság (AU)')
                plt.ylabel('RepMass (M_sun)')
                plt.title(f'Por profil {time_step} éves korban')
                plt.grid(True)
                
                output_filename = f"dust_profile_{time_step:08d}.png"
                plt.savefig(output_filename)
                plt.close()
                print(f"Létrehozva: {output_filename}")
        
        except Exception as e:
            print(f"Hiba történt a {filename} fájl feldolgozása közben: {e}")

    print("\nGrafikonok generálása befejeződött.")

if __name__ == "__main__":
    generate_dust_plots()
