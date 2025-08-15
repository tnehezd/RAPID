import glob
import re
import numpy as np
import matplotlib.pyplot as plt

def generate_dust_plots_with_density(initial_file="output_0005/config/initial_dust_profile.dat"):
    """
    Kiszámítja és ábrázolja a felületi sűrűséget a kezdeti adatokból.

    Args:
        initial_file (str): Az initial_dust_profile.dat fájl elérési útja.
    """
    
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
                        mass_pop2 = float(parts[3])
                        radii.append(radius)
                        rep_mass_pop1.append(mass_pop1)
                        rep_mass_pop2.append(mass_pop2)
                    except (ValueError, IndexError) as e:
                        print(f"Hiba az initial fájl feldolgozásakor: {line.strip()}. Hiba: {e}")
                        continue
        
        # A listákat numpy tömbbé alakítjuk a könnyebb számítások érdekében
        radii = np.array(radii)
        rep_mass_pop1 = np.array(rep_mass_pop1)
        rep_mass_pop2 = np.array(rep_mass_pop2)

        if len(radii) > 1:
            # Gyűrűk határainak kiszámítása
            # A gyűrű belső és külső sugara két szomszédos sugár átlaga
            r_inner = (radii[:-1] + radii[1:]) / 2
            r_outer = (radii[1:] + radii[:-1]) / 2 # Ugyanaz, mint r_inner, de a koncepció kedvéért
            
            # Hozzáadunk egy belső és egy külső határt a teljes tartományhoz
            r_in = np.concatenate(([radii[0] - (radii[1] - radii[0]) / 2], r_inner))
            r_out = np.concatenate((r_outer, [radii[-1] + (radii[-1] - radii[-2]) / 2]))
            
            # A gyűrűk területe (pi * (R_out^2 - R_in^2))
            # Az egységeket AU-ban hagyjuk, mivel a tömeg is relatív (Msun)
            areas = np.pi * (r_out**2 - r_in**2)

            # Felületi sűrűségek kiszámítása
            surface_density_pop1 = rep_mass_pop1 / areas
            surface_density_pop2 = rep_mass_pop2 / areas

            # Grafikon generálása a felületi sűrűséghez
            plt.figure(figsize=(10, 6))
            plt.plot(radii, surface_density_pop1, 'o-', label='RepMass_Pop1 felületi sűrűség')
            plt.plot(radii, surface_density_pop2, 'x-', label='RepMass_Pop2 felületi sűrűség')
            plt.xlabel('Távolság (AU)')
            plt.ylabel(r'Felületi sűrűség ($M_{sun}/AU^2$)')
            plt.title('Kezdeti porprofil felületi sűrűsége')
            plt.legend()
            plt.grid(True)
            plt.savefig("initial_dust_surface_density.png")
            plt.close()
            print("Létrehozva: initial_dust_surface_density.png")
        else:
            print("Nincs elegendő adat a felületi sűrűség kiszámításához.")

    except FileNotFoundError:
        print(f"Hiba: Az initial fájl nem található: {initial_file}")
    except Exception as e:
        print(f"Hiba történt az initial fájl feldolgozása közben: {e}")


if __name__ == "__main__":
    generate_dust_plots_with_density()
