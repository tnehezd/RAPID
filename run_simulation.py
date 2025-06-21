import yaml
import subprocess
import argparse
import os
import datetime

def run_c_program(executable_path, params, arg_mapping, program_name="C Program"):
    """
    Futtat egy C programot a megadott paraméterekkel.
    """
    if not os.path.exists(executable_path):
        print(f"Hiba: A '{program_name}' futtatható fájl nem található itt: {executable_path}. Lefordítottad már?")
        return False, None

    cmd_args = [executable_path]
    for key, value in params.items():
        c_arg_name = arg_mapping.get(key)
        if c_arg_name:
            if c_arg_name == "-i":
                # Csak akkor adjuk át az -i kapcsolót, ha van input_file megadva
                if value is not None and str(value) != "":
                    cmd_args.extend([c_arg_name, str(value)])
            elif isinstance(value, bool): # Ha van még boolean paramétered, ami 0/1-et vár
                cmd_args.extend([c_arg_name, "1" if value else "0"])
            else:
                cmd_args.extend([c_arg_name, str(value)])

    print(f"\n--- Futtatom: {program_name} ---")
    print(f"Parancs: {' '.join(cmd_args)}")

    try:
        process = subprocess.Popen(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
        for line in process.stdout:
            print(line, end='')
        process.wait()

        if process.returncode != 0:
            print(f"{program_name} hibakóddal fejeződött be: {process.returncode}")
            return False, process.returncode
        else:
            print(f"{program_name} sikeresen befejeződött.")
            return True, 0
    except FileNotFoundError:
        print(f"Hiba: A '{executable_path}' parancs nem található. Ellenőrizd az elérési utat és a jogosultságokat.")
        return False, 1
    except Exception as e:
        print(f"Hiba történt a {program_name} futtatása közben: {e}")
        return False, 1


def main():
    parser = argparse.ArgumentParser(description="A szimulációs munkafolyamat futtatása YAML konfigurációval.")
    parser.add_argument("-c", "--config", default="config.yaml",
                        help="A YAML konfigurációs fájl elérési útja.")
    parser.add_argument("-exec", "--executable", default="./bin/simulation", # Feltételezi, hogy 'main' a lefordított program neve
                        help="A fő C program futtatható fájljának elérési útja.")
    args = parser.parse_args()

    config_file = args.config
    main_executable = args.executable

    # 1. Konfiguráció beolvasása
    try:
        with open(config_file, 'r') as f:
            full_config = yaml.safe_load(f) or {}
    except FileNotFoundError:
        print(f"Hiba: A konfigurációs fájl '{config_file}' nem található.")
        return
    except yaml.YAMLError as exc:
        print(f"Hiba a YAML fájl '{config_file}' feldolgozása során: {exc}")
        return

    # Alapértelmezett beállítások (a C program options_t struktúráját tükrözik)
    default_options = {
        "drift": 1.0, "growth": 1.0, "evol": 1.0, "twopop": 1.0,
        "ufrag": 1000.0, "ffrag": 0.37, "ngrid_val": 2000, "input_file": "", "tStep": 0.0,
        "totalTime": 1.0e6, "outputFrequency": 1000.0, "startTime": 0.0,

        # Init_tool által használt paraméterek, most már beolvasztva a fő options-be
        "rmin_val": 1.0, "rmax_val": 100.0, "sigma0_val": 1.0, "sigmap_exp_val": 1.0,
        "alpha_visc_val": 0.01, "star_val": 1.0, "hasp_val": 0.05, "flind_val": 0.5,
        "r_dze_i_val": 0.0, "r_dze_o_val": 0.0, "dr_dze_i_val": 0.0, "dr_dze_o_val": 0.0,
        "a_mod_val": 0.0,
        "md_val": 0.0, "eps_val": 0.0, "ratio_val": 0.0, "mic_val": 0.0, "onesize_val": 0.0,
    }

    all_params = default_options.copy()

    # Szimulációs opciók betöltése
    if "simulation_options" in full_config:
        # A YAML-ből érkező nevek és a C-s `options_t` tagnevek leképezése
        sim_yaml_to_c_mapping = {
            "drift": "drift", "growth": "growth", "evol": "evol", "twopop": "twopop",
            "ufrag": "ufrag", "ffrag": "ffrag", "ngrid": "ngrid_val", "input": "input_file", "tStep": "tStep",
            "totalTime": "totalTime", "outputFrequency": "outputFrequency", "startTime": "startTime",
        }
        for yaml_key, c_key in sim_yaml_to_c_mapping.items():
            if yaml_key in full_config["simulation_options"]:
                all_params[c_key] = full_config["simulation_options"][yaml_key]
    else:
        print("Figyelem: 'simulation_options' szekció hiányzik a YAML-ből. Alapértelmezetteket használunk.")

    # Inicializálási opciók betöltése
    # Ezek a paraméterek akkor is betöltésre kerülnek, ha nincs -i flag,
    # és a C program felhasználja őket az alapértelmezett profil generálásakor.
    if "initialization_options" in full_config:
        init_yaml_to_c_mapping = {
            "n": "ngrid_val", "ri": "rmin_val", "ro": "rmax_val", "sigma0": "sigma0_val",
            "index": "sigmap_exp_val", "rdzei": "r_dze_i_val", "rdzeo": "r_dze_o_val",
            "drdzei": "dr_dze_i_val", "drdzeo": "dr_dze_o_val",
            "alpha": "alpha_visc_val", "amod": "a_mod_val", "h": "hasp_val",
            "flind": "flind_val", "m0": "star_val", "md": "md_val",
            "eps": "eps_val", "ratio": "ratio_val", "mic": "mic_val", "onesize": "onesize_val",
        }
        for yaml_key, c_key in init_yaml_to_c_mapping.items():
            if yaml_key in full_config["initialization_options"]:
                all_params[c_key] = full_config["initialization_options"][yaml_key]
    else:
        print("Figyelem: 'initialization_options' szekció hiányzik a YAML-ből. Az alapértelmezett init paramétereket használjuk.")

    # A C program argumentum leképezése
    # Ezek a kulcsok a Python `all_params` szótárának kulcsai,
    # az értékek pedig a C program parancssori kapcsolói.
    c_arg_mapping = {
        "drift": "-drift", "growth": "-growth", "evol": "-evol", "twopop": "-twopop",
        "ufrag": "-ufrag", "ffrag": "-ffrag", "ngrid_val": "-n", "input_file": "-i", "tStep": "-tStep",
        "totalTime": "-tmax", "outputFrequency": "-outfreq", "startTime": "-curr",
        # NINCS "run_init_only": "-init" többé, teljesen eltávolítva!

        # Init_tool specifikus kapcsolók a C program parse_options függvényének megfelelően
        "rmin_val": "-ri",
        "rmax_val": "-ro",
        "sigma0_val": "-sigma0_init",
        "sigmap_exp_val": "-index_init",
        "r_dze_i_val": "-rdzei",
        "r_dze_o_val": "-rdzeo",
        "dr_dze_i_val": "-drdzei",
        "dr_dze_o_val": "-drdzeo",
        "alpha_visc_val": "-alpha_init",
        "a_mod_val": "-amod",
        "hasp_val": "-h_init",
        "flind_val": "-flind_init",
        "star_val": "-m0_init",
        "md_val": "-md_init",
        "eps_val": "-eps_init",
        "ratio_val": "-ratio_init",
        "mic_val": "-mic_init",
        "onesize_val": "-onesize_init",
    }
    
    # 2. A fő C program futtatása az összes beállított paraméterrel
    success, return_code = run_c_program(main_executable, all_params, c_arg_mapping, "Fő szimuláció program")
    
    if not success:
        print(f"A C program hibával fejeződött be (hibakód: {return_code}).")
        return

    print("\nMinden program sikeresen lefutott.")

if __name__ == "__main__":
    main()