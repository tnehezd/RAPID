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
    for py_key, value in params.items():
        c_arg_name = arg_mapping.get(py_key)
        if c_arg_name:
            # Handle boolean values: convert to "1.0" or "0.0" for C
            if isinstance(value, bool):
                cmd_args.extend([c_arg_name, "1.0" if value else "0.0"])
            # Handle input_file_path specifically: only pass if not empty
            elif c_arg_name == "-i":
                if value is not None and str(value).strip() != "":
                    cmd_args.extend([c_arg_name, str(value)])
            # Handle output_directory_name: ensure it's always passed
            elif c_arg_name == "-o":
                if value is not None and str(value).strip() != "":
                    cmd_args.extend([c_arg_name, str(value)])
                else: # Fallback to a default if somehow empty (shouldn't happen with updated defaults)
                    cmd_args.extend([c_arg_name, "output"])
            # Handle other numeric or string parameters
            else:
                cmd_args.extend([c_arg_name, str(value)])

    print(f"\n--- Futtatom: {program_name} ---")
    print(f"Parancs: {' '.join(cmd_args)}")

    try:
        process = subprocess.Popen(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
        for line in process.stdout:
            print(line, end='') # Print output in real-time
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
    parser.add_argument("-exec", "--executable", default="./bin/simulation",
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
    # Ezeknek **pontosan** meg kell egyezniük a C kód `create_default_options` függvényében beállított alapértékekkel.
    default_options = {
        # Simulation Control Options (Boolean values are handled in run_c_program)
        "drift": 1.0,
        "growth": 1.0,
        "evol": 1.0,
        "twopop": 1.0,
        "ufrag": 1000.0,
        "ffrag": 0.37,

        # Core Disk Parameters (also serve as init_tool defaults if no input file)
        "ngrid_val": 2000,
        "rmin_val": 1.0,
        "rmax_val": 100.0,
        "sigma0_val": 1.0,
        "sigmap_exp_val": 1.0,
        "alpha_visc_val": 0.01,
        "star_val": 1.0,
        "hasp_val": 0.05,
        "flind_val": 0.5,
        "r_dze_i_val": 0.0,
        "r_dze_o_val": 0.0,
        "dr_dze_i_val": 0.0,
        "dr_dze_o_val": 0.0,
        "a_mod_val": 0.0,

        # File Input/Output
        "input_file": "",
        "output_dir_name": "output",

        # Time Parameters
        "tStep": 0.0,
        "totalTime": 1.0e6,
        "outputFrequency": 1000.0,
        "startTime": 0.0,

        # Init tool specific parameters' defaults
        "eps_val": 0.01,
        "ratio_val": 0.85,
        "mic_val": 1e-4,
        "onesize_val": 0.0,
    }

    all_params = default_options.copy()

    # Paraméterek betöltése a YAML-ből
    if "simulation_parameters" in full_config:
        yaml_params = full_config["simulation_parameters"]

        # A YAML kulcsainak leképezése a C-s `options_t` struktúra tagneveire.
        # Minden paraméter egyetlen leképezésben van.
        yaml_to_c_mapping = {
            # Simulation Control Options (new names, mapped to C's internal names)
            "enable_dust_drift": "drift",
            "enable_dust_growth": "growth",
            "enable_gas_evolution": "evol",
            "enable_two_dust_populations": "twopop",
            "fragmentation_velocity": "ufrag",
            "fragmentation_factor": "ffrag",

            # Grid and Disk Initial Parameters (new names)
            "number_of_grid_points": "ngrid_val",
            "inner_radius_au": "rmin_val",
            "outer_radius_au": "rmax_val",
            "initial_gas_sigma0_msun_per_au2": "sigma0_val",
            "sigma_profile_exponent": "sigmap_exp_val",
            "alpha_viscosity": "alpha_visc_val",
            "star_mass_msun": "star_val",
            "aspect_ratio_at_1au": "hasp_val",
            "flaring_index": "flind_val",

            # Dead Zone Parameters (new names)
            "deadzone_inner_radius_au": "r_dze_i_val",
            "deadzone_outer_radius_au": "r_dze_o_val",
            "deadzone_inner_transition_width_mult": "dr_dze_i_val",
            "deadzone_outer_transition_width_mult": "dr_dze_o_val",
            "deadzone_alpha_multiplier": "a_mod_val",

            # Dust Initialization Parameters (new names)
            "initial_dust_to_gas_ratio": "eps_val",
            "population_one_mass_ratio": "ratio_val",
            "micro_particle_size_cm": "mic_val",
            "one_size_particle_value_cm": "onesize_val",

            # File I/O Parameters (new names)
            "input_file_path": "input_file",
            "output_directory_name": "output_dir_name",

            # Time Parameters (new names)
            "fixed_time_step": "tStep",
            "total_simulation_time": "totalTime",
            "output_write_frequency": "outputFrequency",
            "simulation_start_time": "startTime",
        }

        for yaml_key, c_key in yaml_to_c_mapping.items():
            if yaml_key in yaml_params:
                all_params[c_key] = yaml_params[yaml_key]
    else:
        print("Figyelem: 'simulation_parameters' szekció hiányzik a YAML-ből. Alapértelmezetteket használunk.")

    # A C program argumentum leképezése
    # Ezek a kulcsok a Python `all_params` szótárának kulcsai (a C struct tagnevek),
    # az értékek pedig a C program parancssori kapcsolói.
    c_arg_mapping = {
        # Simulation Control (C flags are unchanged)
        "drift": "-drift", "growth": "-growth", "evol": "-evol", "twopop": "-twopop",
        "ufrag": "-ufrag", "ffrag": "-ffrag",

        # Grid and Disk Initial Parameters
        "ngrid_val": "-n", "rmin_val": "-ri", "rmax_val": "-ro",
        "sigma0_val": "-sigma0_init", "sigmap_exp_val": "-index_init",
        "alpha_visc_val": "-alpha_init", "star_val": "-m0_init",
        "hasp_val": "-h_init", "flind_val": "-flind_init",

        # Dead Zone Parameters
        "r_dze_i_val": "-rdzei", "r_dze_o_val": "-rdzeo",
        "dr_dze_i_val": "-drdzei", "dr_dze_o_val": "-drdzeo",
        "a_mod_val": "-amod",

        # Dust Initialization Parameters
        "eps_val": "-eps", "ratio_val": "-ratio", "mic_val": "-mic", "onesize_val": "-onesize",

        # File I/O Parameters
        "input_file": "-i", "output_dir_name": "-o",

        # Time Parameters
        "tStep": "-tStep", "totalTime": "-tmax", "outputFrequency": "-outfreq", "startTime": "-curr",
    }

    # 2. A fő C program futtatása az összes beállított paraméterrel
    success, return_code = run_c_program(main_executable, all_params, c_arg_mapping, "Fő szimuláció program")

    if not success:
        print(f"A C program hibával fejeződött be (hibakód: {return_code}).")
        return

    print("\nMinden program sikeresen lefutott.")

if __name__ == "__main__":
    main()