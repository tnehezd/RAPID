import yaml
import subprocess
import argparse
import os
import datetime

def run_c_program(executable_path, params, arg_mapping, program_name="C Program"):
    """
    Runs a C program with the given parameters.
    """
    if not os.path.exists(executable_path):
        print(f"Error: The '{program_name}' executable was not found at: {executable_path}. Have you compiled it?")
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

    print(f"\n--- Running: {program_name} ---")
    print(f"Command: {' '.join(cmd_args)}")

    # Current environment variables are copied
    current_env = os.environ.copy()
    # Set OMP_NUM_THREADS to 1
    current_env["OMP_NUM_THREADS"] = "1"
    print(f"Setting OMP_NUM_THREADS={current_env['OMP_NUM_THREADS']} for this run.")

    try:
        process = subprocess.Popen(
            cmd_args,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding='cp1252',
            errors='replace',
            bufsize=1,
            env=current_env # Pass the modified environment variables
        )
        for line in process.stdout:
            print(line, end='') # Print output in real-time
        process.wait()

        if process.returncode != 0:
            print(f"{program_name} exited with error code: {process.returncode}")
            return False, process.returncode
        else:
            print(f"{program_name} completed successfully.")
            return True, 0
    except FileNotFoundError:
        print(f"Error: Command '{executable_path}' not found. Check path and permissions.")
        return False, 1
    except Exception as e:
        print(f"An error occurred while running {program_name}: {e}")
        return False, 1


def main():
    parser = argparse.ArgumentParser(description="Runs the simulation workflow with YAML configuration.")
    parser.add_argument("-c", "--config", default="config.yaml",
                        help="Path to the YAML configuration file.")
    parser.add_argument("-exec", "--executable", default="./bin/simulation",
                        help="Path to the main C program executable.")
    args = parser.parse_args()

    config_file = args.config
    main_executable = args.executable

    # 1. Load configuration
    try:
        with open(config_file, 'r') as f:
            full_config = yaml.safe_load(f) or {}
    except FileNotFoundError:
        print(f"Error: Configuration file '{config_file}' not found.")
        return
    except yaml.YAMLError as exc:
        print(f"Error parsing YAML file '{config_file}': {exc}")
        return

    # Default options (reflecting the C program's options_t struct defaults)
    # These should **exactly** match the default values set in the C code's `create_default_options` function.
    default_options = {
        # Simulation Control Options (Boolean values are handled in run_c_program)
        "drift": 1.0,
        "growth": 1.0,
        "evol": 1.0,
        "photoevap": 0.0,
        "twopop": 1.0,
        "ufrag": 1000.0,
        "ffrag": 0.37,

        # Core Disk Parameters (also serve as init_tool defaults if no input file)
        "ngrid_val": 2000,
        "ndust_val": 5000, # <--- MODIFIED: Changed key to 'ndust_val' to match C struct
        "rmin_val": 1.0,
        "rmax_val": 100.0,
        "sigma0_val": 1.0,
        "sigmap_exp_val": 0.5,
        "alpha_visc_val": 0.01,
        "star_val": 1.0,
        "hasp_val": 0.05,
        "flind_val": 0.5,
        "r_dze_i_val": 0.0,
        "r_dze_o_val": 0.0,
        "dr_dze_i_val": 0.0,
        "dr_dze_o_val": 0.0,
        "a_mod_val": 0.0,
        "photoevap_mode": "none",
        "xray_lum": 1e30,
        "cutoff": "none",
        "cutoff_radius": 30.0,
        "cutoff_sharpness": 1.5,

        # File Input/Output
        "input_file": "",
        "output_dir_name": "output",

        # Time Parameters
        "tStep": 0.0,
        "totalTime": 1.0e6,
        "outputFrequency": 1000.0,

        # Init tool specific parameters' defaults
        "eps_val": 0.01,
        "ratio_val": 0.85,
        "mic_val": 1e-4,
        "onesize_val": 0.0,
        "pdensity_val": 1.6,
        "output_format": "ascii",
    }

    all_params = default_options.copy()

    # YAML nested sections
    yaml_sections = [
        "simulation_parameters",
        "disk_parameters",
        "deadzone_parameters",
        "dust_parameters",
        "output_parameters",
        "time_parameters"
    ]

    # YAML → C mapping 
    yaml_to_c_mapping = {
        "enable_dust_drift": "drift",
        "enable_dust_growth": "growth",
        "enable_gas_evolution": "evol",
        "enable_photoevaporation": "photoevap",
        "enable_two_dust_populations": "twopop",
        "fragmentation_velocity": "ufrag",
        "fragmentation_factor": "ffrag",

        "number_of_grid_points": "ngrid_val",
        "number_of_dust_particles": "ndust_val",
        "inner_radius_au": "rmin_val",
        "outer_radius_au": "rmax_val",
        "initial_gas_sigma0_msun_per_au2": "sigma0_val",
        "sigma_profile_exponent": "sigmap_exp_val",
        "alpha_viscosity": "alpha_visc_val",
        "star_mass_msun": "star_val",
        "aspect_ratio_at_1au": "hasp_val",
        "flaring_index": "flind_val",
        "photoevaporation_mode": "photoevap_mode",
        "xray_luminosity_erg_s": "xray_lum",
        "use_cutoff_for_gas": "cutoff",
        "characteristic_cutoff_radius_au": "cutoff_radius",
        "cutoff_sharpness_factor": "cutoff_sharpness",

        "deadzone_inner_radius_au": "r_dze_i_val",
        "deadzone_outer_radius_au": "r_dze_o_val",
        "deadzone_inner_transition_width_mult": "dr_dze_i_val",
        "deadzone_outer_transition_width_mult": "dr_dze_o_val",
        "deadzone_alpha_reduction": "a_mod_val",

        "initial_dust_to_gas_ratio": "eps_val",
        "population_one_mass_ratio": "ratio_val",
        "micro_particle_size_cm": "mic_val",
        "one_size_particle_value_cm": "onesize_val",
        "dust_particle_density_g_cm3": "pdensity_val",

        "input_file_path": "input_file",
        "output_directory_name": "output_dir_name",
        "output_format": "output_format",

        "fixed_time_step": "tStep",
        "total_simulation_time": "totalTime",
        "output_write_frequency": "outputFrequency",
    }

    # Flag: does the file contain one even section
    found_any_section = False

    # Read all the sections
    for section in yaml_sections:
        if section in full_config:
            found_any_section = True
            yaml_params = full_config[section]
            for yaml_key, c_key in yaml_to_c_mapping.items():
                if yaml_key in yaml_params:
                    all_params[c_key] = yaml_params[yaml_key]

    # If none of the above sections found:
    if not found_any_section:
        print("Warning: No recognized parameter sections found in YAML. Using defaults.")


    # C program argument mapping
    # These keys are the keys in the Python `all_params` dictionary (now matching C struct member names),
    # and the values are the C program's command-line flags.
    c_arg_mapping = {
        # Simulation Control (C flags are unchanged)
        "drift": "-drift", "growth": "-growth", "evol": "-evol", "twopop": "-twopop",
        "ufrag": "-ufrag", "ffrag": "-ffrag", "photoevap": "-photoevap",

        # Grid and Disk Initial Parameters
        "ngrid_val": "-n", "ndust_val": "-ndust", # <--- MODIFIED: Key is now 'ndust_val'
        "rmin_val": "-ri", "rmax_val": "-ro",
        "sigma0_val": "-sigma0_init", "sigmap_exp_val": "-index_init",
        "alpha_visc_val": "-alpha_init", "star_val": "-m0_init",
        "hasp_val": "-h_init", "flind_val": "-flind_init", 
        "photoevap_mode": "-photoevap_mode", "xray_lum": "-xray_luminosity",
        "cutoff": "-cutoff", "cutoff_radius": "-cutoff_radius",  "cutoff_sharpness": "-cutoff_sharpness",

        # Dead Zone Parameters
        "r_dze_i_val": "-rdzei", "r_dze_o_val": "-rdzeo",
        "dr_dze_i_val": "-drdzei", "dr_dze_o_val": "-drdzeo",
        "a_mod_val": "-amod",

        # Dust Initialization Parameters
        "eps_val": "-eps", "ratio_val": "-ratio", "mic_val": "-mic", "onesize_val": "-onesize",
        "pdensity_val": "-pdensity",

        # File I/O Parameters
        "input_file": "-i", "output_dir_name": "-o",
        "output_format": "--output-format",


        # Time Parameters
        "tStep": "-tStep", "totalTime": "-tmax", "outputFrequency": "-outfreq"
    }

    # 2. Run the main C program with all configured parameters
    success, return_code = run_c_program(main_executable, all_params, c_arg_mapping, "Main Simulation Program")

    if not success:
        print(f"The C program exited with an error (error code: {return_code}).")
        return

    print("\nAll programs completed successfully.")

if __name__ == "__main__":
    main()