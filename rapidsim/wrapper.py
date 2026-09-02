import argparse
import os
import yaml
import subprocess
import datetime

PREFIX = "[PY-WRAPPER: RAPIDSIM]"

def get_version():
    # 1. Megpróbáljuk kiolvasni a feltelepített csomag metaadatából (ez a legtisztább)
    try:
        import importlib.metadata as importlib_metadata
        return importlib_metadata.version("rapidsim")
    except Exception:
        pass

    # 2. Ha nincs feltelepítve, vagy fejlesztői módban vagyunk, olvassuk ki élőben a pyproject.toml-ből vagy setup.py-ból
    try:
        package_dir = os.path.dirname(os.path.abspath(__file__))
        # Megkeressük a gyökeret (akár a wrapper mappájából, akár egy szinttel feljebb)
        for path_candidate in [
            os.path.join(package_dir, "..", "pyproject.toml"),
            os.path.join(package_dir, "pyproject.toml"),
            os.path.join(package_dir, "..", "setup.py"),
            os.path.join(package_dir, "setup.py")
        ]:
            if os.path.exists(path_candidate):
                with open(path_candidate, "r", encoding="utf-8") as f:
                    content = f.read()
                    match = re.search(r'version\s*=\s*["\']([^"\']+)["\']', content)
                    if match:
                        return match.group(1)
    except Exception:
        pass

    return "unknown"

def main():
    parser = argparse.ArgumentParser(description="Runs the simulation workflow with YAML configuration.")
    parser.add_argument("-c", "--config", default="config.yaml",
                        help="Path to the YAML configuration file.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {get_version()}")
    
    args = parser.parse_args()

    config_file = os.path.abspath(args.config)

    if not os.path.exists(config_file):
        print(f"{PREFIX} Error: Configuration file '{config_file}' not found.")
        return

    try:
        with open(config_file, 'r') as f:
            full_config = yaml.safe_load(f) or {}
    except yaml.YAMLError as exc:
        print(f"{PREFIX} Error parsing YAML file '{config_file}': {exc}")
        return

    all_params = {}

    yaml_sections = [
        "simulation_parameters", "disk_parameters", "boundary_conditions",
        "deadzone_parameters", "dust_parameters", "output_parameters",
        "time_parameters", "log_parameters"
    ]

    yaml_to_c_mapping = {
        "enable_dust_drift": "drift", "enable_dust_growth": "growth",
        "enable_gas_evolution": "evol", "enable_photoevaporation": "photoevap",
        "enable_two_dust_populations": "twopop", "fragmentation_velocity": "ufrag",
        "fragmentation_factor": "ffrag", "inner_boundary_condition": "inner_bc",
        "outer_boundary_condition": "outer_bc", "number_of_grid_points": "ngrid_val",
        "number_of_dust_particles": "ndust_val", "inner_radius_au": "rmin_val",
        "outer_radius_au": "rmax_val", "initial_gas_sigma0_msun_per_au2": "sigma0_val",
        "disk_mass": "disk_mass", "total_disk_mass": "disk_mass",
        "sigma_profile_exponent": "sigmap_exp_val", "alpha_viscosity": "alpha_visc_val",
        "star_mass_msun": "star_val", "aspect_ratio_at_1au": "hasp_val",
        "flaring_index": "flind_val", "photoevaporation_mode": "photoevap_mode",
        "xray_luminosity_erg_s": "xray_lum", "use_cutoff_for_gas": "cutoff",
        "characteristic_cutoff_radius_au": "cutoff_radius",
        "cutoff_sharpness_factor": "cutoff_sharpness", "density_floor": "density_floor",
        "dust_density_floor": "dust_density_floor",
        "gaussian_smoothing_sigma_grid_units": "gaussian_smoothing_sigma_grid_units",
        "gaussian_smoothing_cutoff_sigma": "gaussian_smoothing_cutoff_sigma",
        "deadzone_inner_radius_au": "r_dze_i_val", "deadzone_outer_radius_au": "r_dze_o_val",
        "deadzone_inner_transition_width_mult": "dr_dze_i_val",
        "deadzone_outer_transition_width_mult": "dr_dze_o_val",
        "deadzone_alpha_reduction": "a_mod_val", "initial_dust_to_gas_ratio": "eps_val",
        "population_one_mass_ratio": "ratio_val", "micro_particle_size_cm": "mic_val",
        "one_size_particle_value_cm": "onesize_val", "dust_particle_density_g_cm3": "pdensity_val",
        "input_file_path": "input_file", "output_directory_name": "output_dir_name",
        "output_format": "output_format", "fixed_time_step": "tStep",
        "total_simulation_time": "totalTime", "output_write_frequency": "outputFrequency",
        "dust_smoothing_mode": "dust_smoothing_mode",
    }

    for section in yaml_sections:
        if section in full_config:
            yaml_params = full_config[section]
            if yaml_params:
                for yaml_key, c_key in yaml_to_c_mapping.items():
                    if yaml_key in yaml_params:
                        all_params[c_key] = yaml_params[yaml_key]

    c_arg_mapping = {
        "drift": "-drift", "growth": "-growth", "evol": "-evol", "twopop": "-twopop",
        "ufrag": "-ufrag", "ffrag": "-ffrag", "photoevap": "-photoevap",
        "ngrid_val": "-n", "ndust_val": "-ndust", 
        "rmin_val": "-ri", "rmax_val": "-ro",
        "inner_bc": "-inner_bc", "outer_bc": "-outer_bc",
        "sigma0_val": "-sigma0_init", "sigmap_exp_val": "-index_init",
        "alpha_visc_val": "-alpha_init", "star_val": "-stellar_mass", "disk_mass": "-disk_mass",
        "hasp_val": "-h_init", "flind_val": "-flind_init", 
        "photoevap_mode": "-photoevap_mode", "xray_lum": "-xray_luminosity",
        "cutoff": "-cutoff", "cutoff_radius": "-cutoff_radius",  "cutoff_sharpness": "-cutoff_sharpness",
        "r_dze_i_val": "-rdzei", "r_dze_o_val": "-rdzeo",
        "dr_dze_i_val": "-drdzei", "dr_dze_o_val": "-drdzeo",
        "a_mod_val": "-amod", "density_floor": "-density_floor", "dust_density_floor": "-dust_density_floor",
        "eps_val": "-eps", "ratio_val": "-ratio", "mic_val": "-mic", "onesize_val": "-onesize",
        "pdensity_val": "-pdensity",
        
        # Gaussian smoothing mappings expected by the C parser
        "gaussian_smoothing_sigma_grid_units": "-gaussian_sigma_grid_units",
        "gaussian_smoothing_cutoff_sigma": "-gaussian_cutoff_sigma", 
        
        "input_file": "-i", "output_dir_name": "-o",
        "output_format": "--output-format", "dust_smoothing_mode": "-dust_smoothing",
        "tStep": "-tStep", "totalTime": "-tmax", "outputFrequency": "-outfreq"
    }

    cmd_args = []
    verbosity_level = full_config.get("log_parameters", {}).get("info_level", "none")
    if verbosity_level == "info":
        cmd_args.append("-v")
    elif verbosity_level == "debug":
        cmd_args.append("-vv")

    for py_key, value in all_params.items():
        c_arg_name = c_arg_mapping.get(py_key)
        if c_arg_name:
            if isinstance(value, bool):
                cmd_args.extend([c_arg_name, "1.0" if value else "0.0"])
            elif c_arg_name == "-i":
                if value is not None and str(value).strip() != "":
                    cmd_args.extend([c_arg_name, str(value)])
            elif c_arg_name == "-o":
                if value is not None and str(value).strip() != "":
                    cmd_args.extend([c_arg_name, str(value)])
                else:
                    cmd_args.extend([c_arg_name, "output"])
            else:
                cmd_args.extend([c_arg_name, str(value)])

    # Robust binary location resolver (checking data/ directory first)
    package_dir = os.path.dirname(os.path.abspath(__file__))
    possible_paths = [
        os.path.join(package_dir, "data", "simulation"),          # 1. Inside installed package data/
        os.path.join(package_dir, "bin", "simulation"),          # 2. Legacy / fallback inside package
        os.path.abspath(os.path.join(package_dir, "..", "bin", "simulation")), # 3. Parent relative
        os.path.abspath("./bin/simulation"),                     # 4. Current working directory bin/
        os.path.abspath("../bin/simulation")                     # 5. One level up cwd
    ]

    binary_path = None
    for p in possible_paths:
        if os.path.exists(p):
            binary_path = p
            break

    if not binary_path or not os.path.exists(binary_path):
        print(f"{PREFIX} Error: Could not find 'simulation' binary. Tried paths:\n" + "\n".join(possible_paths))
        return

    full_cmd = [binary_path] + cmd_args

    print(f"\n{PREFIX} --- Running: Main Simulation Program ---")
    print(f"{PREFIX} The current command-line arguments are:\n simulation " + " ".join(cmd_args))

    current_env = os.environ.copy()
    current_env["OMP_NUM_THREADS"] = "1"
    print(f"{PREFIX} Setting OMP_NUM_THREADS={current_env['OMP_NUM_THREADS']} for this run.")
    print(f"{PREFIX} Start running the binary ({binary_path}) at {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}...\n")

    process = None
    try:
        process = subprocess.Popen(
            full_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding='cp1252',
            errors='replace',
            bufsize=1,
            env=current_env
        )

        for line in process.stdout:
            print(line, end='')

        process.wait()

        if process.returncode != 0:
            print(f"\n{PREFIX} Main Simulation Program exited with error code: {process.returncode}")
        else:
            print(f"\n{PREFIX} Main Simulation Program completed successfully.")

    except KeyboardInterrupt:
        print(f"\n{PREFIX} Simulation aborted by user (Ctrl+C). Terminating process...")
        if process:
            process.terminate()
            try:
                process.wait(timeout=2)
            except subprocess.TimeoutExpired:
                process.kill()
        print(f"{PREFIX} The C program exited with an error (error code: 130) with the Python wrapper.")
    except Exception as e:
        print(f"\n{PREFIX} An error occurred: {e}")
        if process:
            process.terminate()

if __name__ == "__main__":
    main()