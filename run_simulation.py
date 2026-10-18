import yaml
import subprocess
import argparse
import os
import datetime
import math

prefix = "[PY-WRAPPER]"


def validate_thermal_config(config: object) -> str:
    if not isinstance(config, dict):
        raise ValueError("Configuration root must be a YAML mapping.")
    if "thermal" in config:
        raise ValueError("Use 'thermodynamics_parameters' instead of 'thermal'.")
    thermal = config.get("thermodynamics_parameters", {})
    if not isinstance(thermal, dict):
        raise ValueError("'thermodynamics_parameters' must be a YAML mapping with a 'model' field.")
    numeric_settings = {
        "mean_molecular_weight": 2.34,
        "gamma": 1.4,
        "background_temperature_k": 10.0,
        "minimum_temperature_k": 1.0,
        "maximum_temperature_k": 1500.0,
        "irradiation_luminosity_lsun": 0.0,
        "irradiation_factor": 0.05,
        "irradiation_surface_height": 4.0,
        "irradiation_stellar_radius_rsun": 2.0,
    }
    for key in thermal:
        if key not in ("model", "irradiation_model", *numeric_settings):
            raise ValueError(f"Unknown thermal setting: {key!r}.")
    if thermal.get("irradiation_model", "prescribed") not in ("prescribed", "geometry"):
        raise ValueError("thermodynamics_parameters.irradiation_model must be 'prescribed' or 'geometry'.")
    model = thermal.get("model", "lociso")
    if model == "locally_isotherm":
        model = "lociso"
    if not isinstance(model, str) or model not in ("lociso", "radiative", "radiative_energy", "polytropic"):
        raise ValueError(
            "thermodynamics_parameters.model must be 'lociso' (or 'locally_isotherm'), "
            "'radiative', 'radiative_energy', or 'polytropic'."
        )
    if model == "polytropic":
        raise ValueError(
            "thermodynamics_parameters.model='polytropic' is not implemented yet. "
            "Use 'lociso', 'locally_isotherm', or 'radiative'."
        )
    for key, default in numeric_settings.items():
        value = thermal.get(key, default)
        irradiation = key in ("irradiation_luminosity_lsun", "irradiation_factor")
        if (isinstance(value, bool) or not isinstance(value, (int, float)) or
                not math.isfinite(value) or (value < 0 if irradiation else value <= 0)):
            requirement = "nonnegative" if irradiation else "positive"
            raise ValueError(f"thermodynamics_parameters.{key} must be a finite {requirement} number.")
        if key == "irradiation_factor" and value > 1:
            raise ValueError("thermodynamics_parameters.irradiation_factor must be in [0, 1].")
        if key == "gamma" and value <= 1:
            raise ValueError("thermodynamics_parameters.gamma must be greater than 1.")
    lower = thermal.get("minimum_temperature_k", 1.0)
    upper = thermal.get("maximum_temperature_k", 1500.0)
    background = thermal.get("background_temperature_k", 10.0)
    if not lower < upper < 1e9 or background > upper:
        raise ValueError("Invalid thermodynamics_parameters temperature bracket or background temperature.")
    return model


def run_c_program(executable_path, params, arg_mapping, verbosity_flag, disable_panels, program_name="C Program", omp_threads=None):
    """
    Runs a C program with the given parameters.
    """
    if not os.path.exists(executable_path):
        print(f"{prefix} Error: The '{program_name}' executable was not found at: {executable_path}.")
        return False, None

    cmd_args = [executable_path]
    
    if verbosity_flag:
        cmd_args.append(verbosity_flag)

    # YAML-based panel disable
    if disable_panels:
        cmd_args.append("-no_panels")
        print(f"{prefix} Terminal panels disabled (YAML setting).")


    for py_key, value in params.items():
        c_arg_name = arg_mapping.get(py_key)
        if c_arg_name:
            if isinstance(value, bool):
                cmd_args.extend([c_arg_name, "1.0" if value else "0.0"])
            elif c_arg_name == "-i" or c_arg_name == "-o" or c_arg_name == "-test":
                if value is not None and str(value).strip() != "":
                    cmd_args.extend([c_arg_name, str(value)])
            else:
                cmd_args.extend([c_arg_name, str(value)])

    print(f"\n{prefix} --- Running: {program_name} ---")
    print(f"{prefix} The current command-line arguments are:\n {' '.join(cmd_args)}")

    current_env = os.environ.copy()
    current_env["OMP_NUM_THREADS"] = (
        str(omp_threads)
        if omp_threads is not None
        else os.environ.get("OMP_NUM_THREADS", "1")
    )
    print(f"{prefix} Setting OMP_NUM_THREADS={current_env['OMP_NUM_THREADS']} for this run.")
    print(f"{prefix} Start running the binary ({executable_path}) at {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}...\n")

    process = None
    try:
        process = subprocess.Popen(
            cmd_args,
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
            print(f"\n{prefix} {program_name} exited with error code: {process.returncode}")
            return False, process.returncode
        else:
            print(f"\n{prefix} {program_name} completed successfully.")
            return True, 0
    except FileNotFoundError:
        print(f"{prefix} Error: Command '{executable_path}' not found.")
        return False, 1
    except KeyboardInterrupt:
        print(f"{prefix} Simulation aborted by user (Ctrl+C). Terminating process...")
        if process:
            process.terminate()
            try:
                process.wait(timeout=2)
            except subprocess.TimeoutExpired:
                process.kill()
        return False, 130
    except Exception as e:
        print(f"{prefix} An error occurred while running {program_name}: {e}")
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
        print(f"{prefix} Error: Configuration file '{config_file}' not found.")
        return
    except yaml.YAMLError as exc:
        print(f"{prefix} Error parsing YAML file '{config_file}': {exc}")
        return

    try:
        thermal_model = validate_thermal_config(full_config)
    except ValueError as exc:
        parser.error(str(exc))
    print(f"{prefix} Thermal model: {thermal_model}.")
    if thermal_model == "lociso" and any(
        key != "model" for key in full_config.get("thermodynamics_parameters", {})
    ):
        print(f"{prefix} Radiative parameters are inactive in locally isothermal mode; "
              "the prescribed aspect ratio and flaring are used.")

    # --- OMP threads from YAML (optional) ---
    omp_threads = full_config.get("simulation_parameters", {}).get("omp_num_threads", None)



    all_params = {}

    # YAML nested sections
    yaml_sections = [
        "simulation_parameters",
        "disk_parameters",
        "boundary_conditions",
        "deadzone_parameters",
        "dust_parameters",
        "output_parameters",
        "time_parameters",
        "log_parameters"
    ]

    # YAML → C mapping 
    yaml_to_c_mapping = {
        "enable_dust_drift": "drift",
        "enable_dust_growth": "growth",
        "enable_gas_evolution": "evol",
        "enable_photoevaporation": "photoevap",
        "enable_selfgravity": "selfgravity",
        "enable_two_dust_populations": "twopop",
        "test_mode": "test_mode",  # <--- Hozzáadva a teszt módszerhez

        "inner_boundary_condition": "inner_bc",
        "outer_boundary_condition": "outer_bc",

        "number_of_grid_points": "ngrid_val",
        "radial_grid_type": "grid_type",
        "number_of_dust_particles": "ndust_val",
        "inner_radius_au": "rmin_val",
        "outer_radius_au": "rmax_val",
        "disk_mass": "disk_mass",                     
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
        "density_floor": "density_floor",
        "dust_density_floor": "dust_density_floor",

        "gaussian_smoothing_sigma_grid_units": "gaussian_smoothing_sigma_grid_units",
        "gaussian_smoothing_cutoff_sigma": "gaussian_smoothing_cutoff_sigma",

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
        "dust_smoothing_mode": "dust_smoothing_mode",
    }

    # Read ONLY what is explicitly written in the YAML file
    for section in yaml_sections:
        if section in full_config:
            yaml_params = full_config[section]
            if yaml_params:
                for yaml_key, c_key in yaml_to_c_mapping.items():
                    if yaml_key in yaml_params:
                        all_params[c_key] = yaml_params[yaml_key]

    # C program argument mapping
    c_arg_mapping = {
        "drift": "-drift", "growth": "-growth", "evol": "-evol", "twopop": "-twopop",
        "photoevap": "-photoevap", "selfgravity": "-selfgravity", "test_mode": "-test",  # <--- Hozzáadva a --test argumentumhoz
        "ngrid_val": "-n", "ndust_val": "-ndust", 
        "grid_type": "-grid_type",
        "rmin_val": "-ri", "rmax_val": "-ro",
        "inner_bc": "-inner_bc", "outer_bc": "-outer_bc",
        "sigma0_val": "-sigma0_init", "sigmap_exp_val": "-index_init",
        "alpha_visc_val": "-alpha_init", "star_val": "-stellar_mass", "disk_mass": "-disk_mass",
        "hasp_val": "-h_init", "flind_val": "-flind_init", 
        "photoevap_mode": "-photoevap_mode", "xray_lum": "-xray_luminosity",
        "cutoff": "-cutoff", "cutoff_radius": "-cutoff_radius", "cutoff_sharpness": "-cutoff_sharpness",
        "r_dze_i_val": "-rdzei", "r_dze_o_val": "-rdzeo",
        "dr_dze_i_val": "-drdzei", "dr_dze_o_val": "-drdzeo",
        "a_mod_val": "-amod", "density_floor": "-density_floor", "dust_density_floor": "-dust_density_floor",
        "eps_val": "-eps", "ratio_val": "-ratio", "mic_val": "-micsize", "onesize_val": "-largesize",
        "pdensity_val": "-pdensity",
        "gaussian_smoothing_sigma_grid_units": "-gaussian_sigma_grid_units",
        "gaussian_smoothing_cutoff_sigma": "-gaussian_cutoff_sigma",
        "input_file": "-i", "output_dir_name": "-o",
        "output_format": "-output_format", "dust_smoothing_mode": "-dust_smoothing",
        "tStep": "-tStep", "totalTime": "-tmax", "outputFrequency": "-outfreq"
    }

    if thermal_model in ("radiative", "radiative_energy"):
        all_params["thermal_model"] = thermal_model
        c_arg_mapping["thermal_model"] = "-thermal_model"
        thermal_flags = {
            "mean_molecular_weight": "-thermal_mu",
            "gamma": "-thermal_gamma",
            "background_temperature_k": "-thermal_background",
            "minimum_temperature_k": "-thermal_tmin",
            "maximum_temperature_k": "-thermal_tmax",
            "irradiation_luminosity_lsun": "-irradiation_luminosity",
            "irradiation_factor": "-irradiation_factor",
            "irradiation_model": "-irradiation_model",
            "irradiation_surface_height": "-irradiation_surface_height",
            "irradiation_stellar_radius_rsun": "-irradiation_stellar_radius",
        }
        thermal = full_config["thermodynamics_parameters"]
        for key, flag in thermal_flags.items():
            if key in thermal:
                all_params[key] = thermal[key]
                c_arg_mapping[key] = flag

    if not all_params:
        print(f"{prefix} Error: No parameters could be parsed from YAML. Exiting.")
        return
    
    # 2. Flag for loggin
    verbosity_level = full_config.get("log_parameters", {}).get("info_level", "none")
    verbosity_flag = ""
    if verbosity_level == "info":
        verbosity_flag = "-v"
    elif verbosity_level == "debug":
        verbosity_flag = "-vv"

    # 3. Panel disable flag (YAML)
    disable_panels = full_config.get("log_parameters", {}).get("disable_terminal_panels", False)

    # 3. Run
    success, return_code = run_c_program(main_executable, all_params, c_arg_mapping, verbosity_flag, disable_panels, "Main Simulation Program", omp_threads)

    if not success:
        print(f"{prefix} The C program exited with an error (error code: {return_code}) with the Python wrapper.")
        raise SystemExit(return_code or 1)

    print(f"{prefix} All programs completed successfully.")

if __name__ == "__main__":
    main()