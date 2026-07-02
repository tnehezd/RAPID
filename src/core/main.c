/**
 * @file main.c
 * @brief Main execution core and simulation orchestrator.
 * @date 2026-06-26
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "config.h"       
#include "init_tool_module.h"
#include "ascii_output.h"        
#include "disk_model.h"      
#include "dust_physics.h"    
#include "simulation_core.h" 
#include "utils.h"           
#include "gas_physics.h"
#include "boundary_conditions.h"
#include "integrator.h"
#include "simulation_types.h"
#include "parser.h"
#include "print_panels.h"    
#include "logger.h" 


extern void initializeDefaultOptions(InitializeDefaultOptions *def);

int main(int argc, const char **argv) {

    time_t wall_start, wall_end;
    time(&wall_start); 

    printSimulationWelcomeBanner();

    char *initial_dir_path = NULL;
    char *kLogFilesDirectory_path = NULL;
    char *cmd_buffer = NULL;
    char *current_inputsig_file = NULL;
    char *current_inputdust_file = NULL;
    char *dens_name_initial = NULL;
    char *actual_run_dir = NULL;

    ParserOptions def;
    createDefaultOptions(&def);

    InitializeDefaultOptions init_tool_params;
    initializeDefaultOptions(&init_tool_params);

    int retCode = parseCLIOptions(argc, argv, &def);

    if (0 != retCode) {
        MSG("************************************************************************");
        LOG_ERROR("Error parsing command-line options. Exiting with code %d.", retCode);
        MSG("************************************************************************");
        return retCode;
    }

    g_log_level = def.verbose;

    LOG_INFO("Command-line options parsed successfully. Logging level set to %d.", g_log_level);
    
    size_t len = strlen(def.output_dir_name);

    if (len > 0 && def.output_dir_name[len - 1] == '/') {
        def.output_dir_name[len - 1] = '\0';
    }

    DiskParameters disk_params; 
    SimulationOptions sim_opts;
    OutputFiles output_files;
    
    output_files.dust_motion_file = NULL;
    output_files.micron_motion_file = NULL;
    output_files.mass_file = NULL;
    output_files.surface_file = NULL;
    output_files.dust_file = NULL;
    output_files.micron_dust_file = NULL;
    output_files.size_file = NULL;

    sim_opts.option_for_evolution = def.option_for_evolution;
    sim_opts.option_for_dust_drift = def.option_for_dust_drift;
    sim_opts.option_for_dust_growth = def.option_for_dust_growth;
    sim_opts.option_for_dust_secondary_population = def.option_for_dust_secondary_population;
    sim_opts.user_defined_time_step = def.user_defined_time_step;
    sim_opts.maximum_simulation_time = def.maximum_simulation_time;
    sim_opts.output_frequency = def.output_frequency;
    sim_opts.number_of_dust_particles = def.number_of_dust_particles;
    sim_opts.output_format = def.output_format;
    sim_opts.dust_smoothing_mode = def.dust_smoothing_mode;
    sim_opts.gaussian_sigma = def.gaussian_sigma;
    sim_opts.gaussian_cutoff     = def.gaussian_cutoff;
    sim_opts.inner_boundary_condition_type = def.inner_boundary_condition_type;
    sim_opts.outer_boundary_condition_type = def.outer_boundary_condition_type;

 

    LOG_DEBUG("def.output_dir_name BEFORE sim_opts population: '%s'", def.output_dir_name);
    LOG_DEBUG("sim_opts.option_for_evolution=%.2f) or drift (sim_opts.option_for_dust_drift=%.2f) is ON. Starting main simulation loop.", sim_opts.option_for_evolution, sim_opts.option_for_dust_drift);

    SnapshotMode mode = determineSnapshotMode(&sim_opts);
    LOG_DEBUG("SnapshotMode = %s", snapshotModeToString(mode));

    disk_params.r_min = def.rmin_val;
    disk_params.r_max = def.rmax_val;
    disk_params.grid_number = def.number_of_grid_points; 
    disk_params.sigma_0 = def.sigma0_val;
    disk_params.sigma_power_law_index = def.sigmap_exp_val;
    disk_params.alpha_parameter = def.alpha_visc_val;
    disk_params.stellar_mass = def.star_val;
    disk_params.h_aspect_ratio = def.hasp_val;
    disk_params.flaring_index = def.flind_val;
    disk_params.r_dze_i = def.r_dze_i_val;
    disk_params.r_dze_o = def.r_dze_o_val;
    disk_params.dr_dze_i = def.dr_dze_i_val;
    disk_params.dr_dze_o = def.dr_dze_o_val;
    disk_params.alpha_parameter_modification = def.a_mod_val;
    disk_params.fragmentation_factor = def.fragmenatation_factor;
    disk_params.fragmentation_velocity = def.fragmenatation_velocity;
    disk_params.drift_factor = 0.55; // set by Birnstiel 2012
    disk_params.particle_density = def.pdensity_val;
    disk_params.total_disk_mass = def.total_disk_mass;    
    disk_params.gaussian_sigma  = sim_opts.gaussian_sigma;
    disk_params.gaussian_cutoff = sim_opts.gaussian_cutoff;

    // --- BOUNDARY CONDITION STRING MAPPING ---
    switch (sim_opts.inner_boundary_condition_type) {
        case 0: strcpy(disk_params.inner_bc_string, "zero"); break;
        case 1: strcpy(disk_params.inner_bc_string, "parabolic"); break;
        case 2: strcpy(disk_params.inner_bc_string, "fixed_flux"); break;
        case 3: strcpy(disk_params.inner_bc_string, "absorbing"); break;
        case 4: strcpy(disk_params.inner_bc_string, "reflecting"); break;
        case 5: strcpy(disk_params.inner_bc_string, "linear"); break;
        case 6: strcpy(disk_params.inner_bc_string, "loggrid"); break;
    }

    switch (sim_opts.outer_boundary_condition_type) {
        case 0: strcpy(disk_params.outer_bc_string, "zero"); break;
        case 1: strcpy(disk_params.outer_bc_string, "parabolic"); break;
        case 2: strcpy(disk_params.outer_bc_string, "fixed_flux"); break;
        case 3: strcpy(disk_params.outer_bc_string, "absorbing"); break;
        case 4: strcpy(disk_params.outer_bc_string, "reflecting"); break;
        case 5: strcpy(disk_params.outer_bc_string, "linear"); break;
        case 6: strcpy(disk_params.outer_bc_string, "loggrid"); break;
    }

    // --- DUST SMOOTHING STRING MAPPING ---
    switch (sim_opts.dust_smoothing_mode) {
        case 0: strcpy(disk_params.dust_smoothing_mode_string, "CIC"); break;
        case 1: strcpy(disk_params.dust_smoothing_mode_string, "NGP"); break;
        case 2: strcpy(disk_params.dust_smoothing_mode_string, "TopHat"); break;
        case 3: strcpy(disk_params.dust_smoothing_mode_string, "Gaussian"); break;
    }




    // --- PHOTOEVAPORATION PARAMETER MAPPING AND WARNING LOGIC ---
    disk_params.enable_photoevaporation = def.enable_photoevaporation;
    disk_params.xray_luminosity = def.xray_luminosity;
    strncpy(disk_params.photoevaporation_mode_string, def.photoevaporation_mode, sizeof(disk_params.photoevaporation_mode_string) - 1);
    disk_params.photoevaporation_mode_string[sizeof(disk_params.photoevaporation_mode_string) - 1] = '\0';

    // Warn the user if a specific model is configured but photoevaporation is globally disabled
    if (!disk_params.enable_photoevaporation) {
        if (strcasecmp(disk_params.photoevaporation_mode_string, "none") != 0 && 
            strcasecmp(disk_params.photoevaporation_mode_string, "") != 0) {
                printWarningPanelForPhotoevaporation(&disk_params);
        }
    }  

    sim_opts.flag_for_deadzone = (disk_params.r_dze_i > 0.0 || disk_params.r_dze_o > 0.0) ? 1.0 : 0.0;

    actual_run_dir = createRunDirectory(def.output_dir_name);

    asprintf(&initial_dir_path, "%s/%s", actual_run_dir, kConfigFilesDirectory);
    createRunDirectory(initial_dir_path);

    asprintf(&kLogFilesDirectory_path, "%s/%s", actual_run_dir, kLogFilesDirectory);
    LOG_DEBUG("kLogFilesDirectory_path assembled as: '%s'", kLogFilesDirectory_path);

    createRunDirectory(kLogFilesDirectory_path);

    strncpy(sim_opts.output_dir_name, actual_run_dir, MAX_PATH_LEN - 1);
    sim_opts.output_dir_name[MAX_PATH_LEN - 1] = '\0'; // Ensure null-termination
    LOG_DEBUG("sim_opts.output_dir_name AFTER population: '%s'", sim_opts.output_dir_name);
    LOG_INFO("Output subdirectories created.");

    int dummy_sys_ret; 
    disk_params.sigma_dot_photoevap = NULL;

    if (def.input_file != NULL && strcmp(def.input_file, "") != 0) {
        strncpy(current_inputsig_file, def.input_file, MAX_PATH_LEN - 1);
        current_inputsig_file[MAX_PATH_LEN - 1] = '\0';
        LOG_DEBUG("Input file specified: '%s'. Attempting to read initial profile.", current_inputsig_file);

        disk_params.grid_number = calculateNumbersOfParticles(current_inputsig_file); 

        if (disk_params.grid_number > 1) {
            disk_params.delta_r = (disk_params.r_max - disk_params.r_min) / (disk_params.grid_number - 1.0);
        } else {
            disk_params.delta_r = 0.0;
        }
        LOG_DEBUG("grid_number set from input file: %d. delta_r calculated as %.4e.", disk_params.grid_number, disk_params.delta_r);

        disk_params.radial_grid = (double *)malloc((disk_params.grid_number + 2) * sizeof(double));
        disk_params.gas_surface_density_vector = (double *)malloc((disk_params.grid_number + 2) * sizeof(double));
        disk_params.gas_pressure_vector = (double *)malloc((disk_params.grid_number + 2) * sizeof(double));
        disk_params.gas_pressure_gradient_vector = (double *)malloc((disk_params.grid_number + 2) * sizeof(double));
        disk_params.gas_velocity_vector = (double *)malloc((disk_params.grid_number + 2) * sizeof(double));

        if (!disk_params.radial_grid || !disk_params.gas_surface_density_vector || !disk_params.gas_pressure_vector || !disk_params.gas_pressure_gradient_vector || !disk_params.gas_velocity_vector) {
            MSG("************************************************************************");
            LOG_ERROR("Failed to allocate memory for disk arrays (input file branch). Exiting.");
            MSG("************************************************************************");
            return 1;
        }
        LOG_DEBUG("Disk profile arrays dynamically allocated with size grid_number+2 = %d (input file branch).", disk_params.grid_number + 2);

        LOG_DEBUG("Calling readDiskParameters to calculate derived disk parameters for main disk_params struct (input file branch).");
        readDiskParameters(&disk_params);
        LOG_DEBUG("readDiskParameters completed (input file branch).");
        asprintf(&cmd_buffer, "cp %s %s/", current_inputsig_file, initial_dir_path);
        dummy_sys_ret = system(cmd_buffer); (void)dummy_sys_ret;
        LOG_DEBUG("Copied initial profile file '%s' to '%s/'.", current_inputsig_file, initial_dir_path);
        asprintf(&cmd_buffer, "cp %s%s %s/", kDiskConfigFile,kFileNamesSuffix, initial_dir_path);
        dummy_sys_ret = system(cmd_buffer); (void)dummy_sys_ret;
        LOG_DEBUG("Copied %s to %s%s/", kDiskConfigFile,kFileNamesSuffix, initial_dir_path);

    } else {
        LOG_WARN("No input file specified (-i flag not used). Generating default grid and profile.");

        // Map cutoff options directly from parser (def) to init_tool_params without touching disk_params struct
        init_tool_params.use_cutoff = def.use_cutoff;
        init_tool_params.r_cutoff = def.r_cutoff;
        init_tool_params.n_for_cutoff = def.n_for_cutoff;

        init_tool_params.n_grid_points = disk_params.grid_number; 
        init_tool_params.r_inner= disk_params.r_min;
        init_tool_params.r_outer = disk_params.r_max;
        init_tool_params.sigma0_gas_au = disk_params.sigma_0;
        init_tool_params.sigma_exponent = disk_params.sigma_power_law_index;
        init_tool_params.deadzone_r_inner = disk_params.r_dze_i;
        init_tool_params.deadzone_r_outer = disk_params.r_dze_o;
        init_tool_params.deadzone_dr_inner = disk_params.dr_dze_i;
        init_tool_params.deadzone_dr_outer = disk_params.dr_dze_o;
        init_tool_params.alpha_viscosity = disk_params.alpha_parameter;
        init_tool_params.deadzone_alpha_mod = disk_params.alpha_parameter_modification;
        init_tool_params.aspect_ratio = disk_params.h_aspect_ratio;
        init_tool_params.flaring_index = disk_params.flaring_index;
        init_tool_params.star_mass = disk_params.stellar_mass;
        init_tool_params.dust_to_gas_ratio = def.eps_val;
        init_tool_params.n_dust_particles = def.number_of_dust_particles; 
        init_tool_params.two_pop_ratio = def.ratio_val;
        init_tool_params.micro_size_cm = def.mic_val;
        init_tool_params.one_size_particle_cm = def.onesize_val;
        init_tool_params.dust_density_g_cm3 = def.pdensity_val;
        init_tool_params.total_disk_mass = disk_params.total_disk_mass;
        init_tool_params.density_floor = disk_params.density_floor;


        LOG_DEBUG("InitializeDefaultOptions (init_tool_params) structure populated for profile generation.");
        LOG_DEBUG("Calling runInitialization(&init_tool_params, &disk_params)...");
        strncpy(init_tool_params.output_base_path, initial_dir_path, MAX_PATH_LEN - 1);
        runInitialization(&init_tool_params, &disk_params,&sim_opts);
        LOG_DEBUG("runInitialization completed. disk_params allocated and populated.");
        asprintf(&current_inputsig_file, "%s/%s%s", initial_dir_path, kInitialGasProfileFileName,kFileNamesSuffix);
        LOG_DEBUG("Generated GAS profile will be loaded from '%s'.", current_inputsig_file);

        disk_params.grid_number = calculateNumbersOfParticles(current_inputsig_file);

        if (disk_params.grid_number > 1) {
            disk_params.delta_r = (disk_params.r_max - disk_params.r_min) / (disk_params.grid_number - 1.0);
        } else {
            disk_params.delta_r = 0.0;
        }
        LOG_DEBUG("grid_number updated from generated file: %d. delta_r calculated as %.4e.", disk_params.grid_number, disk_params.delta_r);
    }

    strncpy(sim_opts.input_filename, current_inputsig_file, MAX_PATH_LEN - 1);
    sim_opts.input_filename[MAX_PATH_LEN - 1] = '\0'; 
    LOG_DEBUG("sim_opts.input_filename set to '%s' for timeIntegrationForTheSystem (gas profile).", sim_opts.input_filename);
    asprintf(&current_inputdust_file, "%s/%s%s", initial_dir_path, kInitialDustProfileFileName, kFileNamesSuffix);
    strncpy(sim_opts.dust_input_filename, current_inputdust_file, MAX_PATH_LEN - 1);
    sim_opts.dust_input_filename[MAX_PATH_LEN - 1] = '\0'; 
    LOG_DEBUG("sim_opts.dust_input_filename set to '%s' for timeIntegrationForTheSystem (dust profile).", sim_opts.dust_input_filename);
    particle_number = calculateNumbersOfParticles(sim_opts.dust_input_filename); 
    LOG_DEBUG("Global particle_number set to %d (from dust input file: %s).", particle_number, sim_opts.dust_input_filename);
    LOG_DEBUG("Initial profile loading for loadGasSurfaceDensityFromFile...");
    loadGasSurfaceDensityFromFile(&disk_params, current_inputsig_file); 
    LOG_DEBUG("loadGasSurfaceDensityFromFile completed. Calling applyBoundaryConditions for disk_params.radial_grid and disk_params.gas_surface_density_vector...");
    sim_opts.current_bc_target = 0;
    applyBoundaryConditions(disk_params.gas_surface_density_vector, &disk_params, &sim_opts);

    // This guarantees that the core integrator always has a valid, zeroed array for photoevaporation
    if (disk_params.sigma_dot_photoevap == NULL) {
        disk_params.sigma_dot_photoevap = (double *)calloc((disk_params.grid_number + 2), sizeof(double));
        if (!disk_params.sigma_dot_photoevap) {
            LOG_ERROR("Failed to allocate memory for sigma_dot_photoevap.");
            return 1;
        }
    }

    LOG_DEBUG("applyBoundaryConditions calls completed for initial profile.");
    LOG_DEBUG("Calling printCurrentInformationAboutRun...");
    printCurrentInformationAboutRun(actual_run_dir, &disk_params, &sim_opts);

    // --- PRINT COMPREHENSIVE RUN STATUS INFO PANEL ---
    printRunConfigurationHeader(&def, &disk_params);
    printBoundaryConditionsStatus(&disk_params);
    printDustSmoothingStatus(&sim_opts, &disk_params);

    if(mode == SnapshotNonevolving) {
        LOG_DEBUG("Evolution (sim_opts.option_for_evolution=%.2f) and drift (sim_opts.option_for_dust_drift=%.2f) are OFF.", sim_opts.option_for_evolution, sim_opts.option_for_dust_drift);

        asprintf(&dens_name_initial, "%s/%s%s", initial_dir_path, kInitialGasProfileFileName,kFileNamesSuffix);
        LOG_DEBUG("Printing initial surface density to %s.", dens_name_initial);

        OutputFiles temp_output_for_initial_print;
        temp_output_for_initial_print.surface_file = fopen(dens_name_initial, "w");
        if (temp_output_for_initial_print.surface_file != NULL) {
            fprintf(temp_output_for_initial_print.surface_file, "# Initial Disk Profile Data (Gas Surface Density, Pressure, Gradient)\n");
            fprintf(temp_output_for_initial_print.surface_file, "# Columns: 1. Radius [AU], 2. Sigma_gas [M_Sun/AU^2],\n");
            fprintf(temp_output_for_initial_print.surface_file, "#   3. Pressure [units], 4. dP/dR [units]\n");
            fprintf(temp_output_for_initial_print.surface_file, "#\n");
            fprintf(temp_output_for_initial_print.surface_file, "# Data generated by Dust Drift Simulation (Initial State)\n");
            fflush(temp_output_for_initial_print.surface_file);
            printGasSurfaceDensityPressurePressureDerivateFile(&disk_params, &temp_output_for_initial_print);
            fclose(temp_output_for_initial_print.surface_file);
            temp_output_for_initial_print.surface_file = NULL;
            LOG_DEBUG("Closed %s.", dens_name_initial);
        } else {
            LOG_ERROR("Could not open %s for initial surface output.", dens_name_initial);
        }

        LOG_DEBUG("printGasSurfaceDensityPressurePressureDerivateFile completed. Program exiting.");
    } else {

        printCenteredText("\nINITIALIZATION COMPLETE.\n");

        printHeader("STARTING MAIN SIMULATION LOOP...");

        LOG_DEBUG("Disk evolution is ON (mode: %s). Starting main simulation loop.", snapshotModeToString(mode));
        LOG_DEBUG("Calling timeIntegrationForTheSystem...");
        timeIntegrationForTheSystem(mode,&disk_params, &sim_opts, &output_files);
        LOG_DEBUG("timeIntegrationForTheSystem completed. Program finished normally.");
    }

    if (disk_params.radial_grid) free(disk_params.radial_grid);
    if (disk_params.gas_surface_density_vector) free(disk_params.gas_surface_density_vector);
    if (disk_params.gas_pressure_vector) free(disk_params.gas_pressure_vector);
    if (disk_params.gas_pressure_gradient_vector) free(disk_params.gas_pressure_gradient_vector);
    if (disk_params.gas_velocity_vector) free(disk_params.gas_velocity_vector);
    if (disk_params.sigma_dot_photoevap) free(disk_params.sigma_dot_photoevap);
    LOG_DEBUG("Dynamically allocated disk arrays freed.");
    free(current_inputsig_file);
    free(current_inputdust_file);
    free(initial_dir_path);
    free(kLogFilesDirectory_path);
    LOG_DEBUG("Program exiting normally.");

    time(&wall_end);
    double elapsed_wall_time = difftime(wall_end, wall_start);

    printFinalSimulationSummary(actual_run_dir, elapsed_wall_time, &sim_opts);

    return 0;
}