#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// Header files
#include "config.h"           // Declarations of global variables and constants
#include "init_tool_module.h" // runInitialization and init_tool_options_t
#include "io_utils.h"         // Functions from io_utils.c
#include "disk_model.h"       // Functions from disk_model.c
#include "dust_physics.h"     // Functions from dust_physics.c
#include "simulation_core.h"  // Functions from simulation_core.c
#include "utils.h"            // Functions from utils.h
#include "gas_physics.h"
#include "boundary_conditions.h"
#include "integrator.h"


// NEW: Include your simulation_types.h and parser.h
#include "simulation_types.h"
#include "parser.h"           // Now includes function declarations for parsing

// Function declaration for default init_tool options, assuming it's in init_tool_module.h
extern void initializeDefaultOptions(init_tool_options_t *def);

// Global variable definition for particle_number (if not defined elsewhere)
// Ensure this is only defined in ONE .c file, and declared as 'extern int particle_number;' in config.h

int main(int argc, const char **argv) {
    // DEBUG: Program start
    // fprintf(stderr, "DEBUG [main]: Program started.\n");

    // --- Declare current_inputsig_file here, at the beginning of main ---
    char current_inputsig_file[MAX_PATH_LEN];
    current_inputsig_file[0] = '\0'; // Initialize to empty string

    // Local structure to store command-line options
    options_t def;
    createDefaultOptions(&def);

    // Local structure to store init_tool parameters (these will be populated from 'def')
    init_tool_options_t init_tool_params;
    initializeDefaultOptions(&init_tool_params);

    // Parse command-line options and populate the 'def' structure
    int retCode = parseCLIOptions(argc, argv, &def);
    if (0 != retCode) {
        fprintf(stderr, "Error parsing command-line options. Exiting with code %d.\n", retCode);
        return retCode;
    }

    // --- Declare instances of the new simulation structs ---
    DiskParameters disk_params; // Main disk parameters struct
    SimulationOptions sim_opts;
    OutputFiles output_files;

    // Initialize output_files pointers to NULL
    output_files.por_motion_file = NULL;
    output_files.micron_motion_file = NULL;
    output_files.mass_file = NULL;
    output_files.surface_file = NULL;
    output_files.dust_file = NULL;
    output_files.micron_dust_file = NULL;
    output_files.size_file = NULL;

    /* Populate the SimulationOptions struct from 'def' (parsed options) */
    sim_opts.evol = def.evol;
    sim_opts.drift = def.drift;
    sim_opts.growth = def.growth;
    sim_opts.twopop = def.twopop;
    sim_opts.DT = def.tStep;
    sim_opts.TMAX = def.totalTime;
    sim_opts.WO = def.outputFrequency;
    sim_opts.TCURR = def.startTime; // Initial current time
    sim_opts.num_dust_particles = def.ndust_val; // NEW: Populate num_dust_particles from parsed options

    // DEBUG: Show def.output_dir_name BEFORE it's used to populate sim_opts.output_dir_name
    fprintf(stderr, "DEBUG [main]: def.output_dir_name BEFORE sim_opts population: '%s'\n", def.output_dir_name);

    fprintf(stderr, "DEBUG [main]: Evolution (sim_opts.evol=%.2f) or drift (sim_opts.drift=%.2f) is ON. Starting main simulation loop.\n", sim_opts.evol, sim_opts.drift);

    // --- Populate DiskParameters with parameters from 'def' ---
    disk_params.r_min = def.rmin_val;
    disk_params.r_max = def.rmax_val;
    disk_params.grid_number = def.ngrid_val; // grid_number (gas grid points) is from parsed options
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
    disk_params.f_frag = def.ffrag;
    disk_params.fragmentation_velocity = def.ufrag;
    disk_params.f_drift = 0.55; // set by Birnstiel 2012
    disk_params.particle_density = def.pdensity_val;

    // Set sim_opts->dzone based on dead zone radii from disk_params
    sim_opts.dzone = (disk_params.r_dze_i > 0.0 || disk_params.r_dze_o > 0.0) ? 1.0 : 0.0;

    // --- Output directory handling ---
    // createRunDirectory contains the numbering logic if the directory already exists.
    createRunDirectory(def.output_dir_name); // Creates the main output folder, potentially with a number.
    fprintf(stderr, "DEBUG [main]: After createRunDirectory (base dir), def.output_dir_name is now: '%s'\n", def.output_dir_name);

    char initial_dir_path[MAX_PATH_LEN];
    char kLogFilesDirectory_path[MAX_PATH_LEN];

    // Create the 'intial' subdirectory using kConfigFilesDirectory
    snprintf(initial_dir_path, sizeof(initial_dir_path), "%s/%s", def.output_dir_name, kConfigFilesDirectory);
    createRunDirectory(initial_dir_path);

    // Create the 'LOGS' subdirectory using kLogFilesDirectory
    snprintf(kLogFilesDirectory_path, sizeof(kLogFilesDirectory_path), "%s/%s", def.output_dir_name, kLogFilesDirectory);
    fprintf(stderr, "DEBUG [main]: kLogFilesDirectory_path assembled as: '%s'\n", kLogFilesDirectory_path);

    createRunDirectory(kLogFilesDirectory_path);

    // CRITICAL: Populate sim_opts.output_dir_name from def.output_dir_name
    strncpy(sim_opts.output_dir_name, def.output_dir_name, MAX_PATH_LEN - 1);
    sim_opts.output_dir_name[MAX_PATH_LEN - 1] = '\0'; // Ensure null-termination
    // DEBUG: Show sim_opts.output_dir_name immediately after population
    fprintf(stderr, "DEBUG [main]: sim_opts.output_dir_name AFTER population: '%s'\n", sim_opts.output_dir_name);

    fprintf(stderr, "Output subdirectories created.\n");

    int dummy_sys_ret; // Dummy variable for system() call return value
    char cmd_buffer[MAX_PATH_LEN * 2]; // Buffer for system commands

    // --- Input file handling logic ---
    // This is where it's decided whether to read from a file or generate a new profile.
    // It's crucial that disk_params members are allocated memory ONLY ONCE.
    if (def.input_file != NULL && strcmp(def.input_file, "") != 0) {
        // If an input file is specified:
        strncpy(current_inputsig_file, def.input_file, MAX_PATH_LEN - 1);
        current_inputsig_file[MAX_PATH_LEN - 1] = '\0';
        fprintf(stderr, "DEBUG [main]: Input file specified: '%s'. Attempting to read initial profile.\n", current_inputsig_file);

        // disk_params.grid_number update from file *before* allocation.
        disk_params.grid_number = calculateNumbersOfParticles(current_inputsig_file); // Update grid_number from file (for GAS grid)

        // Recalculate delta_r based on the updated grid_number
        if (disk_params.grid_number > 1) {
            disk_params.delta_r = (disk_params.r_max - disk_params.r_min) / (disk_params.grid_number - 1.0);
        } else {
            disk_params.delta_r = 0.0;
        }
        fprintf(stderr, "DEBUG [main]: grid_number set from input file: %d. delta_r calculated as %.4e.\n", disk_params.grid_number, disk_params.delta_r);

        // --- Dynamic Memory Allocation for Disk Arrays ---
        // This happens ONLY HERE, because runInitialization is not called!
        disk_params.radial_grid = (double *)malloc((disk_params.grid_number + 2) * sizeof(double));
        disk_params.gas_surface_density_vector = (double *)malloc((disk_params.grid_number + 2) * sizeof(double));
        disk_params.gas_pressure_vector = (double *)malloc((disk_params.grid_number + 2) * sizeof(double));
        disk_params.gas_pressure_gradient_vector = (double *)malloc((disk_params.grid_number + 2) * sizeof(double));
        disk_params.gas_velocity_vector = (double *)malloc((disk_params.grid_number + 2) * sizeof(double));

        if (!disk_params.radial_grid || !disk_params.gas_surface_density_vector || !disk_params.gas_pressure_vector || !disk_params.gas_pressure_gradient_vector || !disk_params.gas_velocity_vector) {
            fprintf(stderr, "ERROR [main]: Failed to allocate memory for disk arrays (input file branch). Exiting.\n");
            return 1;
        }
        fprintf(stderr, "DEBUG [main]: Disk profile arrays dynamically allocated with size grid_number+2 = %d (input file branch).\n", disk_params.grid_number + 2);

        // Call readDiskParameters because runInitialization is not called
        fprintf(stderr, "DEBUG [main]: Calling readDiskParameters to calculate derived disk parameters for main disk_params struct (input file branch).\n");
        readDiskParameters(&disk_params);
        fprintf(stderr, "DEBUG [main]: readDiskParameters completed (input file branch).\n");

        // Copy input profile file to the 'initial' directory
        snprintf(cmd_buffer, sizeof(cmd_buffer), "cp %s %s/", current_inputsig_file, initial_dir_path);
        dummy_sys_ret = system(cmd_buffer); (void)dummy_sys_ret;
        fprintf(stderr, "DEBUG [main]: Copied initial profile file '%s' to '%s/'.\n", current_inputsig_file, initial_dir_path);

        // Copy disk_config.dat (kDiskConfigFile) file to the 'initial' directory
        snprintf(cmd_buffer, sizeof(cmd_buffer), "cp %s%s %s/", kDiskConfigFile,kFileNamesSuffix, initial_dir_path);
        dummy_sys_ret = system(cmd_buffer); (void)dummy_sys_ret;
        fprintf(stderr, "DEBUG [main]: Copied %s to %s%s/\n", kDiskConfigFile,kFileNamesSuffix, initial_dir_path);

    } else {
        // If NO input file is specified:
        fprintf(stderr, "DEBUG [main]: No input file specified (-i flag not used). Generating default grid and profile.\n");

        // Populate init_tool_options_t from 'def' (command-line) values
        init_tool_params.n_grid_points = disk_params.grid_number; // This is the gas grid resolution
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
        init_tool_params.n_dust_particles = def.ndust_val; // Uses the specific dust particle count
        init_tool_params.two_pop_ratio = def.ratio_val;
        init_tool_params.micro_size_cm = def.mic_val;
        init_tool_params.one_size_particle_cm = def.onesize_val;
        init_tool_params.dust_density_g_cm3 = def.pdensity_val;
        fprintf(stderr, "DEBUG [main]: init_tool_options_t (init_tool_params) structure populated for profile generation.\n");

        // --- Generate profile directly into the 'initial' directory ---
        fprintf(stderr, "DEBUG [main]: Calling runInitialization(&init_tool_params, &disk_params)...\n");
        strncpy(init_tool_params.output_base_path, initial_dir_path, MAX_PATH_LEN - 1);
        init_tool_params.output_base_path[MAX_PATH_LEN - 1] = '\0';
        // runInitialization is responsible for allocating and populating disk_params members
        runInitialization(&init_tool_params, &disk_params);
        fprintf(stderr, "DEBUG [main]: runInitialization completed. disk_params allocated and populated.\n");

        // Now current_inputsig_file points to the generated file in initial_dir_path
        // CHANGE HERE: FILENAME_INIT_PROFILE -> kInitialGasProfileFileName
        snprintf(current_inputsig_file, sizeof(current_inputsig_file), "%s/%s%s", initial_dir_path, kInitialGasProfileFileName,kFileNamesSuffix);
        fprintf(stderr, "DEBUG [main]: Generated GAS profile will be loaded from '%s'.\n", current_inputsig_file);

        // --- Update grid_number from the generated file (critical for loadGasSurfaceDensityFromFile sizing) ---
        // Important: Here, the number of lines should be read from kInitialGasProfileFileName,
        // provided that init_tool_module.c creates this file.
        disk_params.grid_number = calculateNumbersOfParticles(current_inputsig_file);

        if (disk_params.grid_number > 1) {
            disk_params.delta_r = (disk_params.r_max - disk_params.r_min) / (disk_params.grid_number - 1.0);
        } else {
            disk_params.delta_r = 0.0;
        }
        fprintf(stderr, "DEBUG [main]: grid_number updated from generated file: %d. delta_r calculated as %.4e.\n", disk_params.grid_number, disk_params.delta_r);

        // No need for 'cp' here for kDiskConfigFile or FILENAME_INIT_PROFILE,
        // since runInitialization created them directly in initial_dir_path.
    }

    // --- CRITICAL STEP: Now that current_inputsig_file is set,
    //     copy it to sim_opts.input_filename for timeIntegrationForTheSystem. ---
    // This is the gas profile filename that loadGasSurfaceDensityFromFile reads.
    strncpy(sim_opts.input_filename, current_inputsig_file, MAX_PATH_LEN - 1);
    sim_opts.input_filename[MAX_PATH_LEN - 1] = '\0'; // Ensure null-termination
    fprintf(stderr, "DEBUG [main]: sim_opts.input_filename set to '%s' for timeIntegrationForTheSystem (gas profile).\n", sim_opts.input_filename);

    // --- NEW PART: Set the dust profile filename in sim_opts.dust_input_filename ---
    // This is the dust profile filename that loadDustParticlesFromFile reads within timeIntegrationForTheSystem.
    char current_inputdust_file[MAX_PATH_LEN];
    snprintf(current_inputdust_file, sizeof(current_inputdust_file), "%s/%s%s", initial_dir_path, kInitialDustProfileFileName, kFileNamesSuffix);
    strncpy(sim_opts.dust_input_filename, current_inputdust_file, MAX_PATH_LEN - 1);
    sim_opts.dust_input_filename[MAX_PATH_LEN - 1] = '\0'; // Ensure null-termination
    fprintf(stderr, "DEBUG [main]: sim_opts.dust_input_filename set to '%s' for timeIntegrationForTheSystem (dust profile).\n", sim_opts.dust_input_filename);

    // --- CRITICAL STEP: Set the global particle_number based on the actual dust particle file. ---
    // This ensures particle_number reflects the *dust* particle count, not the gas grid count.
    particle_number = calculateNumbersOfParticles(sim_opts.dust_input_filename); // Read lines from the dust profile
    fprintf(stderr, "DEBUG [main]: Global particle_number set to %d (from dust input file: %s).\n", particle_number, sim_opts.dust_input_filename);

    // The readDiskParameters call has already occurred in the appropriate branch (if using input file)
    // or was called by runInitialization (if generating).

    fprintf(stderr, "DEBUG [main]: Initial profile loading for loadGasSurfaceDensityFromFile...\n");
    loadGasSurfaceDensityFromFile(&disk_params, current_inputsig_file); // This populates disk_params.gas_surface_density_vector and radial_grid
    fprintf(stderr, "DEBUG [main]: loadGasSurfaceDensityFromFile completed. Calling applyBoundaryConditions for disk_params.radial_grid and disk_params.gas_surface_density_vector...\n");
    applyBoundaryConditions(disk_params.radial_grid, &disk_params);
    applyBoundaryConditions(disk_params.gas_surface_density_vector, &disk_params);
    fprintf(stderr, "DEBUG [main]: applyBoundaryConditions calls completed for initial profile.\n");

    // Print current information
    fprintf(stderr, "DEBUG [main]: Calling printCurrentInformationAboutRun...\n");
    // Here, def.output_dir_name is passed, which already contains the numbered folder name
    printCurrentInformationAboutRun(def.output_dir_name, &disk_params, &sim_opts);

    // Run simulation or exit based on options
    if(sim_opts.evol == 0. && sim_opts.drift == 0.) {
        fprintf(stderr, "DEBUG [main]: Evolution (sim_opts.evol=%.2f) and drift (sim_opts.drift=%.2f) are OFF.\n", sim_opts.evol, sim_opts.drift);

        char dens_name_initial[MAX_PATH_LEN];
        snprintf(dens_name_initial, sizeof(dens_name_initial), "%s/%s%s", initial_dir_path, kInitialGasProfileFileName,kFileNamesSuffix);
        fprintf(stderr, "DEBUG [main]: Printing initial surface density to %s.\n", dens_name_initial);

        // Special handling for printGasSurfaceDensityPressurePressureDerivateFile when only initial output is needed
        OutputFiles temp_output_for_initial_print;
        temp_output_for_initial_print.surface_file = fopen(dens_name_initial, "w");
        if (temp_output_for_initial_print.surface_file != NULL) {
            // Add header to initial surface density file
            fprintf(temp_output_for_initial_print.surface_file, "# Initial Disk Profile Data (Gas Surface Density, Pressure, Gradient)\n");
            fprintf(temp_output_for_initial_print.surface_file, "# Columns: 1. Radius [AU], 2. Sigma_gas [M_Sun/AU^2],\n");
            fprintf(temp_output_for_initial_print.surface_file, "#   3. Pressure [units], 4. dP/dR [units]\n");
            fprintf(temp_output_for_initial_print.surface_file, "#\n");
            fprintf(temp_output_for_initial_print.surface_file, "# Data generated by Dust Drift Simulation (Initial State)\n");
            fflush(temp_output_for_initial_print.surface_file);

            printGasSurfaceDensityPressurePressureDerivateFile(&disk_params, &temp_output_for_initial_print);
            fclose(temp_output_for_initial_print.surface_file);
            temp_output_for_initial_print.surface_file = NULL;
            fprintf(stderr, "DEBUG [main]: Closed %s.\n", dens_name_initial);
        } else {
            fprintf(stderr, "ERROR [main]: Could not open %s for initial surface output.\n", dens_name_initial);
        }

        fprintf(stderr, "DEBUG [main]: printGasSurfaceDensityPressurePressureDerivateFile completed. Program exiting.\n");
    } else {
        fprintf(stderr, "DEBUG [main]: Evolution (sim_opts.evol=%.2f) or drift (sim_opts.drift=%.2f) is ON. Starting main simulation loop.\n", sim_opts.evol, sim_opts.drift);
        fprintf(stderr, "DEBUG [main]: Calling timeIntegrationForTheSystem...\n");
        // Pass sim_opts to timeIntegrationForTheSystem.
        // timeIntegrationForTheSystem must ensure to use the correct (numbered) output_dir_name.
        timeIntegrationForTheSystem(&disk_params, &sim_opts, &output_files);
        fprintf(stderr, "DEBUG [main]: timeIntegrationForTheSystem completed. Program finished normally.\n");
    }

    // --- Free dynamically allocated memory ---
    // Free the memory allocated for disk_params.radial_grid, gas_surface_density_vector, etc.
    // This allocation happened in the 'if' branch if an input file was used,
    // or within runInitialization in the 'else' branch.
    if (disk_params.radial_grid) free(disk_params.radial_grid);
    if (disk_params.gas_surface_density_vector) free(disk_params.gas_surface_density_vector);
    if (disk_params.gas_pressure_vector) free(disk_params.gas_pressure_vector);
    if (disk_params.gas_pressure_gradient_vector) free(disk_params.gas_pressure_gradient_vector);
    if (disk_params.gas_velocity_vector) free(disk_params.gas_velocity_vector);
    fprintf(stderr, "DEBUG [main]: Dynamically allocated disk arrays freed.\n");

    fprintf(stderr, "DEBUG [main]: Program exiting normally.\n");
    return 0;
}
