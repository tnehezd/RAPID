#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// Header files
#include "config.h"           // Declarations of global variables and constants
#include "init_tool_module.h" // run_init_tool and init_tool_options_t
#include "io_utils.h"         // Functions from io_utils.c
#include "disk_model.h"       // Functions from disk_model.c
#include "dust_physics.h"     // Functions from dust_physics.c
#include "simulation_core.h"  // Functions from simulation_core.c
#include "utils.h"            // Functions from utils.h

// NEW: Include your simulation_types.h and parser.h
#include "simulation_types.h"
#include "parser.h"           // Now includes function declarations for parsing

// Function declaration for default init_tool options, assuming it's in init_tool_module.h
extern void create_default_init_tool_options(init_tool_options_t *def);

// Global variable definition for PARTICLE_NUMBER (if not defined elsewhere)
// Ensure this is only defined in ONE .c file, and declared as 'extern int PARTICLE_NUMBER;' in config.h

int main(int argc, const char **argv) {
    // DEBUG: Program start
    // fprintf(stderr, "DEBUG [main]: Program started.\n");

    // --- Declare current_inputsig_file here, at the beginning of main ---
    char current_inputsig_file[MAX_PATH_LEN];
    current_inputsig_file[0] = '\0'; // Initialize to empty string

    // Local structure to store command-line options
    options_t def;
    create_default_options(&def);

    // Local structure to store init_tool parameters (these will be populated from 'def')
    init_tool_options_t init_tool_params;
    create_default_init_tool_options(&init_tool_params);

    // Parse command-line options and populate the 'def' structure
    int retCode = parse_options(argc, argv, &def);
    if (0 != retCode) {
        fprintf(stderr, "Error parsing command-line options. Exiting with code %d.\n", retCode);
        return retCode;
    }

    // --- Declare instances of the new simulation structs ---
    disk_t disk_params; // Main disk parameters struct
    simulation_options_t sim_opts;
    output_files_t output_files;

    // Initialize output_files pointers to NULL
    output_files.por_motion_file = NULL;
    output_files.micron_motion_file = NULL;
    output_files.mass_file = NULL;
    output_files.surface_file = NULL;
    output_files.dust_file = NULL;
    output_files.micron_dust_file = NULL;
    output_files.size_file = NULL;

    /* Populate the simulation_options_t struct from 'def' (parsed options) */
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

    // --- Populate disk_t with parameters from 'def' ---
    disk_params.RMIN = def.rmin_val;
    disk_params.RMAX = def.rmax_val;
    disk_params.NGRID = def.ngrid_val; // NGRID (gas grid points) is from parsed options
    disk_params.SIGMA0 = def.sigma0_val;
    disk_params.SIGMAP_EXP = def.sigmap_exp_val;
    disk_params.alpha_visc = def.alpha_visc_val;
    disk_params.STAR_MASS = def.star_val;
    disk_params.HASP = def.hasp_val;
    disk_params.FLIND = def.flind_val;
    disk_params.r_dze_i = def.r_dze_i_val;
    disk_params.r_dze_o = def.r_dze_o_val;
    disk_params.Dr_dze_i = def.dr_dze_i_val;
    disk_params.Dr_dze_o = def.dr_dze_o_val;
    disk_params.a_mod = def.a_mod_val;
    disk_params.fFrag = def.ffrag;
    disk_params.uFrag = def.ufrag;
    disk_params.fDrift = 0.55; // set by Birnstiel 2012
    disk_params.PDENSITY = def.pdensity_val;

    // Set sim_opts->dzone based on dead zone radii from disk_params
    sim_opts.dzone = (disk_params.r_dze_i > 0.0 || disk_params.r_dze_o > 0.0) ? 1.0 : 0.0;

    // --- Output directory handling ---
    // Mk_Dir contains the numbering logic if the directory already exists.
    Mk_Dir(def.output_dir_name); // Creates the main output folder, potentially with a number.
    fprintf(stderr, "DEBUG [main]: After Mk_Dir (base dir), def.output_dir_name is now: '%s'\n", def.output_dir_name);

    char initial_dir_path[MAX_PATH_LEN];
    char logs_dir_path[MAX_PATH_LEN];

    // Create the 'initial' subdirectory using CONFIG_DIR
    snprintf(initial_dir_path, sizeof(initial_dir_path), "%s/%s", def.output_dir_name, CONFIG_DIR);
    Mk_Dir(initial_dir_path);

    // Create the 'LOGS' subdirectory using LOGS_DIR
    snprintf(logs_dir_path, sizeof(logs_dir_path), "%s/%s", def.output_dir_name, LOGS_DIR);
    fprintf(stderr, "DEBUG [main]: logs_dir_path assembled as: '%s'\n", logs_dir_path);

    Mk_Dir(logs_dir_path);

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

        // disk_params.NGRID update from file *before* allocation.
        disk_params.NGRID = reszecskek_szama(current_inputsig_file); // Update NGRID from file (for GAS grid)

        // Recalculate DD based on the updated NGRID
        if (disk_params.NGRID > 1) {
            disk_params.DD = (disk_params.RMAX - disk_params.RMIN) / (disk_params.NGRID - 1.0);
        } else {
            disk_params.DD = 0.0;
        }
        fprintf(stderr, "DEBUG [main]: NGRID set from input file: %d. DD calculated as %.4e.\n", disk_params.NGRID, disk_params.DD);

        // --- Dynamic Memory Allocation for Disk Arrays ---
        // This happens ONLY HERE, because run_init_tool is not called!
        disk_params.rvec = (double *)malloc((disk_params.NGRID + 2) * sizeof(double));
        disk_params.sigmavec = (double *)malloc((disk_params.NGRID + 2) * sizeof(double));
        disk_params.pressvec = (double *)malloc((disk_params.NGRID + 2) * sizeof(double));
        disk_params.dpressvec = (double *)malloc((disk_params.NGRID + 2) * sizeof(double));
        disk_params.ugvec = (double *)malloc((disk_params.NGRID + 2) * sizeof(double));

        if (!disk_params.rvec || !disk_params.sigmavec || !disk_params.pressvec || !disk_params.dpressvec || !disk_params.ugvec) {
            fprintf(stderr, "ERROR [main]: Failed to allocate memory for disk arrays (input file branch). Exiting.\n");
            return 1;
        }
        fprintf(stderr, "DEBUG [main]: Disk profile arrays dynamically allocated with size NGRID+2 = %d (input file branch).\n", disk_params.NGRID + 2);

        // Call disk_param_be because run_init_tool is not called
        fprintf(stderr, "DEBUG [main]: Calling disk_param_be to calculate derived disk parameters for main disk_params struct (input file branch).\n");
        disk_param_be(&disk_params);
        fprintf(stderr, "DEBUG [main]: disk_param_be completed (input file branch).\n");

        // Copy input profile file to the 'initial' directory
        snprintf(cmd_buffer, sizeof(cmd_buffer), "cp %s %s/", current_inputsig_file, initial_dir_path);
        dummy_sys_ret = system(cmd_buffer); (void)dummy_sys_ret;
        fprintf(stderr, "DEBUG [main]: Copied initial profile file '%s' to '%s/'.\n", current_inputsig_file, initial_dir_path);

        // Copy disk_config.dat (FILENAME_DISK_PARAM) file to the 'initial' directory
        snprintf(cmd_buffer, sizeof(cmd_buffer), "cp %s %s/", FILENAME_DISK_PARAM, initial_dir_path);
        dummy_sys_ret = system(cmd_buffer); (void)dummy_sys_ret;
        fprintf(stderr, "DEBUG [main]: Copied %s to %s/\n", FILENAME_DISK_PARAM, initial_dir_path);

    } else {
        // If NO input file is specified:
        fprintf(stderr, "DEBUG [main]: No input file specified (-i flag not used). Generating default grid and profile.\n");

        // Populate init_tool_options_t from 'def' (command-line) values
        init_tool_params.n_grid_points = disk_params.NGRID; // This is the gas grid resolution
        init_tool_params.r_inner= disk_params.RMIN;
        init_tool_params.r_outer = disk_params.RMAX;
        init_tool_params.sigma0_gas_au = disk_params.SIGMA0;
        init_tool_params.sigma_exponent = disk_params.SIGMAP_EXP;
        init_tool_params.deadzone_r_inner = disk_params.r_dze_i;
        init_tool_params.deadzone_r_outer = disk_params.r_dze_o;
        init_tool_params.deadzone_dr_inner = disk_params.Dr_dze_i;
        init_tool_params.deadzone_dr_outer = disk_params.Dr_dze_o;
        init_tool_params.alpha_viscosity = disk_params.alpha_visc;
        init_tool_params.deadzone_alpha_mod = disk_params.a_mod;
        init_tool_params.aspect_ratio = disk_params.HASP;
        init_tool_params.flaring_index = disk_params.FLIND;
        init_tool_params.star_mass = disk_params.STAR_MASS;
        init_tool_params.dust_to_gas_ratio = def.eps_val;
        init_tool_params.n_dust_particles = def.ndust_val; // Uses the specific dust particle count
        init_tool_params.two_pop_ratio = def.ratio_val;
        init_tool_params.micro_size_cm = def.mic_val;
        init_tool_params.one_size_particle_cm = def.onesize_val;
        init_tool_params.dust_density_g_cm3 = def.pdensity_val;
        fprintf(stderr, "DEBUG [main]: init_tool_options_t (init_tool_params) structure populated for profile generation.\n");

        // --- Generate profile directly into the 'initial' directory ---
        fprintf(stderr, "DEBUG [main]: Calling run_init_tool(&init_tool_params, &disk_params)...\n");
        strncpy(init_tool_params.output_base_path, initial_dir_path, MAX_PATH_LEN - 1);
        init_tool_params.output_base_path[MAX_PATH_LEN - 1] = '\0';
        // run_init_tool is responsible for allocating and populating disk_params members
        run_init_tool(&init_tool_params, &disk_params);
        fprintf(stderr, "DEBUG [main]: run_init_tool completed. disk_params allocated and populated.\n");

        // Now current_inputsig_file points to the generated file in initial_dir_path
        // CHANGE HERE: FILENAME_INIT_PROFILE -> FILENAME_INIT_GAS_PROFILE
        snprintf(current_inputsig_file, sizeof(current_inputsig_file), "%s/%s", initial_dir_path, FILENAME_INIT_GAS_PROFILE);
        fprintf(stderr, "DEBUG [main]: Generated GAS profile will be loaded from '%s'.\n", current_inputsig_file);

        // --- Update NGRID from the generated file (critical for sigIn sizing) ---
        // Important: Here, the number of lines should be read from FILENAME_INIT_GAS_PROFILE,
        // provided that init_tool_module.c creates this file.
        disk_params.NGRID = reszecskek_szama(current_inputsig_file);

        if (disk_params.NGRID > 1) {
            disk_params.DD = (disk_params.RMAX - disk_params.RMIN) / (disk_params.NGRID - 1.0);
        } else {
            disk_params.DD = 0.0;
        }
        fprintf(stderr, "DEBUG [main]: NGRID updated from generated file: %d. DD calculated as %.4e.\n", disk_params.NGRID, disk_params.DD);

        // No need for 'cp' here for FILENAME_DISK_PARAM or FILENAME_INIT_PROFILE,
        // since run_init_tool created them directly in initial_dir_path.
    }

    // --- CRITICAL STEP: Now that current_inputsig_file is set,
    //     copy it to sim_opts.input_filename for tIntegrate. ---
    // This is the gas profile filename that sigIn reads.
    strncpy(sim_opts.input_filename, current_inputsig_file, MAX_PATH_LEN - 1);
    sim_opts.input_filename[MAX_PATH_LEN - 1] = '\0'; // Ensure null-termination
    fprintf(stderr, "DEBUG [main]: sim_opts.input_filename set to '%s' for tIntegrate (gas profile).\n", sim_opts.input_filename);

    // --- NEW PART: Set the dust profile filename in sim_opts.dust_input_filename ---
    // This is the dust profile filename that por_be reads within tIntegrate.
    char current_inputdust_file[MAX_PATH_LEN];
    snprintf(current_inputdust_file, sizeof(current_inputdust_file), "%s/%s", initial_dir_path, FILENAME_INIT_DUST_PROFILE);
    strncpy(sim_opts.dust_input_filename, current_inputdust_file, MAX_PATH_LEN - 1);
    sim_opts.dust_input_filename[MAX_PATH_LEN - 1] = '\0'; // Ensure null-termination
    fprintf(stderr, "DEBUG [main]: sim_opts.dust_input_filename set to '%s' for tIntegrate (dust profile).\n", sim_opts.dust_input_filename);

    // --- CRITICAL STEP: Set the global PARTICLE_NUMBER based on the actual dust particle file. ---
    // This ensures PARTICLE_NUMBER reflects the *dust* particle count, not the gas grid count.
    PARTICLE_NUMBER = reszecskek_szama(sim_opts.dust_input_filename); // Read lines from the dust profile
    fprintf(stderr, "DEBUG [main]: Global PARTICLE_NUMBER set to %d (from dust input file: %s).\n", PARTICLE_NUMBER, sim_opts.dust_input_filename);

    // The disk_param_be call has already occurred in the appropriate branch (if using input file)
    // or was called by run_init_tool (if generating).

    fprintf(stderr, "DEBUG [main]: Initial profile loading for sigIn...\n");
    sigIn(&disk_params, current_inputsig_file); // This populates disk_params.sigmavec and rvec
    fprintf(stderr, "DEBUG [main]: sigIn completed. Calling Perem for disk_params.rvec and disk_params.sigmavec...\n");
    Perem(disk_params.rvec, &disk_params);
    Perem(disk_params.sigmavec, &disk_params);
    fprintf(stderr, "DEBUG [main]: Perem calls completed for initial profile.\n");

    // Print current information
    fprintf(stderr, "DEBUG [main]: Calling infoCurrent...\n");
    // Here, def.output_dir_name is passed, which already contains the numbered folder name
    infoCurrent(def.output_dir_name, &disk_params, &sim_opts);

    // Run simulation or exit based on options
    if(sim_opts.evol == 0. && sim_opts.drift == 0.) {
        fprintf(stderr, "DEBUG [main]: Evolution (sim_opts.evol=%.2f) and drift (sim_opts.drift=%.2f) are OFF.\n", sim_opts.evol, sim_opts.drift);

        char dens_name_initial[MAX_PATH_LEN];
        snprintf(dens_name_initial, sizeof(dens_name_initial), "%s/%s", initial_dir_path, FILENAME_INIT_GAS_PROFILE);
        fprintf(stderr, "DEBUG [main]: Printing initial surface density to %s.\n", dens_name_initial);

        // Special handling for Print_Sigma when only initial output is needed
        output_files_t temp_output_for_initial_print;
        temp_output_for_initial_print.surface_file = fopen(dens_name_initial, "w");
        if (temp_output_for_initial_print.surface_file != NULL) {
            // Add header to initial surface density file
            fprintf(temp_output_for_initial_print.surface_file, "# Initial Disk Profile Data (Gas Surface Density, Pressure, Gradient)\n");
            fprintf(temp_output_for_initial_print.surface_file, "# Columns: 1. Radius [AU], 2. Sigma_gas [M_Sun/AU^2],\n");
            fprintf(temp_output_for_initial_print.surface_file, "#   3. Pressure [units], 4. dP/dR [units]\n");
            fprintf(temp_output_for_initial_print.surface_file, "#\n");
            fprintf(temp_output_for_initial_print.surface_file, "# Data generated by Dust Drift Simulation (Initial State)\n");
            fflush(temp_output_for_initial_print.surface_file);

            Print_Sigma(&disk_params, &temp_output_for_initial_print);
            fclose(temp_output_for_initial_print.surface_file);
            temp_output_for_initial_print.surface_file = NULL;
            fprintf(stderr, "DEBUG [main]: Closed %s.\n", dens_name_initial);
        } else {
            fprintf(stderr, "ERROR [main]: Could not open %s for initial surface output.\n", dens_name_initial);
        }

        fprintf(stderr, "DEBUG [main]: Print_Sigma completed. Program exiting.\n");
    } else {
        fprintf(stderr, "DEBUG [main]: Evolution (sim_opts.evol=%.2f) or drift (sim_opts.drift=%.2f) is ON. Starting main simulation loop.\n", sim_opts.evol, sim_opts.drift);
        fprintf(stderr, "DEBUG [main]: Calling tIntegrate...\n");
        // Pass sim_opts to tIntegrate.
        // tIntegrate must ensure to use the correct (numbered) output_dir_name.
        tIntegrate(&disk_params, &sim_opts, &output_files);
        fprintf(stderr, "DEBUG [main]: tIntegrate completed. Program finished normally.\n");
    }

    // --- Free dynamically allocated memory ---
    // Free the memory allocated for disk_params.rvec, sigmavec, etc.
    // This allocation happened in the 'if' branch if an input file was used,
    // or within run_init_tool in the 'else' branch.
    if (disk_params.rvec) free(disk_params.rvec);
    if (disk_params.sigmavec) free(disk_params.sigmavec);
    if (disk_params.pressvec) free(disk_params.pressvec);
    if (disk_params.dpressvec) free(disk_params.dpressvec);
    if (disk_params.ugvec) free(disk_params.ugvec);
    fprintf(stderr, "DEBUG [main]: Dynamically allocated disk arrays freed.\n");

    fprintf(stderr, "DEBUG [main]: Program exiting normally.\n");
    return 0;
}
