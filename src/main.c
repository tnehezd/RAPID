#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> // For M_PI, if not defined in config.h
#include <time.h> // Not directly used in main now, but good to keep if needed

// Header files
#include "config.h"           // Declarations of global variables and constants
#include "io_utils.h"         // Functions from io_utils.c (needs refactoring for structs)
#include "disk_model.h"       // Functions from disk_model.c (ALREADY refactored for structs)
#include "dust_physics.h"     // Functions from dust_physics.c (needs refactoring for structs)
#include "simulation_core.h"  // Functions from simulation_core.c (needs refactoring for structs)
#include "utils.h"            // Functions from utils.c (needs refactoring for structs)
#include "init_tool_module.h" // run_init_tool and init_tool_options_t (ALREADY refactored for structs)

// --- NEW: Include your simulation_types.h and parser.h ---
#include "simulation_types.h"
#include "parser.h" // Now includes function declarations for parsing

// Function declaration for default init_tool options, assuming it's in init_tool_module.h
extern void create_default_init_tool_options(init_tool_options_t *def);

// Assuming FILENAME_INIT_PROFILE and FILENAME_DISK_PARAM are external from config.c/h
// (as defined in your config.c)

int main(int argc, const char **argv) {
    // DEBUG: Program start
    fprintf(stderr, "DEBUG [main]: Program started.\n");

    // --- Declare current_inputsig_file here, at the beginning of main ---
    const char *current_inputsig_file = NULL;

    // Local structure to store command-line options
    options_t def;
    // Set default values
    create_default_options(&def);
    fprintf(stderr, "DEBUG [main]: default_options.output_dir_name = '%s'\n", def.output_dir_name);

    // Local structure to store init_tool parameters (these will be populated from 'def')
    init_tool_options_t init_tool_params;
    // Set default values for init_tool. This is useful if 'def' doesn't override all of them.
    create_default_init_tool_options(&init_tool_params);
    fprintf(stderr, "DEBUG [main]: Default init_tool options created.\n");

    // Parse command-line options and populate the 'def' structure
    int retCode = parse_options(argc, argv, &def);

    if (0 != retCode) {
        // Exit on error (usage printed by parse_options)
        fprintf(stderr, "DEBUG [main]: Error parsing command-line options. Exiting with code %d.\n", retCode);
        return retCode;
    }
    fprintf(stderr, "DEBUG [main]: Command-line options parsed successfully.\n");
    // DEBUG: After parsing, check output directory name from def
    fprintf(stderr, "DEBUG [main]: after_parse_options.output_dir_name = '%s'\n", def.output_dir_name);


    // --- Declare instances of the new simulation structs ---
    disk_t disk_params;
    simulation_options_t sim_opts;
    output_files_t output_files;

    // Initialize output_files pointers to NULL
    output_files.por_motion_file = NULL;
    output_files.micron_motion_file = NULL;
    output_files.mass_file = NULL;
    output_files.surface_file = NULL;
    output_files.dust_file = NULL;
    output_files.micron_dust_file = NULL;
    output_files.size_file = NULL; // Initialize this one too
    fprintf(stderr, "DEBUG [main]: New simulation structs declared and output_files_t initialized to NULL.\n");

    /* Populate the simulation_options_t struct from 'def' (parsed options) */
    sim_opts.evol = def.evol;
    sim_opts.drift = def.drift;
    sim_opts.growth = def.growth;
    sim_opts.twopop = def.twopop;
    sim_opts.DT = def.tStep;
    sim_opts.TMAX = def.totalTime;
    sim_opts.WO = def.outputFrequency;
    sim_opts.TCURR = def.startTime; // Initial current time

    // --- CRITICAL: Populate sim_opts.output_dir_name from def.output_dir_name ---
    strncpy(sim_opts.output_dir_name, def.output_dir_name, MAX_PATH_LEN - 1);
    sim_opts.output_dir_name[MAX_PATH_LEN - 1] = '\0'; // Ensure null-termination
    fprintf(stderr, "DEBUG [main]: sim_opts.output_dir_name set to '%s' from def.output_dir_name.\n", sim_opts.output_dir_name);
    // --- END CRITICAL ADDITION ---

    fprintf(stderr, "DEBUG [main]: sim_opts structure populated (including output_dir_name).\n");
    fprintf(stderr, "DEBUG [main]:   evol=%.2f, drift=%.2f, growth=%.2f, twopop=%.2f, \n",
           sim_opts.evol, sim_opts.drift, sim_opts.growth, sim_opts.twopop);
    fprintf(stderr, "DEBUG [main]:   Time parameters set: TMAX=%.2e, WO=%.2e, TCURR=%.2e, DT=%.2e\n",
           sim_opts.TMAX, sim_opts.WO, sim_opts.TCURR, sim_opts.DT);


    // --- Populate disk_t with parameters from 'def' ---
    disk_params.RMIN = def.rmin_val;
    disk_params.RMAX = def.rmax_val;
    disk_params.NGRID = def.ngrid_val; // Initial NGRID, might be updated by input file reading
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
    disk_params.fDrift = 0.55;       // set by Birnstiel 2012

    fprintf(stderr, "DEBUG [main]: disk_params structure populated with initial parameters.\n");
    fprintf(stderr, "DEBUG [main]:   RMIN=%.2f, RMAX=%.2f, NGRID=%d, SIGMA0=%.2e, SIGMAP_EXP=%.2f\n",
           disk_params.RMIN, disk_params.RMAX, disk_params.NGRID, disk_params.SIGMA0, disk_params.SIGMAP_EXP);
    fprintf(stderr, "DEBUG [main]:   alpha_visc=%.2e, STAR_MASS=%.2f, HASP=%.2f, FLIND=%.2f\n",
           disk_params.alpha_visc, disk_params.STAR_MASS, disk_params.HASP, disk_params.FLIND);
    fprintf(stderr, "DEBUG [main]:   r_dze_i=%.2f, r_dze_o=%.2f, Dr_dze_i=%.2f, Dr_dze_o=%.2f, a_mod=%.2f, fFrag=%.2f, uFrag=%.2f\n",
           disk_params.r_dze_i, disk_params.r_dze_o, disk_params.Dr_dze_i, disk_params.Dr_dze_o, disk_params.a_mod, disk_params.fFrag, disk_params.uFrag);


    // Set sim_opts->dzone based on dead zone radii from disk_params
    sim_opts.dzone = (disk_params.r_dze_i > 0.0 || disk_params.r_dze_o > 0.0) ? 1.0 : 0.0;

    fprintf(stderr, "DEBUG [main]: Dead zone parameter (sim_opts.dzone) set.\n");


    // --- Input file handling logic ---
    if (def.input_file != NULL && strcmp(def.input_file, "") != 0) {
        current_inputsig_file = def.input_file;
        fprintf(stderr, "DEBUG [main]: Input file specified: '%s'. Attempting to read initial profile.\n", current_inputsig_file);

        // disk_params.NGRID needs to be updated by the file content *before* allocation.
        // reszecskek_szama takes the filename to count particles.
        disk_params.NGRID = reszecskek_szama(current_inputsig_file); // Update NGRID from file

        // Re-calculate DD based on updated NGRID
        if (disk_params.NGRID > 1) { // Avoid division by zero if NGRID is 1
            disk_params.DD = (disk_params.RMAX - disk_params.RMIN) / (disk_params.NGRID - 1.0);
        } else {
            disk_params.DD = 0.0; // Or handle as an error if NGRID <= 1 is invalid for your simulation
        }

        fprintf(stderr, "DEBUG [main]: NGRID set from input file: %d. DD calculated as %.4e.\n", disk_params.NGRID, disk_params.DD);
    } else {
        fprintf(stderr, "DEBUG [main]: No input file specified (-i flag not used). Generating default grid and profile.\n");
        // NGRID is already set from def.ngrid_val

        // --- Generate profile and set current_inputsig_file to init_data.dat ---
        fprintf(stderr, "DEBUG [main]: Calling run_init_tool(&init_tool_params)...\n");
        run_init_tool(&init_tool_params);
        fprintf(stderr, "DEBUG [main]: run_init_tool completed.\n");

        current_inputsig_file = "init_data.dat"; // This is the default output of init_tool
        fprintf(stderr, "DEBUG [main]: Generated profile will be loaded from '%s'.\n", current_inputsig_file);

        // --- Update NGRID from the generated file (crucial for accurate size for sigIn) ---
        disk_params.NGRID = reszecskek_szama(current_inputsig_file);
        // --- END ADDITION ---

        if (disk_params.NGRID > 1) {
            disk_params.DD = (disk_params.RMAX - disk_params.RMIN) / (disk_params.NGRID - 1.0);
        } else {
            disk_params.DD = 0.0;
        }
        fprintf(stderr, "DEBUG [main]: NGRID set to default/command-line value: %d. DD calculated as %.4e.\n", disk_params.NGRID, disk_params.DD);
    }

    // --- CRITICAL STEP: Now that current_inputsig_file is definitively set,
    //                  copy it to sim_opts.input_filename for tIntegrate. ---
    fprintf(stderr, "DEBUG [main]: Preparing to set sim_opts.input_filename. current_inputsig_file = '%s' (address %p).\n",
                current_inputsig_file ? current_inputsig_file : "(NULL or empty)", (void*)current_inputsig_file);
    if (current_inputsig_file != NULL) {
        strncpy(sim_opts.input_filename, current_inputsig_file, MAX_PATH_LEN - 1); // Use MAX_PATH_LEN - 1 to ensure null-termination
        sim_opts.input_filename[MAX_PATH_LEN - 1] = '\0'; // Explicitly ensure null-termination
    } else {
        sim_opts.input_filename[0] = '\0'; // Set to empty string if NULL, though should not happen here.
    }
    fprintf(stderr, "DEBUG [main]: sim_opts.input_filename set to '%s'.\n", sim_opts.input_filename);


    // --- CRITICAL STEP: Ensure global PARTICLE_NUMBER matches disk_params.NGRID ---
    PARTICLE_NUMBER = disk_params.NGRID;
    fprintf(stderr, "DEBUG [main]: Global PARTICLE_NUMBER set to %d (from disk_params.NGRID).\n", PARTICLE_NUMBER);


    // --- Dynamic Memory Allocation for Disk Arrays ---
    // Allocate NGRID+2 elements for boundary conditions
    disk_params.rvec = (double *)malloc((disk_params.NGRID + 2) * sizeof(double));
    disk_params.sigmavec = (double *)malloc((disk_params.NGRID + 2) * sizeof(double));
    disk_params.pressvec = (double *)malloc((disk_params.NGRID + 2) * sizeof(double));
    disk_params.dpressvec = (double *)malloc((disk_params.NGRID + 2) * sizeof(double));
    disk_params.ugvec = (double *)malloc((disk_params.NGRID + 2) * sizeof(double));

    if (!disk_params.rvec || !disk_params.sigmavec || !disk_params.pressvec || !disk_params.dpressvec || !disk_params.ugvec) {
        fprintf(stderr, "ERROR [main]: Failed to allocate memory for disk arrays. Exiting.\n");
        return 1; // Exit on memory allocation failure
    }
    fprintf(stderr, "DEBUG [main]: Disk profile arrays dynamically allocated with size NGRID+2 = %d.\n", disk_params.NGRID + 2);

    // Call disk_param_be to set PDENSITY and PDENSITYDIMLESS within disk_params
    fprintf(stderr, "DEBUG [main]: Calling disk_param_be to calculate derived disk parameters.\n");
    disk_param_be(&disk_params);
    // After disk_param_be, the PDENSITY and PDENSITYDIMLESS are now correctly set in disk_params.
    fprintf(stderr, "DEBUG [main]: disk_param_be completed. disk_params.PDENSITY=%.2e, disk_params.PDENSITYDIMLESS=%.2e.\n", disk_params.PDENSITY, disk_params.PDENSITYDIMLESS);


    // --- Populate init_tool_options_t from 'def' (command line) values ---
    init_tool_params.n_grid_points = disk_params.NGRID; // Use the NGRID from disk_params, which might be updated by input file
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
    init_tool_params.two_pop_ratio = def.ratio_val;
    init_tool_params.micro_size_cm = def.mic_val;
    init_tool_params.one_size_particle_cm = def.onesize_val;
    fprintf(stderr, "DEBUG [main]: init_tool_options_t (init_tool_params) structure populated for profile generation.\n");
    fprintf(stderr, "DEBUG [main]:   init_tool_params.n=%d, init_tool_params.ri=%.2f, init_tool_params.ro=%.2f, init_tool_params.sigma0=%.2e\n",
           init_tool_params.n_grid_points, init_tool_params.r_inner, init_tool_params.r_outer, init_tool_params.sigma0_gas_au);


    fprintf(stderr, "DEBUG [main]: Time parameters are now directly in sim_opts. No separate timePar() call needed if it just sets globals.\n");


    // --- Initial profile loading for sigIn ---
    // Note: current_inputsig_file is already correctly set to either def.input_file or "init_data.dat"
    // and disk_params.NGRID is already updated to reflect the file's content.
    fprintf(stderr, "DEBUG [main]: Loading initial profile from '%s'.\n", current_inputsig_file);
    fprintf(stderr, "DEBUG [main]: Calling sigIn(disk_params.sigmavec, disk_params.rvec, &disk_params, '%s')...\n", current_inputsig_file);
    sigIn(disk_params.sigmavec, disk_params.rvec, &disk_params, current_inputsig_file);
    fprintf(stderr, "DEBUG [main]: sigIn completed. Calling Perem for disk_params.rvec and disk_params.sigmavec...\n");
    Perem(disk_params.rvec,&disk_params);
    Perem(disk_params.sigmavec,&disk_params);
    fprintf(stderr, "DEBUG [main]: Perem calls completed for initial profile.\n");

    // --- End of initial profile handling ---

    // Initialize pressure and gas velocity profiles
    fprintf(stderr, "DEBUG [main]: Initializing pressure and gas velocity profiles...\n");
    Initial_Press(&disk_params);
    Initial_dPress(&disk_params);
    Initial_Ugas(&disk_params);
    fprintf(stderr, "DEBUG [main]: Pressure and gas velocity profiles initialized.\n");

    // --- Output directory handling ---
    fprintf(stderr, "DEBUG [main]: Creating main output directory: '%s'.\n", def.output_dir_name);
    Mk_Dir(def.output_dir_name); // Creates the main output folder

    char initial_dir_path[MAX_PATH_LEN]; // Using MAX_PATH_LEN for consistency
    char logs_dir_path[MAX_PATH_LEN];    // Using MAX_PATH_LEN for consistency

    // Create the 'initial' subdirectory
    snprintf(initial_dir_path, sizeof(initial_dir_path), "%s/initial", def.output_dir_name);
    fprintf(stderr, "DEBUG [main]: Creating initial files directory: '%s'.\n", initial_dir_path);
    Mk_Dir(initial_dir_path);

    // Create the 'LOGS' subdirectory
    snprintf(logs_dir_path, sizeof(logs_dir_path), "%s/LOGS", def.output_dir_name);
    fprintf(stderr, "DEBUG [main]: Creating LOGS directory: '%s'.\n", logs_dir_path);
    Mk_Dir(logs_dir_path);

    fprintf(stderr, "DEBUG [main]: Output subdirectories created.\n");

    int dummy_sys_ret; // Dummy variable for system() call return value
    char cmd_buffer[MAX_PATH_LEN * 2]; // Buffer for system commands, needs to accommodate two paths and "cp "

    // Copy config and initial profile files to the 'initial' directory
    fprintf(stderr, "DEBUG [main]: Copying config and initial profile files to '%s'.\n", initial_dir_path);

    // Copy disk_param.dat
    snprintf(cmd_buffer, sizeof(cmd_buffer), "cp %s %s/", FILENAME_DISK_PARAM, initial_dir_path);
    dummy_sys_ret = system(cmd_buffer); (void)dummy_sys_ret; // Cast to void to suppress unused result warning
    fprintf(stderr, "DEBUG [main]: Copied %s to %s/\n", FILENAME_DISK_PARAM, initial_dir_path);

    // Copy input profile file (or generated init_data.dat)
    if(def.input_file != NULL && strcmp(def.input_file, "") != 0) {
        snprintf(cmd_buffer, sizeof(cmd_buffer), "cp %s %s/", def.input_file, initial_dir_path);
        dummy_sys_ret = system(cmd_buffer); (void)dummy_sys_ret;
        fprintf(stderr, "DEBUG [main]: Copied initial profile file '%s' to '%s/'.\n", def.input_file, initial_dir_path);
    } else {
        // Use current_inputsig_file which is correctly set to "init_data.dat" in this branch
        snprintf(cmd_buffer, sizeof(cmd_buffer), "cp %s %s/", current_inputsig_file, initial_dir_path);
        dummy_sys_ret = system(cmd_buffer); (void)dummy_sys_ret;
        fprintf(stderr, "DEBUG [main]: Copied generated profile '%s' to '%s/'.\n", current_inputsig_file, initial_dir_path);
    }

    // Print current information
    fprintf(stderr, "DEBUG [main]: Calling infoCurrent...\n");
    // This infoCurrent is likely for the main output directory, not 'initial' specifically.
    infoCurrent(def.output_dir_name, &disk_params, &sim_opts);
    fprintf(stderr, "DEBUG [main]: infoCurrent completed.\n");

    // --- Final Debug Check before main simulation loop ---
    fprintf(stderr, "DEBUG [main]: Final check before tIntegrate:\n");
    fprintf(stderr, "DEBUG [main]:   sim_opts.output_dir_name = '%s'\n", sim_opts.output_dir_name);
    fprintf(stderr, "DEBUG [main]:   disk_params.RMIN = %.2f, RMAX = %.2f, NGRID = %d, DD = %.2e\n",
            disk_params.RMIN, disk_params.RMAX, disk_params.NGRID, disk_params.DD);
    fprintf(stderr, "DEBUG [main]:   disk_params.SIGMA0 = %.2e, SIGMAP_EXP = %.2f, alpha_visc = %.2e\n",
            disk_params.SIGMA0, disk_params.SIGMAP_EXP, disk_params.alpha_visc);
    fprintf(stderr, "DEBUG [main]:   disk_params.STAR_MASS = %.2f, HASP = %.2f, FLIND = %.2f\n",
            disk_params.STAR_MASS, disk_params.HASP, disk_params.FLIND);
    fprintf(stderr, "DEBUG [main]:   disk_params.PDENSITY = %.2e, PDENSITYDIMLESS = %.2e\n",
            disk_params.PDENSITY, disk_params.PDENSITYDIMLESS);
    fprintf(stderr, "DEBUG [main]:   disk_params.rvec[0] = %.2e (sample first rvec element)\n",
            disk_params.rvec ? disk_params.rvec[0] : -999.9);
    // --- END Final Debug Check ---


    // Run simulation or exit based on options
    if(sim_opts.evol == 0. && sim_opts.drift == 0.) {
        fprintf(stderr, "DEBUG [main]: Evolution (sim_opts.evol=%.2f) and drift (sim_opts.drift=%.2f) are OFF.\n", sim_opts.evol, sim_opts.drift);
        printf("A megadott opciok szerint nem szamol sem sigmat, sem driftet, ezert a program kilep!\n\nA kezdeti file-ok a %s mappaban talalhatoak!\n", initial_dir_path);

        char dens_name_initial[MAX_PATH_LEN]; // Using MAX_PATH_LEN for consistency
        // Corrected typo: removed 'da' from filename
        snprintf(dens_name_initial, sizeof(dens_name_initial), "%s/%s", initial_dir_path,INITIAL_SURFACE_DENSITY_FILE);
        fprintf(stderr, "DEBUG [main]: Printing initial surface density to %s.\n", dens_name_initial);

        // --- Special handling for Print_Sigma when only initial output is needed ---
        // Open the file directly here as it's a one-off write to the 'initial' folder.
        // The output_files_t struct is usually for files opened in tIntegrate.
        output_files_t temp_output_for_initial_print;
        temp_output_for_initial_print.surface_file = fopen(dens_name_initial, "w");
        if (temp_output_for_initial_print.surface_file != NULL) {
            Print_Sigma(&disk_params, &temp_output_for_initial_print);
            fclose(temp_output_for_initial_print.surface_file);
            temp_output_for_initial_print.surface_file = NULL; // Reset pointer
            fprintf(stderr, "DEBUG [main]: Closed %s.\n", dens_name_initial);
        } else {
            fprintf(stderr, "ERROR [main]: Could not open %s for initial surface output.\n", dens_name_initial);
        }

        fprintf(stderr, "DEBUG [main]: Print_Sigma completed. Program exiting.\n");
    } else {
        fprintf(stderr, "DEBUG [main]: Evolution (sim_opts.evol=%.2f) or drift (sim_opts.drift=%.2f) is ON. Starting main simulation loop.\n", sim_opts.evol, sim_opts.drift);
        fprintf(stderr, "DEBUG [main]: Calling tIntegrate...\n");
        // Pass the LOGS directory path to tIntegrate
        tIntegrate(&disk_params, &sim_opts, &output_files);
        fprintf(stderr, "DEBUG [main]: tIntegrate completed. Program finished normally.\n");
    }

    // --- Free dynamically allocated memory ---
    free(disk_params.rvec);
    free(disk_params.sigmavec);
    free(disk_params.pressvec);
    free(disk_params.dpressvec);
    free(disk_params.ugvec);
    fprintf(stderr, "DEBUG [main]: Dynamically allocated disk arrays freed.\n");

    fprintf(stderr, "DEBUG [main]: Program exiting normally.\n");
    return 0;
}