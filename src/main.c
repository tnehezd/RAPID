#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// Header files
#include "config.h"           // Declarations of global variables and constants (to be reduced)
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


int main(int argc, const char **argv) {
    // DEBUG: Program start
    printf("DEBUG [main]: Program started.\n");

    // Local structure to store command-line options
    options_t def;
    // Set default values
    create_default_options(&def);
    printf("DEBUG [main]: Default options created.\n");

    // Local structure to store init_tool parameters (these will be populated from 'def')
    init_tool_options_t init_tool_params;
    // Set default values for init_tool. This is useful if 'def' doesn't override all of them.
    create_default_init_tool_options(&init_tool_params);
    printf("DEBUG [main]: Default init_tool options created.\n");

    // Parse command-line options and populate the 'def' structure
    int retCode = parse_options(argc, argv, &def);

    if (0 != retCode) {
        // Exit on error (usage printed by parse_options)
        printf("DEBUG [main]: Error parsing command-line options. Exiting with code %d.\n", retCode);
        return retCode;
    }
    printf("DEBUG [main]: Command-line options parsed successfully.\n");

    // --- Declare instances of the new simulation structs ---
    disk_t my_disk;
    simulation_options_t sim_opts;
    // (dust_population_t por_dust, micron_dust, output_files_t output_files, if needed later)

    printf("DEBUG [main]: New simulation structs declared.\n");

    /* Populate the simulation_options_t struct from 'def' (parsed options) */
    sim_opts.evol = def.evol;
    sim_opts.drift = def.drift;
    sim_opts.growth = def.growth;
    sim_opts.twopop = def.twopop;
    sim_opts.fFrag = def.ffrag;
    sim_opts.uFrag = def.ufrag;
    sim_opts.DT = def.tStep;
    sim_opts.TMAX = def.totalTime;
    sim_opts.WO = def.outputFrequency;
    sim_opts.TCURR = def.startTime; // Initial current time
    printf("DEBUG [main]: sim_opts structure populated.\n");
    printf("DEBUG [main]:   evol=%.2f, drift=%.2f, growth=%.2f, twopop=%.2f, fFrag=%.2f, uFrag=%.2f\n",
           sim_opts.evol, sim_opts.drift, sim_opts.growth, sim_opts.twopop, sim_opts.fFrag, sim_opts.uFrag);
    printf("DEBUG [main]:   Time parameters set: TMAX=%.2e, WO=%.2e, TCURR=%.2e, DT=%.2e\n",
           sim_opts.TMAX, sim_opts.WO, sim_opts.TCURR, sim_opts.DT);


    // --- Populate disk_t with parameters from 'def' ---
    my_disk.RMIN = def.rmin_val;
    my_disk.RMAX = def.rmax_val;
    my_disk.NGRID = def.ngrid_val; // Initial NGRID, might be updated by input file reading
    my_disk.SIGMA0 = def.sigma0_val;
    my_disk.SIGMAP_EXP = def.sigmap_exp_val;
    my_disk.alpha_visc = def.alpha_visc_val;
    my_disk.STAR_MASS = def.star_val;
    my_disk.HASP = def.hasp_val;
    my_disk.FLIND = def.flind_val;
    my_disk.r_dze_i = def.r_dze_i_val;
    my_disk.r_dze_o = def.r_dze_o_val;
    my_disk.Dr_dze_i = def.dr_dze_i_val;
    my_disk.Dr_dze_o = def.dr_dze_o_val;
    my_disk.a_mod = def.a_mod_val;

    printf("DEBUG [main]: my_disk structure populated with initial parameters.\n");
    printf("DEBUG [main]:   RMIN=%.2f, RMAX=%.2f, NGRID=%d, SIGMA0=%.2e, SIGMAP_EXP=%.2f\n",
           my_disk.RMIN, my_disk.RMAX, my_disk.NGRID, my_disk.SIGMA0, my_disk.SIGMAP_EXP);
    printf("DEBUG [main]:   alpha_visc=%.2e, STAR_MASS=%.2f, HASP=%.2f, FLIND=%.2f\n",
           my_disk.alpha_visc, my_disk.STAR_MASS, my_disk.HASP, my_disk.FLIND);
    printf("DEBUG [main]:   r_dze_i=%.2f, r_dze_o=%.2f, Dr_dze_i=%.2f, Dr_dze_o=%.2f, a_mod=%.2f\n",
           my_disk.r_dze_i, my_disk.r_dze_o, my_disk.Dr_dze_i, my_disk.Dr_dze_o, my_disk.a_mod);


    // --- CRITICAL STEP: TRANSFER PARSED VALUES FROM STRUCTS TO GLOBAL VARIABLES ---
    // This is the missing link that causes "optdr is OFF"
    printf("DEBUG [main]: Transferring parsed values from structs to GLOBAL variables.\n");
    optdr     = sim_opts.drift;
    optgr     = sim_opts.growth;
    optev     = sim_opts.evol;
    opttwopop = sim_opts.twopop;

    // Set optdze based on dead zone radii from my_disk
    optdze = (my_disk.r_dze_i > 0.0 || my_disk.r_dze_o > 0.0) ? 1.0 : 0.0;
    printf("DEBUG [main]: Global optdr set to %.2f (from sim_opts.drift).\n", optdr);
    printf("DEBUG [main]: Global optgr set to %.2f (from sim_opts.growth).\n", optgr);
    printf("DEBUG [main]: Global optev set to %.2f (from sim_opts.evol).\n", optev);
    printf("DEBUG [main]: Global opttwopop set to %.2f (from sim_opts.twopop).\n", opttwopop);
    printf("DEBUG [main]: Global optdze set to %.2f (based on my_disk.r_dze_i/o).\n", optdze);

    // Also transfer other global disk parameters that are still used by legacy functions
    RMIN = my_disk.RMIN;
    RMAX = my_disk.RMAX;
    NGRID = my_disk.NGRID;
    SIGMA0 = my_disk.SIGMA0;
    SIGMAP_EXP = my_disk.SIGMAP_EXP;
    alpha_visc = my_disk.alpha_visc;
    STAR = my_disk.STAR_MASS;
    HASP = my_disk.HASP;
    FLIND = my_disk.FLIND;
    r_dze_i = my_disk.r_dze_i;
    r_dze_o = my_disk.r_dze_o;
    Dr_dze_i = my_disk.Dr_dze_i;
    Dr_dze_o = my_disk.Dr_dze_o;
    a_mod = my_disk.a_mod;
    DD = my_disk.DD; // DD is calculated after NGRID is finalized

    TMAX = sim_opts.TMAX;
    WO = sim_opts.WO;
    TCURR = sim_opts.TCURR;
    // --- END OF CRITICAL GLOBAL UPDATES ---


    // --- Input file handling logic ---
    const char *current_inputsig_file = NULL;

    if (def.input_file != NULL && strcmp(def.input_file, "") != 0) {
        current_inputsig_file = def.input_file;
        printf("DEBUG [main]: Input file specified: '%s'. Attempting to read initial profile.\n", current_inputsig_file);

        // Update NGRID based on file content. reszecskek_szama needs refactoring to just take filename.
        // It should *not* rely on global 'lout'.
        // Assuming reszecskek_szama is available and returns the correct NGRID from the file.
        // It's critical that this NGRID accurately reflects the file, otherwise malloc will be wrong.
        // Pass the output_files for its file pointers if needed, or remove 'lout' entirely.
        int dummy_lout = 0; // Temporarily keeping dummy_lout for reszecskek_szama
        my_disk.NGRID = reszecskek_szama(dummy_lout, current_inputsig_file); // Update NGRID from file

        // Re-calculate DD based on updated NGRID
        if (my_disk.NGRID > 1) { // Avoid division by zero if NGRID is 1
            my_disk.DD = (my_disk.RMAX - my_disk.RMIN) / (my_disk.NGRID - 1.0);
        } else {
            my_disk.DD = 0.0; // Or handle as an error if NGRID <= 1 is invalid for your simulation
        }
        // Update global DD after recalculation
        DD = my_disk.DD;

        printf("DEBUG [main]: NGRID set from input file: %d. DD calculated as %.4e.\n", my_disk.NGRID, my_disk.DD);
    } else {
        printf("DEBUG [main]: No input file specified (-i flag not used). Generating default grid and profile.\n");
        // NGRID is already set from def.ngrid_val
        if (my_disk.NGRID > 1) {
            my_disk.DD = (my_disk.RMAX - my_disk.RMIN) / (my_disk.NGRID - 1.0);
        } else {
            my_disk.DD = 0.0;
        }
        // Update global DD after calculation
        DD = my_disk.DD;
        printf("DEBUG [main]: NGRID set to default/command-line value: %d. DD calculated as %.4e.\n", my_disk.NGRID, my_disk.DD);
    }

     // --- CRITICAL STEP: Ensure global PARTICLE_NUMBER matches my_disk.NGRID ---
    PARTICLE_NUMBER = my_disk.NGRID;
    printf("DEBUG [main]: Global PARTICLE_NUMBER set to %d (from my_disk.NGRID).\n", PARTICLE_NUMBER);


    // --- Dynamic Memory Allocation for Disk Arrays ---
    // Allocate NGRID+2 elements for boundary conditions
    my_disk.rvec = (double *)malloc((my_disk.NGRID + 2) * sizeof(double));
    my_disk.sigmavec = (double *)malloc((my_disk.NGRID + 2) * sizeof(double));
    my_disk.pressvec = (double *)malloc((my_disk.NGRID + 2) * sizeof(double));
    my_disk.dpressvec = (double *)malloc((my_disk.NGRID + 2) * sizeof(double));
    my_disk.ugvec = (double *)malloc((my_disk.NGRID + 2) * sizeof(double));

    if (!my_disk.rvec || !my_disk.sigmavec || !my_disk.pressvec || !my_disk.dpressvec || !my_disk.ugvec) {
        fprintf(stderr, "ERROR [main]: Failed to allocate memory for disk arrays. Exiting.\n");
        return 1; // Exit on memory allocation failure
    }
    printf("DEBUG [main]: Disk profile arrays dynamically allocated with size NGRID+2 = %d.\n", my_disk.NGRID + 2);

    // Call disk_param_be to set PDENSITY and PDENSITYDIMLESS within my_disk
    // This function is already refactored to use `disk_t *`
    printf("DEBUG [main]: Calling disk_param_be to calculate derived disk parameters.\n");
    disk_param_be(&my_disk);
    // After disk_param_be, the PDENSITY and PDENSITYDIMLESS are now correctly set in my_disk.
    // Also, update global PDENSITY/PDENSITYDIMLESS for legacy functions.
    PDENSITY = my_disk.PDENSITY;
    PDENSITYDIMLESS = my_disk.PDENSITYDIMLESS;
    printf("DEBUG [main]: disk_param_be completed. my_disk.PDENSITY=%.2e, my_disk.PDENSITYDIMLESS=%.2e.\n", my_disk.PDENSITY, my_disk.PDENSITYDIMLESS);


    // --- Populate init_tool_options_t from 'def' (command line) values ---
    // (since these are the finalized parameters for profile generation)
    // Also ensuring that NGRID used by init_tool matches my_disk.NGRID.
    init_tool_params.n_grid_points = my_disk.NGRID; // Use the NGRID from my_disk, which might be updated by input file
    init_tool_params.r_inner= my_disk.RMIN;
    init_tool_params.r_outer = my_disk.RMAX;
    init_tool_params.sigma0_gas_au = my_disk.SIGMA0;
    //init_tool_params.sigma0cgs = my_disk.SIGMA0 / SDCONV; // SDCONV needs to be a constant/macro or part of sim_opts
    init_tool_params.sigma_exponent = my_disk.SIGMAP_EXP;
    init_tool_params.deadzone_r_inner = my_disk.r_dze_i;
    init_tool_params.deadzone_r_outer = my_disk.r_dze_o;
    init_tool_params.deadzone_dr_inner = my_disk.Dr_dze_i;
    init_tool_params.deadzone_dr_outer = my_disk.Dr_dze_o;
    init_tool_params.alpha_viscosity = my_disk.alpha_visc;
    init_tool_params.deadzone_alpha_mod = my_disk.a_mod;
    init_tool_params.aspect_ratio = my_disk.HASP;
    init_tool_params.flaring_index = my_disk.FLIND;
    init_tool_params.star_mass = my_disk.STAR_MASS;
    // The following init_tool specific options come from 'def' (command line)
    // Make sure 'md_val' from options_t is correctly mapped if it exists,
    // otherwise use the default or a derived value.
    // init_tool_params.disk_mass_dust = def.md_val; // If 'md_val' is a direct input for init_tool
    init_tool_params.dust_to_gas_ratio = def.eps_val;
    init_tool_params.two_pop_ratio = def.ratio_val;
    init_tool_params.micro_size_cm = def.mic_val; // Corrected from micro_size_cm
    init_tool_params.one_size_particle_cm = def.onesize_val;
    printf("DEBUG [main]: init_tool_options_t (init_tool_params) structure populated for profile generation.\n");
    printf("DEBUG [main]:   init_tool_params.n=%d, init_tool_params.ri=%.2f, init_tool_params.ro=%.2f, init_tool_params.sigma0=%.2e\n",
           init_tool_params.n_grid_points, init_tool_params.r_inner, init_tool_params.r_outer, init_tool_params.sigma0_gas_au);


    // Time parameters: No need to call timePar if it just updates globals.
    // Instead, just pass sim_opts to relevant functions.
    printf("DEBUG [main]: Time parameters are now directly in sim_opts. No separate timePar() call needed if it just sets globals.\n");


    // --- Initial profile loading or generation ---
    if(current_inputsig_file != NULL) {
        printf("DEBUG [main]: current_inputsig_file is NOT NULL. Loading initial profile from '%s'.\n", current_inputsig_file);
        // sigIn must be refactored to accept disk_t*
        printf("DEBUG [main]: Calling sigIn(my_disk.sigmavec, my_disk.rvec, &my_disk, '%s')...\n", current_inputsig_file);
        // If sigIn needs the filename, it must be passed! Assuming a refactored sigIn like:
        sigIn(my_disk.sigmavec, my_disk.rvec, &my_disk, current_inputsig_file);
        printf("DEBUG [main]: sigIn completed. Calling Perem for my_disk.rvec and my_disk.sigmavec...\n");
        Perem(my_disk.rvec);    // Perem conditions for rvec (needs refactoring to take NGRID/DD from my_disk)
        Perem(my_disk.sigmavec); // Perem conditions for sigmavec (needs refactoring)
        printf("DEBUG [main]: Perem calls completed for initial profile.\n");
    } else {
        printf("DEBUG [main]: current_inputsig_file is NULL. Generating initial profile with init_tool...\n");
        printf("DEBUG [main]: Calling run_init_tool(&init_tool_params)...\n");
        run_init_tool(&init_tool_params); // This is refactored.
        printf("DEBUG [main]: run_init_tool completed.\n");

        current_inputsig_file = "init_data.dat"; // This is the default output of init_tool
        printf("DEBUG [main]: Generated profile will be loaded from '%s'.\n", current_inputsig_file);

        // NGRID might have been reset or confirmed by init_tool and init_tool_params.n.
        // If init_tool could change NGRID in a way that wasn't already updated in my_disk.NGRID
        // you would need to recalculate my_disk.DD here and potentially realloc.
        // For now, assuming init_tool_params.n directly controls my_disk.NGRID.
        printf("DEBUG [main]: Calling sigIn(my_disk.sigmavec, my_disk.rvec, &my_disk, '%s') for generated profile...\n", current_inputsig_file);
        sigIn(my_disk.sigmavec, my_disk.rvec, &my_disk, current_inputsig_file); // sigIn must be refactored to take filename
        printf("DEBUG [main]: sigIn completed for generated profile. Calling Perem...\n");
        Perem(my_disk.rvec);
        Perem(my_disk.sigmavec);
        printf("DEBUG [main]: Perem calls completed for generated profile.\n");
    }
    // --- End of initial profile handling ---

    // Initialize pressure and gas velocity profiles
    printf("DEBUG [main]: Initializing pressure and gas velocity profiles...\n");
    // These functions (Initial_Press, Initial_dPress, Initial_Ugas) need to be refactored
    // to take `disk_t *` as argument to get rvec, sigmavec, pressvec, dpressvec, ugvec, NGRID, DD etc.
    Initial_Press(my_disk.pressvec, my_disk.sigmavec, my_disk.rvec);
    Initial_dPress(my_disk.dpressvec, my_disk.pressvec);
    Initial_Ugas(my_disk.sigmavec, my_disk.rvec, my_disk.ugvec);
    printf("DEBUG [main]: Pressure and gas velocity profiles initialized.\n");

    // Dead Zone option (optdze)
    // This should ideally be moved into sim_opts or disk_t.
    // It's already handled in the "CRITICAL GLOBAL UPDATES" section above.
    // Remove this duplicate logic here:
    /*
    optdze = 0; // Default to inactive (global)
    if(my_disk.r_dze_i != 0.0 || my_disk.r_dze_o != 0.0) {
        optdze = 1.;
        printf("DEBUG [main]: Dead Zone active (r_dze_i != 0.0 or r_dze_o != 0.0).\n");
    } else {
        printf("DEBUG [main]: Dead Zone inactive (both r_dze_i and r_dze_o are 0.0).\n");
    }
    */
    // The debug print for optdze is already done above.


    // Output directory handling
    printf("DEBUG [main]: Creating output directory: '%s'.\n", def.output_dir_name);
    // Mk_Dir should just create the directory and not rely on or set global 'nev'.
    // If it currently sets 'nev', you'll need to refactor Mk_Dir itself.
    // For now, assuming it handles its internal state and 'def.output_dir_name' is the source of truth.
    Mk_Dir(def.output_dir_name);
    printf("DEBUG [main]: Output directory created: %s\n", def.output_dir_name);

    int dummy_sys_ret; // Dummy variable for system() call return value

    // Copy config and time files to the output directory
    printf("DEBUG [main]: Copying config and time files to output directory.\n");
    char cmd_buffer[4096];
    // filenev2 ("disk_param.dat") is still global.
    // These need to be refactored into 'output_files_t' or passed as string literals.
    snprintf(cmd_buffer, sizeof(cmd_buffer), "cp %s %s/", filenev2, def.output_dir_name);
    dummy_sys_ret = system(cmd_buffer); (void)dummy_sys_ret; // Cast to void to suppress unused result warning

    // Copy input profile file to the output directory
    if(def.input_file != NULL && strcmp(def.input_file, "") != 0) {
        printf("DEBUG [main]: Copying initial profile file '%s' to output directory.\n", def.input_file);
        snprintf(cmd_buffer, sizeof(cmd_buffer), "cp %s %s/", def.input_file, def.output_dir_name);
        dummy_sys_ret = system(cmd_buffer); (void)dummy_sys_ret;
        printf("DEBUG [main]: Initial profile file copied.\n");
    } else {
        // If no input file was given, init_tool generated 'init_data.dat' (filenev1).
        // Copy the generated file.
        printf("DEBUG [main]: No specific input file to copy. Copying generated profile '%s' to output directory.\n", filenev1);
        snprintf(cmd_buffer, sizeof(cmd_buffer), "cp %s %s/", filenev1, def.output_dir_name);
        dummy_sys_ret = system(cmd_buffer); (void)dummy_sys_ret;
        printf("DEBUG [main]: Generated profile file copied.\n");
    }

    // Print current information
    printf("DEBUG [main]: Calling infoCurrent...\n");
    // infoCurrent needs to be refactored to take `simulation_options_t *` and `disk_t *`
    // and the output directory name, instead of relying on global 'nev'.
    infoCurrent(def.output_dir_name); // Assuming it uses the directory name passed to it.
    printf("DEBUG [main]: infoCurrent completed.\n");

    // Run simulation or exit based on options
    if(sim_opts.evol == 0. && sim_opts.drift == 0.) {
        printf("DEBUG [main]: Evolution (sim_opts.evol=%.2f) and drift (sim_opts.drift=%.2f) are OFF.\n", sim_opts.evol, sim_opts.drift);
        printf("A megadott opciok szerint nem szamol sem sigmat, sem driftet, ezert a program kilep!\n\nA kezdeti file-ok a %s mappaban talalhatoak!\n", def.output_dir_name);
        char dens_name[4096];
        snprintf(dens_name, sizeof(dens_name), "%s/surface.dat", def.output_dir_name);
        printf("DEBUG [main]: Printing initial surface density to %s.\n", dens_name);
        // Print_Sigma needs to be refactored to take `disk_t *` and `output_files_t *`
        // For now, it will use global variables it depends on.
        Print_Sigma(dens_name, my_disk.rvec, my_disk.sigmavec, my_disk.pressvec, my_disk.dpressvec);
        printf("DEBUG [main]: Print_Sigma completed. Program exiting.\n");
    } else {
        printf("DEBUG [main]: Evolution (sim_opts.evol=%.2f) or drift (sim_opts.drift=%.2f) is ON. Starting main simulation loop.\n", sim_opts.evol, sim_opts.drift);
        printf("DEBUG [main]: Calling tIntegrate...\n");
        // tIntegrate needs to be refactored to take `disk_t *`, `simulation_options_t *`,
        // `dust_population_t *` (for both populations), and `output_files_t *`.
        // For now, it relies on the global variables updated above.
        tIntegrate(def.output_dir_name, &my_disk);
        printf("DEBUG [main]: tIntegrate completed. Program finished normally.\n");
    }

    // --- Free dynamically allocated memory ---
    free(my_disk.rvec);
    free(my_disk.sigmavec);
    free(my_disk.pressvec);
    free(my_disk.dpressvec);
    free(my_disk.ugvec);
    printf("DEBUG [main]: Dynamically allocated disk arrays freed.\n");

    printf("DEBUG [main]: Program exiting normally.\n");
    return 0;
}