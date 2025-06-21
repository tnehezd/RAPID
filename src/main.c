#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <string.h> // For strcpy


// Header files
#include "config.h"           // Declarations of global variables and constants (we will reduce reliance on these)
#include "io_utils.h"         // Functions from io_utils.c (needs refactoring for structs)
#include "disk_model.h"       // Functions from disk_model.c (ALREADY refactored for structs)
#include "dust_physics.h"     // Functions from dust_physics.c (needs refactoring for structs)
#include "simulation_core.h"  // Functions from simulation_core.c (needs refactoring for structs)
#include "utils.h"            // Functions from utils.c (needs refactoring for structs)
#include "init_tool_module.h" // run_init_tool and init_tool_options_t (ALREADY refactored for structs)

// --- NEW: Include your simulation_types.h ---
#include "simulation_types.h"

// options_t structure: Gathers all command-line options.
// This structure also contains parameters used by init_tool
// for initial profile generation.
// IMPORTANT: This 'options_t' is primarily for parsing command-line args.
// Its values are then copied into the more specific 'disk_t' and 'simulation_options_t' structs.
typedef struct options {
    double drift;
    double growth;
    double evol;
    double twopop;
    double ufrag;
    double ffrag;

    // Core simulation and init_tool parameters
    int ngrid_val;             // NGRID
    double rmin_val;          // RMIN
    double rmax_val;          // RMAX
    double sigma0_val;        // SIGMA0
    double sigmap_exp_val;    // SIGMAP_EXP (Note: this is the positive exponent, init_tool will negate if needed)
    double alpha_visc_val;    // alpha_visc
    double star_val;          // STAR
    double hasp_val;          // HASP
    double flind_val;         // FLIND
    double r_dze_i_val;       // r_dze_i
    double r_dze_o_val;       // r_dze_o
    double dr_dze_i_val;      // Dr_dze_i (This is the multiplier for the transition width)
    double dr_dze_o_val;      // Dr_dze_o (This is the multiplier for the transition width)
    double a_mod_val;         // a_mod

    const char *input_file;   // Filename for initialization/loading. NULL if not specified.
    char output_dir_name[256]; // Name of the output directory (allocate fixed space)


    double tStep;             // Fixed time step (DT in sim_opts)
    double totalTime;         // TMAX in sim_opts
    double outputFrequency;   // WO in sim_opts
    double startTime;         // TCURR in sim_opts

    // Init tool specific options needed for profile generation
    double md_val;
    long double eps_val;
    long double ratio_val;
    long double mic_val;
    long double onesize_val;

} options_t;

// Function declarations
void create_default_options(options_t *def);
int parse_options(int argc, const char **argv, options_t *def);


// --- Global variables: ---
// IMPORTANT: These are still here because other (unrefactored) functions might use them.
// The goal is to remove all of them eventually by passing structs instead.
// For now, we update them from our structs to maintain compatibility with legacy code.
// Consider moving these to config.h if they aren't already there.
extern double TMAX, WO, TCURR;
extern double RMIN, RMAX, SIGMA0, SIGMAP_EXP, alpha_visc, STAR, HASP, FLIND,
              r_dze_i, r_dze_o, Dr_dze_i, Dr_dze_o, a_mod;
extern int NGRID;
extern double DD;
extern double PDENSITY, PDENSITYDIMLESS; // Set by disk_param_be
extern double optdze; // Dead Zone option
extern char filenev1[1024], filenev2[1024], filenev3[1024], nev[1024];


int main(int argc, const char **argv) {
    // DEBUG: Program start
    printf("DEBUG [main]: Program started.\n");

    // Local structure to store command-line options
    options_t def;
    // Set default values
    create_default_options(&def);
    printf("DEBUG [main]: Default options created.\n");

    // Local structure to store init_tool parameters
    init_tool_options_t init_def;
    // Set default values for init_tool
    create_default_init_tool_options(&init_def);
    printf("DEBUG [main]: Default init_tool options created.\n");

    // Parse command-line options and populate the 'def' structure
    int retCode = parse_options(argc, argv, &def);

    if (0 != retCode) {
        // Exit on error
        printf("DEBUG [main]: Error parsing command-line options. Exiting with code %d.\n", retCode);
        return retCode;
    }
    printf("DEBUG [main]: Command-line options parsed successfully.\n");

    // --- Declare instances of the new simulation structs ---
    // These will hold all the simulation state and parameters
    disk_t my_disk;
    simulation_options_t sim_opts;
//    dust_population_t por_dust;    // Primary (cm-sized) dust population
//    dust_population_t micron_dust; // Secondary (micron-sized) dust population (if twopop is active)
//    output_files_t output_files;   // For organizing file pointers

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

        printf("DEBUG [main]: NGRID set from input file: %d. DD calculated as %.4e.\n", my_disk.NGRID, my_disk.DD);
    } else {
        printf("DEBUG [main]: No input file specified (-i flag not used). Generating default grid and profile.\n");
        // NGRID is already set from def.ngrid_val
        if (my_disk.NGRID > 1) {
            my_disk.DD = (my_disk.RMAX - my_disk.RMIN) / (my_disk.NGRID - 1.0);
        } else {
            my_disk.DD = 0.0;
        }
        printf("DEBUG [main]: NGRID set to default/command-line value: %d. DD calculated as %.4e.\n", my_disk.NGRID, my_disk.DD);
    }

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

    // --- TEMPORARY: Update global variables for functions that *still* rely on them ---
    // This section is a compatibility layer. It should be removed once ALL functions
    // (sigIn, Perem, Initial_Press, Initial_dPress, Initial_Ugas, Print_Sigma, tIntegrate, Mk_Dir, infoCurrent, etc.)
    // are refactored to take `disk_t *` and `simulation_options_t *` etc.
    printf("DEBUG [main]: Temporarily updating global disk and time parameters from structs for legacy functions.\n");
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
    DD = my_disk.DD;

    TMAX = sim_opts.TMAX;
    WO = sim_opts.WO;
    TCURR = sim_opts.TCURR;
    // --- END OF TEMPORARY GLOBAL UPDATES ---


    // Call disk_param_be to set PDENSITY and PDENSITYDIMLESS within my_disk
    // This function is already refactored to use `disk_t *`
    printf("DEBUG [main]: Calling disk_param_be to calculate derived disk parameters.\n");
    disk_param_be(&my_disk);
    // After disk_param_be, the PDENSITY and PDENSITYDIMLESS are now correctly set in my_disk.
    // Also, update global PDENSITY/PDENSITYDIMLESS for legacy functions.
    PDENSITY = my_disk.PDENSITY;
    PDENSITYDIMLESS = my_disk.PDENSITYDIMLESS;
    printf("DEBUG [main]: disk_param_be completed. my_disk.PDENSITY=%.2e, my_disk.PDENSITYDIMLESS=%.2e.\n", my_disk.PDENSITY, my_disk.PDENSITYDIMLESS);


    // --- Populate init_tool_options_t from 'my_disk' values ---
    // (since these are the finalized parameters for profile generation)
    // Also ensuring that NGRID used by init_tool matches my_disk.NGRID.
    init_def.n = my_disk.NGRID;
    init_def.ri = my_disk.RMIN;
    init_def.ro = my_disk.RMAX;
    init_def.sigma0 = my_disk.SIGMA0;
    init_def.sigma0cgs = my_disk.SIGMA0 / SDCONV; // SDCONV needs to be a constant/macro or part of sim_opts
    init_def.index = my_disk.SIGMAP_EXP;
    init_def.rdze_i = my_disk.r_dze_i;
    init_def.rdze_o = my_disk.r_dze_o;
    init_def.drdze_i = my_disk.Dr_dze_i;
    init_def.drdze_o = my_disk.Dr_dze_o;
    init_def.alphaParam = my_disk.alpha_visc;
    init_def.amod = my_disk.a_mod;
    init_def.h = my_disk.HASP;
    init_def.flind = my_disk.FLIND;
    init_def.m0 = my_disk.STAR_MASS;
    // The following init_tool specific options come from 'def' (command line)
    init_def.md = def.md_val;
    init_def.eps = def.eps_val;
    init_def.ratio = def.ratio_val;
    init_def.mic = def.mic_val;
    init_def.onesize = def.onesize_val;
    printf("DEBUG [main]: init_tool_options_t (init_def) structure populated for profile generation.\n");
    printf("DEBUG [main]:   init_def.n=%d, init_def.ri=%.2f, init_def.ro=%.2f, init_def.sigma0=%.2e\n",
           init_def.n, init_def.ri, init_def.ro, init_def.sigma0);


    // Time parameters: No need to call timePar if it just updates globals.
    // Instead, just pass sim_opts to relevant functions.
    // If timePar *does* more than just setting globals (e.g., prints a report),
    // then uncomment and refactor its signature to `timePar(&sim_opts);`
    printf("DEBUG [main]: Time parameters are now directly in sim_opts. No separate timePar() call needed if it just sets globals.\n");


    // --- Initial profile loading or generation ---
    if(current_inputsig_file != NULL) {
        printf("DEBUG [main]: current_inputsig_file is NOT NULL. Loading initial profile from '%s'.\n", current_inputsig_file);
        // sigIn must be refactored to accept disk_t*
        printf("DEBUG [main]: Calling sigIn(my_disk.sigmavec, my_disk.rvec, &my_disk, '%s')...\n", current_inputsig_file);
        // If sigIn needs the filename, it must be passed! Assuming a refactored sigIn like:
        sigIn(my_disk.sigmavec, my_disk.rvec, &my_disk, current_inputsig_file);
        printf("DEBUG [main]: sigIn completed. Calling Perem for my_disk.rvec and my_disk.sigmavec...\n");
        Perem(my_disk.rvec);     // Perem conditions for rvec (needs refactoring to take NGRID/DD from my_disk)
        Perem(my_disk.sigmavec); // Perem conditions for sigmavec (needs refactoring)
        printf("DEBUG [main]: Perem calls completed for initial profile.\n");
    } else {
        printf("DEBUG [main]: current_inputsig_file is NULL. Generating initial profile with init_tool...\n");
        printf("DEBUG [main]: Calling run_init_tool(&init_def)...\n");
        run_init_tool(&init_def); // This is refactored.
        printf("DEBUG [main]: run_init_tool completed.\n");

        current_inputsig_file = "init_data.dat"; // This is the default output of init_tool
        printf("DEBUG [main]: Generated profile will be loaded from '%s'.\n", current_inputsig_file);

        // NGRID might have been reset or confirmed by init_tool and init_def.n.
        // If init_tool could change NGRID in a way that wasn't already updated in my_disk.NGRID
        // you would need to recalculate my_disk.DD here and potentially realloc.
        // For now, assuming init_def.n directly controls my_disk.NGRID.
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
    optdze = 0; // Default to inactive (global)
    if(my_disk.r_dze_i != 0.0 || my_disk.r_dze_o != 0.0) {
        optdze = 1.;
        printf("DEBUG [main]: Dead Zone active (r_dze_i != 0.0 or r_dze_o != 0.0).\n");
    } else {
        printf("DEBUG [main]: Dead Zone inactive (both r_dze_i and r_dze_o are 0.0).\n");
    }

    // Output directory handling
    // The `Mk_Dir` function needs to be refactored to simply create the directory
    // based on the provided name and not rely on the global `nev`.
    // It should also return the path or confirm creation.
    printf("DEBUG [main]: Creating output directory: '%s'.\n", def.output_dir_name);
    Mk_Dir(def.output_dir_name); // Mk_Dir should just create the directory. It might still set global 'nev'.
    // Assuming Mk_Dir successfully creates the directory and 'def.output_dir_name' now holds the final path.
    // If 'Mk_Dir' populates 'nev' global, then 'nev' will hold the path.
    // It's best to remove 'nev' entirely.
    // For the `cp` commands below, we'll use `def.output_dir_name`.
    printf("DEBUG [main]: Output directory created: %s\n", def.output_dir_name);

    int dummy_sys_ret; // Dummy variable for system() call return value

    // Copy config and time files to the output directory
    printf("DEBUG [main]: Copying config and time files to output directory.\n");
    char cmd_buffer[4096];
    // filenev2 ("disk_param.dat") and filenev3 ("time.dat") are still global.
    // These need to be refactored into 'output_files_t' or passed as string literals.
    // Assuming filenev2 and filenev3 are defined (e.g., in config.h).
    snprintf(cmd_buffer, sizeof(cmd_buffer), "cp %s %s/", filenev2, def.output_dir_name);
    dummy_sys_ret = system(cmd_buffer); (void)dummy_sys_ret; // Cast to void to suppress unused result warning
    snprintf(cmd_buffer, sizeof(cmd_buffer), "cp %s %s/", filenev3, def.output_dir_name);
    dummy_sys_ret = system(cmd_buffer); (void)dummy_sys_ret;
    printf("DEBUG [main]: disk_param.dat and time.dat copied.\n");

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

/* --- create_default_options: Sets default values for the 'options_t' structure --- */
void create_default_options(options_t *opt) {
    printf("DEBUG [create_default_options]: Setting default values for options_t.\n");
    // Simulation control options
    opt->drift           = 1.;
    opt->growth          = 1.;
    opt->evol            = 1.;
    opt->twopop          = 1.;
    opt->ufrag           = 1000.0;
    opt->ffrag           = 0.37;

    // Core disk parameters (also serve as init_tool defaults)
    opt->ngrid_val       = 2000;
    opt->rmin_val        = 1.0;
    opt->rmax_val        = 100.0;
    opt->sigma0_val      = 1.0;
    opt->sigmap_exp_val  = 1.0;    // e.g., 1.0 for r^-1 profile. This is the positive exponent.
    opt->alpha_visc_val  = 0.01;
    opt->star_val        = 1.0;
    opt->hasp_val        = 0.05;
    opt->flind_val       = 0.5;
    opt->r_dze_i_val     = 0.0;
    opt->r_dze_o_val     = 0.0;
    opt->dr_dze_i_val    = 0.0; // Multiplier for transition width (input parameter for init_tool)
    opt->dr_dze_o_val    = 0.0; // Multiplier for transition width (input parameter for init_tool)
    opt->a_mod_val       = 0.0;

    // File input/output
    opt->input_file      = ""; // Default to empty string if no -i flag
    strcpy(opt->output_dir_name, "output");
    
    // Time parameters
    opt->tStep           = 0.; // This is the numerical time step, might be calculated later
    opt->totalTime       = 1.0e6;
    opt->outputFrequency = 1000.0;
    opt->startTime       = 0.0;

    // Init tool specific parameters' defaults
    opt->md_val = 0.01;
    opt->eps_val = 0.01;
    opt->ratio_val = 0.85;
    opt->mic_val = 1e-4;
    opt->onesize_val = 1.0;
    printf("DEBUG [create_default_options]: Default options setting complete.\n");
}

/* --- parse_options: Parses command-line arguments and fills the 'options_t' struct --- */
int parse_options(int argc, const char **argv, options_t *opt){
    printf("DEBUG [parse_options]: Parsing command-line arguments (%d total).\n", argc);
    int i = 1;

    while (i < argc) {
        printf("DEBUG [parse_options]: Processing argument %d: %s\n", i, argv[i]);
        if(strcmp(argv[i], "-drift") == 0) {
            i++;
            if (i < argc) opt->drift = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -drift.\n"); return 1; }
            printf("DEBUG [parse_options]:   -drift set to %.2f\n", opt->drift);
        }
        else if (strcmp(argv[i], "-growth") == 0) {
            i++;
            if (i < argc) opt->growth = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -growth.\n"); return 1; }
            printf("DEBUG [parse_options]:   -growth set to %.2f\n", opt->growth);
        }
        else if (strcmp(argv[i], "-evol") == 0) {
            i++;
            if (i < argc) opt->evol = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -evol.\n"); return 1; }
            printf("DEBUG [parse_options]:   -evol set to %.2f\n", opt->evol);
        }
        else if (strcmp(argv[i], "-twopop") == 0) {
            i++;
            if (i < argc) opt->twopop = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -twopop.\n"); return 1; }
            printf("DEBUG [parse_options]:   -twopop set to %.2f\n", opt->twopop);
        }
        else if (strcmp(argv[i], "-ufrag") == 0) {
            i++;
            if (i < argc) opt->ufrag = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -ufrag.\n"); return 1; }
            printf("DEBUG [parse_options]:   -ufrag set to %.2f\n", opt->ufrag);
        }
        else if (strcmp(argv[i], "-ffrag") == 0) {
            i++;
            if (i < argc) opt->ffrag = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -ffrag.\n"); return 1; }
            printf("DEBUG [parse_options]:   -ffrag set to %.2f\n", opt->ffrag);
        }
        else if (strcmp(argv[i], "-tStep") == 0) {
            i++;
            if (i < argc) opt->tStep = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -tStep.\n"); return 1; }
            printf("DEBUG [parse_options]:   -tStep set to %.2e\n", opt->tStep);
        }
        else if (strcmp(argv[i], "-n") == 0) { // Main simulation NGRID
            i++;
            if (i < argc) opt->ngrid_val = atoi(argv[i]); else { fprintf(stderr, "Error: Missing value for -n.\n"); return 1; }
            printf("DEBUG [parse_options]:   -n (NGRID) set to %d\n", opt->ngrid_val);
        }
        else if (strcmp(argv[i], "-i") == 0) {
            i++;
            if (i < argc) opt->input_file = argv[i]; else { fprintf(stderr, "Error: Missing value for -i.\n"); return 1; }
            printf("DEBUG [parse_options]:   -i (input_file) set to '%s'\n", opt->input_file);
        }
        else if (strcmp(argv[i], "-o") == 0) { // Output directory name
            i++;
            if (i < argc) { // Check if there is an argument after -o
                strcpy(opt->output_dir_name, argv[i]); // <--- CORRECTED: Use strcpy
            } else { // No argument after -o
            fprintf(stderr, "Error: Missing value for -o.\n");
            return 1;
            }
        }
        else if (strcmp(argv[i], "-tmax") == 0) {
            i++;
            if (i < argc) opt->totalTime = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -tmax.\n"); return 1; }
            printf("DEBUG [parse_options]:   -tmax set to %.2e\n", opt->totalTime);
        }
        else if (strcmp(argv[i], "-outfreq") == 0) {
            i++;
            if (i < argc) opt->outputFrequency = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -outfreq.\n"); return 1; }
            printf("DEBUG [parse_options]:   -outfreq set to %.2e\n", opt->outputFrequency);
        }
        else if (strcmp(argv[i], "-curr") == 0) {
            i++;
            if (i < argc) opt->startTime = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -curr.\n"); return 1; }
            printf("DEBUG [parse_options]:   -curr (startTime) set to %.2e\n", opt->startTime);
        }
        // --- Init_tool specific options processed in the main parser ---
        else if (strcmp(argv[i], "-n_init") == 0) { i++; if (i < argc) opt->ngrid_val = atoi(argv[i]); else { fprintf(stderr, "Error: Missing value for -n_init.\n"); return 1; } printf("DEBUG [parse_options]:   -n_init (NGRID for init) set to %d\n", opt->ngrid_val); }
        else if (strcmp(argv[i], "-ri") == 0) { i++; if (i < argc) opt->rmin_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -ri.\n"); return 1; } printf("DEBUG [parse_options]:   -ri (RMIN for init) set to %.2f\n", opt->rmin_val); }
        else if (strcmp(argv[i], "-ro") == 0) { i++; if (i < argc) opt->rmax_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -ro.\n"); return 1; } printf("DEBUG [parse_options]:   -ro (RMAX for init) set to %.2f\n", opt->rmax_val); }
        else if (strcmp(argv[i], "-sigma0_init") == 0) { i++; if (i < argc) opt->sigma0_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -sigma0_init.\n"); return 1; } printf("DEBUG [parse_options]:   -sigma0_init set to %.2e\n", opt->sigma0_val); }
        else if (strcmp(argv[i], "-index_init") == 0) { i++; if (i < argc) opt->sigmap_exp_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -index_init.\n"); return 1; } printf("DEBUG [parse_options]:   -index_init set to %.2f\n", opt->sigmap_exp_val); }
        else if (strcmp(argv[i], "-rdzei") == 0) { i++; if (i < argc) opt->r_dze_i_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -rdzei.\n"); return 1; } printf("DEBUG [parse_options]:   -rdzei set to %.2f\n", opt->r_dze_i_val); }
        else if (strcmp(argv[i], "-rdzeo") == 0) { i++; if (i < argc) opt->r_dze_o_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -rdzeo.\n"); return 1; } printf("DEBUG [parse_options]:   -rdzeo set to %.2f\n", opt->r_dze_o_val); }
        else if (strcmp(argv[i], "-drdzei") == 0) { i++; if (i < argc) opt->dr_dze_i_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -drdzei.\n"); return 1; } printf("DEBUG [parse_options]:   -drdzei set to %.2f\n", opt->dr_dze_i_val); }
        else if (strcmp(argv[i], "-drdzeo") == 0) { i++; if (i < argc) opt->dr_dze_o_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -drdzeo.\n"); return 1; } printf("DEBUG [parse_options]:   -drdzeo set to %.2f\n", opt->dr_dze_o_val); }
        else if (strcmp(argv[i], "-alpha_init") == 0) { i++; if (i < argc) opt->alpha_visc_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -alpha_init.\n"); return 1; } printf("DEBUG [parse_options]:   -alpha_init set to %.2e\n", opt->alpha_visc_val); }
        else if (strcmp(argv[i], "-amod") == 0) { i++; if (i < argc) opt->a_mod_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -amod.\n"); return 1; } printf("DEBUG [parse_options]:   -amod set to %.2f\n", opt->a_mod_val); }
        else if (strcmp(argv[i], "-h_init") == 0) { i++; if (i < argc) opt->hasp_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -h_init.\n"); return 1; } printf("DEBUG [parse_options]:   -h_init set to %.2f\n", opt->hasp_val); }
        else if (strcmp(argv[i], "-flind_init") == 0) { i++; if (i < argc) opt->flind_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -flind_init.\n"); return 1; } printf("DEBUG [parse_options]:   -flind_init set to %.2f\n", opt->flind_val); }
        else if (strcmp(argv[i], "-m0_init") == 0) { i++; if (i < argc) opt->star_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -m0_init.\n"); return 1; } printf("DEBUG [parse_options]:   -m0_init set to %.2f\n", opt->star_val); }
        else if (strcmp(argv[i], "-md_init") == 0) { i++; if (i < argc) opt->md_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -md_init.\n"); return 1; } printf("DEBUG [parse_options]:   -md_init set to %.2e\n", opt->md_val); }
        else if (strcmp(argv[i], "-eps_init") == 0) { i++; if (i < argc) opt->eps_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -eps_init.\n"); return 1; } printf("DEBUG [parse_options]:   -eps_init set to %Le\n", opt->eps_val); }
        else if (strcmp(argv[i], "-ratio_init") == 0) { i++; if (i < argc) opt->ratio_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -ratio_init.\n"); return 1; } printf("DEBUG [parse_options]:   -ratio_init set to %Le\n", opt->ratio_val); }
        else if (strcmp(argv[i], "-onesize_init") == 0) { i++; if (i < argc) opt->onesize_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -onesize_init.\n"); return 1; } printf("DEBUG [parse_options]:   -onesize_init set to %Le\n", opt->onesize_val); }
        else if (strcmp(argv[i], "-mic_init") == 0) { i++; if (i < argc) opt->mic_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -mic_init.\n"); return 1; } printf("DEBUG [parse_options]:   -mic_init set to %Le\n", opt->mic_val); }
        // --- End of init_tool option processing ---
        else {
            fprintf(stderr, "ERROR [parse_options]: Invalid switch on command-line: %s!\n", argv[i]);
            printf("\n\n**** Invalid switch on command-line: %s! ****\n\n", argv[i]);
            printf("**** Try following parameters: ****\n\n-drift <val>\n-growth <val>\n-evol <val>\n-twopop <val>\n-ufrag <val>\n-ffrag <val>\n-tStep <val>\n-n <val>\n-i <filename>\n-o <output_directory_name>\n-tmax <val>\n-outfreq <val>\n-curr <val>\n");
            printf("\n**** Initial profile parameters (used if -i is not provided): ****\n-n_init <val>\n-ri <val>\n-ro <val>\n-sigma0_init <val>\n-index_init <val>\n-rdzei <val>\n-rdzeo <val>\n-drdzei <val>\n-drdzeo <val>\n-alpha_init <val>\n-amod <val>\n-h_init <val>\n-flind_init <val>\n-m0_init <val>\n-md_init <val>\n-eps_init <val>\n-ratio_init <val>\n-onesize_init <val>\n-mic_init <val>\n");
            return 1;
        }
        i++;
    }

    printf("DEBUG [parse_options]: Command-line parsing complete.\n");
    return 0;
}