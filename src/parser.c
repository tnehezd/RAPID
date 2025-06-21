// parser.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h> // For true/false

#include "parser.h" // Include its own header
#include "simulation_types.h" // In case simulation_types.h defines something else needed here

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
    // ADJUST THESE DEFAULTS FOR REALISTIC VALUES, AS DISCUSSED PREVIOUSLY!
    opt->sigma0_val      = 1.0; // Consider changing this default (e.g., 1.0e3)
    opt->sigmap_exp_val  = 1.0; // Consider changing this default (e.g., 1.5 for r^-1.5)
    opt->alpha_visc_val  = 0.01;
    opt->star_val        = 1.0;
    opt->hasp_val        = 0.05;
    opt->flind_val       = 0.5;
    opt->r_dze_i_val     = 0.0;
    opt->r_dze_o_val     = 0.0;
    opt->dr_dze_i_val    = 0.0;
    opt->dr_dze_o_val    = 0.0;
    opt->a_mod_val       = 0.0;

    // File input/output
    opt->input_file      = NULL; // Default to NULL if no -i flag, easier to check
    strcpy(opt->output_dir_name, "output");
    
    // Time parameters
    opt->tStep           = 0.;
    opt->totalTime       = 1.0e6;
    opt->outputFrequency = 1000.0;
    opt->startTime       = 0.0;

    // Init tool specific parameters' defaults
    // ADJUST THESE DEFAULTS FOR REALISTIC VALUES, AS DISCUSSED PREVIOUSLY!
    opt->eps_val = 0.01; // e.g., 0.01 for 1% dust-to-gas ratio
    opt->ratio_val = 0.85; // e.g., 0.85 for 85% Pop1
    opt->mic_val = 1e-4; // e.g., 1e-4 for 100 micron (0.01 cm)
    opt->onesize_val = 0.0; // 0.0 for size distribution, 1.0 for single size (e.g., mic_val)
    printf("DEBUG [create_default_options]: Default options setting complete.\n");
}

/* --- print_usage: Prints a help message for command-line arguments --- */
void print_usage() {
    fprintf(stderr, "Usage: simulation [OPTIONS]\n");
    fprintf(stderr, "Simulation control options:\n");
    fprintf(stderr, "  -drift <val>   Enable/disable drift (0. or 1., default: 1.0)\n");
    fprintf(stderr, "  -growth <val>  Enable/disable growth (0. or 1., default: 1.0)\n");
    fprintf(stderr, "  -evol <val>    Enable/disable evolution (0. or 1., default: 1.0)\n");
    fprintf(stderr, "  -twopop <val>  Enable/disable two dust populations (0. or 1., default: 1.0)\n");
    fprintf(stderr, "  -ufrag <val>   Fragmentation velocity (default: 1000.0)\n");
    fprintf(stderr, "  -ffrag <val>   Fragmentation factor (default: 0.37)\n");
    fprintf(stderr, "Time parameters:\n");
    fprintf(stderr, "  -tStep <val>   Fixed time step (default: 0.0)\n");
    fprintf(stderr, "  -tmax <val>    Total simulation time (default: 1.0e6)\n");
    fprintf(stderr, "  -outfreq <val> Output frequency (default: 1000.0)\n");
    fprintf(stderr, "  -curr <val>    Current start time (default: 0.0)\n");
    fprintf(stderr, "File I/O:\n");
    fprintf(stderr, "  -i <file>      Input profile file (e.g., init_data.dat)\n");
    fprintf(stderr, "  -o <dir>       Output directory name (default: 'output')\n");
    fprintf(stderr, "Initial profile generation options (used if -i is not provided):\n");
    fprintf(stderr, "  -n <val>       Number of grid points (default: 2000)\n"); // This is common for sim and init
    fprintf(stderr, "  -ri <val>      Inner radius (AU, default: 1.0)\n");
    fprintf(stderr, "  -ro <val>      Outer radius (AU, default: 100.0)\n");
    fprintf(stderr, "  -sigma0_init <val> Initial gas surface density at 1 AU (M_sun/AU^2, default: 1.0)\n");
    fprintf(stderr, "  -index_init <val> Exponent of surface density profile (positive value, default: 1.0 for r^-1)\n");
    fprintf(stderr, "  -alpha_init <val> Alpha viscosity (default: 0.01)\n");
    fprintf(stderr, "  -m0_init <val> Star mass (M_sun, default: 1.0)\n");
    fprintf(stderr, "  -h_init <val>  Aspect ratio at 1 AU (H/R, default: 0.05)\n");
    fprintf(stderr, "  -flind_init <val> Flaring index (default: 0.5)\n");
    fprintf(stderr, "  -rdzei <val>   Inner dead zone radius (AU, default: 0.0)\n");
    fprintf(stderr, "  -rdzeo <val>   Outer dead zone radius (AU, default: 0.0)\n");
    fprintf(stderr, "  -drdzei <val>  Inner dead zone transition width multiplier (default: 0.0)\n");
    fprintf(stderr, "  -drdzeo <val>  Outer dead zone transition width multiplier (default: 0.0)\n");
    fprintf(stderr, "  -amod <val>    Alpha viscosity multiplier in dead zone (default: 0.0)\n");
    fprintf(stderr, "  -eps <val>     Dust-to-gas ratio (default: 0.01)\n");
    fprintf(stderr, "  -ratio <val>   Ratio of Pop1 dust mass to total dust mass (default: 0.85)\n");
    fprintf(stderr, "  -mic <val>     Micro-sized particle radius (cm, default: 1e-4)\n");
    fprintf(stderr, "  -onesize <val> Use one size particles (0 for distribution, 1 for mic_val, default: 0.0)\n");
    fprintf(stderr, "Other:\n");
    fprintf(stderr, "  -h or --help   Display this help message\n");
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
        else if (strcmp(argv[i], "-n") == 0) { // Main simulation NGRID (and also init_tool's NGRID)
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
            if (i < argc) {
                // Ensure output_dir_name has enough space and use strncpy for safety
                strncpy(opt->output_dir_name, argv[i], sizeof(opt->output_dir_name) - 1);
                opt->output_dir_name[sizeof(opt->output_dir_name) - 1] = '\0'; // Ensure null termination
                printf("DEBUG [parse_options]:   -o (output_dir_name) set to '%s'\n", opt->output_dir_name);
            } else {
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
        // Note: -n is already handled above for both sim and init
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
        else if (strcmp(argv[i], "-eps") == 0) { i++; if (i < argc) opt->eps_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -eps.\n"); return 1; } printf("DEBUG [parse_options]:   -eps set to %.2e\n", opt->eps_val); }
        else if (strcmp(argv[i], "-ratio") == 0) { i++; if (i < argc) opt->ratio_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -ratio.\n"); return 1; } printf("DEBUG [parse_options]:   -ratio set to %.2e\n", opt->ratio_val); }
        else if (strcmp(argv[i], "-mic") == 0) { i++; if (i < argc) opt->mic_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -mic.\n"); return 1; } printf("DEBUG [parse_options]:   -mic set to %.2e\n", opt->mic_val); }
        else if (strcmp(argv[i], "-onesize") == 0) { i++; if (i < argc) opt->onesize_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -onesize.\n"); return 1; } printf("DEBUG [parse_options]:   -onesize set to %.2e\n", opt->onesize_val); }
        else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage();
            return 1; // Exit after printing usage
        }
        else {
            fprintf(stderr, "ERROR [parse_options]: Invalid switch on command-line: %s!\n", argv[i]);
            print_usage(); // Print usage on unknown option
            return 1;
        }
        i++;
    }

    printf("DEBUG [parse_options]: Command-line parsing complete.\n");
    return 0;
}