// parser.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h> // For true/false

#include "parser.h" // Include its own header
#include "simulation_types.h" // In case simulation_types.h defines something else needed here

/* --- createDefaultOptions: Sets default values for the 'ParserOptions' structure --- */
void createDefaultOptions(ParserOptions *opt) {
    fprintf(stderr, "DEBUG [createDefaultOptions]: Setting default values for ParserOptions.\n");
    // Simulation control options
    opt->option_for_dust_drift           = 1.;
    opt->option_for_dust_growth          = 1.;
    opt->option_for_evolution            = 1.;
    opt->option_for_dust_secondary_population          = 1.;
    opt->ufrag           = 1000.0;
    opt->ffrag           = 0.37;

    // Core disk parameters (also serve as init_tool defaults)
    opt->ngrid_val       = 2000;
    opt->ndust_val       = 5000;
    opt->rmin_val        = 1.0;
    opt->rmax_val        = 100.0;
    // ADJUST THESE DEFAULTS FOR REALISTIC VALUES, AS DISCUSSED PREVIOUSLY!
    opt->sigma0_val      = 1.0; // Consider changing this default (e.g., 1.0e3)
    opt->sigmap_exp_val  = 0.5; // KORRIGÁLVA: Pozitív értékre, pl. 0.5 a Sigma ~ r^-0.5 esetén
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
    strncpy(opt->output_dir_name, "output", sizeof(opt->output_dir_name) - 1);
    opt->output_dir_name[sizeof(opt->output_dir_name) - 1] = '\0'; // Ensure null termination

    // Time parameters
    opt->user_defined_time_step           = 0.;
    opt->maximum_simulation_time       = 1.0e6;
    opt->output_frequency = 1000.0;

    // Init tool specific parameters' defaults
    // ADJUST THESE DEFAULTS FOR REALISTIC VALUES, AS DISCUSSED PREVIOUSLY!
    opt->eps_val = 0.01; // e.g., 0.01 for 1% dust-to-gas ratio
    opt->ratio_val = 0.85; // e.g., 0.85 for 85% Pop1
    opt->mic_val = 1e-4; // e.g., 1e-4 for 100 micron (0.01 cm)
    opt->onesize_val = 0.0; // 0.0 for size distribution, 1.0 for single size (e.g., mic_val)

    // NEW: particle_density (dust particle density) default value
    opt->pdensity_val = 1.6; // Por sűrűsége (g/cm^3) - ez az alapértelmezett 1.6
    fprintf(stderr, "Default options setting complete.\n");
}

/* --- printUsageToTerminal: Prints a help message for command-line arguments --- */
void printUsageToTerminal() {
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
    fprintf(stderr, "File I/O:\n");
    fprintf(stderr, "  -i <file>      Input profile file (e.g., init_data.dat)\n");
    fprintf(stderr, "  -o <dir>       Output directory name (default: 'output')\n");
    fprintf(stderr, "Initial profile generation options (used if -i is not provided):\n");
    fprintf(stderr, "  -n <val>       Number of grid points (default: 2000)\n"); // This is common for sim and init
    fprintf(stderr, "  -ri <val>      Inner radius (AU, default: 1.0)\n");
    fprintf(stderr, "  -ro <val>      Outer radius (AU, default: 100.0)\n");
    fprintf(stderr, "  -sigma0_init <val> Initial gas surface density at 1 AU (M_sun/AU^2, default: 1.0)\n");
    fprintf(stderr, "  -index_init <val> Exponent of surface density profile (positive value, default: 0.5 for r^-0.5)\n"); // KORRIGÁLVA a súgóban
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
    fprintf(stderr, "  -ndust <val>   The number of the particles, default: 5000)\n");
    // NEW: particle_density usage message
    fprintf(stderr, "Other:\n");
    fprintf(stderr, "  -pdensity <val> Dust particle density (g/cm^3, default: 1.6)\n"); // NEW usage line
    fprintf(stderr, "  -h or --help   Display this help message\n");
}


/* --- parseCLIOptions: Parses command-line arguments and fills the 'ParserOptions' struct --- */
int parseCLIOptions(int argc, const char **argv, ParserOptions *opt){
    fprintf(stderr, "DEBUG [parseCLIOptions]: Parsing command-line arguments (%d total).\n", argc);
    int i = 1;

    while (i < argc) {
//        fprintf(stderr, "DEBUG [parseCLIOptions]: Processing argument %d: %s\n", i, argv[i]);
        if(strcmp(argv[i], "-drift") == 0) {
            i++;
            if (i < argc) opt->option_for_dust_drift = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -drift.\n"); return 1; }
//            fprintf(stderr, "DEBUG [parseCLIOptions]:   -drift set to %.2f\n", opt->option_for_dust_drift);
        }
        else if (strcmp(argv[i], "-growth") == 0) {
            i++;
            if (i < argc) opt->option_for_dust_growth = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -growth.\n"); return 1; }
//            fprintf(stderr, "DEBUG [parseCLIOptions]:   -growth set to %.2f\n", opt->option_for_dust_growth);
        }
        else if (strcmp(argv[i], "-evol") == 0) {
            i++;
            if (i < argc) opt->option_for_evolution = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -evol.\n"); return 1; }
//            fprintf(stderr, "DEBUG [parseCLIOptions]:   -evol set to %.2f\n", opt->option_for_evolution);
        }
        else if (strcmp(argv[i], "-twopop") == 0) {
            i++;
            if (i < argc) opt->option_for_dust_secondary_population = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -twopop.\n"); return 1; }
//            fprintf(stderr, "DEBUG [parseCLIOptions]:   -twopop set to %.2f\n", opt->option_for_dust_secondary_population);
        }
        else if (strcmp(argv[i], "-ufrag") == 0) {
            i++;
            if (i < argc) opt->ufrag = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -ufrag.\n"); return 1; }
//            fprintf(stderr, "DEBUG [parseCLIOptions]:   -ufrag set to %.2f\n", opt->ufrag);
        }
        else if (strcmp(argv[i], "-ffrag") == 0) {
            i++;
            if (i < argc) opt->ffrag = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -ffrag.\n"); return 1; }
//            fprintf(stderr, "DEBUG [parseCLIOptions]:   -ffrag set to %.2f\n", opt->ffrag);
        }
        else if (strcmp(argv[i], "-tStep") == 0) {
            i++;
            if (i < argc) opt->user_defined_time_step = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -tStep.\n"); return 1; }
//            fprintf(stderr, "DEBUG [parseCLIOptions]:   -tStep set to %.2e\n", opt->user_defined_time_step);
        }
        else if (strcmp(argv[i], "-n") == 0) { // Main simulation NGRID (and also init_tool's NGRID)
            i++;
            if (i < argc) opt->ngrid_val = atoi(argv[i]); else { fprintf(stderr, "Error: Missing value for -n.\n"); return 1; }
//            fprintf(stderr, "DEBUG [parseCLIOptions]:   -n (NGRID) set to %d\n", opt->ngrid_val);
        }
        else if (strcmp(argv[i], "-ndust") == 0) { 
            i++;
            if (i < argc) opt->ndust_val = atoi(argv[i]); else { fprintf(stderr, "Error: Missing value for -ndust.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-i") == 0) {
            i++;
            if (i < argc) opt->input_file = argv[i]; else { fprintf(stderr, "Error: Missing value for -i.\n"); return 1; }
//            fprintf(stderr, "DEBUG [parseCLIOptions]:   -i (input_file) set to '%s'\n", opt->input_file);
        }
        else if (strcmp(argv[i], "-o") == 0) { // Output directory name
            i++;
            if (i < argc) {
                // Ensure output_dir_name has enough space and use strncpy for safety
                strncpy(opt->output_dir_name, argv[i], sizeof(opt->output_dir_name) - 1);
                opt->output_dir_name[sizeof(opt->output_dir_name) - 1] = '\0'; // Ensure null termination
//                fprintf(stderr, "DEBUG [parseCLIOptions]:   -o (output_dir_name) set to '%s'\n", opt->output_dir_name);
            } else {
                fprintf(stderr, "Error: Missing value for -o.\n");
                return 1;
            }
        }
        else if (strcmp(argv[i], "-tmax") == 0) {
            i++;
            if (i < argc) opt->maximum_simulation_time = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -tmax.\n"); return 1; }
//            fprintf(stderr, "DEBUG [parseCLIOptions]:   -tmax set to %.2e\n", opt->maximum_simulation_time);
        }
        else if (strcmp(argv[i], "-outfreq") == 0) {
            i++;
            if (i < argc) opt->output_frequency = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -outfreq.\n"); return 1; }
//            fprintf(stderr, "DEBUG [parseCLIOptions]:   -outfreq set to %.2e\n", opt->output_frequency);
        }
        // --- Init_tool specific options processed in the main parser ---
        // Note: -n is already handled above for both sim and init
        else if (strcmp(argv[i], "-ri") == 0) { i++; if (i < argc) opt->rmin_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -ri.\n"); return 1; }; } //fprintf(stderr, "DEBUG [parseCLIOptions]:   -ri (RMIN for init) set to %.2f\n", opt->rmin_val); }
        else if (strcmp(argv[i], "-ro") == 0) { i++; if (i < argc) opt->rmax_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -ro.\n"); return 1; }; } // fprintf(stderr, "DEBUG [parseCLIOptions]:   -ro (RMAX for init) set to %.2f\n", opt->rmax_val); }
        else if (strcmp(argv[i], "-sigma0_init") == 0) { i++; if (i < argc) opt->sigma0_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -sigma0_init.\n"); return 1; }; } // fprintf(stderr, "DEBUG [parseCLIOptions]:   -sigma0_init set to %.2e\n", opt->sigma0_val); }
        else if (strcmp(argv[i], "-index_init") == 0) { i++; if (i < argc) opt->sigmap_exp_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -index_init.\n"); return 1; }; } // fprintf(stderr, "DEBUG [parseCLIOptions]:   -index_init set to %.2f\n", opt->sigmap_exp_val); }
        else if (strcmp(argv[i], "-rdzei") == 0) { i++; if (i < argc) opt->r_dze_i_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -rdzei.\n"); return 1; }; } // fprintf(stderr, "DEBUG [parseCLIOptions]:   -rdzei set to %.2f\n", opt->r_dze_i_val); }
        else if (strcmp(argv[i], "-rdzeo") == 0) { i++; if (i < argc) opt->r_dze_o_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -rdzeo.\n"); return 1; }; } // fprintf(stderr, "DEBUG [parseCLIOptions]:   -rdzeo set to %.2f\n", opt->r_dze_o_val); }
        else if (strcmp(argv[i], "-drdzei") == 0) { i++; if (i < argc) opt->dr_dze_i_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -drdzei.\n"); return 1; }; } // fprintf(stderr, "DEBUG [parseCLIOptions]:   -drdzei set to %.2f\n", opt->dr_dze_i_val); }
        else if (strcmp(argv[i], "-drdzeo") == 0) { i++; if (i < argc) opt->dr_dze_o_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -drdzeo.\n"); return 1; }; } // fprintf(stderr, "DEBUG [parseCLIOptions]:   -drdzeo set to %.2f\n", opt->dr_dze_o_val); }
        else if (strcmp(argv[i], "-alpha_init") == 0) { i++; if (i < argc) opt->alpha_visc_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -alpha_init.\n"); return 1; }; } // fprintf(stderr, "DEBUG [parseCLIOptions]:   -alpha_init set to %.2e\n", opt->alpha_visc_val); }
        else if (strcmp(argv[i], "-amod") == 0) { i++; if (i < argc) opt->a_mod_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -amod.\n"); return 1; }; } // fprintf(stderr, "DEBUG [parseCLIOptions]:   -amod set to %.2f\n", opt->a_mod_val); }
        else if (strcmp(argv[i], "-h_init") == 0) { i++; if (i < argc) opt->hasp_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -h_init.\n"); return 1; }; } // fprintf(stderr, "DEBUG [parseCLIOptions]:   -h_init set to %.2f\n", opt->hasp_val); }
        else if (strcmp(argv[i], "-flind_init") == 0) { i++; if (i < argc) opt->flind_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -flind_init.\n"); return 1; }; } // fprintf(stderr, "DEBUG [parseCLIOptions]:   -flind_init set to %.2f\n", opt->flind_val); }
        else if (strcmp(argv[i], "-m0_init") == 0) { i++; if (i < argc) opt->star_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -m0_init.\n"); return 1; }; } // fprintf(stderr, "DEBUG [parseCLIOptions]:   -m0_init set to %.2f\n", opt->star_val); }
        else if (strcmp(argv[i], "-eps") == 0) { i++; if (i < argc) opt->eps_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -eps.\n"); return 1; } ; } //fprintf(stderr, "DEBUG [parseCLIOptions]:   -eps set to %.2e\n", opt->eps_val); }
        else if (strcmp(argv[i], "-ratio") == 0) { i++; if (i < argc) opt->ratio_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -ratio.\n"); return 1; }; } // fprintf(stderr, "DEBUG [parseCLIOptions]:   -ratio set to %.2e\n", opt->ratio_val); }
        else if (strcmp(argv[i], "-mic") == 0) { i++; if (i < argc) opt->mic_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -mic.\n"); return 1; }; } // fprintf(stderr, "DEBUG [parseCLIOptions]:   -mic set to %.2e\n", opt->mic_val); }
        else if (strcmp(argv[i], "-onesize") == 0) { i++; if (i < argc) opt->onesize_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -onesize.\n"); return 1; }; } // fprintf(stderr, "DEBUG [parseCLIOptions]:   -onesize set to %.2e\n", opt->onesize_val); }
        // NEW: Handle -pdensity flag
        else if (strcmp(argv[i], "-pdensity") == 0) {
            i++;
            if (i < argc) opt->pdensity_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -pdensity.\n"); return 1; }
//            fprintf(stderr, "DEBUG [parseCLIOptions]:   -pdensity set to %.2f\n", opt->pdensity_val);
        }
        else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
//            printUsageToTerminal();
            return 1; // Exit after printing usage
        }
        else {
            fprintf(stderr, "ERROR [parseCLIOptions]: Invalid switch on command-line: %s!\n", argv[i]);
  //          printUsageToTerminal(); // Print usage on unknown option
            return 1;
        }
        i++;
    }

    fprintf(stderr, "DEBUG [parseCLIOptions]: Command-line parsing complete.\n");
    return 0;
}