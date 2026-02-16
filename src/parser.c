// src/parser.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h> 
#include "parser.h"
#include "simulation_types.h" 

void createDefaultOptions(ParserOptions *opt) {
    fprintf(stderr, "DEBUG [createDefaultOptions]: Setting default values for ParserOptions.\n");
    // Simulation control options
    opt->option_for_dust_drift                      = 1.;
    opt->option_for_dust_growth                     = 1.;
    opt->option_for_evolution                       = 1.;
    opt->option_for_dust_secondary_population       = 1.;
    opt->fragmenatation_velocity                    = 1000.0;
    opt->fragmenatation_factor                      = 0.37;
    opt->number_of_grid_points                      = 2000;
    opt->number_of_vertical_grid_points             = 100;
    opt->number_of_dust_particles                   = 5000;
    opt->rmin_val                                   = 1.0;
    opt->rmax_val                                   = 100.0;
    opt->zmax_val                                   = 4.0;
    opt->sigma0_val                                 = 1.0; 
    opt->sigmap_exp_val                             = 0.5; 
    opt->alpha_visc_val                             = 0.01;
    opt->star_val                                   = 1.0;
    opt->hasp_val                                   = 0.05;
    opt->flind_val                                  = 0.5;
    opt->r_dze_i_val                                = 0.0;
    opt->r_dze_o_val                                = 0.0;
    opt->dr_dze_i_val                               = 0.0;
    opt->dr_dze_o_val                               = 0.0;
    opt->a_mod_val                                  = 0.0;
    opt->input_file                                 = NULL;

    strncpy(opt->output_dir_name, "output", sizeof(opt->output_dir_name) - 1);
    opt->output_dir_name[sizeof(opt->output_dir_name) - 1] = '\0'; 

    opt->user_defined_time_step                     = 0.;
    opt->maximum_simulation_time                    = 1.0e6;
    opt->output_frequency                           = 1000.0;
    opt->eps_val                                    = 0.01; 
    opt->ratio_val                                  = 0.85;
    opt->mic_val                                    = 1e-4; 
    opt->onesize_val                                = 0.0; 
    opt->pdensity_val                               = 1.6; 
    opt->output_format                              = OUTPUT_ASCII;
    opt->disk_dimension = 1;   // default = 1D

    fprintf(stderr, "Default options setting complete.\n");
}

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
    fprintf(stderr, "  -n <val>       Number of radial grid points (default: 2000)\n");
    fprintf(stderr, "  -nz <val>      Number of vertical grid points (default: 100)\n");
    fprintf(stderr, "  -ri <val>      Inner radius (AU, default: 1.0)\n");
    fprintf(stderr, "  -ro <val>      Outer radius (AU, default: 100.0)\n");
    fprintf(stderr, "  -sigma0_init <val> Initial gas surface density at 1 AU (M_sun/AU^2, default: 1.0)\n");
    fprintf(stderr, "  -index_init <val> Exponent of surface density profile (positive value, default: 0.5 for r^-0.5)\n"); 
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
    fprintf(stderr, "Other:\n");
    fprintf(stderr, "  -pdensity <val> Dust particle density (g/cm^3, default: 1.6)\n"); 
    fprintf(stderr, "  -h or --help   Display this help message\n");
    fprintf(stderr, "Output format options:\n");
    fprintf(stderr, "  --output-format <ascii|hdf5>  Select output format (default: ascii)\n");
    fprintf(stderr, "  -ascii                        Shortcut for --output-format ascii\n");
    fprintf(stderr, "  -hdf5                       Shortcut for --output-format hdf5\n");
    fprintf(stderr, "Geometry:\n");
    fprintf(stderr, "  -dim <1|2>    Disk dimensionality (1=radial, 2=radial+vertical, default: 1)\n");
}

int parseCLIOptions(int argc, const char **argv, ParserOptions *opt){

    fprintf(stderr, "DEBUG [parseCLIOptions]: Parsing command-line arguments (%d total).\n", argc);
    int i = 1;

    while (i < argc) {
        if(strcmp(argv[i], "-drift") == 0) {
            i++;
            if (i < argc) opt->option_for_dust_drift = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -drift.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-growth") == 0) {
            i++;
            if (i < argc) opt->option_for_dust_growth = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -growth.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-evol") == 0) {
            i++;
            if (i < argc) opt->option_for_evolution = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -evol.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-twopop") == 0) {
            i++;
            if (i < argc) opt->option_for_dust_secondary_population = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -twopop.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-ufrag") == 0) {
            i++;
            if (i < argc) opt->fragmenatation_velocity = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -ufrag.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-ffrag") == 0) {
            i++;
            if (i < argc) opt->fragmenatation_factor = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -ffrag.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-tStep") == 0) {
            i++;
            if (i < argc) opt->user_defined_time_step = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -tStep.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-n") == 0) {
            i++;
            if (i < argc) opt->number_of_grid_points = atoi(argv[i]); else { fprintf(stderr, "Error: Missing value for -n.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-nz") == 0) {
            i++;
            if (i < argc) opt->number_of_vertical_grid_points = atoi(argv[i]); else { fprintf(stderr, "Error: Missing value for -nz.\n"); return 1; }
        }        
        else if (strcmp(argv[i], "-ndust") == 0) { 
            i++;
            if (i < argc) opt->number_of_dust_particles = atoi(argv[i]); else { fprintf(stderr, "Error: Missing value for -ndust.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-i") == 0) {
            i++;
            if (i < argc) opt->input_file = argv[i]; else { fprintf(stderr, "Error: Missing value for -i.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-o") == 0) { // Output directory name
            i++;
            if (i < argc) {
                strncpy(opt->output_dir_name, argv[i], sizeof(opt->output_dir_name) - 1);
                opt->output_dir_name[sizeof(opt->output_dir_name) - 1] = '\0';
            } else {
                fprintf(stderr, "Error: Missing value for -o.\n");
                return 1;
            }
        }
        else if (strcmp(argv[i], "-tmax") == 0) {
            i++;
            if (i < argc) opt->maximum_simulation_time = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -tmax.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-outfreq") == 0) {
            i++;
            if (i < argc) opt->output_frequency = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -outfreq.\n"); return 1; }
        }

        else if (strcmp(argv[i], "-ri") == 0) { 
            i++; 
            if (i < argc) opt->rmin_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -ri.\n"); return 1; }; 
        } 
        else if (strcmp(argv[i], "-ro") == 0) { 
            i++; 
            if (i < argc) opt->rmax_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -ro.\n"); return 1; }; 
        } 
        else if (strcmp(argv[i], "-zmax") == 0) { 
            i++; 
            if (i < argc) opt->zmax_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -zmax.\n"); return 1; }; 
        } 
        else if (strcmp(argv[i], "-sigma0_init") == 0) {
            i++; 
            if (i < argc) opt->sigma0_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -sigma0_init.\n"); return 1; }; 
        }
        else if (strcmp(argv[i], "-index_init") == 0) { 
            i++; 
            if (i < argc) opt->sigmap_exp_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -index_init.\n"); return 1; }; 
        } 
        else if (strcmp(argv[i], "-rdzei") == 0) { 
            i++; 
            if (i < argc) opt->r_dze_i_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -rdzei.\n"); return 1; }; 
        }
        else if (strcmp(argv[i], "-rdzeo") == 0) { 
            i++; 
            if (i < argc) opt->r_dze_o_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -rdzeo.\n"); return 1; }; 
        }
        else if (strcmp(argv[i], "-drdzei") == 0) { 
            i++; 
            if (i < argc) opt->dr_dze_i_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -drdzei.\n"); return 1; }; 
        }
        else if (strcmp(argv[i], "-drdzeo") == 0) { 
            i++; 
            if (i < argc) opt->dr_dze_o_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -drdzeo.\n"); return 1; }; 
        }
        else if (strcmp(argv[i], "-alpha_init") == 0) { 
            i++; 
            if (i < argc) opt->alpha_visc_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -alpha_init.\n"); return 1; }; 
        } 
        else if (strcmp(argv[i], "-amod") == 0) { 
            i++; 
            if (i < argc) opt->a_mod_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -amod.\n"); return 1; }; 
        }
        else if (strcmp(argv[i], "-h_init") == 0) { 
            i++; 
            if (i < argc) opt->hasp_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -h_init.\n"); return 1; }; 
        } 
        else if (strcmp(argv[i], "-flind_init") == 0) { 
            i++; 
            if (i < argc) opt->flind_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -flind_init.\n"); return 1; }; 
        }
        else if (strcmp(argv[i], "-m0_init") == 0) { 
            i++; 
            if (i < argc) opt->star_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -m0_init.\n"); return 1; }; 
        } 
        else if (strcmp(argv[i], "-eps") == 0) { 
            i++; 
            if (i < argc) opt->eps_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -eps.\n"); return 1; }; 
        } 
        else if (strcmp(argv[i], "-ratio") == 0) { 
            i++; 
            if (i < argc) opt->ratio_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -ratio.\n"); return 1; }; 
        } 
        else if (strcmp(argv[i], "-mic") == 0) { 
            i++; 
            if (i < argc) opt->mic_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -mic.\n"); return 1; }; 
        } 
        else if (strcmp(argv[i], "-onesize") == 0) { 
            i++; 
            if (i < argc) opt->onesize_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -onesize.\n"); return 1; }; 
        }
        else if (strcmp(argv[i], "-pdensity") == 0) {
            i++;
            if (i < argc) opt->pdensity_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -pdensity.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
//            printUsageToTerminal();
            return 1;
        }

        else if (strcmp(argv[i], "--output-format") == 0) {
            i++;
            if (i < argc) {
                if (strcmp(argv[i], "ascii") == 0)
                    opt->output_format = OUTPUT_ASCII;
                else if (strcmp(argv[i], "hdf5") == 0)
                    opt->output_format = OUTPUT_HDF5;
                else {
                    fprintf(stderr, "Error: Unknown output format '%s'. Use ascii or hdf5.\n", argv[i]);
                    return 1;
                }
            } else {
                fprintf(stderr, "Error: Missing value for --output-format.\n");
                return 1;
            }
        }
        else if (strcmp(argv[i], "-ascii") == 0) {
            opt->output_format = OUTPUT_ASCII;
        }
        else if (strcmp(argv[i], "-hdf5") == 0) {
            opt->output_format = OUTPUT_HDF5;
        }
        else if (strcmp(argv[i], "-dim") == 0) {
            i++;
            if (i < argc) {
                opt->disk_dimension = atoi(argv[i]);
                if (opt->disk_dimension != 1 && opt->disk_dimension != 2) {
                    fprintf(stderr, "Error: -dim must be 1 or 2.\n");
                    return 1;
                }
            } else {
                fprintf(stderr, "Error: Missing value for -dim.\n");
                return 1;
            }
        }
        else {
            fprintf(stderr, "ERROR [parseCLIOptions]: Invalid switch on command-line: %s!\n", argv[i]);
  //          printUsageToTerminal(); 
            return 1;
        }
        i++;
    }

    fprintf(stderr, "DEBUG [parseCLIOptions]: Command-line parsing complete.\n");
    return 0;
}