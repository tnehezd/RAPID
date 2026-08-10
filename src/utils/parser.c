// src/parser.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h> 
#include "parser.h"
#include "simulation_types.h" 
#include "logger.h"
#include "print_panels.h"

void createDefaultOptions(ParserOptions *opt) {
    LOG_DEBUG("Setting default values for ParserOptions.\n");
    // Simulation control options
    opt->option_for_dust_drift                      = 1.;
    opt->option_for_dust_growth                     = 1.;
    opt->option_for_evolution                       = 1.;
    opt->option_for_dust_secondary_population       = 1.;
    opt->fragmenatation_velocity                    = 1000.0;
    opt->fragmenatation_factor                      = 0.37;
    opt->number_of_grid_points                      = 2000;
    opt->number_of_dust_particles                   = 5000;
    opt->rmin_val                                   = 1.0;
    opt->rmax_val                                   = 100.0;
    opt->sigma0_val                                 = 0.0; 
    opt->total_disk_mass                            = 0.0; 
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
    opt->dust_smoothing_mode                        = 0;
    opt->density_floor                              = 1e-12;  
    opt->dust_density_floor                         = 1e-15;
    
    opt->gaussian_sigma                             = 1.0;      /**< default Gaussian kernel width */
    opt->gaussian_cutoff                            = 3.0;     /**< default cutoff radius (3σ) */

    opt->verbose                                    = 0; // Default: Only errors and warnings

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

    opt->enable_photoevaporation                    = false;       // off by default
    strncpy(opt->photoevaporation_mode, "None", sizeof(opt->photoevaporation_mode) - 1);
    opt->photoevaporation_mode[sizeof(opt->photoevaporation_mode) - 1] = '\0';
    opt->xray_luminosity                            = 1.0e30;      // Standard X-ray luminosity [erg/s]
    
    opt->enable_cutoff                          = 0.0; 
    opt->use_cutoff                             = false; // Default: Normal profile
    opt->r_cutoff                               = 30.0;  
    opt->n_for_cutoff                           = 1.5;

    opt->inner_boundary_condition_type          = 0;   // default: zero-gradient inner
    opt->outer_boundary_condition_type          = 1;   // default: parabolic outer


    opt->dust_smoothing_mode = SMOOTHING_CIC; 
    opt->test_enabled = false;
    strcpy(opt->test_type, "");
    opt->ring_center = 10.0;
    opt->ring_width = 2.0;
    opt->fixed_ring_viscosity = 1e-5; // Default viscosity for ring test
    opt->ring_mass = 1.0; // Default gas surface density at ring center


    LOG_DEBUG("Default options setting complete.\n");
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
    fprintf(stderr, "  -n <val>       Number of grid points (default: 2000)\n"); // This is common for sim and init
    fprintf(stderr, "  -ri <val>      Inner radius (AU, default: 1.0)\n");
    fprintf(stderr, "  -ro <val>      Outer radius (AU, default: 100.0)\n");
    fprintf(stderr, "  -disk_mass <val>   Total gas disk mass in Solar Masses (Alternative to -sigma0_init)\n");
    fprintf(stderr, "  -sigma0_init <val> Initial gas surface density at 1 AU (Alternative to -disk_mass)\n");
    fprintf(stderr, "  -index_init <val> Exponent of surface density profile (positive value, default: 0.5 for r^-0.5)\n"); 
    fprintf(stderr, "  -alpha_init <val> Alpha viscosity (default: 0.01)\n");
    fprintf(stderr, "  -stellar_mass <val> Star mass (M_sun, default: 1.0)\n");
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
    fprintf(stderr, "  -density_floor <val>      Minimum allowed gas surface density (default: 1e-12)\n");  
    fprintf(stderr, "  -dust_density_floor <val> Minimum allowed dust surface density (default: 1e-15)\n");    
    fprintf(stderr, "  -dust_smoothing <cic|ngp|tophat|gaussian>  Select dust smoothing method\n");
    fprintf(stderr, "                                             for mapping Lagrangian particles to the Euler grid\n");
    fprintf(stderr, "  -gaussian_sigma_grid <val>   Width of the Gaussian kernel in grid-cell units (Δr).\n");
    fprintf(stderr, "                               Controls how far dust mass spreads around a particle.\n");
    fprintf(stderr, "  -gaussian_cutoff <val>       Kernel cutoff expressed in multiples of sigma.\n");
    fprintf(stderr, "                               Contributions beyond cutoff*sigma are ignored.\n");
    fprintf(stderr, "Output format options:\n");
    fprintf(stderr, "  --output-format <ascii|hdf5>  Select output format (default: ascii)\n");
    fprintf(stderr, "  -ascii                        Shortcut for --output-format ascii\n");
    fprintf(stderr, "  -hdf5                       Shortcut for --output-format hdf5\n");
    fprintf(stderr, "Photoevaporation options:\n");
    fprintf(stderr, "  -photoevap <0|1>  Enable/disable photoevaporation globally (default: 0 = false)\n");
    fprintf(stderr, "  -photoevap_mode <string>  Model selection: 'Owen' or 'Picogna' (default: 'None')\n");
    fprintf(stderr, "  -xray_lumosity <val>         Stellar X-ray luminosity in erg/s (default: 1.0e30)\n");
    fprintf(stderr,"Initial Condition Profiles:\n");
    fprintf(stderr,"  -cutoff <0|1>               Enable Anna's self-similar profile with exponential taper\n");
    fprintf(stderr,"  -cutoff_radius <double>     Tapering radius in AU for cutoff style (default: 30.0)\n");
    fprintf(stderr,"  -cutoff_sharpness <double>     Sharpness factor for exponential cutoff (default: 1.5)\n\n");
    fprintf(stderr, "Boundary condition options:\n");
    fprintf(stderr, "  -inner_bc <name>   Inner boundary type. Available names:\n");
    fprintf(stderr, "                     zero        = Zero-gradient (Neumann)\n");
    fprintf(stderr, "                     parabolic   = Smooth parabolic extrapolation\n");
    fprintf(stderr, "                     fixed_flux  = Lynden–Bell viscous flux condition\n");
    fprintf(stderr, "                     absorbing   = Dust sink (ghost cell = 0)\n");
    fprintf(stderr, "                     reflecting  = Mirror boundary (non-physical)\n");
    fprintf(stderr, "                     linear      = Linear extrapolation outflow\n");
    fprintf(stderr, "                     loggrid     = Logarithmic-grid extrapolation\n\n");

    fprintf(stderr, "  -outer_bc <name>   Outer boundary type. Same names as inner_bc.\n");
    fprintf(stderr, "                     Recommended: parabolic or linear for smooth outflow.\n\n");

    fprintf(stderr, "Log control:\n");
    fprintf(stderr, "  -v             Enable INFO logging level\n");
    fprintf(stderr, "  -vv            Enable DEBUG logging level (very verbose)\n\n");
    fprintf(stderr, "Test options:\n");
    fprintf(stderr, "  -test <test_name>   Enable a specific test. Available tests:\n");
    fprintf(stderr, "                      ring  = Viscous ring test\n");
    fprintf(stderr, "  -ring_center <val>  Center of the ring for the viscous ring test (default: 10.0 AU)\n");
    fprintf(stderr, "  -ring_width <val>   Width of the ring for the viscous ring test (default: 2.0 AU)\n");
    fprintf(stderr, "  -ring_visc <val>    Fixed viscosity for the viscous ring test (default: 1e-5)\n");
    fprintf(stderr, "  -fixed_ring_visc <val>  Fixed viscosity for the viscous ring test (default: 1e-5)\n");   
    fprintf(stderr, "  -ring_mass <val>  Mass of the ring for the viscous ring test (default: 1.0)\n");

}

int parseCLIOptions(int argc, const char **argv, ParserOptions *opt){

    LOG_DEBUG("Parsing command-line arguments (%d total).\n", argc);
    int i = 1;

    bool flag_has_disk_mass = false;
    bool flag_has_sigma0    = false;

    while (i < argc) {
        if(strcmp(argv[i], "-drift") == 0) {
            i++;
            if (i < argc) opt->option_for_dust_drift = atof(argv[i]); else { LOG_ERROR("Missing value for -drift.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-growth") == 0) {
            i++;
            if (i < argc) opt->option_for_dust_growth = atof(argv[i]); else { LOG_ERROR("Missing value for -growth.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-evol") == 0) {
            i++;
            if (i < argc) opt->option_for_evolution = atof(argv[i]); else { LOG_ERROR("Missing value for -evol.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-twopop") == 0) {
            i++;
            if (i < argc) opt->option_for_dust_secondary_population = atof(argv[i]); else { LOG_ERROR("Missing value for -twopop.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-ufrag") == 0) {
            i++;
            if (i < argc) opt->fragmenatation_velocity = atof(argv[i]); else { LOG_ERROR("Missing value for -ufrag.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-ffrag") == 0) {
            i++;
            if (i < argc) opt->fragmenatation_factor = atof(argv[i]); else { LOG_ERROR("Missing value for -ffrag.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-tStep") == 0) {
            i++;
            if (i < argc) opt->user_defined_time_step = atof(argv[i]); else { LOG_ERROR("Missing value for -tStep.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-n") == 0) {
            i++;
            if (i < argc) opt->number_of_grid_points = atoi(argv[i]); else { LOG_ERROR("Missing value for -n.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-ndust") == 0) { 
            i++;
            if (i < argc) opt->number_of_dust_particles = atoi(argv[i]); else { LOG_ERROR("Missing value for -ndust.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-i") == 0) {
            i++;
            if (i < argc) opt->input_file = argv[i]; else { LOG_ERROR("Missing value for -i.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-o") == 0) { // Output directory name
            i++;
            if (i < argc) {
                strncpy(opt->output_dir_name, argv[i], sizeof(opt->output_dir_name) - 1);
                opt->output_dir_name[sizeof(opt->output_dir_name) - 1] = '\0';
            } else {
                LOG_ERROR("Missing value for -o.\n");
                return 1;
            }
        }
        else if (strcmp(argv[i], "-tmax") == 0) {
            i++;
            if (i < argc) opt->maximum_simulation_time = atof(argv[i]); else { LOG_ERROR("Missing value for -tmax.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-outfreq") == 0) {
            i++;
            if (i < argc) opt->output_frequency = atof(argv[i]); else { LOG_ERROR("Missing value for -outfreq.\n"); return 1; }
        }

        else if (strcmp(argv[i], "-ri") == 0) { 
            i++; 
            if (i < argc) opt->rmin_val = atof(argv[i]); else { LOG_ERROR("Missing value for -ri.\n"); return 1; }; 
        } 
        else if (strcmp(argv[i], "-ro") == 0) { 
            i++; 
            if (i < argc) opt->rmax_val = atof(argv[i]); else { LOG_ERROR("Missing value for -ro.\n"); return 1; }; 
        } 
        // --- PARSING THE NEW MASS OPTIONS ---
        else if (strcmp(argv[i], "-disk_mass") == 0) {
            i++;
            if (i < argc) {
                opt->total_disk_mass = atof(argv[i]);
                flag_has_disk_mass = true;
            } else { LOG_ERROR("Missing value for -disk_mass.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-sigma0_init") == 0) {
            i++; 
            if (i < argc) {
                opt->sigma0_val = atof(argv[i]);
                flag_has_sigma0 = true;
            } else { LOG_ERROR("Missing value for -sigma0_init.\n"); return 1; }; 
        }
        else if (strcmp(argv[i], "-index_init") == 0) { 
            i++; 
            if (i < argc) opt->sigmap_exp_val = atof(argv[i]); else { LOG_ERROR("Missing value for -index_init.\n"); return 1; }; 
        } 
        else if (strcmp(argv[i], "-rdzei") == 0) { 
            i++; 
            if (i < argc) opt->r_dze_i_val = atof(argv[i]); else { LOG_ERROR("Missing value for -rdzei.\n"); return 1; }; 
        }
        else if (strcmp(argv[i], "-rdzeo") == 0) { 
            i++; 
            if (i < argc) opt->r_dze_o_val = atof(argv[i]); else { LOG_ERROR("Missing value for -rdzeo.\n"); return 1; }; 
        }
        else if (strcmp(argv[i], "-drdzei") == 0) { 
            i++; 
            if (i < argc) opt->dr_dze_i_val = atof(argv[i]); else { LOG_ERROR("Missing value for -drdzei.\n"); return 1; }; 
        }
        else if (strcmp(argv[i], "-drdzeo") == 0) { 
            i++; 
            if (i < argc) opt->dr_dze_o_val = atof(argv[i]); else { LOG_ERROR("Missing value for -drdzeo.\n"); return 1; }; 
        }
        else if (strcmp(argv[i], "-alpha_init") == 0) { 
            i++; 
            if (i < argc) opt->alpha_visc_val = atof(argv[i]); else { LOG_ERROR("Missing value for -alpha_init.\n"); return 1; }; 
        } 
        else if (strcmp(argv[i], "-amod") == 0) { 
            i++; 
            if (i < argc) opt->a_mod_val = atof(argv[i]); else { LOG_ERROR("Missing value for -amod.\n"); return 1; }; 
        }
        else if (strcmp(argv[i], "-h_init") == 0) { 
            i++; 
            if (i < argc) opt->hasp_val = atof(argv[i]); else { LOG_ERROR("Missing value for -h_init.\n"); return 1; }; 
        } 
        else if (strcmp(argv[i], "-flind_init") == 0) { 
            i++; 
            if (i < argc) opt->flind_val = atof(argv[i]); else { LOG_ERROR("Missing value for -flind_init.\n"); return 1; }; 
        }
        else if (strcmp(argv[i], "-stellar_mass") == 0) { 
            i++; 
            if (i < argc) opt->star_val = atof(argv[i]); else { LOG_ERROR("Missing value for -stellar_mass.\n"); return 1; }; 
        } 
        else if (strcmp(argv[i], "-eps") == 0) { 
            i++; 
            if (i < argc) opt->eps_val = atof(argv[i]); else { LOG_ERROR("Missing value for -eps.\n"); return 1; }; 
        } 
        else if (strcmp(argv[i], "-ratio") == 0) { 
            i++; 
            if (i < argc) opt->ratio_val = atof(argv[i]); else { LOG_ERROR("Missing value for -ratio.\n"); return 1; }; 
        } 
        else if (strcmp(argv[i], "-mic") == 0) { 
            i++; 
            if (i < argc) opt->mic_val = atof(argv[i]); else { LOG_ERROR("Missing value for -mic.\n"); return 1; }; 
        } 
        else if (strcmp(argv[i], "-onesize") == 0) { 
            i++; 
            if (i < argc) opt->onesize_val = atof(argv[i]); else { LOG_ERROR("Missing value for -onesize.\n"); return 1; }; 
        }
        else if (strcmp(argv[i], "-pdensity") == 0) {
            i++;
            if (i < argc) opt->pdensity_val = atof(argv[i]); else { LOG_ERROR("Missing value for -pdensity.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printUsageToTerminal();
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
                    LOG_ERROR("Unknown output format '%s'. Use ascii or hdf5.\n", argv[i]);
                    return 1;
                }
            } else {
                LOG_ERROR("Missing value for --output-format.\n");
                return 1;
            }
        }

        else if (strcmp(argv[i], "-ring_test") == 0) {
            opt->test_enabled = true;
            strcpy(opt->test_type, "ring");
        }
        else if (strcmp(argv[i], "-ring_center") == 0) {
            i++;
            if (i < argc) opt->ring_center = atof(argv[i]);
            else { LOG_ERROR("Missing value for -ring_center.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-ring_width") == 0) {
            i++;
            if (i < argc) opt->ring_width = atof(argv[i]);
            else { LOG_ERROR("Missing value for -ring_width.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-fixed_ring_visc") == 0) {
            i++;
            if (i < argc) opt->fixed_ring_viscosity = atof(argv[i]);
            else { LOG_ERROR("Missing value for -fixed_ring_visc.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-ring_mass") == 0) {
            i++;
            if (i < argc) opt->ring_mass = atof(argv[i]);
            else { LOG_ERROR("Missing value for -ring_mass.\n"); return 1; }
        }


        else if (strcmp(argv[i], "-photoevap") == 0) {
            i++;
            if (i < argc) {
                int val = atoi(argv[i]);
                opt->enable_photoevaporation = (val != 0);
            } else { LOG_ERROR("Missing value for -photoevap.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-photoevap_mode") == 0) {
            i++;
            if (i < argc) {
                strncpy(opt->photoevaporation_mode, argv[i], sizeof(opt->photoevaporation_mode) - 1);
                opt->photoevaporation_mode[sizeof(opt->photoevaporation_mode) - 1] = '\0';
                
                if (strcasecmp(opt->photoevaporation_mode, "none") != 0) {
                    opt->enable_photoevaporation = true;
                }
            } else { LOG_ERROR("Missing value for -photoevap_mode.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-xray_luminosity") == 0) {
            i++;
            if (i < argc) opt->xray_luminosity = atof(argv[i]); else { LOG_ERROR("Missing value for -xray_luminosity.\n"); return 1; }
        }

        else if (strcmp(argv[i], "-cutoff") == 0) {
            opt->enable_cutoff = atof(argv[++i]); // Beolvassa a mögötte lévő 1.0-át!
            if (opt->enable_cutoff == 1.0) {
                opt->use_cutoff = true;
            }
        }
        else if (strcmp(argv[i], "-cutoff_radius") == 0) {
            if (i + 1 < argc) opt->r_cutoff = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "-cutoff_sharpness") == 0) {
            if (i + 1 < argc) opt->n_for_cutoff = atof(argv[++i]);
        }        

        else if (strcmp(argv[i], "-ascii") == 0) {
            opt->output_format = OUTPUT_ASCII;
        }
        else if (strcmp(argv[i], "-hdf5") == 0) {
            opt->output_format = OUTPUT_HDF5;
        }
        else if (strcmp(argv[i], "-gaussian_sigma_grid_units") == 0) {
            i++;
            if (i < argc) opt->gaussian_sigma = atof(argv[i]);
        }

        else if (strcmp(argv[i], "-gaussian_cutoff_sigma") == 0) {
            i++;
            if (i < argc) opt->gaussian_cutoff = atof(argv[i]);
        }
        else if (strcmp(argv[i], "-density_floor") == 0) {
            i++;
            if (i < argc) opt->density_floor = atof(argv[i]);
            else { LOG_ERROR("Missing value for -density_floor.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-dust_density_floor") == 0) {
            i++;
            if (i < argc) opt->dust_density_floor = atof(argv[i]);
            else { LOG_ERROR("Missing value for -dust_density_floor.\n"); return 1; }
        }
        else if (strcmp(argv[i], "-dust_smoothing") == 0) {
            i++;
            if (i < argc) {
                if (strcmp(argv[i], "cic") == 0)
                    opt->dust_smoothing_mode = 0;
                else if (strcmp(argv[i], "ngp") == 0)
                    opt->dust_smoothing_mode = 1;
                else if (strcmp(argv[i], "tophat") == 0)
                    opt->dust_smoothing_mode = 2;
                else if (strcmp(argv[i], "gaussian") == 0)
                    opt->dust_smoothing_mode = 3;
                else {
                    LOG_ERROR("Unknown dust smoothing mode '%s'. Use cic|ngp|tophat|gaussian.\n", argv[i]);
                    return 1;
                }
            } else {
                LOG_ERROR("Missing value for -dust_smoothing.\n");
                return 1;
            }
        }

        else if (strcmp(argv[i], "-inner_bc") == 0) {
            i++;
            if (i < argc) {
                if      (strcmp(argv[i], "zero_gradient") == 0)         opt->inner_boundary_condition_type = 0;
                else if (strcmp(argv[i], "parabolic") == 0)             opt->inner_boundary_condition_type = 1;
                else if (strcmp(argv[i], "fixed_flux") == 0)            opt->inner_boundary_condition_type = 2;
                else if (strcmp(argv[i], "absorbing") == 0)             opt->inner_boundary_condition_type = 3;
                else if (strcmp(argv[i], "reflecting") == 0)            opt->inner_boundary_condition_type = 4;
                else if (strcmp(argv[i], "linear_extrapolation") == 0)  opt->inner_boundary_condition_type = 5;
                else if (strcmp(argv[i], "loggrid") == 0)               opt->inner_boundary_condition_type = 6;
                else {
                    LOG_ERROR("Unknown inner BC '%s'.\n", argv[i]);
                    return 1;
                }
            }
        }
        else if (strcmp(argv[i], "-outer_bc") == 0) {
            i++;
            if (i < argc) {
                if      (strcmp(argv[i], "zero_gradient") == 0)         opt->outer_boundary_condition_type = 0;
                else if (strcmp(argv[i], "parabolic") == 0)             opt->outer_boundary_condition_type = 1;
                else if (strcmp(argv[i], "fixed_flux") == 0)            opt->outer_boundary_condition_type = 2;
                else if (strcmp(argv[i], "absorbing") == 0)             opt->outer_boundary_condition_type = 3;
                else if (strcmp(argv[i], "reflecting") == 0)            opt->outer_boundary_condition_type = 4;
                else if (strcmp(argv[i], "linear_extrapolation") == 0)  opt->outer_boundary_condition_type = 5;
                else if (strcmp(argv[i], "loggrid") == 0)               opt->outer_boundary_condition_type = 6;
                else {
                    LOG_ERROR("Unknown outer BC '%s'.\n", argv[i]);
                    return 1;
                }
            }
        }

        else if (strcmp(argv[i], "-v") == 0) {
            if (opt->verbose == 2) {
                LOG_WARN("You provided both -v and -vv. Ignoring -v and using -vv (DEBUG level).\n");
            } else {
                opt->verbose = 1; // Info szint
            }
        }
        else if (strcmp(argv[i], "-vv") == 0) {
            if (opt->verbose == 1) {
                LOG_WARN("You provided both -v and -vv. Using -vv (DEBUG level) instead of -v.\n");
            }
            opt->verbose = 2; // Debug szint
        }


        else {
            LOG_ERROR("ERROR [parseCLIOptions]: Invalid switch on command-line: %s!\n", argv[i]);
  //          printUsageToTerminal(); 
            return 1;
        }
        i++;
    }


        // --- CRITICAL CONFLICT RESOLUTION GUARD ---
        if (flag_has_disk_mass && flag_has_sigma0) {
            printFatalErrorMessageForDiskMass();
            return 1;
        }

    // Fallback fall-throughs if neither was explicitly set by the user
    if (!flag_has_disk_mass && !flag_has_sigma0) {
        opt->sigma0_val = 1.0; // Apply the safe default fallback to Sigma0
    }


    LOG_DEBUG("Command-line parsing complete.\n");
    return 0;
}