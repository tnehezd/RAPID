/**
 * @file parser.h
 * @brief Command-line option parsing for the simulation.
 *
 * Provides functions for initializing default parameters, parsing user
 * input from the command line, and printing usage information.
 * The parsed values populate a ParserOptions structure used to configure
 * the simulation.
 */

#ifndef PARSER_H
#define PARSER_H

#include <stdbool.h>  
#include "config.h"   
#include "simulation_types.h" 

#define MAX_OUTPUT_DIR_LEN MAX_PATH_LEN

/**
 * @brief Structure holding all command-line parser options for the simulation.
 *
 * This structure stores all user-configurable parameters parsed from the
 * command line or configuration files. These values are later used to
 * initialize SimulationOptions and DiskParameters.
 */
typedef struct {
    double option_for_dust_drift;                       /**< Enable/disable dust drift module. */
    double option_for_dust_growth;                      /**< Enable/disable dust growth module. */
    double option_for_evolution;                        /**< Enable/disable disk evolution. */
    double option_for_dust_secondary_population;        /**< Enable/disable secondary dust population. */
    double fragmenatation_velocity;                     /**< Fragmentation velocity threshold. */
    double fragmenatation_factor;                       /**< Fragmentation mass reduction factor. */
    int    number_of_grid_points;                       /**< Number of radial grid points. */
    int    number_of_dust_particles;                    /**< Number of dust particles to simulate. */
    double rmin_val;                                    /**< Inner disk radius. */
    double rmax_val;                                    /**< Outer disk radius. */
    double sigma0_val;                                  /**< Gas surface density normalization. */
    double total_disk_mass;                             /**< Total mass of the disk in MSun unit */
    double sigmap_exp_val;                              /**< Power-law exponent of surface density. */
    double alpha_visc_val;                              /**< Disk viscosity parameter. */
    double star_val;                                    /**< Stellar mass. */
    double hasp_val;                                    /**< Disk aspect ratio parameter. */
    double flind_val;                                   /**< Disk flaring index. */
    double r_dze_i_val;                                 /**< Inner dead-zone edge radius. */
    double r_dze_o_val;                                 /**< Outer dead-zone edge radius. */
    double dr_dze_i_val;                                /**< Inner dead-zone width. */
    double dr_dze_o_val;                                /**< Outer dead-zone width. */
    double a_mod_val;                                   /**< Alpha modification factor. */
    const char *input_file;                             /**< Path to input configuration file. */
    char output_dir_name[MAX_OUTPUT_DIR_LEN];           /**< Output directory name. */
    double user_defined_time_step;                      /**< User-specified simulation time step. */
    double maximum_simulation_time;                     /**< Maximum allowed simulation time. */
    double output_frequency;                            /**< Frequency of output snapshots. */
    double eps_val;                                     /**< Numerical epsilon parameter. */
    double ratio_val;                                   /**< Ratio parameter for dust/gas. */
    double mic_val;                                     /**< Micron dust fraction. */
    double onesize_val;                                 /**< Single-size dust mode parameter. */
    double pdensity_val;                                /**< Dust particle density. */
    int output_format;                                  /**< Output format selector (0 = ASCII, 1 = HDF5). */
    int dust_smoothing_mode;                            /**< Smoothing method for mapping Lagrangian dust surface density to Eulerian grid. */
    double gaussian_sigma;                              /**< Width of the Gaussian kernel in grid-cell units (Δr). Controls how far dust mass spreads around a particle. */
    double gaussian_cutoff;                             /**< Kernel cutoff in multiples of sigma. Contributions beyond cutoff*sigma are ignored. */
    bool enable_photoevaporation;                       /**< Option for enabling photoevaporation */
    char photoevaporation_mode[32];                     /**< Switch between Picogna, or Owen model */
    double xray_luminosity;                             /**< X-ray luminosity of the central star */
    double density_floor;                               /**< Minimum allowed gas surface density */
    double enable_cutoff;                               /**< Enable cutoff radius */
    bool   use_cutoff;                                  /**< Profile style: false = Normal Power-law, true = Cutoff */
    double r_cutoff;                                    /**< Scaling/tapering radius for cutoff profile [AU] */
    double n_for_cutoff;                                /**< Shape/sharpness parameter for exponential cutoff */
    int inner_boundary_condition_type;                  /**< Inner boundary condition */
    int outer_boundary_condition_type;                  /**< Outer boundary condition */

} ParserOptions;

/**
 * @brief Initializes a ParserOptions structure with default values.
 *
 * @param opt Pointer to the ParserOptions structure to initialize.
 */
void createDefaultOptions(ParserOptions *opt);

/**
 * @brief Parses command-line arguments and fills a ParserOptions structure.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of argument strings.
 * @param opt Pointer to the ParserOptions structure to populate.
 * @return 0 on success, non-zero on parsing error.
 */
int parseCLIOptions(int argc, const char **argv, ParserOptions *opt);

/**
 * @brief Prints usage information and available command-line options to the terminal.
 */
void printUsageToTerminal();

#endif // PARSER_H