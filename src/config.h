#ifndef CONFIG_H
#define CONFIG_H

#include <stdio.h> 
#include <math.h>  
#include "globals.h"
#include "simulation_types.h" 

/**
 * @file config.h
 * @brief Global configuration parameters, file pointers, and constants for the simulation.
 * @details This header defines a central repository for various parameters for the simulation,
 * including the total number of particles and file handles for I/O operations. All parameters
 * are declared with `extern`, signifying that their definitions (memory allocation)
 * are in the `config.c` file.
 */

// --- Global Variable Declarations (extern) ---

// Dust parameter
/**
 * @brief Global variable storing the total number of particles in the simulation.
 * @details This variable defines the number of individual dust particles. This
 * value is set during the first step and is fixed for the whole simulation. 
 * Declared as `extern` here, defined in `config.c`.
 */
extern int PARTICLE_NUMBER;


// --- Global File Pointer Declarations (extern) ---


/**
 * @brief File pointer for output related to simulation timescales.
 * @details This file handle is specifically intended for logging and storing
 * data related to various characteristic timescales computed during the simulation,
 * such as dynamical, thermal, or accretion timescales. It complements other
 * general output files by providing a dedicated stream for temporal analysis.
 * Declared as `extern` here, defined and initialized (to `NULL`) in `config.c`.
 */
extern FILE *timescale_output_file;

/**
 * @brief File pointer to store a summary of the current simulation run's configuration and status.
 * @details This file is typically opened at the beginning of a simulation run and
 * logs critical input parameters, timestamp, and other relevant metadata for
 * reproducibility and quick overview. Declared as `extern` here,
 * defined and initialized (to `NULL`) in `config.c`.
 */
extern FILE *info_current_file;

/**
 * @brief First input file pointer for general input data.
 * @details Used to read primary input data streams for the simulation.
 * Declared as `extern` here, defined and initialized (to `NULL`) in `config.c`.
 */
extern FILE *file_in;


// --- Global Filename Declarations (extern const char * const) ---
// Define distinct names for gas and dust initial profiles
/**
 * @brief Critical: Defines the filename for the initial gas surface density profile.
 * @details This constant string specifies the name of the file expected to contain
 * the initial radial distribution of gas within the simulated environment.
 * Declared as `extern const char * const` here, defined in `config.c`.
 */
extern const char * const FILENAME_INIT_GAS_PROFILE;

/**
 * @brief Critical: Defines the filename for the initial dust surface density profile.
 * @details This constant string specifies the name of the file expected to contain
 * the initial radial distribution of dust within the simulated environment.
 * Declared as `extern const char * const` here, defined in `config.c`.
 */
extern const char * const FILENAME_INIT_DUST_PROFILE;

/**
 * @brief Prefix for output filenames containing time-dependent density profiles.
 * @details This constant string defines the base name for output files that store
 * the evolving density profiles over time. Actual filenames will typically be formed
 * by appending a timestamp or iteration number (e.g., `density_profile_TIME.dat`).
 * Declared as `extern const char * const` here, defined in `config.c`.
 */
extern const char * const FILE_DENS_PREFIX;

/**
 * @brief Filename for the output file tracking mass accumulation at the DZE edge.
 * @details This constant string specifies the name of the file used to log gas and/or
 * dust mass flux and accumulation rates at the inner and outer edges of the
 * Dead Zone Edge (DZE) or other relevant boundaries.
 * Declared as `extern const char * const` here, defined in `config.c`.
 */
extern const char * const FILE_MASS_ACCUMULATE;

/**
 * @brief Filename for the output file detailing the evolution of individual dust particles.
 * @details This constant string specifies the name of the file intended to contain
 * detailed, time-series data for each simulated particle, such as positions,
 * velocities, sizes, or other relevant properties.
 * Declared as `extern const char * const` here, defined in `config.c`.
 */
extern const char * const FILE_DUST_EVOLUTION;

/**
 * @brief Filename for the configuration file containing global disk parameters.
 * @details This constant string specifies the name of the primary configuration file
 * for initializing the simulation. It typically defines properties like disk mass,
 * radius, viscosity, and other global physical parameters.
 * Declared as `extern const char * const` here, defined in `config.c`.
 */
extern const char * const FILENAME_DISK_PARAM;
/**
 * @brief Constant string defining the name of the directory where all simulation log files are stored.
 * @details This constant helps in organizing simulation output by centralizing
 * all log, diagnostic, and summary files into a dedicated directory.
 * Declared as `extern const char * const` here, defined in `config.c`.
 */

extern const char * const LOGS_DIR;

/**
 * @brief Filename for the output file detailing simulation timescales.
 * @details This constant string specifies the name of the file used to log
 * various characteristic timescales of the simulation, such as dynamical
 * timescales, accretion timescales, or other relevant temporal parameters.
 * It helps in analyzing the evolution of different physical processes over time.
 * Declared as `extern const char * const` here, its definition (the actual string literal)
 * should reside in `config.c`.
 */
extern const char * const FILE_TIMESCALE;

/**
 * @brief Constant string defining the name of the directory where configuration files are located.
 * @details This provides a standard, relative path for locating input configuration
 * files and other setup-related resources.
 * Declared as `extern const char * const` here, defined in `config.c`.
 */
extern const char * const CONFIG_DIR;

/**
 * @brief Global pointer to a character string specifying the input sigma profile filename.
 * @details This string holds the name of the file from which the initial sigma profile
 * (e.g., surface density distribution) for the disk is read. Its value is typically
 * determined during configuration parsing at runtime.
 * Declared as `extern` here, defined (initialized to `NULL`) in `config.c`.
 */
extern const char *inputsig; // Parameter controls


#endif // CONFIG_H