#ifndef CONFIG_H
#define CONFIG_H

#include <stdio.h> 
#include <math.h>  
#include "globals.h"
#include "simulation_types.h" 

/**
 * @file config.h
 * @brief Global configuration parameters and file pointers.
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
 * @brief File pointer for drift timescales output.
 * @details This file (`drift_timescale.dat`) contains the precalculated values of the dust depletion timescale
 * at a given radius within your predefined `output` dir. The timescale values are presented in years, consistent with
 * the simulation's astrophysical unit system. In this system, the gravitational
 * constant $G=1$, and one orbital period at 1 Astronomical Unit (AU) is defined as
 * $2\pi$ simulation time units, which corresponds to one year.
 *
 */
extern FILE *timescale_output_file;

/**
 * @brief File pointer to store timestamp and a summary log of the current simulation run's configuration and status.
 * @details This file is opened at the beginning of a simulation run and
 * logs critical input parameters of the simulation (¨summary.dat`) within your predefined `output` dir. 
 * Declared as `extern` here, defined and initialized (to `NULL`) in `config.c`.
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
 * @details This constant string specifies the name of the file (`initial_gas_profile.dat`) that contains
 * the initial radial distribution of gas in the `config` directory of the predefined `output` dir.
 * Declared as `extern const char * const` here, defined in `config.c`.
 */
extern const char * const FILENAME_INIT_GAS_PROFILE;

/**
 * @brief Critical: Defines the filename for the initial dust surface density profile.
 * @details This constant string specifies the name of the file (`initial_dust_profile.dat`) that contains
 * the initial radial distribution of dus in the `config` directory of the predefined `output` dir.
 * Declared as `extern const char * const` here, defined in `config.c`.
 */
extern const char * const FILENAME_INIT_DUST_PROFILE;

/**
 * @brief Prefix for output filenames containing time-dependent gas surface density profiles.
 * @details This constant string defines the base name for output files that store
 * the evolving density profiles over time. Actual filenames will be formed
 * by appending a zero-padded number representing the simulation time
 * (e.g., `density_profile_00000010.dat`) within `LOGS` in the predefined `output` directory.
 * Declared as `extern const char * const` here, defined in `config.c`.
 */
extern const char * const FILE_DENS_PREFIX;

/**
 * @brief Prefix for output filenames containing time-dependent dust density profiles.
 * @details This constant string defines the base name for output files that store
 * the evolving density profiles over time. Actual filenames will be formed
 * by appending a zero-padded number representing the simulation time
 * (e.g., `dust_profile_00000010.dat`) within `LOGS` in the predefined `output` directory.
 * Declared as `extern const char * const` here, defined in `config.c`.
 */
extern const char * const FILE_DUST_PREFIX;

/**
 * @brief Filename for the output file tracking mass accumulation at the DZE edge.
 * @details This constant string specifies the name of the file used to log 
 * dust accumulation at the inner and outer edges of the
 * Dead Zone Edge (DZE) within `LOGS` in the predefined `output` directory.
 * Declared as `extern const char * const` here, defined in `config.c`.
 */
extern const char * const FILE_MASS_ACCUMULATE;

/**
 * @brief Filename for the output file detailing the evolution of individual dust particles.
 * @details This constant string specifies the name of the file intended to contain
 * detailed, time-series data for each simulated particle, such as positions,
 * velocities, sizes, or other relevant properties within `LOGS` in the predefined `output` directory.
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
 * @brief Filename for the output file detailing simulation timescales.
 * @details This constant string specifies the name of the file used to log
 * various characteristic timescales of the simulation, such as dynamical
 * timescales, accretion timescales, or other relevant temporal parameters.
 * It helps in analyzing the evolution of different physical processes over time.
 * Declared as `extern const char * const` here, its definition (the actual string literal)
 * should reside in `config.c`.
 */
extern const char * const FILE_SUMMARY;

/**
 * @brief Constant string defining the name of the directory where all simulation log files are stored.
 * @details This constant helps in organizing simulation output by centralizing
 * all log, diagnostic, and summary files into a dedicated directory.
 * Declared as `extern const char * const` here, defined in `config.c`.
 */

extern const char * const LOGS_DIR;


/**
 * @brief Constant string defining the name of the directory where configuration files are located.
 * @details This provides a standard, relative path for locating input configuration
 * files and other setup-related resources.
 * Declared as `extern const char * const` here, defined in `config.c`.
 */
extern const char * const CONFIG_DIR;

/**
 * @brief Global pointer to a character string specifying the input surface density profile filename.
 * @details This string holds the name of the file from which the initial surfae density profile for the disk is read. 
 * Its value is typically determined during configuration parsing at runtime.
 * Declared as `extern` here, defined (initialized to `NULL`) in `config.c`.
 * @note CHECK THE FUNCTIONALITY!!!
 */
extern const char *inputsig; // Parameter controls


#endif // CONFIG_H