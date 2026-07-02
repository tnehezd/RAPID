/**
 * @file config.h
 * @brief Global configuration, physical constants, and file handles.
 * * This file contains global variable declarations, preprocessor macros for 
 * physical constants, unit conversion factors, and external file pointers 
 * used for simulation logging and data management.
 */

#ifndef CONFIG_H
#define CONFIG_H

#include <stdio.h> 
#include <math.h>  

/** @brief Total number of representative dust particles in the simulation. */
extern int particle_number;

/** @name Simulation Output File Pointers
 * External file handles used for real-time tracking and data export.
 */
/**@{*/
extern FILE *drift_timescale_file;      /**< Pointer to the file tracking dust radial drift timescales. */
extern FILE *current_info_file;         /**< Pointer to the file for general simulation status and logging. */	
extern FILE *load_dust_particles_file;  /**< Pointer to the input file used for loading initial particle states. */
/**@}*/

#include "simulation_types.h"

/** @name Physical Constants and Unit Conversions
 * Macros defining fundamental constants and factors to maintain the simulation's unit system.
 */
/**@{*/
/** @brief Factor to convert surface density to internal simulation units. */
#define SURFACE_DENSITY_CONVERSION_FACTOR       1.12521e-7   

/** @brief Gravitational constant in dimensionless form (standardized to 1.0). */
#define G_DIMENSIONLESS                         1.0          

/** @brief Conversion factor from \f$cm/s\f$ to \f$AU/year\f$, normalized by \f$2\pi\f$. */
#define CM_PER_SEC_TO_AU_PER_YEAR_2PI           3.35725e-07  

/** @brief Mass of the Sun in grams (CGS units). */
#define SOLAR_MASS_IN_GRAMS                     1.989e33     

/** @brief One Astronomical Unit in centimeters (CGS units). */
#define AU_IN_CM                                1.496e13     

/** @brief Factor to convert days to years. */
#define DAYS_PER_YEAR_CONVERSION_FACTOR         2.737850787e-3

/** @brief Factor to convert years to days. */
#define YEARS_PER_DAY_CONVERSION_FACTOR         365.242199

/** @brief Safety maximum value for timestep calculations. */
#define TIMESTEP_MAX_SAFETY_LIMIT               -10000.0

/** @brief Numerical safety factor used for rounding or thresholding operations. */
#define ROUNDING_FACTOR                         1.0          

/** @brief Current version string of the RAPID simulation core. */
#define SIM_VERSION                             "2.2.0"
/**@}*/

/** @name File and Directory Name Constants
 * Global string constants defining the naming convention for input/output files.
 */
/**@{*/
extern const char * const kInitialGasProfileFileName;   			/**< Filename for the initial gas profile. */
extern const char * const kInitialDustProfileFileName;  			/**< Filename for the initial dust profile. */
extern const char * const kGasDensityProfileFilePrefix; 			/**< Prefix for time-dependent gas density snapshots. */
extern const char * const kDustAccumulationFileName;    			/**< Filename for tracking accumulated dust mass. */
extern const char * const kDustParticleEvolutionFile;   			/**< Filename for Lagrangian particle trajectory data. */
extern const char * const kDriftTimescaleFileName;      			/**< Filename for storing drift timescale analysis. */
extern const char * const kDustDensityProfileFilePrefix; 			/**< Prefix for Eulerian dust density snapshots. */
extern const char * const kMicronDustDensityProfileFilePrefix; 		/**< Prefix for small (micron) dust snapshots. */
extern const char * const kDustMicronParticleEvolutionFile; 		/**< Evolution file for small dust populations. */
extern const char * const kDiskConfigFile;              			/**< Filename for the primary disk configuration input. */
extern const char * const kLogFilesDirectory;           			/**< Directory path where log files are stored. */
extern const char * const kConfigFilesDirectory;        			/**< Directory path where configuration files reside. */
extern const char * const kFileNamesSuffix;             			/**< Common suffix for output files for ASCII files (e.g., .dat or .txt). */
extern const char * const kCurrentInfoFile;             			/**< Filename for general runtime info. */
extern const char * const kDustParticleSizeFileName;    			/**< Filename for particle size distribution data. */
extern const char * const kCurrentRuntimeInfoFile;      			/**< Filename for high-frequency runtime diagnostics. */
extern const char * const kSnapshotOutputFileNamePrefix;			/**< Filename for HDF5 snapshot output. */
extern const char * const kTimeSeriesForMassAccumulatinFileName; 	/**< Filename for HDF mass accumulation file. */
extern const char * const kFileNamesHDF5Suffix;						/**< Common suffix for output files for HDF5 files (e.g., .dat or .txt). */

extern const int kTerminalWidth;                                    /**< Width of the terminal for formatted output. */
/**@}*/

#endif // CONFIG_H