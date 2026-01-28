#ifndef CONFIG_H
#define CONFIG_H

#include <stdio.h> // Required for FILE* type
#include <math.h>  // Required for M_PI (for TWOPI macro)

// --- Global Variable Declarations (extern) ---
// Dust parameters
extern int particle_number;


// --- Global File Pointer Declarations (extern) ---
extern FILE *drift_timescale_file;
extern FILE *current_info_file;
extern FILE *load_dust_particles_file;

#include "simulation_types.h"

// --- Physical Constants (Macros) ---
// These constants have been moved here from init_tool.c for global access,
// to prevent duplication or static limitations.

#define SURFACE_DENSITY_CONVERSION_FACTOR            1.12521e-7          // Surface density conversion factor
#define ICE_LINE_DUST_ENHANCEMENT_FACTOR         3.0                 // Factor for dust density beyond snowline
#define SNOWLINE_RADIUS_AU          2.7                 // Snowline radius in AU
#define G_DIMENSIONLESS 1.0 // Gravitációs Konstans (dimenziótlan G=1 rendszerben)
// Ha az init_tool_module.c is ezt a G-t akarja használni, akkor itt lehet egy alias
#define CM_PER_SEC_TO_AU_PER_YEAR_2PI 3.35725e-07         // cm/sec to AU/(yr/2pi) conversion

#define SOLAR_MASS_IN_GRAMS 1.989e33            // Solar Mass in grams (M_solar -> g) - PLEASE VERIFY THIS VALUE!
#define AU_IN_CM 1.496e13            // Astronomical Unit in centimeters (AU -> cm) - PLEASE VERIFY THIS VALUE!
#define ROUNDING_FACTOR 1.0                 // Make sure this value is correct for your physics model!


// --- Global Filename Declarations (extern) ---
// Define distinct names for gas and dust initial profiles
extern const char * const kInitialGasProfileFileName;   // NEW: For initial gas profile
extern const char * const kInitialDustProfileFileName;  // NEW: For initial dust profile
// You can remove or repurpose FILENAME_INIT_PROFILE if it's no longer generic.
// For clarity, I recommend using the more specific names above.

extern const char * const kGasDensityProfileFilePrefix;
extern const char * const kDustAccumulationFileName;
extern const char * const kDustParticleEvolutionFile;
extern const char * const kDriftTimescaleFileName;
extern const char * const kDustMicronParticleEvolutionFile;
extern const char * const kDiskConfigFile;
extern const char * const kLogFilesDirectory;
extern const char * const kConfigFilesDirectory;
extern const char * const kFileNamesSuffix;
extern const char * const kCurrentInfoFile;


#endif // CONFIG_H