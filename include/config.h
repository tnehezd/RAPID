#ifndef CONFIG_H
#define CONFIG_H

#include <stdio.h> 
#include <math.h>  

extern int particle_number;
extern FILE *drift_timescale_file;
extern FILE *current_info_file;
extern FILE *load_dust_particles_file;

#include "simulation_types.h"

#define SURFACE_DENSITY_CONVERSION_FACTOR       1.12521e-7   
#define G_DIMENSIONLESS 						1.0 		
#define CM_PER_SEC_TO_AU_PER_YEAR_2PI			3.35725e-07 
#define SOLAR_MASS_IN_GRAMS 					1.989e33    
#define AU_IN_CM 								1.496e13    
#define ROUNDING_FACTOR 						1.0         
#define SIM_VERSION 							"2.1-beta"

extern const char * const kInitialGasProfileFileName;   
extern const char * const kInitialDustProfileFileName;  
extern const char * const kGasDensityProfileFilePrefix;
extern const char * const kDustAccumulationFileName;
extern const char * const kDustParticleEvolutionFile;
extern const char * const kDriftTimescaleFileName;
extern const char * const kDustDensityProfileFilePrefix;
extern const char * const kMicronDustDensityProfileFilePrefix;
extern const char * const kDustMicronParticleEvolutionFile;
extern const char * const kDiskConfigFile;
extern const char * const kLogFilesDirectory;
extern const char * const kConfigFilesDirectory;
extern const char * const kFileNamesSuffix;
extern const char * const kCurrentInfoFile;
extern const char * const kDustParticleSizeFileName;
extern const char * const kCurrentRuntimeInfoFile;

#endif // CONFIG_H