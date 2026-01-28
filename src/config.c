// src/config.c
#include "config.h"
#include "simulation_types.h"

int PARTICLE_NUMBER;

// --- Global File Pointer Definitions ---
FILE *drift_timescale_file = NULL;
FILE *current_info_file = NULL;
FILE *load_dust_particles_file = NULL;
const char *inputsig = NULL;

// --- Global Filename Definitions (Constant Pointers) ---
// Using const char * const for truly constant string literals.
// These cannot be modified at runtime.

const char * const kInitialGasProfileFileName = "initial_gas_profile";   
const char * const kInitialDustProfileFileName = "initial_dust_profile"; 
const char * const kGasDensityProfileFilePrefix = "density_profile";
const char * const kDriftTimescaleFileName = "drift_timecale";
const char * const kDustAccumulationFileName = "mass_accumulation_dze_edge";
const char * const kDustParticleEvolutionFile = "dust_particle_evolution";
const char * const kDustMicronParticleEvolutionFile = "micron_dust_particle_evolution";
const char * const kDiskConfigFile = "disk_config";
const char * const kFileNamesSuffix = ".dat";
const char * const kLogFilesDirectory = "LOGS";
const char * const kConfigFilesDirectory = "config";
const char * const kCurrentInfoFile = "current_run_info";