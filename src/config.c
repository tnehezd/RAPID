// src/config.c
#include "config.h"
#include "simulation_types.h"

int PARTICLE_NUMBER;

// --- Global File Pointer Definitions ---
FILE *fmo = NULL;
FILE *fout = NULL;
FILE *fout2 = NULL;
FILE *fout3 = NULL;
FILE *foutmicr = NULL;
FILE *massfil = NULL;
FILE *jelfut = NULL;
FILE *fin1 = NULL;
FILE *fin2 = NULL;
FILE *fil = NULL;
const char *inputsig = NULL;

// --- Global Filename Definitions (Constant Pointers) ---
// Using const char * const for truly constant string literals.
// These cannot be modified at runtime.

const char * const kInitialGasProfileFileName = "initial_gas_profile.dat";   
const char * const kInitialDustProfileFileName = "initial_dust_profile.dat"; 
const char * const kGasDensityProfileFilePrefix = "density_profile";
const char * const kDustAccumulationFileName = "mass_accumulation_dze_edge.dat";
const char * const kDustParticleEvolutionFile = "dust_particle_evolution.dat";
const char * const kDiskConfigFile = "disk_config.dat";
const char * const kLogFilesDirectory = "LOGS";
const char * const kConfigFilesDirectory = "config";