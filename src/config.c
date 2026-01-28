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

// CRITICAL CHANGE: Define separate names for gas and dust profiles
const char * const kInitialGasProfileFileName = "initial_gas_profile.dat";   // This is the one you need for loadGasSurfaceDensityFromFile!
const char * const kInitialDustProfileFileName = "initial_dust_profile.dat"; // This is your existing dust file

// You can now remove 'FILENAME_INIT_PROFILE' as it's redundant/conflicting
// with the more specific names above, unless it serves another generic purpose.
// If you keep it, make sure its usage is unambiguous.

const char * const kGasDensityProfileFilePrefix = "density_profile";
const char * const kDustAccumulationFileName = "mass_accumulation_dze_edge.dat";
const char * const kDustParticleEvolutionFile = "dust_particle_evolution.dat";



const char * const kDiskConfigFile = "disk_config.dat";



const char * const kLogFilesDirectory = "LOGS";
const char * const kConfigFilesDirectory = "config";