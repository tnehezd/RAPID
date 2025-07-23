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
const char * const FILENAME_INIT_GAS_PROFILE = "initial_gas_profile.dat";   // This is the one you need for sigIn!
const char * const FILENAME_INIT_DUST_PROFILE = "initial_dust_profile.dat"; // This is your existing dust file

// You can now remove 'FILENAME_INIT_PROFILE' as it's redundant/conflicting
// with the more specific names above, unless it serves another generic purpose.
// If you keep it, make sure its usage is unambiguous.

const char * const FILE_DENS_PREFIX = "density_profile";
const char * const FILE_MASS_ACCUMULATE = "mass_accumulation_dze_edge.dat";
const char * const FILE_DUST_EVOLUTION = "dust_particle_evolution.dat";



const char * const FILENAME_DISK_PARAM = "disk_config.dat";



const char * const LOGS_DIR = "LOGS";
const char * const CONFIG_DIR = "config";