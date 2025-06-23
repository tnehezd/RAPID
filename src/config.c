// src/config.c
#include "config.h"

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
const char * const FILENAME_INIT_PROFILE = "initial_dust_profile.dat";
const char * const FILENAME_DISK_PARAM = "disk_config.dat";
const char * const INITIAL_SURFACE_DENSITY_FILE = "initial_gas_surface_density.dat";



