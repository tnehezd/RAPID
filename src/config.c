// src/config.c

/**
 * @file config.c
 * @brief This file contains the definitions (memory allocation and initialization)
 * for global configuration variables, file pointers, and constants used throughout the simulation.
 * @details The corresponding declarations for these global elements are found in `globals.h`.
 */

#include "config.h"
#include "simulation_types.h"
#include "globals.h"

int PARTICLE_NUMBER;

// --- Global File Pointer Definitions ---
FILE *timescale_output_file = NULL;
FILE *info_current_file = NULL;
FILE *file_in = NULL;

const char *inputsig = NULL;

// --- Global Filename Definitions (Constant Pointers) ---
// These cannot be modified at runtime.
const char * const FILENAME_INIT_GAS_PROFILE = "initial_gas_profile.dat";   
const char * const FILENAME_INIT_DUST_PROFILE = "initial_dust_profile.dat"; 
const char * const FILE_DENS_PREFIX = "density_profile";
const char * const FILE_DUST_PREFIX = "dust_profile";
const char * const FILE_MASS_ACCUMULATE = "mass_accumulation_dze_edge.dat";
const char * const FILE_DUST_EVOLUTION = "dust_particle_evolution.dat";
const char * const FILENAME_DISK_PARAM = "disk_config.dat";
const char * const FILE_TIMESCALE = "drift_timescale.dat";
const char * const FILE_SUMMARY = "summary.dat";
const char * const LOGS_DIR = "LOGS";
const char * const CONFIG_DIR = "config";