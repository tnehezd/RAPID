#ifndef INIT_TOOL_MODULE_H
#define INIT_TOOL_MODULE_H

// Required includes if not already present globally
#include <stdbool.h> // For 'true'/'false'

#include "simulation_types.h"


// Structure to hold initialization options
typedef struct init_tool_options {
    // Grid and Physical Parameters
    int     n_grid_points;
    int     n_dust_particles;  // NEW: Number of initial dust particles written to file header

    double  r_inner;
    double  r_outer;
    double  sigma0_gas_au;      // Gas surface density at 1 AU [M_Sun/AU^2]
    double  sigma_exponent;     // Positive exponent for surface density profile (Sigma ~ r^(-index))
    double  alpha_viscosity;    // Alpha viscosity parameter
    double  star_mass;          // Central star mass [M_Sun]
    double  aspect_ratio;       // Disk aspect ratio (H/r)
    double  flaring_index;      // Flaring index for disk height (H ~ r^(1+flind))

    // Dead Zone Parameters
    double  deadzone_r_inner;
    double  deadzone_r_outer;
    double  deadzone_dr_inner;  // Transition width multiplier (e.g., in units of H)
    double  deadzone_dr_outer;  // Transition width multiplier
    double  deadzone_alpha_mod; // Alpha reduction factor in dead zone (e.g., 0.01 for 1% of original alpha)

    // Dust Parameters
    double  dust_to_gas_ratio;      // Initial dust-to-gas ratio (epsilon)
    double  disk_mass_dust;         // Total dust disk mass [M_Sun] - if > 0, Sigma0 is calculated from this
    double  one_size_particle_cm;   // If > 0, particles are fixed to this size, overrides two_pop_ratio
    double  two_pop_ratio;          // Ratio of mass in larger particles for two-population model (0.0 - 1.0)
    double  micro_size_cm;          // Size of micron-sized particles for two-population model
    double  f_drift;                // Factor for drift-limited size (default value, adjust as needed)
    double  f_frag;                 // Factor for fragmentation-limited size (default value, adjust as needed)

    char output_base_path[MAX_PATH_LEN]; // NEW: Base path where init_tool should create its output files
    double dust_density_g_cm3; // NEW: Por szemcse sűrűsége (g/cm^3)


} init_tool_options_t;

// Function prototypes
void initializeDefaultOptions(init_tool_options_t *opt);
int runInitialization(init_tool_options_t *opts, DiskParameters *output_disk_params);

#endif // INIT_TOOL_MODULE_H