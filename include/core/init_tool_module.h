/**
 * @file init_tool_module.h
 * @brief Initialization tools and default option management.
 *
 * This module defines the setup structures and functions used to bridge
 * configuration files with the internal disk simulation parameters.
 */

#ifndef INIT_TOOL_MODULE_H
#define INIT_TOOL_MODULE_H

#include <stdbool.h>
#include "simulation_types.h"

/**
 * @struct InitializeDefaultOptions
 * @brief Configuration structure for setting up the initial disk model.
 * * This structure holds all the raw input parameters required to define
 * the disk's physical and numerical setup before the simulation starts.
 */
typedef struct {
    bool    use_cutoff;                     /**< Initial condition profile type: false = Normal Power-Law, true = Exponential Cutoff. */
    int     n_grid_points;                  /**< Number of radial grid cells. */
    int     n_dust_particles;               /**< Total number of Lagrangian dust particles. */
    double  r_inner;                        /**< Inner disk boundary in AU. */
    double  r_outer;                        /**< Outer disk boundary in AU. */
    double  sigma0_gas_au;                  /**< Gas surface density normalization at 1 AU. */
    double  sigma_exponent;                 /**< Power-law exponent for the surface density profile. */
    double  r_cutoff;                       /**< Characteristic exponential scaling/tapering radius in AU (used if ic_type == 2). */
    double  n_for_cutoff;                   /**< Shape/sharpness parameter for the exponential cutoff (typically 2.0 - gamma). */
    double  alpha_viscosity;                /**< Background Shakura-Sunyaev alpha viscosity. */
    double  star_mass;                      /**< Central star mass in solar units. */
    double  aspect_ratio;                   /**< Disk aspect ratio, \f$h=H/r\f$, at 1 AU. */
    double  flaring_index;                  /**< Power-law index for disk flaring. */
    double  density_floor;                  /**< Minimum allowed gas surface density (Σ_min) to prevent numerical underflow or negative densities. */
    double  deadzone_r_inner;               /**< Inner radius of the Dead Zone in AU. */
    double  deadzone_r_outer;               /**< Outer radius of the Dead Zone in AU. */
    double  deadzone_dr_inner;              /**< Smoothing width for the inner Dead Zone edge in AU. */
    double  deadzone_dr_outer;              /**< Smoothing width for the outer Dead Zone edge in AU. */
    double  deadzone_alpha_mod;             /**< Reduced alpha viscosity within the Dead Zone. */
    double  dust_to_gas_ratio;              /**< Global dust-to-gas mass ratio, \f$\epsilon\f$. */
    double  total_disk_mass;                /**< Total mass of the two components of the disk. */
    double  one_size_particle_cm;           /**< Size of the large dust population in \f$cm\f$. */
    double  two_pop_ratio;                  /**< Mass ratio between the two dust populations. */
    double  micro_size_cm;                  /**< Size of the small (micron) dust population in \f$cm\f$. */
    double  drift_factor;                   /**< Scaling factor for radial drift velocity. */
    double  fragmentation_factor;           /**< Efficiency factor for dust fragmentation. */
    char output_base_path[MAX_PATH_LEN];    /**< Base directory path for simulation output. */
    double dust_density_g_cm3;              /**< Intrinsic material density of dust grains (e.g., \f 1.6g/cm^3\f$). */
} InitializeDefaultOptions;

/**
 * @brief populates the InitializeDefaultOptions structure with hardcoded defaults.
 * * Useful for ensuring all parameters have safe initial values before 
 * reading from a configuration file.
 * * @param[out] Pointer to the options structure to be initialized.
 */
void initializeDefaultOptions(InitializeDefaultOptions *opt);

/**
 * @brief Executes the full initialization sequence.
 * * Converts the raw options into a physical disk state, populating the 
 * internal arrays for density, pressure, and velocity.
 * * @param[in]  opts                Pointer to the configured initialization options.
 * @param[out] output_disk_params  Pointer to the disk parameters structure to be populated.
 * @return                     Returns 0 on success, non-zero on failure.
 */
int runInitialization(InitializeDefaultOptions *opts, DiskParameters *output_disk_params, SimulationOptions *sim_opts);

#endif // INIT_TOOL_MODULE_H