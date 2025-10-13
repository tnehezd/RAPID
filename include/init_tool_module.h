#ifndef INIT_TOOL_MODULE_H
#define INIT_TOOL_MODULE_H

#include <stdbool.h>
#include "simulation_types.h"


/**
 * @struct init_tool_options
 * @brief  Structure to hold all configuration options for the initialization tool.
 *
 * This structure is used to pass all necessary parameters from the main program
 * to the `run_init_tool` function, which sets up the initial state of the
 * protoplanetary disk and dust particles.
 */
typedef struct init_tool_options {
    // Grid and Physical Parameters
    /**
     * @brief Number of radial grid points for the gas disk.
     *
     * This determines the resolution of the gas surface density and pressure profiles.
     */
    int     n_grid_points;

    /**
     * @brief Number of initial dust particles to generate.
     *
     * This defines how many individual particles are created, each with its own
     * properties (position, size, mass).
     */
    int     n_dust_particles;

    /**
     * @brief Inner radius of the disk in astronomical units (AU).
     */
    double  r_inner;

    /**
     * @brief Outer radius of the disk in astronomical units (AU).
     */
    double  r_outer;

    /**
     * @brief Gas surface density at 1 AU [M_Sun/AU^2].
     *
     * This value sets the overall mass scale of the gas disk.
     * It is used if `disk_mass_dust` is not specified.
     */
    double  sigma0_gas_au;

    /**
     * @brief Power-law exponent for the gas surface density profile.
     *
     * The profile follows the form Sigma(r) ~ r^(-sigma_exponent).
     * This value should be positive.
     */
    double  sigma_exponent;

    /**
     * @brief The alpha viscosity parameter.
     *
     * This value quantifies the turbulent viscosity of the gas disk,
     * which affects gas dynamics and particle-gas coupling.
     */
    double  alpha_viscosity;

    /**
     * @brief Mass of the central star in solar masses [M_Sun].
     */
    double  star_mass;

    /**
     * @brief The disk aspect ratio (H/r).
     *
     * This is the ratio of the disk's scale height (H) to the radial distance (r),
     * which is a measure of the disk's vertical thickness.
     */
    double  aspect_ratio;

    /**
     * @brief The flaring index for the disk scale height.
     *
     * The scale height H(r) scales with radius as H ~ r^(1 + flaring_index).
     */
    double  flaring_index;

    // Dead Zone Parameters
    /**
     * @brief Inner radius of the dead zone in AU.
     *
     * In this region, turbulence (and thus alpha viscosity) is significantly reduced.
     */
    double  deadzone_r_inner;

    /**
     * @brief Outer radius of the dead zone in AU.
     */
    double  deadzone_r_outer;

    /**
     * @brief Transition width multiplier for the inner dead zone edge.
     *
     * This value, typically in units of H/r, determines how smoothly the
     * alpha viscosity transitions at the inner dead zone boundary.
     */
    double  deadzone_dr_inner;

    /**
     * @brief Transition width multiplier for the outer dead zone edge.
     */
    double  deadzone_dr_outer;

    /**
     * @brief Alpha viscosity reduction factor within the dead zone.
     *
     * The alpha value is multiplied by this factor inside the dead zone.
     */
    double  deadzone_alpha_mod;

    // Dust Parameters
    /**
     * @brief The initial dust-to-gas ratio (epsilon).
     *
     * This determines the dust surface density relative to the gas.
     */
    double  dust_to_gas_ratio;

    /**
     * @brief The total mass of dust in the disk [M_Sun].
     *
     * If this value is greater than a default threshold, it is used to
     * calculate `sigma0_gas_au`, overriding any user-provided value.
     */
    double  disk_mass_dust;

    /**
     * @brief Fixed particle size in centimeters.
     *
     * If greater than zero, all dust particles will have this size, and
     * the two-population model will be disabled.
     */
    double  one_size_particle_cm;

    /**
     * @brief Mass ratio of larger particles in a two-population model (0.0 - 1.0).
     *
     * This specifies the fraction of the total dust mass that is assigned
     * to the larger particle size population.
     */
    double  two_pop_ratio;

    /**
     * @brief The size of the micron-sized particle population in centimeters.
     */
    double  micro_size_cm;

    /**
     * @brief Factor for the drift-limited particle size.
     *
     * This is a scaling factor for the maximum particle size determined
     * by radial drift  defined by Birnstiel et al 2012.
     */
    double  f_drift;

    /**
     * @brief Factor for the fragmentation-limited particle size.
     *
     * This is a scaling factor for the maximum particle size determined
     * by collisional fragmentation defined by Birnstiel et al 2012.
     */
    double  f_frag;

    /**
     * @brief Base path for output files.
     *
     * The `init_tool` will create its output files (e.g., dust profiles,
     * gas profiles) in this directory.
     */
    char output_base_path[MAX_PATH_LEN];

    /**
     * @brief Dust particle material density in g/cm^3.
     */
    double dust_density_g_cm3;

} init_tool_options_t;

/**
 * @brief Initializes an `init_tool_options_t` struct with default values.
 *
 * This function provides a baseline configuration, which can then be
 * overridden by user-specified values or command-line arguments.
 *
 * @param opt Pointer to the `init_tool_options_t` struct to be initialized.
 */
void create_default_init_tool_options(init_tool_options_t *opt);

/**
 * @brief Main function to run the initialization tool.
 *
 * This function performs the core logic of the tool: validating input options,
 * calculating the physical properties of the gas and dust disks, and writing
 * the initial conditions to output files.
 *
 * @param opts Pointer to the `init_tool_options_t` struct containing the
 * configuration options.
 * @param output_disk_params Pointer to a `disk_t` struct where the calculated
 * gas grid properties will be stored.
 * @return Returns 0 on success, or a non-zero value on error.
 */
int run_init_tool(init_tool_options_t *opts, disk_t *output_disk_params, simulation_options_t *simopts);


#endif // INIT_TOOL_MODULE_H
