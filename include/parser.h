#ifndef PARSER_H
#define PARSER_H

#include <stdbool.h>      // For bool type
#include "config.h"       // Keep this for other config definitions if needed
#include "simulation_types.h" // CRITICAL: This line provides MAX_PATH_LEN

// Define a maximum length for the output directory name
// Use MAX_PATH_LEN from simulation_types.h for consistency
#define MAX_OUTPUT_DIR_LEN MAX_PATH_LEN

/*
 * @brief Structure to hold all parsed command-line options.
 * This replaces many individual global variables for options.
 */
typedef struct {
    // Simulation control options
    double drift;
    double growth;
    double evol;
    double twopop;
    double ufrag;
    double ffrag;

    // Core disk parameters (also serve as init_tool defaults if no input file)
    int    ngrid_val;       // Number of grid points
    int    ndust_val;
    double rmin_val;        // Inner radius (AU)
    double rmax_val;        // Outer radius (AU)
    double sigma0_val;      // Initial gas surface density at 1 AU (M_sun/AU^2)
    double sigmap_exp_val;  // Exponent of surface density profile (positive value, e.g., 1.0 for r^-1)
    double alpha_visc_val;  // Alpha viscosity
    double star_val;        // Star mass (M_sun)
    double hasp_val;        // Aspect ratio at 1 AU (H/R)
    double flind_val;       // Flaring index
    double r_dze_i_val;     // Inner dead zone radius (AU)
    double r_dze_o_val;     // Outer dead zone radius (AU)
    double dr_dze_i_val;    // Inner dead zone transition width multiplier
    double dr_dze_o_val;    // Outer dead zone transition width multiplier
    double a_mod_val;       // Alpha viscosity multiplier in dead zone

    // File input/output
    const char *input_file;      // Path to the input profile file, NULL if not provided
    char output_dir_name[MAX_OUTPUT_DIR_LEN]; // Name of the output directory

    // Time parameters
    double tStep;           // Fixed time step
    double totalTime;       // Total simulation time
    double outputFrequency; // Output frequency
    double startTime;       // Current start time (for resuming or initial time)

    // Init tool specific parameters (if initial profile is generated)
    // These could also be merged into the general disk parameters if they always overlap
    double eps_val;         // Dust-to-gas ratio
    double ratio_val;       // Ratio of Pop1 dust mass to total dust mass
    double mic_val;         // Micro-sized particle radius (cm)
    double onesize_val;     // Flag/value for one size particles (0.0 for distribution, 1.0 for mic_val)

    
    // NEW: Add PDENSITY (dust particle density) parameter
    double pdensity_val;    // Dust particle density [g/cm^3]

} options_t;

/*
 * @brief Initializes the options_t structure with default values.
 * @param opt Pointer to the options_t structure to be initialized.
 */
void create_default_options(options_t *opt);

/*
 * @brief Parses command-line arguments and populates the options_t structure.
 * @param argc The number of command-line arguments.
 * @param argv An array of strings containing the command-line arguments.
 * @param opt Pointer to the options_t structure where parsed values will be stored.
 * @return 0 on successful parsing, 1 on error (e.g., missing value, unknown option).
 */
int parse_options(int argc, const char **argv, options_t *opt);

/*
 * @brief Prints the command-line usage information to stderr.
 */
void print_usage();

#endif // PARSER_H