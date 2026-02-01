#include "init_tool_module.h"
#include "config.h" 
#include "disk_model.h" 
#include "dust_physics.h" 
#include "utils.h" 
#include "io_utils.h" 
#include "gas_physics.h"
#include "boundary_conditions.h"
#include "simulation_types.h"


#include <stdio.h>
#include <stdlib.h>
#include <string.h> // For strcpy, snprintf
#include <math.h>   // For M_PI, pow, fabs, sqrt, tanh, log

// --- External declaration for the existing linearInterpolation function ---
// This function is assumed to be implemented in utils.c (or similar)
// and its prototype should be in utils.h.
//extern void linearInterpolation(double *invec, double *radial_grid, double pos, double *out, double rd, int opt, const DiskParameters *disk_params);


void initializeDefaultOptions(InitializeDefaultOptions *def) {
    // Set default values for the InitializeDefaultOptions struct
    def->n_grid_points = 1000; // Number of radial gas grid points
    def->n_dust_particles = 2000; // NEW: Number of initial dust particles
    def->r_inner = 0.1;
    def->r_outer = 5.0;
    def->sigma0_gas_au = 0.01; // Gas surface density at 1 AU [M_Sun/AU^2] (increased from 0.0001)
    def->sigma_exponent = 0.5;   // CORRECTED: Positive exponent for Sigma ~ r^(-index) form
    def->alpha_viscosity = 1.0e-2; // Alpha viscosity parameter
    def->star_mass = 1.0;    // Central star mass [M_Sun]
    def->aspect_ratio = 5.0e-2; // Disk aspect ratio (H/r)
    def->flaring_index = 0.0;    // Flaring index for disk height (H ~ r^(1+flind))

    // Dead Zone Parameters (0.0 implies inactive)
    def->deadzone_r_inner = 0.0;
    def->deadzone_r_outer = 0.0;
    def->deadzone_dr_inner = 0.0; // Transition width multiplier
    def->deadzone_dr_outer = 0.0; // Transition width multiplier
    def->deadzone_alpha_mod = 0.01; // Alpha reduction factor in dead zone

    // Dust Parameters
    def->dust_to_gas_ratio = 0.01; // Initial dust-to-gas ratio (epsilon)
    def->disk_mass_dust = 0.0100000001; // Total dust disk mass [M_Sun]
    def->one_size_particle_cm = 1.0; // If > 0, particles are fixed to this size
    def->two_pop_ratio = 0.85; // Ratio of mass in larger particles for two-population model
    def->micro_size_cm = 1e-4; // Size of micron-sized particles for two-population model
    def->f_drift = 1.0;  // Factor for drift-limited size (default value, adjust as needed)
    def->f_frag = 1.0;   // Factor for fragmentation-limited size (default value, adjust as needed)

    def->output_base_path[0] = '\0'; // This will be set by main.c
    def->dust_density_g_cm3 = 1.6; // NEW: Dust particle density (g/cm^3) - default value
}

// Calculates the gas surface density normalization constant (Sigma0)
// based on total dust disk mass (Md), in M_Sun / AU / AU.
static long double calculateSigm0FromDiskMass(InitializeDefaultOptions *init_opts) {
    // Sigma ~ r^(-index), so integral is r^(2-index)
    double exponent_for_integral = -init_opts->sigma_exponent + 2.0;

    double denominator;
    if (fabs(exponent_for_integral) < 1e-9) { // Handle case where exponent is close to 0 (logarithmic integral)
        // This case implies Sigma ~ r^(-2), integral is ln(r)
        denominator = log(init_opts->r_outer) - log(init_opts->r_inner);
    } else {
        denominator = (pow(init_opts->r_outer, exponent_for_integral) - pow(init_opts->r_inner, exponent_for_integral)) / exponent_for_integral;
    }

    if (fabs(denominator) < 1e-12) {
        fprintf(stderr, "Error: Denominator is zero or too small in Sigma0 calculation! Check r_max, r_min, and SIGMA_EXP values.\n");
        return 0.0;
    }
    // Md = 2 * PI * Sigma0 * epsilon * Integral(r^(1-index) dr) from Ri to Ro
    // So, Sigma0 = Md / (2 * PI * epsilon * Integral(...))
    // Integral(r^(1-index) dr) = [r^(2-index)] / (2-index) from Ri to Ro
    // = (Ro^(2-index) - Ri^(2-index)) / (2-index)
    // Hence, Sigma0 = Md * (2-index) / (2 * PI * epsilon * (Ro^(2-index) - Ri^(2-index)))
    // Or rewritten: Sigma0 = Md * (2-index) / (2 * PI * epsilon * (Ro^(2-index) - Ri^(2-index)))
    return (long double)init_opts->disk_mass_dust /
           (2.0 * M_PI * init_opts->dust_to_gas_ratio * denominator);
}

// Calculates gas surface density at radial position r [M_Sun / AU / AU].
static long double calculateGasSurfaceDensityInitTool(double r_au, InitializeDefaultOptions *init_opts, long double current_sigma0) {

    return current_sigma0 * pow(r_au, init_opts->sigma_exponent);
}

// Calculates dust surface density at radial position r [M_Sun / AU / AU].
// Includes handling for snowline/ice factor if enabled (currently commented out as in original).
static long double calculateDustSurfaceDensityInitTool(double r_au, InitializeDefaultOptions *init_opts, long double current_sigma0) {

    long double sigma_dust = calculateGasSurfaceDensityInitTool(r_au, init_opts, current_sigma0) * init_opts->dust_to_gas_ratio;

    // --- Snowline and Ice Factor Handling (Uncomment and adjust if needed) ---
    // if (r_au >= SNOWLINE_RADIUS_AU) {
    //     sigma_dust *= ICE_LINE_DUST_ENHANCEMENT_FACTOR;
    // }
    // ---------------------------------------------------------------------

    return sigma_dust;
}

// Finds the minimum of three double values.
static double findMinimumForThreeNumbersInitTool(double value1, double value2, double value3) {

        return fmin(value1, fmin(value2, value3));
}


static int validateInitializationInputs(const InitializeDefaultOptions *opts) {

    if (opts->r_inner <= 0.0) {
        fprintf(stderr,
            "ERROR [runInitialization]: Inner radius (r_inner) must be positive. Current value: %lg\n",
            opts->r_inner);
        return 1;
    }

    if (opts->r_outer <= opts->r_inner) {
        fprintf(stderr,
            "ERROR [runInitialization]: Outer radius (r_outer) must be greater than inner radius (r_inner). "
            "R_inner: %lg, R_outer: %lg\n",
            opts->r_inner, opts->r_outer);
        return 1;
    }

    if (opts->n_grid_points <= 0) {
        fprintf(stderr,
            "ERROR [runInitialization]: Number of gas grid points (n_grid_points) must be positive. "
            "Current value: %d\n",
            opts->n_grid_points);
        return 1;
    }

    if (opts->n_dust_particles <= 0) {
        fprintf(stderr,
            "ERROR [runInitialization]: Number of dust particles (n_dust_particles) must be positive. "
            "Current value: %d\n",
            opts->n_dust_particles);
        return 1;
    }

    // Dead zone INNER: if radius is set but width is missing → error
    if (opts->deadzone_r_inner > 0.0 && opts->deadzone_dr_inner <= 0.0) {
        fprintf(stderr,
            "ERROR [runInitialization]: deadzone_r_inner is set (%lg), "
            "but deadzone_dr_inner is zero or missing. Please provide a non-zero width.\n",
            opts->deadzone_r_inner);
        return 1;
    }

    // Dead zone OUTER: if radius is set but width is missing → error
    if (opts->deadzone_r_outer > 0.0 && opts->deadzone_dr_outer <= 0.0) {
        fprintf(stderr,
            "ERROR [runInitialization]: deadzone_r_outer is set (%lg), "
            "but deadzone_dr_outer is zero or missing. Please provide a non-zero width.\n",
            opts->deadzone_r_outer);
        return 1;
    }


    return 0;   // All good
}



// --- Main Init Tool Function ---

int runInitialization(InitializeDefaultOptions *default_options, DiskParameters *disk_params) {

    FILE *dust_ouputput_file = NULL; // For dust particle profile
    FILE *disk_parameters_output_file = NULL; // For disk parameters
    FILE *gas_parameters_output_file = NULL; // For gas density profile

    long double current_sigma0_gas;

    // --- 0. Input validation ---
    if (validateInitializationInputs(default_options) != 0) { 
        return 1; 
    }

    // --- 1. Dead zone width calculations (now correctly placed) ---
    double drdze_inner_calculated =
        pow(default_options->deadzone_r_inner, 1.0 + default_options->flaring_index)
        * default_options->aspect_ratio
        * default_options->deadzone_dr_inner;

    double drdze_outer_calculated =
        pow(default_options->deadzone_r_outer, 1.0 + default_options->flaring_index)
        * default_options->aspect_ratio
        * default_options->deadzone_dr_outer;


    const double DEFAULT_disk_mass_DUST = 0.01;
    if (fabs(default_options->disk_mass_dust - DEFAULT_disk_mass_DUST) > 1e-9) {
        current_sigma0_gas = calculateSigm0FromDiskMass(default_options);
        fprintf(stderr,"Sigma0 calculated from total dust disk mass (Md): %Lg M_Sun/AU^2\n", current_sigma0_gas);
    } else {
        current_sigma0_gas = default_options->sigma0_gas_au;
        fprintf(stderr,"Using explicit Sigma0 (gas surface density at 1 AU): %Lg M_Sun/AU^2\n", current_sigma0_gas);
    }

    const double DEFAULT_ONE_SIZE = 1.0;
    if (fabs(default_options->one_size_particle_cm - DEFAULT_ONE_SIZE) > 1e-9 && default_options->one_size_particle_cm > 0) {
        default_options->two_pop_ratio = 1.0;
    }

    char full_init_dust_profile_path[MAX_PATH_LEN];
    char full_disk_param_path[MAX_PATH_LEN];
    char full_init_density_path[MAX_PATH_LEN];

    snprintf(full_init_dust_profile_path, sizeof(full_init_dust_profile_path), "%s/%s%s", default_options->output_base_path, kInitialDustProfileFileName,kFileNamesSuffix);
    snprintf(full_disk_param_path, sizeof(full_disk_param_path), "%s/%s%s", default_options->output_base_path, kDiskConfigFile,kFileNamesSuffix);
    snprintf(full_init_density_path, sizeof(full_init_density_path), "%s/%s%s", default_options->output_base_path, kInitialGasProfileFileName,kFileNamesSuffix);

    // --- Open Output Files with full paths ---
    dust_ouputput_file = fopen(full_init_dust_profile_path, "w");
    if (dust_ouputput_file == NULL) {
        perror("Error opening initial dust profile file in init_tool_module");
        return 1;
    }

    disk_parameters_output_file = fopen(full_disk_param_path, "w");
    if (disk_parameters_output_file == NULL) {
        perror("Error opening disk parameters file in init_tool_module");
        fclose(dust_ouputput_file);
        return 1;
    }

    gas_parameters_output_file = fopen(full_init_density_path, "w");
    if (gas_parameters_output_file == NULL) {
        perror("Error opening initial density profile file in init_tool_module");
        fclose(dust_ouputput_file);
        fclose(disk_parameters_output_file);
        return 1;
    }

    fprintf(stderr,"\n--- Simulation Parameters ---\n");
    fprintf(stderr,"Total dust disk mass (Solar Mass): %lg\n", default_options->disk_mass_dust);
    fprintf(stderr,"Inner disk edge (AU): %lg\n", default_options->r_inner);
    fprintf(stderr,"Outer disk edge (AU): %lg\n", default_options->r_outer);
    fprintf(stderr,"Surface density profile exponent: %lg\n", -default_options->sigma_exponent);
//    fprintf(stderr,"Snowline position (AU): %lg\n", SNOWLINE_RADIUS_AU);
//    fprintf(stderr,"Ice factor beyond snowline: %lg\n", ICE_LINE_DUST_ENHANCEMENT_FACTOR);
    fprintf(stderr,"Gas surface density at 1 AU (Solar Mass/AU^2): %Lg\n", current_sigma0_gas);
    fprintf(stderr,"Dust to gas ratio: %lg\n", default_options->dust_to_gas_ratio);
    fprintf(stderr,"Number of gas grid points: %d\n", default_options->n_grid_points); // Clarified for gas grid
    fprintf(stderr,"Number of dust particles to generate: %d\n", default_options->n_dust_particles); // NEW: Corrected for dust particles
    fprintf(stderr,"Dust particle density (g/cm^3): %lg\n", default_options->dust_density_g_cm3);
    fprintf(stderr,"------------------------------\n\n");

    // --- Prepare HeaderData for initial files ---
    HeaderData initial_header_data;
    initial_header_data.current_time = 0.0;
    initial_header_data.is_initial_data = 1;
    initial_header_data.R_in = default_options->r_inner;
    initial_header_data.R_out = default_options->r_outer;
    initial_header_data.sigma_exponent = default_options->sigma_exponent;
    initial_header_data.sigma0_gas_au = current_sigma0_gas;
    initial_header_data.grav_const = G_DIMENSIONLESS;
    initial_header_data.dz_r_inner = default_options->deadzone_r_inner;
    initial_header_data.dz_r_outer = default_options->deadzone_r_outer;
    initial_header_data.dz_dr_inner_calc = drdze_inner_calculated;
    initial_header_data.dz_dr_outer_calc = drdze_outer_calculated;
    initial_header_data.dz_alpha_mod = default_options->deadzone_alpha_mod;
    initial_header_data.dust_density_g_cm3 = default_options->dust_density_g_cm3;
    initial_header_data.alpha_viscosity = default_options->alpha_viscosity;
    initial_header_data.star_mass = default_options->star_mass;
    initial_header_data.flaring_index = default_options->flaring_index;
    initial_header_data.n_grid_points = default_options->n_grid_points;
    // initial_header_data.n_dust_particles = default_options->n_dust_particles; // This line was removed in the previous fix

    // --- Write Headers using the new function ---
    if (dust_ouputput_file != NULL) {
        printFileHeader(dust_ouputput_file, FILE_TYPE_PARTICLE_SIZE, &initial_header_data);
    }
    if (gas_parameters_output_file != NULL) {
        printFileHeader(gas_parameters_output_file, FILE_TYPE_GAS_DENSITY, &initial_header_data);
    }
    if (disk_parameters_output_file != NULL) {
        printFileHeader(disk_parameters_output_file, FILE_TYPE_DISK_PARAM, &initial_header_data);
    }

    // Physical Constants
    double u_frag_cm_s = 1000.0;
    double u_frag_au_yr2pi = u_frag_cm_s * CM_PER_SEC_TO_AU_PER_YEAR_2PI;
    double u_frag_sq_au_yr2pi_sq = u_frag_au_yr2pi * u_frag_au_yr2pi;

    // Populate disk_params structure and allocate its arrays
    disk_params->grid_number = default_options->n_grid_points; // Gas grid resolution
    disk_params->r_min = default_options->r_inner;
    disk_params->r_max = default_options->r_outer;
    disk_params->sigma_0 = current_sigma0_gas;
    disk_params->sigma_power_law_index = default_options->sigma_exponent;
    disk_params->r_dze_i = default_options->deadzone_r_inner;
    disk_params->r_dze_o = default_options->deadzone_r_outer;
    disk_params->dr_dze_i = default_options->deadzone_dr_inner;
    disk_params->dr_dze_o = default_options->deadzone_dr_outer;
    disk_params->alpha_parameter = default_options->alpha_viscosity;
    disk_params->alpha_parameter_modification = default_options->deadzone_alpha_mod;
    disk_params->h_aspect_ratio = default_options->aspect_ratio;
    disk_params->flaring_index = default_options->flaring_index;
    disk_params->stellar_mass = default_options->star_mass;
    disk_params->particle_density = default_options->dust_density_g_cm3;
    if (disk_params->grid_number > 1) {
        disk_params->delta_r = (disk_params->r_max - disk_params->r_min) / ((double)disk_params->grid_number - 1.0);
    } else {
        disk_params->delta_r = 0.0;
    }
    // Allocate arrays for the GAS grid (based on grid_number)
    disk_params->radial_grid = (double *)malloc((disk_params->grid_number + 2) * sizeof(double));
    disk_params->gas_surface_density_vector = (double *)malloc((disk_params->grid_number + 2) * sizeof(double));
    disk_params->gas_pressure_vector = (double *)malloc((disk_params->grid_number + 2) * sizeof(double));
    disk_params->gas_pressure_gradient_vector = (double *)malloc((disk_params->grid_number + 2) * sizeof(double));
    disk_params->gas_velocity_vector = (double *)malloc((disk_params->grid_number + 2) * sizeof(double));

    if (!disk_params->radial_grid || !disk_params->gas_surface_density_vector || !disk_params->gas_pressure_vector || !disk_params->gas_pressure_gradient_vector || !disk_params->gas_velocity_vector) {
        fprintf(stderr, "ERROR [runInitialization]: Failed to allocate disk arrays. Exiting.\n");
        if (dust_ouputput_file) fclose(dust_ouputput_file);
        if (disk_parameters_output_file) fclose(disk_parameters_output_file);
        if (gas_parameters_output_file) fclose(gas_parameters_output_file);
        // Free already allocated memory before exiting
        if (disk_params->radial_grid) free(disk_params->radial_grid);
        if (disk_params->gas_surface_density_vector) free(disk_params->gas_surface_density_vector);
        if (disk_params->gas_pressure_vector) free(disk_params->gas_pressure_vector);
        if (disk_params->gas_pressure_gradient_vector) free(disk_params->gas_pressure_gradient_vector);
        if (disk_params->gas_velocity_vector) free(disk_params->gas_velocity_vector);
        return 1;
    }

    // Initialize gas disk properties on the grid_number grid
    readDiskParameters(disk_params);
    createRadialGrid(disk_params);
    createInitialGasSurfaceDensity(disk_params);
    createInitialGasPressure(disk_params);
    createInitialGasPressureGradient(disk_params);
    createInitialGasVelocity(disk_params);

    // --- NEW SECTION: Write Gas Density Profile (to gas_parameters_output_file) ---
    // This loop iterates over the gas grid points (default_options->n_grid_points)
    for (int i_loop = 0; i_loop < default_options->n_grid_points; i_loop++) {
        double r_gas_grid_au = disk_params->radial_grid[i_loop + 1]; // Use 1-indexed for physical grid

        if (r_gas_grid_au <= 0) {
            fprintf(stderr, "ERROR: Calculated gas grid radial position is non-positive at index %d (%lg AU). Skipping this point.\n", i_loop, r_gas_grid_au);
            continue;
        }

        // Get gas properties directly from the disk_params arrays for the gas grid
        long double sigma_gas_local_val = disk_params->gas_surface_density_vector[i_loop + 1];
        double pressure_local_val = disk_params->gas_pressure_vector[i_loop + 1];
        double dPdr_local_val = disk_params->gas_pressure_gradient_vector[i_loop + 1];

        fprintf(gas_parameters_output_file, "%-15.6e %-15.6Lg %-15.6e %-15.6e\n",
            r_gas_grid_au,
            sigma_gas_local_val,
            pressure_local_val,
            dPdr_local_val);
    }
    fflush(gas_parameters_output_file); // Flush after writing gas data

    // --- EXISTING SECTION: Loop through the DESIRED NUMBER OF DUST PARTICLES to calculate and write their data ---
    // This loop iterates based on default_options->n_dust_particles
    for (int i_loop = 0; i_loop < default_options->n_dust_particles; i_loop++) {
        // Calculate the radial position for the current dust particle
        double r_dust_particle_au;
        if (default_options->n_dust_particles > 1) {
            r_dust_particle_au = default_options->r_inner + (default_options->r_outer - default_options->r_inner) * i_loop / ((double)default_options->n_dust_particles - 1.0);
        } else { // Handle case for single particle
            r_dust_particle_au = default_options->r_inner; // Or (r_inner + r_outer) / 2.0; depending on desired behavior
        }

        if (r_dust_particle_au <= 0) {
            fprintf(stderr, "ERROR: Calculated radial position is non-positive at dust particle index %d (%lg AU). Skipping this point.\n", i_loop, r_dust_particle_au);
            continue;
        }

        // --- Interpolate gas disk properties at the dust particle's radial position using the existing 'linearInterpolation' function ---
        double temp_sigma, temp_pressure, temp_dPdr;

        linearInterpolation(disk_params->gas_surface_density_vector, disk_params->radial_grid, r_dust_particle_au, &temp_sigma, disk_params->delta_r, 0, disk_params);
        long double sigma_gas_local = temp_sigma;

        linearInterpolation(disk_params->gas_pressure_vector, disk_params->radial_grid, r_dust_particle_au, &temp_pressure, disk_params->delta_r, 0, disk_params);
        double pressure_local = temp_pressure;

        linearInterpolation(disk_params->gas_pressure_gradient_vector, disk_params->radial_grid, r_dust_particle_au, &temp_dPdr, disk_params->delta_r, 0, disk_params);
        double dPdr_local = temp_dPdr;


        double s_max_cm;

        if (fabs(default_options->one_size_particle_cm - DEFAULT_ONE_SIZE) > 1e-9 && default_options->one_size_particle_cm > 0) {
            s_max_cm = default_options->one_size_particle_cm;
        } else {
            // Use linearInterpolationated gas properties for calculations
            // double H_au = calculatePressureScaleHeight(r_dust_particle_au, disk_params); // Removed: unused
            double calculateKeplerianVelocity_au_yr2pi = calculateKeplerianVelocity(r_dust_particle_au, disk_params);
            // double omega_yr2pi = calculateKeplerianFrequency(r_dust_particle_au, disk_params); // Removed: unused
            double sound_speed_au_yr2pi = calculateLocalSoundSpeed(r_dust_particle_au, disk_params);
            double sound_speed_sq = sound_speed_au_yr2pi * sound_speed_au_yr2pi;

            long double sigma_dust_local = calculateDustSurfaceDensityInitTool(r_dust_particle_au, default_options, current_sigma0_gas);
            long double sigma_dust_local_cgs = sigma_dust_local / SURFACE_DENSITY_CONVERSION_FACTOR;
            double sigma_gas_local_cgs = (double)sigma_gas_local / SURFACE_DENSITY_CONVERSION_FACTOR;


            double dlnPdlnr_local;
            if (fabs(pressure_local) < 1e-12) {
                fprintf(stderr, "Error: Pressure is near zero in dlnPdlnr calculation at r = %lg. Check input parameters. Setting dlnPdlnr = 0.\n", r_dust_particle_au);
                dlnPdlnr_local = 0.0;
            } else {
                dlnPdlnr_local = r_dust_particle_au / pressure_local * dPdr_local;
            }

            double s_drift = default_options->f_drift * 2.0 / M_PI * sigma_dust_local_cgs / default_options->dust_density_g_cm3 *
                             (calculateKeplerianVelocity_au_yr2pi * calculateKeplerianVelocity_au_yr2pi) / sound_speed_sq * fabs(1.0 / dlnPdlnr_local);

            double s_frag = default_options->f_frag * 2.0 / (3.0 * M_PI) * sigma_gas_local_cgs /
                            (default_options->dust_density_g_cm3 * calculateTurbulentAlpha(r_dust_particle_au, disk_params)) *
                            u_frag_sq_au_yr2pi_sq / sound_speed_sq;

            double dlnPdlnr_abs_cs2_half = fabs(dlnPdlnr_local * sound_speed_sq * 0.5);
            double s_df;
            if (dlnPdlnr_abs_cs2_half < 1e-12) {
                fprintf(stderr, "Error: Denominator is near zero in s_df calculation at r = %lg. Check dlnPdlnr value. Setting s_df to a large value.\n", r_dust_particle_au);
                s_df = 1e99;
            } else {
                s_df = u_frag_au_yr2pi * calculateKeplerianVelocity_au_yr2pi / dlnPdlnr_abs_cs2_half * 2.0 * sigma_gas_local_cgs / (M_PI * default_options->dust_density_g_cm3);
            }

            s_max_cm = findMinimumForThreeNumbersInitTool(s_drift, s_frag, s_df);
        }

        if (s_max_cm <= 0) {
            fprintf(stderr, "Warning: s_max_cm <= 0 at r = %lg. This might indicate problematic physical parameters. Setting to a small positive value.\n", r_dust_particle_au);
            s_max_cm = 1e-10;
        }

        long double representative_mass_total_in_cell = 2.0 * M_PI * r_dust_particle_au *
                                                        ((default_options->r_outer - default_options->r_inner) / ((double)default_options->n_dust_particles - 1.0)) * // Use dust particle spacing
                                                        calculateDustSurfaceDensityInitTool(r_dust_particle_au, default_options, current_sigma0_gas);

        long double repr_mass_pop1 = representative_mass_total_in_cell * default_options->two_pop_ratio;
        long double repr_mass_pop2 = representative_mass_total_in_cell * (1.0 - default_options->two_pop_ratio);

        // Write data to kInitialDustProfileFileName (dust_ouputput_file)
        fprintf(dust_ouputput_file, "%-5d %-15.6e %-20.12Lg %-20.12Lg %-15.6e %-15.6e\n",
                i_loop, r_dust_particle_au, // Use r_dust_particle_au
                repr_mass_pop1,
                repr_mass_pop2,
                s_max_cm, default_options->micro_size_cm);
    }

    fflush(dust_ouputput_file);
    fclose(dust_ouputput_file);
    dust_ouputput_file = NULL;

    // gas_parameters_output_file is already flushed above, just close it here
    fclose(gas_parameters_output_file);
    gas_parameters_output_file = NULL;

    fprintf(stderr,"Particle data file created (%s). Writing disk parameters file!\n\n", full_init_dust_profile_path);

    // Write parameters to kDiskConfigFile (disk_parameters_output_file)
    // CORRECTION: Negate SigmaExp value so that the actual negative exponent is written to the file.
    // Here, we are not calling printFileHeader, because data writing is a separate line, not part of the header.
    // The header was already written by the printFileHeader(disk_parameters_output_file, FILE_TYPE_DISK_PARAM, &initial_header_data); call above.
    fprintf(disk_parameters_output_file, "%-15.6e %-15.6e %-10d %-15.6e %-20.12Lg %-15.6e %-15.6e %-15.6e %-20.12e %-20.12e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e\n",
            default_options->r_inner, default_options->r_outer, default_options->n_grid_points, -default_options->sigma_exponent, current_sigma0_gas,
            G_DIMENSIONLESS, default_options->deadzone_r_inner, default_options->deadzone_r_outer,
            drdze_inner_calculated, drdze_outer_calculated,
            default_options->deadzone_alpha_mod, default_options->dust_density_g_cm3, default_options->alpha_viscosity, default_options->star_mass, default_options->flaring_index);

    fflush(disk_parameters_output_file);
    fclose(disk_parameters_output_file);
    disk_parameters_output_file = NULL;

    fprintf(stderr,"Disk parameters file created (%s).\n\n", full_disk_param_path);

    return 0;
}