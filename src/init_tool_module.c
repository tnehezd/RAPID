#include "init_tool_module.h"
#include "config.h" // For FILENAME_INIT_DUST_PROFILE, FILENAME_DISK_PARAM, SDCONV, G_GRAV_CONST, SNOWLINE, ICEFACTOR, CMPSECTOAUPYRP2PI
#include "disk_model.h" // Contains declarations for load_R, Initial_Profile, Initial_Press, Initial_dPress, Initial_Ugas, disk_param_be, scale_height, v_kep, kep_freq, c_sound, press, rho_mp
#include "dust_physics.h" // May contain GetMass, etc.
#include "utils.h" // For find_max, find_min, etc. and for interpolation functions
#include "io_utils.h" // For Mk_Dir (if used internally here)

#include "simulation_types.h"


#include <stdio.h>
#include <stdlib.h>
#include <string.h> // For strcpy, snprintf
#include <math.h>   // For M_PI, pow, fabs, sqrt, tanh, log

// --- External declaration for the existing interpol function ---
// This function is assumed to be implemented in utils.c (or similar)
// and its prototype should be in utils.h.
//extern void interpol(double *invec, double *rvec, double pos, double *out, double rd, int opt, const disk_t *disk_params);


void create_default_init_tool_options(init_tool_options_t *def) {
    // Set default values for the init_tool_options_t struct
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
    def->disk_mass_dust = 0.01; // Total dust disk mass [M_Sun]
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
static long double calculate_sigma0_from_disk_mass(init_tool_options_t *init_opts) {
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
        fprintf(stderr, "Error: Denominator is zero or too small in Sigma0 calculation! Check RMAX, RMIN, and SIGMA_EXP values.\n");
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
static long double calculate_gas_surface_density(double r_au, init_tool_options_t *init_opts, long double current_sigma0) {
    return current_sigma0 * pow(r_au, init_opts->sigma_exponent);
}

// Calculates dust surface density at radial position r [M_Sun / AU / AU].
// Includes handling for snowline/ice factor if enabled (currently commented out as in original).
static long double calculate_dust_surface_density(double r_au, init_tool_options_t *init_opts, long double current_sigma0) {
    long double sigma_dust = calculate_gas_surface_density(r_au, init_opts, current_sigma0) * init_opts->dust_to_gas_ratio;

    // --- Snowline and Ice Factor Handling (Uncomment and adjust if needed) ---
    // if (r_au >= SNOWLINE) {
    //     sigma_dust *= ICEFACTOR;
    // }
    // ---------------------------------------------------------------------

    return sigma_dust;
}

// Finds the minimum of three double values.
static double find_minimum_double(double s1, double s2, double s3) {
    double min_val = s1;
    if (s2 < min_val) {
        min_val = s2;
    }
    if (s3 < min_val) {
        min_val = s3;
    }
    return min_val;
}

// --- Main Init Tool Function ---

int run_init_tool(init_tool_options_t *opts, disk_t *disk_params) {
    FILE *fout_data = NULL; // For dust particle profile
    FILE *fout_params = NULL; // For disk parameters
    FILE *fout_dens = NULL; // For gas density profile

    long double current_sigma0_gas;

    // --- Input parameter validation ---
    if (opts->r_inner <= 0.0) {
        fprintf(stderr, "ERROR [run_init_tool]: Inner radius (r_inner) must be positive. Current value: %lg\n", opts->r_inner);
        return 1;
    }
    if (opts->r_outer <= opts->r_inner) {
        fprintf(stderr, "ERROR [run_init_tool]: Outer radius (r_outer) must be greater than inner radius (r_inner). R_inner: %lg, R_outer: %lg\n", opts->r_inner, opts->r_outer);
        return 1;
    }
    if (opts->n_grid_points <= 0) {
        fprintf(stderr, "ERROR [run_init_tool]: Number of gas grid points (n_grid_points) must be positive. Current value: %d\n", opts->n_grid_points);
        return 1;
    }
    if (opts->n_dust_particles <= 0) { // NEW: Validate n_dust_particles
        fprintf(stderr, "ERROR [run_init_tool]: Number of dust particles (n_dust_particles) must be positive. Current value: %d\n", opts->n_dust_particles);
        return 1;
    }
    // --- End of input parameter validation ---


    double drdze_inner_calculated = pow(opts->deadzone_r_inner, 1.0 + opts->flaring_index) * opts->aspect_ratio *
                                    ((opts->deadzone_dr_inner == 0.0) ? 1e-6 : opts->deadzone_dr_inner);
    double drdze_outer_calculated = pow(opts->deadzone_r_outer, 1.0 + opts->flaring_index) * opts->aspect_ratio *
                                    ((opts->deadzone_dr_outer == 0.0) ? 1e-6 : opts->deadzone_dr_outer);

    const double DEFAULT_DISK_MASS_DUST = 0.01;
    if (fabs(opts->disk_mass_dust - DEFAULT_DISK_MASS_DUST) > 1e-9) {
        current_sigma0_gas = calculate_sigma0_from_disk_mass(opts);
        fprintf(stderr,"Sigma0 calculated from total dust disk mass (Md): %Lg M_Sun/AU^2\n", current_sigma0_gas);
    } else {
        current_sigma0_gas = opts->sigma0_gas_au;
        fprintf(stderr,"Using explicit Sigma0 (gas surface density at 1 AU): %Lg M_Sun/AU^2\n", current_sigma0_gas);
    }

    const double DEFAULT_ONE_SIZE = 1.0;
    if (fabs(opts->one_size_particle_cm - DEFAULT_ONE_SIZE) > 1e-9 && opts->one_size_particle_cm > 0) {
        opts->two_pop_ratio = 1.0;
    }

    char full_init_dust_profile_path[MAX_PATH_LEN];
    char full_disk_param_path[MAX_PATH_LEN];
    char full_init_density_path[MAX_PATH_LEN];

    snprintf(full_init_dust_profile_path, sizeof(full_init_dust_profile_path), "%s/%s", opts->output_base_path, FILENAME_INIT_DUST_PROFILE);
    snprintf(full_disk_param_path, sizeof(full_disk_param_path), "%s/%s", opts->output_base_path, FILENAME_DISK_PARAM);
    snprintf(full_init_density_path, sizeof(full_init_density_path), "%s/%s", opts->output_base_path, FILENAME_INIT_GAS_PROFILE);

    // --- Open Output Files with full paths ---
    fout_data = fopen(full_init_dust_profile_path, "w");
    if (fout_data == NULL) {
        perror("Error opening initial dust profile file in init_tool_module");
        return 1;
    }

    fout_params = fopen(full_disk_param_path, "w");
    if (fout_params == NULL) {
        perror("Error opening disk parameters file in init_tool_module");
        fclose(fout_data);
        return 1;
    }

    fout_dens = fopen(full_init_density_path, "w");
    if (fout_dens == NULL) {
        perror("Error opening initial density profile file in init_tool_module");
        fclose(fout_data);
        fclose(fout_params);
        return 1;
    }

    fprintf(stderr,"\n--- Simulation Parameters ---\n");
    fprintf(stderr,"Total dust disk mass (Solar Mass): %lg\n", opts->disk_mass_dust);
    fprintf(stderr,"Inner disk edge (AU): %lg\n", opts->r_inner);
    fprintf(stderr,"Outer disk edge (AU): %lg\n", opts->r_outer);
    fprintf(stderr,"Surface density profile exponent: %lg\n", -opts->sigma_exponent);
    fprintf(stderr,"Snowline position (AU): %lg\n", SNOWLINE);
    fprintf(stderr,"Ice factor beyond snowline: %lg\n", ICEFACTOR);
    fprintf(stderr,"Gas surface density at 1 AU (Solar Mass/AU^2): %Lg\n", current_sigma0_gas);
    fprintf(stderr,"Dust to gas ratio: %lg\n", opts->dust_to_gas_ratio);
    fprintf(stderr,"Number of gas grid points: %d\n", opts->n_grid_points); // Clarified for gas grid
    fprintf(stderr,"Number of dust particles to generate: %d\n", opts->n_dust_particles); // NEW: Corrected for dust particles
    fprintf(stderr,"Dust particle density (g/cm^3): %lg\n", opts->dust_density_g_cm3);
    fprintf(stderr,"------------------------------\n\n");

    // --- Prepare HeaderData_t for initial files ---
    HeaderData_t initial_header_data;
    initial_header_data.current_time = 0.0;
    initial_header_data.is_initial_data = 1;
    initial_header_data.R_in = opts->r_inner;
    initial_header_data.R_out = opts->r_outer;
    initial_header_data.sigma_exponent = opts->sigma_exponent;
    initial_header_data.sigma0_gas_au = current_sigma0_gas;
    initial_header_data.grav_const = G_GRAV_CONST;
    initial_header_data.dz_r_inner = opts->deadzone_r_inner;
    initial_header_data.dz_r_outer = opts->deadzone_r_outer;
    initial_header_data.dz_dr_inner_calc = drdze_inner_calculated;
    initial_header_data.dz_dr_outer_calc = drdze_outer_calculated;
    initial_header_data.dz_alpha_mod = opts->deadzone_alpha_mod;
    initial_header_data.dust_density_g_cm3 = opts->dust_density_g_cm3;
    initial_header_data.alpha_viscosity = opts->alpha_viscosity;
    initial_header_data.star_mass = opts->star_mass;
    initial_header_data.flaring_index = opts->flaring_index;
    initial_header_data.n_grid_points = opts->n_grid_points;
    // initial_header_data.n_dust_particles = opts->n_dust_particles; // This line was removed in the previous fix

    // --- Write Headers using the new function ---
    if (fout_data != NULL) {
        print_file_header(fout_data, FILE_TYPE_PARTICLE_SIZE, &initial_header_data);
    }
    if (fout_dens != NULL) {
        print_file_header(fout_dens, FILE_TYPE_GAS_DENSITY, &initial_header_data);
    }
    if (fout_params != NULL) {
        print_file_header(fout_params, FILE_TYPE_DISK_PARAM, &initial_header_data);
    }

    // Physical Constants
    double u_frag_cm_s = 1000.0;
    double u_frag_au_yr2pi = u_frag_cm_s * CMPSECTOAUPYRP2PI;
    double u_frag_sq_au_yr2pi_sq = u_frag_au_yr2pi * u_frag_au_yr2pi;

    // Populate disk_params structure and allocate its arrays
    disk_params->NGRID = opts->n_grid_points; // Gas grid resolution
    disk_params->RMIN = opts->r_inner;
    disk_params->RMAX = opts->r_outer;
    disk_params->SIGMA0 = current_sigma0_gas;
    disk_params->SIGMAP_EXP = opts->sigma_exponent;
    disk_params->r_dze_i = opts->deadzone_r_inner;
    disk_params->r_dze_o = opts->deadzone_r_outer;
    disk_params->Dr_dze_i = opts->deadzone_dr_inner;
    disk_params->Dr_dze_o = opts->deadzone_dr_outer;
    disk_params->alpha_visc = opts->alpha_viscosity;
    disk_params->a_mod = opts->deadzone_alpha_mod;
    disk_params->HASP = opts->aspect_ratio;
    disk_params->FLIND = opts->flaring_index;
    disk_params->STAR_MASS = opts->star_mass;
    disk_params->PDENSITY = opts->dust_density_g_cm3;
    if (disk_params->NGRID > 1) {
        disk_params->DD = (disk_params->RMAX - disk_params->RMIN) / ((double)disk_params->NGRID - 1.0);
    } else {
        disk_params->DD = 0.0;
    }
    // Allocate arrays for the GAS grid (based on NGRID)
    disk_params->rvec = (double *)malloc((disk_params->NGRID + 2) * sizeof(double));
    disk_params->sigmavec = (double *)malloc((disk_params->NGRID + 2) * sizeof(double));
    disk_params->pressvec = (double *)malloc((disk_params->NGRID + 2) * sizeof(double));
    disk_params->dpressvec = (double *)malloc((disk_params->NGRID + 2) * sizeof(double));
    disk_params->ugvec = (double *)malloc((disk_params->NGRID + 2) * sizeof(double));

    if (!disk_params->rvec || !disk_params->sigmavec || !disk_params->pressvec || !disk_params->dpressvec || !disk_params->ugvec) {
        fprintf(stderr, "ERROR [run_init_tool]: Failed to allocate disk arrays. Exiting.\n");
        if (fout_data) fclose(fout_data);
        if (fout_params) fclose(fout_params);
        if (fout_dens) fclose(fout_dens);
        // Free already allocated memory before exiting
        if (disk_params->rvec) free(disk_params->rvec);
        if (disk_params->sigmavec) free(disk_params->sigmavec);
        if (disk_params->pressvec) free(disk_params->pressvec);
        if (disk_params->dpressvec) free(disk_params->dpressvec);
        if (disk_params->ugvec) free(disk_params->ugvec);
        return 1;
    }

    // Initialize gas disk properties on the NGRID grid
    disk_param_be(disk_params);
    load_R(disk_params);
    Initial_Profile(disk_params);
    Initial_Press(disk_params);
    Initial_dPress(disk_params);
    Initial_Ugas(disk_params);

    // --- NEW SECTION: Write Gas Density Profile (to fout_dens) ---
    // This loop iterates over the gas grid points (opts->n_grid_points)
    for (int i_loop = 0; i_loop < opts->n_grid_points; i_loop++) {
        double r_gas_grid_au = disk_params->rvec[i_loop + 1]; // Use 1-indexed for physical grid

        if (r_gas_grid_au <= 0) {
            fprintf(stderr, "ERROR: Calculated gas grid radial position is non-positive at index %d (%lg AU). Skipping this point.\n", i_loop, r_gas_grid_au);
            continue;
        }

        // Get gas properties directly from the disk_params arrays for the gas grid
        long double sigma_gas_local_val = disk_params->sigmavec[i_loop + 1];
        double pressure_local_val = disk_params->pressvec[i_loop + 1];
        double dPdr_local_val = disk_params->dpressvec[i_loop + 1];

        fprintf(fout_dens, "%-15.6e %-15.6Lg %-15.6e %-15.6e\n",
            r_gas_grid_au,
            sigma_gas_local_val,
            pressure_local_val,
            dPdr_local_val);
    }
    fflush(fout_dens); // Flush after writing gas data

    // --- EXISTING SECTION: Loop through the DESIRED NUMBER OF DUST PARTICLES to calculate and write their data ---
    // This loop iterates based on opts->n_dust_particles
    for (int i_loop = 0; i_loop < opts->n_dust_particles; i_loop++) {
        // Calculate the radial position for the current dust particle
        double r_dust_particle_au;
        if (opts->n_dust_particles > 1) {
            r_dust_particle_au = opts->r_inner + (opts->r_outer - opts->r_inner) * i_loop / ((double)opts->n_dust_particles - 1.0);
        } else { // Handle case for single particle
            r_dust_particle_au = opts->r_inner; // Or (r_inner + r_outer) / 2.0; depending on desired behavior
        }

        if (r_dust_particle_au <= 0) {
            fprintf(stderr, "ERROR: Calculated radial position is non-positive at dust particle index %d (%lg AU). Skipping this point.\n", i_loop, r_dust_particle_au);
            continue;
        }

        // --- Interpolate gas disk properties at the dust particle's radial position using the existing 'interpol' function ---
        double temp_sigma, temp_pressure, temp_dPdr;

        interpol(disk_params->sigmavec, disk_params->rvec, r_dust_particle_au, &temp_sigma, disk_params->DD, 0, disk_params);
        long double sigma_gas_local = temp_sigma;

        interpol(disk_params->pressvec, disk_params->rvec, r_dust_particle_au, &temp_pressure, disk_params->DD, 0, disk_params);
        double pressure_local = temp_pressure;

        interpol(disk_params->dpressvec, disk_params->rvec, r_dust_particle_au, &temp_dPdr, disk_params->DD, 0, disk_params);
        double dPdr_local = temp_dPdr;


        double s_max_cm;

        if (fabs(opts->one_size_particle_cm - DEFAULT_ONE_SIZE) > 1e-9 && opts->one_size_particle_cm > 0) {
            s_max_cm = opts->one_size_particle_cm;
        } else {
            // Use interpolated gas properties for calculations
            // double H_au = scale_height(r_dust_particle_au, disk_params); // Removed: unused
            double v_kep_au_yr2pi = v_kep(r_dust_particle_au, disk_params);
            // double omega_yr2pi = kep_freq(r_dust_particle_au, disk_params); // Removed: unused
            double sound_speed_au_yr2pi = c_sound(r_dust_particle_au, disk_params);
            double sound_speed_sq = sound_speed_au_yr2pi * sound_speed_au_yr2pi;

            long double sigma_dust_local = calculate_dust_surface_density(r_dust_particle_au, opts, current_sigma0_gas);
            long double sigma_dust_local_cgs = sigma_dust_local / SDCONV;
            double sigma_gas_local_cgs = (double)sigma_gas_local / SDCONV;


            double dlnPdlnr_local;
            if (fabs(pressure_local) < 1e-12) {
                fprintf(stderr, "Error: Pressure is near zero in dlnPdlnr calculation at r = %lg. Check input parameters. Setting dlnPdlnr = 0.\n", r_dust_particle_au);
                dlnPdlnr_local = 0.0;
            } else {
                dlnPdlnr_local = r_dust_particle_au / pressure_local * dPdr_local;
            }

            double s_drift = opts->f_drift * 2.0 / M_PI * sigma_dust_local_cgs / opts->dust_density_g_cm3 *
                             (v_kep_au_yr2pi * v_kep_au_yr2pi) / sound_speed_sq * fabs(1.0 / dlnPdlnr_local);

            double s_frag = opts->f_frag * 2.0 / (3.0 * M_PI) * sigma_gas_local_cgs /
                            (opts->dust_density_g_cm3 * calculate_turbulent_alpha(r_dust_particle_au, disk_params)) *
                            u_frag_sq_au_yr2pi_sq / sound_speed_sq;

            double dlnPdlnr_abs_cs2_half = fabs(dlnPdlnr_local * sound_speed_sq * 0.5);
            double s_df;
            if (dlnPdlnr_abs_cs2_half < 1e-12) {
                fprintf(stderr, "Error: Denominator is near zero in s_df calculation at r = %lg. Check dlnPdlnr value. Setting s_df to a large value.\n", r_dust_particle_au);
                s_df = 1e99;
            } else {
                s_df = u_frag_au_yr2pi * v_kep_au_yr2pi / dlnPdlnr_abs_cs2_half * 2.0 * sigma_gas_local_cgs / (M_PI * opts->dust_density_g_cm3);
            }

            s_max_cm = find_minimum_double(s_drift, s_frag, s_df);
        }

        if (s_max_cm <= 0) {
            fprintf(stderr, "Warning: s_max_cm <= 0 at r = %lg. This might indicate problematic physical parameters. Setting to a small positive value.\n", r_dust_particle_au);
            s_max_cm = 1e-10;
        }

        long double representative_mass_total_in_cell = 2.0 * M_PI * r_dust_particle_au *
                                                        ((opts->r_outer - opts->r_inner) / ((double)opts->n_dust_particles - 1.0)) * // Use dust particle spacing
                                                        calculate_dust_surface_density(r_dust_particle_au, opts, current_sigma0_gas);

        long double repr_mass_pop1 = representative_mass_total_in_cell * opts->two_pop_ratio;
        long double repr_mass_pop2 = representative_mass_total_in_cell * (1.0 - opts->two_pop_ratio);

        // Write data to FILENAME_INIT_DUST_PROFILE (fout_data)
        fprintf(fout_data, "%-5d %-15.6e %-20.12Lg %-20.12Lg %-15.6e %-15.6e\n",
                i_loop, r_dust_particle_au, // Use r_dust_particle_au
                repr_mass_pop1,
                repr_mass_pop2,
                s_max_cm, opts->micro_size_cm);
    }

    fflush(fout_data);
    fclose(fout_data);
    fout_data = NULL;

    // fout_dens is already flushed above, just close it here
    fclose(fout_dens);
    fout_dens = NULL;

    fprintf(stderr,"Particle data file created (%s). Writing disk parameters file!\n\n", full_init_dust_profile_path);

    // Write parameters to FILENAME_DISK_PARAM (fout_params)
    // CORRECTION: Negate SigmaExp value so that the actual negative exponent is written to the file.
    // Here, we are not calling print_file_header, because data writing is a separate line, not part of the header.
    // The header was already written by the print_file_header(fout_params, FILE_TYPE_DISK_PARAM, &initial_header_data); call above.
    fprintf(fout_params, "%-15.6e %-15.6e %-10d %-15.6e %-20.12Lg %-15.6e %-15.6e %-15.6e %-20.12e %-20.12e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e\n",
            opts->r_inner, opts->r_outer, opts->n_grid_points, -opts->sigma_exponent, current_sigma0_gas,
            G_GRAV_CONST, opts->deadzone_r_inner, opts->deadzone_r_outer,
            drdze_inner_calculated, drdze_outer_calculated,
            opts->deadzone_alpha_mod, opts->dust_density_g_cm3, opts->alpha_viscosity, opts->star_mass, opts->flaring_index);

    fflush(fout_params);
    fclose(fout_params);
    fout_params = NULL;

    fprintf(stderr,"Disk parameters file created (%s).\n\n", full_disk_param_path);

    return 0;
}