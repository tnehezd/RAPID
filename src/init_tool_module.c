#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h> // For 'true'/'false'

#include "init_tool_module.h"
#include "config.h" // For SDCONV, G_GRAV_CONST, SNOWLINE, ICEFACTOR, CMPSECTOAUPYRP2PI

// --- Helper Functions (Static for internal use) ---

// Sets default options for the initialization tool.
void create_default_init_tool_options(init_tool_options_t *opt) {
    // Grid and Physical Parameters
    opt->n_grid_points   = 1000;
    opt->r_inner         = 0.1;
    opt->r_outer         = 5.0;
    // JAVÍTVA: Megnövelt sigma0 alapértelmezett érték
    opt->sigma0_gas_au   = 0.01; // Gas surface density at 1 AU [M_Sun/AU^2] (increased from 0.0001)
    opt->sigma_exponent  = 0.5;    // Positive exponent for surface density profile (Sigma ~ r^(-index))
    opt->alpha_viscosity = 1.0e-2; // Alpha viscosity parameter
    opt->star_mass       = 1.0;    // Central star mass [M_Sun]
    opt->aspect_ratio    = 5.0e-2; // Disk aspect ratio (H/r)
    opt->flaring_index   = 0.0;    // Flaring index for disk height (H ~ r^(1+flind))

    // Dead Zone Parameters (0.0 implies inactive)
    opt->deadzone_r_inner    = 0.0;
    opt->deadzone_r_outer    = 0.0;
    opt->deadzone_dr_inner   = 0.0; // Transition width multiplier
    opt->deadzone_dr_outer   = 0.0; // Transition width multiplier
    opt->deadzone_alpha_mod  = 0.01; // Alpha reduction factor in dead zone

    // Dust Parameters
    opt->dust_to_gas_ratio = 0.01; // Initial dust-to-gas ratio (epsilon)
    opt->disk_mass_dust    = 0.01; // Total dust disk mass [M_Sun]
    opt->one_size_particle_cm = 1.0; // If > 0, particles are fixed to this size
    opt->two_pop_ratio     = 0.85; // Ratio of mass in larger particles for two-population model
    opt->micro_size_cm     = 1e-4; // Size of micron-sized particles for two-population model
    opt->f_drift           = 1.0;  // Factor for drift-limited size (default value, adjust as needed)
    opt->f_frag            = 1.0;  // Factor for fragmentation-limited size (default value, adjust as needed)
}

// Calculates the turbulent alpha parameter, considering dead zones.
// Uses parameters directly from init_opts.
static double calculate_turbulent_alpha(double r_au, init_tool_options_t *init_opts) {
    // Calculate effective transition widths, ensuring non-zero denominator for tanh.
    // If drdze_i/o is 0, use a small default to prevent division by zero in tanh argument.
    // The original code used 1e-6; a more physically based non-zero default might be init_opts->h * r_au
    // or a fixed small fraction of the grid cell size, but for now we retain 1e-6.
    double drdze_inner_eff = pow(init_opts->deadzone_r_inner, 1.0 + init_opts->flaring_index) * init_opts->aspect_ratio *
                             ((init_opts->deadzone_dr_inner == 0.0) ? 1e-6 : init_opts->deadzone_dr_inner);
    double drdze_outer_eff = pow(init_opts->deadzone_r_outer, 1.0 + init_opts->flaring_index) * init_opts->aspect_ratio *
                             ((init_opts->deadzone_dr_outer == 0.0) ? 1e-6 : init_opts->deadzone_dr_outer);

    double alpha_reduction_factor = 1.0 - 0.5 * (1.0 - init_opts->deadzone_alpha_mod) *
                                   (tanh((r_au - init_opts->deadzone_r_inner) / drdze_inner_eff) +
                                    tanh((init_opts->deadzone_r_outer - r_au) / drdze_outer_eff));

    return alpha_reduction_factor * init_opts->alpha_viscosity;
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
    // Hence, Sigma0 = Md / (2 * PI * epsilon * (Ro^(2-index) - Ri^(2-index)) / (2-index))
    // Or rewritten: Sigma0 = Md * (2-index) / (2 * PI * epsilon * (Ro^(2-index) - Ri^(2-index)))
    return (long double)init_opts->disk_mass_dust /
           (2.0 * M_PI * init_opts->dust_to_gas_ratio * denominator);
}

// Calculates gas surface density at radial position r [M_Sun / AU / AU].
static long double calculate_gas_surface_density(double r_au, init_tool_options_t *init_opts, long double current_sigma0) {
    return current_sigma0 * pow(r_au, -init_opts->sigma_exponent);
}

// Calculates dust surface density at radial position r [M_Sun / AU / AU].
// Includes handling for snowline/ice factor if enabled (currently commented out as in original).
static long double calculate_dust_surface_density(double r_au, init_tool_options_t *init_opts, long double current_sigma0) {
    long double sigma_dust = calculate_gas_surface_density(r_au, init_opts, current_sigma0) * init_opts->dust_to_gas_ratio;

    // --- Snowline and Ice Factor Handling (Uncomment and adjust if needed) ---
    // if (r_au >= SNOWLINE) {
    //      sigma_dust *= ICEFACTOR;
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

int run_init_tool(init_tool_options_t *init_opts) {
    FILE *fout_data = NULL;
    FILE *fout_params = NULL;

    long double current_sigma0_gas; // Local variable for the determined gas sigma0

    // Calculated dead zone transition widths (for output in disk_param.dat)
    double drdze_inner_calculated = pow(init_opts->deadzone_r_inner, 1.0 + init_opts->flaring_index) * init_opts->aspect_ratio *
                                    ((init_opts->deadzone_dr_inner == 0.0) ? 1e-6 : init_opts->deadzone_dr_inner);
    double drdze_outer_calculated = pow(init_opts->deadzone_r_outer, 1.0 + init_opts->flaring_index) * init_opts->aspect_ratio *
                                    ((init_opts->deadzone_dr_outer == 0.0) ? 1e-6 : init_opts->deadzone_dr_outer);

    // Determine current_sigma0_gas: if disk_mass_dust is explicitly set (i.e., not default 0.01),
    // calculate Sigma0 from it. Otherwise, use the explicit sigma0_gas_au value.
    const double DEFAULT_DISK_MASS_DUST = 0.01; // Define this as a constant if not already
    
    // JAVÍTVA: pontosabb kiírás a sigma0 értékének forrásáról
    if (fabs(init_opts->disk_mass_dust - DEFAULT_DISK_MASS_DUST) > 1e-9) {
        current_sigma0_gas = calculate_sigma0_from_disk_mass(init_opts);
        printf("Sigma0 calculated from total dust disk mass (Md): %Lg M_Sun/AU^2\n", current_sigma0_gas);
    } else {
        current_sigma0_gas = init_opts->sigma0_gas_au;
        printf("Using explicit Sigma0 (gas surface density at 1 AU): %Lg M_Sun/AU^2\n", current_sigma0_gas);
    }

    // If one_size_particle_cm is set (i.e., not default 1.0), override two_pop_ratio to 1.0 (single population).
    const double DEFAULT_ONE_SIZE = 1.0; // Define this as a constant if not already
    if (fabs(init_opts->one_size_particle_cm - DEFAULT_ONE_SIZE) > 1e-9) {
        init_opts->two_pop_ratio = 1.0;
    }

    printf("Surface density profile exponent: %lg\n", -init_opts->sigma_exponent);

    // --- Open Output Files ---
    fout_data = fopen("init_data.dat", "w");
    if (fout_data == NULL) {
        perror("Error opening init_data.dat in init_tool_module");
        return 1;
    }

    fout_params = fopen("disk_param.dat", "w");
    if (fout_params == NULL) {
        perror("Error opening disk_param.dat in init_tool_module");
        fclose(fout_data); // Close already opened file
        return 1;
    }

    // Calculate grid spacing (delta R for linear grid)
    double delta_r_grid = (init_opts->r_outer - init_opts->r_inner) / ((double)init_opts->n_grid_points - 1.0); // N grid points means N-1 intervals

    // --- Print Simulation Parameters to Console ---
    printf("\n--- Simulation Parameters ---\n");
    printf("Total dust disk mass (Solar Mass): %lg\n", init_opts->disk_mass_dust);
    printf("Inner disk edge (AU): %lg\n", init_opts->r_inner);
    printf("Outer disk edge (AU): %lg\n", init_opts->r_outer);
    printf("Surface density profile exponent: %lg\n", -init_opts->sigma_exponent);
    printf("Snowline position (AU): %lg\n", SNOWLINE); // Assuming SNOWLINE is a global constant
    printf("Ice factor beyond snowline: %lg\n", ICEFACTOR); // Assuming ICEFACTOR is a global constant
    printf("Gas surface density at 1 AU (Solar Mass/AU^2): %Lg\n", current_sigma0_gas);
    printf("Dust to gas ratio: %lg\n", init_opts->dust_to_gas_ratio);
    printf("Number of representative particles: %d\n", init_opts->n_grid_points);
    printf("------------------------------\n\n");

    // --- Write Header to init_data.dat ---
    fprintf(fout_data, "# Initial Particle Data\n");
    fprintf(fout_data, "# Generated by init_tool_module (Date: %s %s)\n", __DATE__, __TIME__);
    fprintf(fout_data, "#--------------------------------------------------------------------------\n");
    fprintf(fout_data, "# %-5s %-15s %-20s %-20s %-15s %-15s\n",
              "Index", "Radius_AU", "RepMass_Pop1_Msun", "RepMass_Pop2_Msun", "MaxPartSize_cm", "MicroSize_cm");
    fprintf(fout_data, "#--------------------------------------------------------------------------\n");

    // --- Physical Constants (often defined globally or in a separate constants file) ---
    double dust_particle_density_g_cm3 = 1.6; // Por sűrűsége (g/cm^3)
    double u_frag_cm_s = 1000.0;             // Ütközési sebesség (cm/s)
    double u_frag_au_yr2pi = u_frag_cm_s * CMPSECTOAUPYRP2PI; // Convert to AU / (yr/2pi)
    double u_frag_sq_au_yr2pi_sq = u_frag_au_yr2pi * u_frag_au_yr2pi;

    // --- Loop through grid points to calculate and write particle data ---
    for (int i = 0; i < init_opts->n_grid_points; i++) {
        // Calculate radial position at cell center.
        double r_cell_center_au = init_opts->r_inner + i * delta_r_grid + delta_r_grid / 2.0;
        double r_cell_center_sq_au = r_cell_center_au * r_cell_center_au;

        // Representative particle mass for the current cell
        long double representative_mass = 2.0 * M_PI * r_cell_center_au * delta_r_grid *
                                         calculate_dust_surface_density(r_cell_center_au, init_opts, current_sigma0_gas);

        double s_max_cm; // Maximum particle size in cm

        // Handle one_size particle case
        if (fabs(init_opts->one_size_particle_cm - DEFAULT_ONE_SIZE) > 1e-9) { // Check if one_size is explicitly set
            s_max_cm = init_opts->one_size_particle_cm;
        } else {
            // Calculate disk height (H) and Keplerian velocity (v_kep)
            double H_au = pow(r_cell_center_au, 1.0 + init_opts->flaring_index) * init_opts->aspect_ratio;
            double v_kep_au_yr2pi = sqrt(G_GRAV_CONST * init_opts->star_mass / r_cell_center_au); // G_GRAV_CONST should be in AU^3 / (M_Sun * (yr/2pi)^2)
            double omega_yr2pi = v_kep_au_yr2pi / r_cell_center_au; // Omega = v_kep / r

            double sound_speed_au_yr2pi = omega_yr2pi * H_au;
            double sound_speed_sq = sound_speed_au_yr2pi * sound_speed_au_yr2pi;

            // Gas surface density and convert to CGS
            long double sigma_gas_local = calculate_gas_surface_density(r_cell_center_au, init_opts, current_sigma0_gas);
            double sigma_gas_local_cgs = (double)sigma_gas_local / SDCONV;

            // Dust surface density and convert to CGS
            long double sigma_dust_local = calculate_dust_surface_density(r_cell_center_au, init_opts, current_sigma0_gas);
            long double sigma_dust_local_cgs = sigma_dust_local / SDCONV;

            // Midplane gas density and pressure
            double rho_midplane_gas_local = 1.0 / sqrt(2.0 * M_PI) * sigma_gas_local / H_au;
            double pressure_local = rho_midplane_gas_local * sound_speed_sq;

            // Calculate dP/dr (radial pressure gradient)
            // This assumes P ~ r^(flind - index - 2) as in the original code's dPdr.
            // Original: (FLIND + SIGMAP_EXP - 2) * pow(reval,(FLIND + SIGMAP_EXP - 3.0)) * HASP * G2 * STAR * SIGMA0 / sqrt(2.0 * M_PI);
            // Where SIGMAP_EXP is negative. In modern code, sigma_exponent is positive, so -sigma_exponent.
            // P ~ rho * cs^2 ~ (Sigma/H) * (Omega * H)^2 ~ Sigma * Omega^2 * H
            // Sigma ~ r^(-sigma_exponent)
            // H ~ r^(1+flaring_index)
            // Omega^2 ~ r^(-3)
            // P ~ r^(-sigma_exponent) * r^(-3) * r^(1+flaring_index) = r^(-sigma_exponent - 2 + flaring_index)
            // Derivative of r^N is N * r^(N-1)
            double exponent_for_pressure = -init_opts->sigma_exponent + init_opts->flaring_index - 2.0;

            // Reconstructing dPdr based on P ~ r^(exponent_for_pressure)
            // P = C * r^(exponent_for_pressure) where C is the constant part
            // C = (1.0 / sqrt(2.0 * M_PI)) * current_sigma0_gas * G_GRAV_CONST * init_opts->star_mass * init_opts->aspect_ratio;
            // The original dPdr expression has G2, HASP, SIGMA0 in it, suggesting it's derived from a direct formula for pressure.
            // Let's stick to the modern code's pressure calculation (rho_mp * sound_speed_sq) for consistency,
            // and use numerical differentiation if analytic derivation is complex, or
            // directly use the original's structure for dPdr, ensuring consistent exponents.
            // Original dPdr term: (FLIND + SIGMAP_EXP - 2)
            // SIGMAP_EXP in original code is already negative.
            // In init_tool_options_t, sigma_exponent is POSITIVE. So we use -init_opts->sigma_exponent.
            // (init_opts->flaring_index + (-init_opts->sigma_exponent) - 2.0)
            double dPdr_local = (init_opts->flaring_index - init_opts->sigma_exponent - 2.0) *
                                 pow(r_cell_center_au, (init_opts->flaring_index - init_opts->sigma_exponent - 3.0)) *
                                 init_opts->aspect_ratio * G_GRAV_CONST * init_opts->star_mass * current_sigma0_gas / sqrt(2.0 * M_PI);


            double dlnPdlnr_local;
            if (fabs(pressure_local) < 1e-12) { // Check for near-zero pressure
                fprintf(stderr, "Error: Pressure is near zero in dlnPdlnr calculation at r = %lg. Check input parameters. Setting dlnPdlnr = 0.\n", r_cell_center_au);
                dlnPdlnr_local = 0.0;
            } else {
                dlnPdlnr_local = r_cell_center_au / pressure_local * dPdr_local;
            }

            // Calculate limiting particle sizes based on physics
            // NOW USING init_opts->f_drift and init_opts->f_frag
            double s_drift = init_opts->f_drift * 2.0 / M_PI * sigma_dust_local_cgs / dust_particle_density_g_cm3 *
                             (v_kep_au_yr2pi * v_kep_au_yr2pi) / sound_speed_sq * fabs(1.0 / dlnPdlnr_local);
            
            // Replaced alpha_turb_it with direct call to calculate_turbulent_alpha
            double s_frag = init_opts->f_frag * 2.0 / (3.0 * M_PI) * sigma_gas_local_cgs /
                            (dust_particle_density_g_cm3 * calculate_turbulent_alpha(r_cell_center_au, init_opts)) *
                            u_frag_sq_au_yr2pi_sq / sound_speed_sq;

            double dlnPdlnr_abs_cs2_half = fabs(dlnPdlnr_local * sound_speed_sq * 0.5);
            double s_df; // Shear fragmentation (or related) size
            if (dlnPdlnr_abs_cs2_half < 1e-12) {
                fprintf(stderr, "Error: Denominator is near zero in s_df calculation at r = %lg. Check dlnPdlnr value. Setting s_df to a large value.\n", r_cell_center_au);
                s_df = 1e99; // Set to a very large value so it doesn't limit s_max
            } else {
                s_df = u_frag_au_yr2pi * v_kep_au_yr2pi / dlnPdlnr_abs_cs2_half * 2.0 * sigma_gas_local_cgs / (M_PI * dust_particle_density_g_cm3);
            }

            s_max_cm = find_minimum_double(s_drift, s_frag, s_df);
        }

        // Ensure s_max_cm is positive
        if (s_max_cm <= 0) {
            fprintf(stderr, "Warning: s_max_cm <= 0 at r = %lg. This might indicate problematic physical parameters. Setting to a small positive value.\n", r_cell_center_au);
            s_max_cm = 1e-10;
        }

        // Write data to init_data.dat
        fprintf(fout_data, "%-5d %-15.6e %-20.12Lg %-20.12Lg %-15.6e %-15.6e\n",
                i, r_cell_center_au,
                representative_mass * init_opts->two_pop_ratio,
                representative_mass * (1.0 - init_opts->two_pop_ratio),
                s_max_cm, init_opts->micro_size_cm);
    }

    fclose(fout_data);

    printf("Particle data file created (init_data.dat). Writing disk parameters file!\n\n");
    // Removed getchar(); to allow non-interactive runs

    // --- Write Header to disk_param.dat ---
    fprintf(fout_params, "# Disk Parameters\n");
    fprintf(fout_params, "# Generated by init_tool_module (Date: %s %s)\n", __DATE__, __TIME__);
    fprintf(fout_params, "#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
    fprintf(fout_params, "# %-15s %-15s %-10s %-15s %-20s %-15s %-15s %-15s %-20s %-20s %-15s %-15s %-15s %-15s %-15s\n",
              "R_Min_AU", "R_Max_AU", "N_Grid", "SigmaExp", "Sigma0_gas_Msun_AU2",
              "G_GravConst", "DzR_Inner_AU", "DzR_Outer_AU", "DzDr_Inner_Calc_AU", "DzDr_Outer_Calc_AU",
              "DzAlphaMod", "DustDensity_g_cm3", "AlphaViscosity", "StarMass_Msun", "FlaringIndex");
    fprintf(fout_params, "#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");


    // Write parameters to disk_param.dat
    fprintf(fout_params, "%-15.6e %-15.6e %-10d %-15.6e %-20.12Lg %-15.6e %-15.6e %-15.6e %-20.12e %-20.12e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e\n",
              init_opts->r_inner, init_opts->r_outer, init_opts->n_grid_points, init_opts->sigma_exponent, current_sigma0_gas,
              G_GRAV_CONST, init_opts->deadzone_r_inner, init_opts->deadzone_r_outer,
              drdze_inner_calculated, drdze_outer_calculated,
              init_opts->deadzone_alpha_mod, dust_particle_density_g_cm3, init_opts->alpha_viscosity, init_opts->star_mass, init_opts->flaring_index);

    fclose(fout_params);

    printf("Disk parameters file created (disk_param.dat).\n\n");

    return 0;
}