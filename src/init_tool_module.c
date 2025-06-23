#include "init_tool_module.h"
#include "config.h" // For FILENAME_INIT_PROFILE, FILENAME_DISK_PARAM, SDCONV, G_GRAV_CONST, SNOWLINE, ICEFACTOR, CMPSECTOAUPYRP2PI
#include "disk_model.h" // Contains declarations for load_R, Initial_Profile, Initial_Press, Initial_dPress, Initial_Ugas, disk_param_be, scale_height, v_kep, kep_freq, c_sound, press, rho_mp
#include "dust_physics.h" // May contain GetMass, etc.
#include "utils.h" // For find_max, find_min, etc.
#include "io_utils.h" // For Mk_Dir (if used internally here)

#include "simulation_types.h"


#include <stdio.h>
#include <stdlib.h>
#include <string.h> // For strcpy, snprintf
#include <math.h>   // For M_PI, pow, fabs, sqrt, tanh, log

// Removed redundant extern declarations assuming they are in disk_model.h or dust_physics.h
// and are correctly included by the #include directives above.

void create_default_init_tool_options(init_tool_options_t *def) {
    // Set default values for the init_tool_options_t struct
    def->n_grid_points = 1000;
    def->r_inner = 0.1;
    def->r_outer = 5.0;
    def->sigma0_gas_au = 0.01; // Gas surface density at 1 AU [M_Sun/AU^2] (increased from 0.0001)
    def->sigma_exponent = 0.5;   // JAVÍTVA: Pozitív exponent a Sigma ~ r^(-index) formához
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
    def->dust_density_g_cm3 = 1.6; // NEW: Por szemcse sűrűsége (g/cm^3) - alapértelmezett érték
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
    // Hence, Sigma0 = Md * (2-index) / (2 * PI * epsilon * (Ro^(2-index) - Ri^(2-index)))
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

int run_init_tool(init_tool_options_t *opts) {
    FILE *fout_data = NULL;
    FILE *fout_params = NULL;

    long double current_sigma0_gas; // Local variable for the determined gas sigma0

    // Calculated dead zone transition widths (for output in disk_param.dat)
    double drdze_inner_calculated = pow(opts->deadzone_r_inner, 1.0 + opts->flaring_index) * opts->aspect_ratio *
                                    ((opts->deadzone_dr_inner == 0.0) ? 1e-6 : opts->deadzone_dr_inner);
    double drdze_outer_calculated = pow(opts->deadzone_r_outer, 1.0 + opts->flaring_index) * opts->aspect_ratio *
                                    ((opts->deadzone_dr_outer == 0.0) ? 1e-6 : opts->deadzone_dr_outer);

    // Determine current_sigma0_gas: if disk_mass_dust is explicitly set (i.e., not default 0.01),
    // calculate Sigma0 from it. Otherwise, use the explicit sigma0_gas_au value.
    const double DEFAULT_DISK_MASS_DUST = 0.01; // Define this as a constant if not already
    
    // JAVÍTVA: pontosabb kiírás a sigma0 értékének forrásáról
    if (fabs(opts->disk_mass_dust - DEFAULT_DISK_MASS_DUST) > 1e-9) {
        current_sigma0_gas = calculate_sigma0_from_disk_mass(opts);
        printf("Sigma0 calculated from total dust disk mass (Md): %Lg M_Sun/AU^2\n", current_sigma0_gas);
    } else {
        current_sigma0_gas = opts->sigma0_gas_au;
        printf("Using explicit Sigma0 (gas surface density at 1 AU): %Lg M_Sun/AU^2\n", current_sigma0_gas);
    }

    // If one_size_particle_cm is set (i.e., not default 1.0), override two_pop_ratio to 1.0 (single population).
    const double DEFAULT_ONE_SIZE = 1.0; // Define this as a constant if not already
    if (fabs(opts->one_size_particle_cm - DEFAULT_ONE_SIZE) > 1e-9 && opts->one_size_particle_cm > 0) { // Check if one_size is explicitly set and positive
        opts->two_pop_ratio = 1.0;
    }

    printf("Surface density profile exponent: %lg\n", -opts->sigma_exponent); // Itt továbbra is a negált értéket írjuk ki a konzolra, hogy tükrözze a Sigma ~ r^(-exponent) formát.

    // --- NEW: Construct full file paths ---
    char full_init_profile_path[MAX_PATH_LEN];
    char full_disk_param_path[MAX_PATH_LEN];

    snprintf(full_init_profile_path, sizeof(full_init_profile_path), "%s/%s", opts->output_base_path, FILENAME_INIT_PROFILE);
    snprintf(full_disk_param_path, sizeof(full_disk_param_path), "%s/%s", opts->output_base_path, FILENAME_DISK_PARAM);


    // --- Open Output Files with full paths ---
    fout_data = fopen(full_init_profile_path, "w");
    if (fout_data == NULL) {
        perror("Error opening initial dust profile file in init_tool_module");
        return 1;
    }

    fout_params = fopen(full_disk_param_path, "w");
    if (fout_params == NULL) {
        perror("Error opening disk parameters file in init_tool_module");
        fclose(fout_data); // Close already opened file
        return 1;
    }

    // Calculate grid spacing (delta R for linear grid)
    double delta_r_grid = (opts->r_outer - opts->r_inner) / ((double)opts->n_grid_points - 1.0); // N grid points means N-1 intervals

    // --- Print Simulation Parameters to Console ---
    printf("\n--- Simulation Parameters ---\n");
    printf("Total dust disk mass (Solar Mass): %lg\n", opts->disk_mass_dust);
    printf("Inner disk edge (AU): %lg\n", opts->r_inner);
    printf("Outer disk edge (AU): %lg\n", opts->r_outer);
    printf("Surface density profile exponent: %lg\n", -opts->sigma_exponent); // Negált érték kiírása
    printf("Snowline position (AU): %lg\n", SNOWLINE); // Assuming SNOWLINE is a global constant
    printf("Ice factor beyond snowline: %lg\n", ICEFACTOR); // Assuming ICEFACTOR is a global constant
    printf("Gas surface density at 1 AU (Solar Mass/AU^2): %Lg\n", current_sigma0_gas);
    printf("Dust to gas ratio: %lg\n", opts->dust_to_gas_ratio);
    printf("Number of representative particles: %d\n", opts->n_grid_points);
    printf("Dust particle density (g/cm^3): %lg\n", opts->dust_density_g_cm3); // NEW: PDENSITY kiírás
    printf("------------------------------\n\n");

    // --- Write Header to FILENAME_INIT_PROFILE ---
    fprintf(fout_data, "# Initial Particle Data\n");
    fprintf(fout_data, "# Generated by init_tool_module (Date: %s %s)\n", __DATE__, __TIME__);
    fprintf(fout_data, "#--------------------------------------------------------------------------\n");
    fprintf(fout_data, "# %-5s %-15s %-20s %-20s %-15s %-15s\n",
             "Index", "Radius_AU", "RepMass_Pop1_Msun", "RepMass_Pop2_Msun", "MaxPartSize_cm", "MicroSize_cm");
    fprintf(fout_data, "#--------------------------------------------------------------------------\n");

    // --- Physical Constants ---
    // Használd a opts->dust_density_g_cm3-at a hardkódolt 1.6 helyett!
    // double dust_particle_density_g_cm3 = 1.6; // Por sűrűsége (g/cm^3) - ELTÁVOLÍTVA
    double u_frag_cm_s = 1000.0;             // Ütközési sebesség (cm/s)
    double u_frag_au_yr2pi = u_frag_cm_s * CMPSECTOAUPYRP2PI; // Convert to AU / (yr/2pi)
    double u_frag_sq_au_yr2pi_sq = u_frag_au_yr2pi * u_frag_au_yr2pi;

    // --- Temporary disk_t structure to pass to disk_model functions ---
    disk_t temp_disk_params;
    temp_disk_params.NGRID = opts->n_grid_points;
    temp_disk_params.RMIN = opts->r_inner;
    temp_disk_params.RMAX = opts->r_outer;
    temp_disk_params.SIGMA0 = current_sigma0_gas; // Use the *determined* Sigma0
    temp_disk_params.SIGMAP_EXP = opts->sigma_exponent; // Itt a pozitív értéket tároljuk, ha ez a konvenció
    temp_disk_params.r_dze_i = opts->deadzone_r_inner;
    temp_disk_params.r_dze_o = opts->deadzone_r_outer;
    temp_disk_params.Dr_dze_i = opts->deadzone_dr_inner;
    temp_disk_params.Dr_dze_o = opts->deadzone_dr_outer;
    temp_disk_params.alpha_visc = opts->alpha_viscosity;
    temp_disk_params.a_mod = opts->deadzone_alpha_mod;
    temp_disk_params.HASP = opts->aspect_ratio;
    temp_disk_params.FLIND = opts->flaring_index;
    temp_disk_params.STAR_MASS = opts->star_mass;
    temp_disk_params.PDENSITY = opts->dust_density_g_cm3; // NEW: PDENSITY átadása temp_disk_params-nek
    // Calculate DD based on the current NGRID, RMIN, RMAX
    if (temp_disk_params.NGRID > 1) {
        temp_disk_params.DD = (temp_disk_params.RMAX - temp_disk_params.RMIN) / (temp_disk_params.NGRID - 1.0);
    } else {
        temp_disk_params.DD = 0.0;
    }
    // Allocate temporary arrays for this run
    temp_disk_params.rvec = (double *)malloc((temp_disk_params.NGRID + 2) * sizeof(double));
    temp_disk_params.sigmavec = (double *)malloc((temp_disk_params.NGRID + 2) * sizeof(double));
    temp_disk_params.pressvec = (double *)malloc((temp_disk_params.NGRID + 2) * sizeof(double));
    temp_disk_params.dpressvec = (double *)malloc((temp_disk_params.NGRID + 2) * sizeof(double));
    temp_disk_params.ugvec = (double *)malloc((temp_disk_params.NGRID + 2) * sizeof(double));

    if (!temp_disk_params.rvec || !temp_disk_params.sigmavec || !temp_disk_params.pressvec || !temp_disk_params.dpressvec || !temp_disk_params.ugvec) {
        fprintf(stderr, "ERROR [run_init_tool]: Failed to allocate temporary disk arrays. Exiting.\n");
        // Clean up already opened files before exit
        if (fout_data) fclose(fout_data);
        if (fout_params) fclose(fout_params);
        return 1;
    }

    // Populate remaining derived disk_params fields (PDENSITY and PDENSITYDIMLESS)
    // Megjegyzés: a disk_param_be most már feltételezi, hogy a PDENSITY a disk_t struktúra PDENSITY tagjába van töltve.
    disk_param_be(&temp_disk_params);


    // Populate rvec and sigmavec for temp_disk_params, needed for pressure/velocity calcs
    load_R(&temp_disk_params); // This function should populate temp_disk_params.rvec
    Initial_Profile(&temp_disk_params); // This function should populate temp_disk_params.sigmavec

    // KULCSFONTOSSÁGÚ JAVÍTÁS:
    // Hívjuk meg a nyomás és nyomásgradiens inicializáló függvényeket,
    // mielőtt a dpressvec-et használnánk az s_max_cm számításokhoz!
    Initial_Press(&temp_disk_params);
    Initial_dPress(&temp_disk_params);
    Initial_Ugas(&temp_disk_params); // Ez is ide illik logikailag

    // Ensure boundary conditions are applied if load_R/Initial_Profile don't handle them for the full range
    // Perem(temp_disk_params.rvec, &temp_disk_params); // Assuming Perem handles disk_t*
    // Perem(temp_disk_params.sigmavec, &temp_disk_params);
    // Valószínűleg ezeket a Perem hívásokat nem itt kellene megismételni,
    // mivel az Initial_Profile már meghívja a sigmavec-re, és a load_R is.


    // --- Loop through grid points to calculate and write particle data ---
    for (int i = 0; i < opts->n_grid_points; i++) {
        // Calculate radial position at cell center.
        // For 1-based indexing in disk_params.rvec used by other functions, use i+1
        double r_cell_center_au = temp_disk_params.rvec[i+1];
        // Ensure r_cell_center_au is within valid range to avoid issues with pow(negative, non-integer)
        if (r_cell_center_au <= 0) {
             fprintf(stderr, "ERROR: r_cell_center_au is non-positive at index %d (%lg). Skipping this point.\n", i, r_cell_center_au);
             continue; // Skip to next iteration
        }


        double s_max_cm; // Maximum particle size in cm

        // Handle one_size particle case
        if (fabs(opts->one_size_particle_cm - DEFAULT_ONE_SIZE) > 1e-9 && opts->one_size_particle_cm > 0) { // Check if one_size is explicitly set and positive
            s_max_cm = opts->one_size_particle_cm;
        } else {
            // Calculate disk height (H) and Keplerian velocity (v_kep) using temp_disk_params
            // Assuming these functions take disk_t* and use its relevant members
            double H_au = scale_height(r_cell_center_au, &temp_disk_params);
            double v_kep_au_yr2pi = v_kep(r_cell_center_au, &temp_disk_params);
            double omega_yr2pi = kep_freq(r_cell_center_au, &temp_disk_params);

            double sound_speed_au_yr2pi = c_sound(r_cell_center_au, &temp_disk_params);
            double sound_speed_sq = sound_speed_au_yr2pi * sound_speed_au_yr2pi;

            // Gas surface density and convert to CGS
            long double sigma_gas_local = calculate_gas_surface_density(r_cell_center_au, opts, current_sigma0_gas);
            double sigma_gas_local_cgs = (double)sigma_gas_local / SDCONV;

            // Dust surface density and convert to CGS
            long double sigma_dust_local = calculate_dust_surface_density(r_cell_center_au, opts, current_sigma0_gas);
            long double sigma_dust_local_cgs = sigma_dust_local / SDCONV;

            // Midplane gas density and pressure (using functions from dust_physics.h/c)
            double rho_midplane_gas_local = rho_mp((double)sigma_gas_local, r_cell_center_au, &temp_disk_params);
            double pressure_local = press((double)sigma_gas_local, r_cell_center_au, &temp_disk_params);

            // Calculate dP/dr (radial pressure gradient)
            // Most megbízhatóan a temp_disk_params.dpressvec-ből vesszük, feltételezve, hogy az Initial_dPress
            // már meghívásra került és feltöltötte.
            double dPdr_local;
            if (i + 1 >= 0 && i + 1 <= temp_disk_params.NGRID + 1) {
                 dPdr_local = temp_disk_params.dpressvec[i+1];
            } else {
                 // Fallback - this should ideally not be hit if Initial_dPress works correctly
                 fprintf(stderr, "WARNING: dpressvec index out of bounds for r=%lg. Estimating dPdr numerically.\n", r_cell_center_au);
                 if (i > 0 && i < opts->n_grid_points - 1) {
                     double r_minus = temp_disk_params.rvec[i];
                     double r_plus = temp_disk_params.rvec[i+2];
                     double p_minus = press((double)calculate_gas_surface_density(r_minus, opts, current_sigma0_gas), r_minus, &temp_disk_params);
                     double p_plus = press((double)calculate_gas_surface_density(r_plus, opts, current_sigma0_gas), r_plus, &temp_disk_params);
                     dPdr_local = (p_plus - p_minus) / (r_plus - r_minus);
                 } else { // Edge case or single point
                     dPdr_local = 0.0;
                 }
            }


            double dlnPdlnr_local;
            if (fabs(pressure_local) < 1e-12) { // Check for near-zero pressure
                fprintf(stderr, "Error: Pressure is near zero in dlnPdlnr calculation at r = %lg. Check input parameters. Setting dlnPdlnr = 0.\n", r_cell_center_au);
                dlnPdlnr_local = 0.0;
            } else {
                dlnPdlnr_local = r_cell_center_au / pressure_local * dPdr_local;
            }

            // Calculate limiting particle sizes based on physics
            // NOW USING opts->f_drift and opts->f_frag
            // Használjuk az opts->dust_density_g_cm3-at!
            double s_drift = opts->f_drift * 2.0 / M_PI * sigma_dust_local_cgs / opts->dust_density_g_cm3 *
                             (v_kep_au_yr2pi * v_kep_au_yr2pi) / sound_speed_sq * fabs(1.0 / dlnPdlnr_local);
            
            // Replaced alpha_turb_it with direct call to calculate_turbulent_alpha
            // Használjuk az opts->dust_density_g_cm3-at!
            double s_frag = opts->f_frag * 2.0 / (3.0 * M_PI) * sigma_gas_local_cgs /
                            (opts->dust_density_g_cm3 * calculate_turbulent_alpha(r_cell_center_au, opts)) *
                            u_frag_sq_au_yr2pi_sq / sound_speed_sq;

            double dlnPdlnr_abs_cs2_half = fabs(dlnPdlnr_local * sound_speed_sq * 0.5);
            double s_df; // Shear fragmentation (or related) size
            if (dlnPdlnr_abs_cs2_half < 1e-12) {
                fprintf(stderr, "Error: Denominator is near zero in s_df calculation at r = %lg. Check dlnPdlnr value. Setting s_df to a large value.\n", r_cell_center_au);
                s_df = 1e99; // Set to a very large value so it doesn't limit s_max
            } else {
                // Használjuk az opts->dust_density_g_cm3-at!
                s_df = u_frag_au_yr2pi * v_kep_au_yr2pi / dlnPdlnr_abs_cs2_half * 2.0 * sigma_gas_local_cgs / (M_PI * opts->dust_density_g_cm3);
            }

            s_max_cm = find_minimum_double(s_drift, s_frag, s_df);
        }

        // Ensure s_max_cm is positive
        if (s_max_cm <= 0) {
            fprintf(stderr, "Warning: s_max_cm <= 0 at r = %lg. This might indicate problematic physical parameters. Setting to a small positive value.\n", r_cell_center_au);
            s_max_cm = 1e-10;
        }

        // Calculate representative mass for each population
        long double representative_mass_total_in_cell = 2.0 * M_PI * r_cell_center_au * delta_r_grid *
                                                       calculate_dust_surface_density(r_cell_center_au, opts, current_sigma0_gas);

        long double repr_mass_pop1 = representative_mass_total_in_cell * opts->two_pop_ratio;
        long double repr_mass_pop2 = representative_mass_total_in_cell * (1.0 - opts->two_pop_ratio);

        // Write data to FILENAME_INIT_PROFILE
        fprintf(fout_data, "%-5d %-15.6e %-20.12Lg %-20.12Lg %-15.6e %-15.6e\n",
                i, r_cell_center_au,
                repr_mass_pop1,
                repr_mass_pop2,
                s_max_cm, opts->micro_size_cm);
    }

    fflush(fout_data);
    fclose(fout_data);
    fout_data = NULL; // Reset pointer after closing


    printf("Particle data file created (%s). Writing disk parameters file!\n\n", full_init_profile_path);
    // Removed getchar(); to allow non-interactive runs

    // --- Write Header to FILENAME_DISK_PARAM ---
    fprintf(fout_params, "# Disk Parameters\n");
    fprintf(fout_params, "# Generated by init_tool_module (Date: %s %s)\n", __DATE__, __TIME__);
    fprintf(fout_params, "#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
    fprintf(fout_params, "# %-15s %-15s %-10s %-15s %-20s %-15s %-15s %-15s %-20s %-20s %-15s %-15s %-15s %-15s %-15s\n",
              "R_Min_AU", "R_Max_AU", "N_Grid", "SigmaExp", "Sigma0_gas_Msun_AU2",
              "G_GravConst", "DzR_Inner_AU", "DzR_Outer_AU", "DzDr_Inner_Calc_AU", "DzDr_Outer_Calc_AU",
              "DzAlphaMod", "DustDensity_g_cm3", "AlphaViscosity", "StarMass_Msun", "FlaringIndex");
    fprintf(fout_params, "#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");


    // Write parameters to FILENAME_DISK_PARAM
    // JAVÍTÁS: A SigmaExp értékét negáljuk, hogy a fájlba a valódi, negatív kitevő kerüljön.
    fprintf(fout_params, "%-15.6e %-15.6e %-10d %-15.6e %-20.12Lg %-15.6e %-15.6e %-15.6e %-20.12e %-20.12e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e\n",
              opts->r_inner, opts->r_outer, opts->n_grid_points, -opts->sigma_exponent, current_sigma0_gas, // ITT A JAVÍTÁS
              G_GRAV_CONST, opts->deadzone_r_inner, opts->deadzone_r_outer,
              drdze_inner_calculated, drdze_outer_calculated,
              opts->deadzone_alpha_mod, opts->dust_density_g_cm3, opts->alpha_viscosity, opts->star_mass, opts->flaring_index); // Itt is opts->dust_density_g_cm3

    fflush(fout_params);
    fclose(fout_params);
    fout_params = NULL; // Reset pointer after closing

    // Free temporary disk_t arrays
    free(temp_disk_params.rvec);
    free(temp_disk_params.sigmavec);
    free(temp_disk_params.pressvec);
    free(temp_disk_params.dpressvec);
    free(temp_disk_params.ugvec);


    printf("Particle data file created (%s). Writing disk parameters file!\n\n", full_init_profile_path);
    printf("Disk parameters file created (%s).\n\n", full_disk_param_path);

    return 0; // Success
}