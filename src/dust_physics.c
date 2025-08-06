// src/dust_physics.c

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <omp.h>

#include "disk_model.h"
#include "dust_physics.h"
#include "config.h"
#include "simulation_types.h"
#include "globals.h"
#include "io_utils.h"
#include "simulation_core.h"
#include "utils.h"
#include "particle_data.h"

// Calculating dust surface density (from Lagrangian to Eulerian)
void calculate_dust_surface_density_profile(
    double *output_sigma_d_grid,      // Output: The calculated dust surface density (per grid cell)
    double *output_r_grid_centers,    // Output: Radial positions of the grid cell centers
    const double radin[][2],          // Input: Particle radial positions (const)
    const double *massin,             // Input: Particle masses (const)
    int n_particles,                  // Input: Number of particles
    int n_grid_cells,                 // Input: Number of grid cells for density calculation
    const disk_t *disk_params         // Input: The disk parameters (const)
) {
    // Memory allocation and initialization

    // Temporary array to store mass collected in bins.
    // Important: due to the CIC method, a size of NGRID + 2 is needed (one extra bin at each edge)
    double *total_mass_in_bins = (double *)calloc(n_grid_cells + 2, sizeof(double));
    if (total_mass_in_bins == NULL) {
        fprintf(stderr, "ERROR: Memory allocation failed for total_mass_in_bins in calculate_dust_surface_density_profile.\n");
        exit(EXIT_FAILURE);
    }

    // Array to store the exact area of the bins.
    double *bin_areas = (double *)calloc(n_grid_cells + 2, sizeof(double));
    if (bin_areas == NULL) {
        fprintf(stderr, "ERROR: Memory allocation failed for bin_areas in calculate_dust_surface_density_profile.\n");
        free(total_mass_in_bins); // Free the first allocated memory block in case of an error
        exit(EXIT_FAILURE);
    }

    // The width of the grid cells (delta r).
    double dr_cell_width = (disk_params->RMAX - disk_params->RMIN) / (double)n_grid_cells;

    // 1. step: Calculate the area of the grid cells and populate the centers
    for (int j = 0; j < n_grid_cells; j++) {
        double r_inner = disk_params->RMIN + j * dr_cell_width;
        double r_outer = disk_params->RMIN + (j + 1) * dr_cell_width;

        // The area of the ring: pi * (R_outer^2 - R_inner^2)
        bin_areas[j] = M_PI * (r_outer * r_outer - r_inner * r_inner);

        // Fill the output grid points with the radii (cell centers)
        output_r_grid_centers[j] = r_inner + 0.5 * dr_cell_width;

        // Initialize the output density array with 0.
        output_sigma_d_grid[j] = 0.0;
    }

    // 2. step: Distribute particle masses among bins (Cloud-in-Cell - CIC)
    // Iterate through each dust particle. Can be parallelized.
    #pragma omp parallel for
    for (int i = 0; i < n_particles; i++) {
        double current_r = radin[i][0]; // The particle's current radial position
        double current_mass = massin[i]; // The particle's mass

        // Only consider particles that are within the simulation domain.
        // If a particle moved outside RMIN, it's discarded and does not contribute to the density.
        if (current_r < disk_params->RMIN || current_r >= disk_params->RMAX) {
            continue; // Jump to the next particle
        }

        // Calculate the particle's "normalized" radial position relative to the grid.
        double r_normalized = (current_r - disk_params->RMIN) / dr_cell_width;

        // Find the index of the lower bin where the particle is located.
        int lower_bin_idx = (int)floor(r_normalized);

        // Calculate how "inside" the particle is within the lower bin (a value between 0 and 1).
        double fraction_in_lower_bin = r_normalized - lower_bin_idx;

        // Add a portion of the mass to the lower bin.
        // A critical section is necessary because multiple threads can write to the same total_mass_in_bins element.
        if (lower_bin_idx >= 0 && lower_bin_idx < n_grid_cells) {
            #pragma omp atomic
            total_mass_in_bins[lower_bin_idx] += current_mass * (1.0 - fraction_in_lower_bin);
        }

        // Add the other portion of the mass to the upper bin (if it exists).
        // Important: the size of the total_mass_in_bins array is n_grid_cells + 2, so the indices
        // are valid from 0 to n_grid_cells + 1. lower_bin_idx + 1 can be at most n_grid_cells.
        // The condition must ensure that lower_bin_idx + 1 does not exceed the end of the array.
        if ((lower_bin_idx + 1) < (n_grid_cells + 2)) { // Corrected upper bound for array access
            #pragma omp atomic
            total_mass_in_bins[lower_bin_idx + 1] += current_mass * fraction_in_lower_bin;
        }
    }

    // 3. step: Calculate density from the collected masses and fill the output array
    // Iterate through each grid cell. Can be parallelized.
    #pragma omp parallel for
    for (int j = 0; j < n_grid_cells; j++) {
        // Check that bin_areas[j] is not zero to avoid floating-point errors.
        if (bin_areas[j] > 1e-18) { // Compare to a very small positive value
            output_sigma_d_grid[j] = total_mass_in_bins[j] / bin_areas[j];
        } else {
            output_sigma_d_grid[j] = 0.0; // If the area is zero, the density is also zero
        }
    }

    // Freeing memory
    free(total_mass_in_bins);
    free(bin_areas);
}

// Calculating the turbulent alpha parameter --> reducing alpha with alpha_r
double calculate_turbulent_alpha(double r, const disk_t *disk_params) {
    double alpha_r;
    alpha_r = 1.0 - 0.5 * (1.0 - disk_params->a_mod) * (tanh((r - disk_params->r_dze_i) / disk_params->Dr_dze_i) + tanh((disk_params->r_dze_o - r) / disk_params->Dr_dze_o));
    return alpha_r * disk_params->alpha_visc;
}



// Calculate Stokes Number
double Stokes_Number(double pradius_au, double sigma_msun_au2, double r_au, const disk_t *disk_params) {

    // --- Conversions to CGS unit system for calculation ---

    // 1. Particle radius: AU -> cm
    double pradius_cm = pradius_au * AU_TO_CM; // AU_TO_CM: 1.495978707e13

    // 2. Radial distance: AU -> cm
    double r_cm = r_au * AU_TO_CM;

    // 3. Gas surface density: M_Sun/AU^2 -> g/cm^2
    //    1 M_Sun = SUN_MASS_TO_GRAMS g (1.98847e33 g)
    //    1 AU^2 = (AU_TO_CM)^2 cm^2
    double sigma_g_cm2 = sigma_msun_au2 * SUN_MASS_TO_GRAMS / (AU_TO_CM * AU_TO_CM);

    // 4. Particle internal density: g/cm^3 (assuming disk_params->PDENSITY is already in this unit)
    double rho_s_g_cm3 = disk_params->PDENSITY; // Check this value! E.g., 3.0 g/cm^3

    // 5. Keplerian velocity (or angular velocity): AU/time_unit -> cm/s
    //    calculate_keplerian_angular_velocity(r_au, disk_params) returns a value in 1/time_unit (e.g., 2pi/year if G=1, M=1, R=1 AU is the period)
    //    This needs to be converted to rad/s.
    double omega_k_internal_units = calculate_keplerian_angular_velocity(r_au, disk_params);
    double omega_k_rad_s = omega_k_internal_units / INTERNAL_TIME_TO_SEC; // INTERNAL_TIME_TO_SEC: 31557600.0 / (2.0 * PI)

    // 6. Sound speed: AU/time_unit -> cm/s
    //    calculate_local_sound_speed(r_au, disk_params) return a value in AU/time_unit
    //    This needs to be converted to cm/s.
    double csound_au_per_timeunit = calculate_local_sound_speed(r_au, disk_params);
    double csound_cm_s = csound_au_per_timeunit * AU_TO_CM / INTERNAL_TIME_TO_SEC;

    // --- Check before division by zero ---
    if (sigma_g_cm2 <= 0.0 || csound_cm_s <= 0.0 || omega_k_rad_s <= 0.0) {
        fprintf(stderr, "ERROR [Stokes_Number]: Invalid CGS parameters for calculation (sigma_g_cm2: %.2e, csound_cm_s: %.2e, omega_k_rad_s: %.2e, r_au: %.2e).\n",
                 sigma_g_cm2, csound_cm_s, omega_k_rad_s, r_au);
        return 0.0;
    }

    // --- The Stokes Number Formula (Epstein drag) in CGS units ---
    // St = (rho_s * a / Sigma_g) * (Omega_K * r / c_s)
    double Stokes = (rho_s_g_cm3 * pradius_cm / sigma_g_cm2) * (omega_k_rad_s * r_cm / csound_cm_s);

    return Stokes;
}


// Particle growth (Birnstiel et al. 2012)
// Determining the initial size of the representative particle
// 1. Maximum size determined by radial drift --> output in cm!
double a_drift(double sigmad, double r, double p, double dp, double rho_p, const disk_t *disk_params) {

    double Sigmad_cgs = sigmad / GAS_SD_CONV_RATE;

    double vkep = calculate_keplerian_velocity(r,disk_params);
    double vkep2 = vkep * vkep;
    double c_s = calculate_local_sound_speed(r,disk_params);
    double c_s2 = c_s * c_s;
    double dlnPdlnr = r / p * dp;
    double s_drift =  disk_params->fDrift * 2.0 / M_PI * Sigmad_cgs / rho_p * vkep2 / c_s2 * fabs(1.0 / dlnPdlnr);
    return s_drift;
}

// 2. Maximum size according to fragmentation caused by small-scale turbulence --> output in cm!
double a_turb(double sigma, double r, double rho_p, const disk_t *disk_params) {

    double s_frag, u_frag, u_frag2, Sigma_cgs, c_s, c_s2;

    u_frag = disk_params->uFrag * CM_PER_SEC_TO_AU_PER_YEAR_OVER_2PI; // cm/sec --> AU / (yr/2pi) 
    u_frag2 = u_frag * u_frag;
    Sigma_cgs = sigma / GAS_SD_CONV_RATE;
    c_s = calculate_local_sound_speed(r,disk_params);
    c_s2 = c_s * c_s;

    s_frag = disk_params->fFrag * 2.0 / (3.0 * M_PI) * Sigma_cgs / (rho_p * calculate_turbulent_alpha(r,disk_params)) * u_frag2 / c_s2;

    return s_frag;
}

// 3. Maximum size according to fragmentation caused by radial drift --> output in cm!
double a_df(double sigma, double r, double p, double dp, double rho_p, const disk_t *disk_params) {

    double u_frag, vkep, dlnPdlnr, c_s, c_s2, s_df, Sigma_cgs;

    u_frag = disk_params->uFrag * CM_PER_SEC_TO_AU_PER_YEAR_OVER_2PI; // cm/sec --> AU / (yr/2pi)
    Sigma_cgs = sigma / GAS_SD_CONV_RATE;
    c_s = calculate_local_sound_speed(r,disk_params);
    c_s2 = c_s * c_s;
    dlnPdlnr = r / p * dp;
    vkep = calculate_keplerian_velocity(r,disk_params);

    s_df = u_frag * vkep / fabs(dlnPdlnr * c_s2 * 0.5) * 2.0 * Sigma_cgs / (M_PI * rho_p);

    return s_df;
}

// timescale for particle growth
double tauGr(double r, double eps,const disk_t *disk_params) {
    double omega = calculate_keplerian_angular_velocity(r,disk_params);
    double taugr = eps / omega;
    return taugr;
}

// calculates the particle size at a given location --> Birnstiel et al 2012
double getSize(double prad, double pdens, double sigma, double sigmad, double y, double p, double dpress_val, double dt, const disk_t *disk_params) {

    double sturb = a_turb(sigma, y, pdens, disk_params);      // in cm
    double sdf = a_df(sigma, y, p, dpress_val, pdens,disk_params); // in cm
    double srdf = a_drift(sigmad, y, p, dpress_val, pdens, disk_params); // in cm
    double smin = find_min(sturb, sdf, srdf);         // in cm -- gives the smaller size from the two particle limits above (this is the upper limit for particle growth)
    double eps = sigmad / sigma; 
    double tau_gr = tauGr(y, eps, disk_params);
    double rt = 0.0;

    smin = smin / AU_TO_CM; // in AU

    // calculates whether the above smin or the size resulting from the growth timescale limits the particle size
    if (prad < smin) {
        rt = find_min(prad * exp(dt / tau_gr), smin, HUGE_VAL);
    } else { // prad >= smin
        rt = smin;
    }

    rt = rt;
    return rt;
}

// Modified implementation of the Get_Sigmad function
void Get_Sigmad(const ParticleData_t *p_data, disk_t *disk_params, const simulation_options_t *sim_opts) {
    
    double *sigma_d_temp = (double *)calloc(disk_params->NGRID, sizeof(double));
    double *sigma_dm_temp = (double *)calloc(disk_params->NGRID, sizeof(double));

    if (sigma_d_temp == NULL || sigma_dm_temp == NULL) {
        fprintf(stderr, "ERROR: Memory allocation failed for sigma_d_temp or sigma_dm_temp in Get_Sigmad.\n");
        free(sigma_d_temp);
        free(sigma_dm_temp);
        exit(EXIT_FAILURE);
    }

    if (p_data->num_particles_pop1 > 0) {
        double (*rad_pop1_2d)[2] = (double (*)[2])calloc(p_data->num_particles_pop1, sizeof(double[2]));
        double *mass_pop1_1d = (double *)calloc(p_data->num_particles_pop1, sizeof(double));

        if (rad_pop1_2d == NULL || mass_pop1_1d == NULL) {
            fprintf(stderr, "ERROR: Memory allocation failed for Pop1 temp arrays in Get_Sigmad.\n");
            free(sigma_d_temp); free(sigma_dm_temp);
            free(rad_pop1_2d); free(mass_pop1_1d);
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < p_data->num_particles_pop1; i++) {
            rad_pop1_2d[i][0] = p_data->particles_pop1[i].distance_au;
            rad_pop1_2d[i][1] = p_data->particles_pop1[i].current_size_au;
            mass_pop1_1d[i] = p_data->particles_pop1[i].initial_mass_msun;
        }

        calculate_dust_surface_density_profile(sigma_d_temp, disk_params->rvec,
                                               rad_pop1_2d, mass_pop1_1d,
                                               p_data->num_particles_pop1, disk_params->NGRID, disk_params);

        free(rad_pop1_2d);
        free(mass_pop1_1d);
    }

    if (sim_opts->twopop == 1.0 && p_data->num_particles_pop2 > 0) {
        double (*rad_pop2_2d)[2] = (double (*)[2])calloc(p_data->num_particles_pop2, sizeof(double[2]));
        double *mass_pop2_1d = (double *)calloc(p_data->num_particles_pop2, sizeof(double));

        if (rad_pop2_2d == NULL || mass_pop2_1d == NULL) {
            fprintf(stderr, "ERROR: Memory allocation failed for Pop2 temp arrays in Get_Sigmad.\n");
            free(sigma_d_temp); free(sigma_dm_temp);
            free(rad_pop2_2d); free(mass_pop2_1d);
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < p_data->num_particles_pop2; i++) {
            rad_pop2_2d[i][0] = p_data->particles_pop2[i].distance_au;
            rad_pop2_2d[i][1] = p_data->particles_pop2[i].current_size_au;
            mass_pop2_1d[i] = p_data->particles_pop2[i].initial_mass_msun;
        }

        calculate_dust_surface_density_profile(sigma_dm_temp, disk_params->rvec,
                                               rad_pop2_2d, mass_pop2_1d,
                                               p_data->num_particles_pop2, disk_params->NGRID, disk_params);

        free(rad_pop2_2d);
        free(mass_pop2_1d);
    }

    for (int i = 0; i < disk_params->NGRID; i++) {
        disk_params->sigmadustvec[i] = sigma_d_temp[i];
        // disk_params->sigmadustmicrvec[i] = sigma_dm_temp[i]; // If such an array exists
    }

    free(sigma_d_temp);
    free(sigma_dm_temp);
    
}


void Get_Radius(dust_particle_t *particles_array, int num_particles, double deltat, double t,
                 const simulation_options_t *sim_opts, const disk_t *disk_params) {

    if (disk_params == NULL) {
        fprintf(stderr, "ERROR [Get_Radius]: disk_params pointer IS NULL on entry!\n");
        exit(EXIT_FAILURE); // Exit if disk_params is NULL
    }

    if (particles_array == NULL || num_particles <= 0) {
        return;
    }

    const double *aggregated_sigmad_ptr = disk_params->sigmadustvec;
    const double *aggregated_rdvec_ptr = disk_params->rvec;
    int num_grid_points = disk_params->NGRID;
    double grid_dd = disk_params->DD;
    double grid_rmin = disk_params->RMIN;
    double grid_rmax = disk_params->RMAX;

    #pragma omp parallel for
    for (int i = 0; i < num_particles; ++i) {

        // --- 1. Early check: if the particle size is 0.0.
        // This "kill" flag was set in the previous step if it went outside RMIN/RMAX.
        if (particles_array[i].current_size_au <= 0.0) {
            // Ensure its distance and velocity are also 0.0.
            particles_array[i].distance_au = 0.0;
            particles_array[i].size_reciprocal = 0.0;
            particles_array[i].drdt = 0.0; 
            continue; // Jump to the next particle, do not process it further.
        }

        double current_r = particles_array[i].distance_au;
        double current_particle_size = particles_array[i].current_size_au;

        double new_r = 0.0; // int_step will update this
        double new_size = current_particle_size; // int_step will update this (if growth=1)

        // Calling int_step with current data and freshly calculated sigmad values
        int_step(t,
                  current_particle_size,
                  current_r,
                  deltat,
                  &new_r, // pointer to new radial position
                  &new_size, // pointer to new particle size
                  aggregated_sigmad_ptr,
                  aggregated_rdvec_ptr,
                  num_grid_points,
                  grid_dd,
                  grid_rmin,
                  disk_params,
                  sim_opts);

        // Update particle's position and size
        particles_array[i].distance_au = new_r;

        // Only update size if growth is enabled
        if (sim_opts->growth == 1.0) {
            particles_array[i].current_size_au = new_size;
            // Update reciprocal value too
            if (particles_array[i].current_size_au > 0.0) {
                particles_array[i].size_reciprocal = 1.0 / particles_array[i].current_size_au;
            } else {
                particles_array[i].size_reciprocal = 0.0;
            }
        }
        
        // If growth is not enabled, particle size remains unchanged.
        
        // --- 2. Clamping AND "kill" process AFTER int_step RUNS ---
        // If the new position falls below RMIN or above RMAX, "kill" the particle.
        if (particles_array[i].distance_au <= grid_rmin || particles_array[i].distance_au >= grid_rmax) {
//            fprintf(stderr, "DEBUG [Get_Radius]: Particle %d (R=%.4e AU) went out of bounds (RMIN=%.4e, RMAX=%.4e). Killing it.\n",
//                            i, particles_array[i].distance_au, grid_rmin, grid_rmax);
            particles_array[i].distance_au = 0.0;
            particles_array[i].current_size_au = 0.0; // This is the main "kill" flag for the next step
            particles_array[i].size_reciprocal = 0.0;
            particles_array[i].drdt = 0.0; // Important: reset the velocity as well
        }
    }
}