#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <unistd.h> // For access()

#include "config.h"
#include "io_utils.h"
#include "disk_model.h"
#include "dust_physics.h"
#include "utils.h"
#include "simulation_core.h"
#include "particle_data.h"
#include "globals.h"

// A pici érték (epsilon) definíciója, a nullával való osztás
// elkerülésére. Ezt a fájl elejére kell beilleszteni.
#define EPSILON_SIZE 1e-30

// NOTE: PARTICLE_NUMBER is assumed to be defined as a global variable (e.g., in config.h)
// and represents the number of particles in a *single* population.
// If twopop is enabled, both populations will have PARTICLE_NUMBER particles.
extern int PARTICLE_NUMBER;

/* Kiszámolja az 1D-s driftet */
void eqrhs(double prad, double dp, double sigma, double ug, double r, double *drdt, const disk_t *disk_params) {
    double P, H, dPdr, St, csound;

    *drdt = 0.0;

    if (isnan(prad) || isinf(prad) || prad < EPSILON_SIZE) {
        fprintf(stderr, "ERROR [eqrhs]: Invalid prad (%.10lg) for r=%.10lg. Setting drdt to 0.\n", prad, r);
        return;
    }
    if (isnan(sigma) || isinf(sigma) || sigma < EPSILON_SIZE) {
        fprintf(stderr, "ERROR [eqrhs]: Invalid sigma (%.10lg) for r=%.10lg. Setting drdt to 0.\n", sigma, r);
        return;
    }
    if (isnan(r) || isinf(r) || r < EPSILON_SIZE || r > disk_params->RMAX) {
        if (r > disk_params->RMAX) {
            fprintf(stderr, "INFO [eqrhs]: Particle is outside RMAX (r=%.10lg, RMAX=%.10lg). Setting drdt to 0.\n", r, disk_params->RMAX);
        } else {
            fprintf(stderr, "ERROR [eqrhs]: Invalid r (%.10lg). Setting drdt to 0.\n", r);
        }
        return;
    }
    if (isnan(ug) || isinf(ug)) {
        fprintf(stderr, "ERROR [eqrhs]: Invalid ug (%.10lg) for r=%.10lg. Setting drdt to 0.\n", ug, r);
        return;
    }
    
    if (isnan(dp) || isinf(dp)) {
        fprintf(stderr, "ERROR [eqrhs]: Invalid dp (%.10lg) for r=%.10lg. Setting drdt to 0.\n", dp, r);
        return;
    }

    St = calculate_stokes_number(prad, sigma, r, disk_params);
    H = calculate_scale_height(r, disk_params);
    P = calculate_gas_pressure(sigma, r, disk_params);
    dPdr = dp;
    csound = calculate_local_sound_speed(r, disk_params);

    if (isnan(St) || isinf(St) || isnan(H) || isinf(H) || isnan(P) || isinf(P) || isnan(csound) || isinf(csound)) {
        fprintf(stderr, "ERROR [eqrhs]: One of the intermediate values is invalid. St: %.10lg, H: %.10lg, P: %.10lg, csound: %.10lg. Setting drdt to 0.\n", St, H, P, csound);
        return;
    }

    double denominator = 1.0 + St * St;
    if (fabs(denominator) < EPSILON_SIZE) {
        fprintf(stderr, "ERROR [eqrhs]: Denominator (1+St*St) is too close to zero (%.10lg) for r=%.10lg. Setting drdt to 0.\n", denominator, r);
        return;
    }
    
    double term1 = ug / denominator;
    double term2 = St / denominator * H / P * dPdr * csound;
    
    *drdt = term1 + term2;

    if (isnan(*drdt) || isinf(*drdt)) {
        fprintf(stderr, "CRITICAL ERROR [eqrhs]: Final drdt BECAME NaN/INF (%.10lg) for particle at r=%.10lg!\n", *drdt, r);
        fprintf(stderr, "DEBUG [eqrhs]: Components: ug=%.10lg, St=%.10lg, H=%.10lg, P=%.10lg, dPdr=%.10lg, csound=%.10lg\n", ug, St, H, P, dPdr, csound);
        *drdt = 0.0;
    }
}

/* for solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*hu/(2*r)*dsigma/dr   --> 3*nu = Coeff_1  */
double Coeff_1(double r, const disk_t *disk_params){              
    double A;
    A = 3.0 * calculate_gas_viscosity(r, disk_params);
//        fprintf(stderr, "DEBUG: Coeff_1 result: %e\n", A);

    return A;
}

/* for solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*hu/(2*r)*dsigma/dr   --> 9*nu /(2*r) = Coeff_2    */
double Coeff_2(double r, const disk_t *disk_params){                    
    double B;
    B = 9.0 * calculate_gas_viscosity(r, disk_params) / (2.0 * r);
//    fprintf(stderr, "DEBUG: Coeff_2 result: %e\n", B);

    return B;
}

double time_step(const disk_t *disk_params) {
    double A_max, stepping;
    int i;

    A_max = -10000.0;
    
    for(i = 0; i < disk_params->NGRID; i++) {
        if(Coeff_1(disk_params->rvec[i], disk_params) > A_max) {
            A_max = Coeff_1(disk_params->rvec[i], disk_params);
        }
    }
    stepping = disk_params->DD * disk_params->DD / (2.0 * A_max);
    return stepping;
}

/* Runge-Kutta4 integrator */
void int_step(double t,
              double psize,
              double current_r,
              double step,
              double *new_prad_ptr,
              double *new_size_ptr,
              const double *aggregated_sigmad_ptr,
              const double *aggregated_rdvec_ptr,
              int num_grid_points,
              double grid_dd,
              double grid_rmin,
              const disk_t *disk_params,
              const simulation_options_t *sim_opts) {
    
    if (psize < EPSILON_SIZE) {
        psize = EPSILON_SIZE;
    }
    
    double ytemp, dpress_temp, sigma_gas_temp, ugas_temp, sigma_dust_temp;
    double dy1, dy2, dy3, dy4;

    // Runge-Kutta 4th-order method
    // 1. lépés
    interpol(disk_params->sigmavec, disk_params->rvec, current_r, &sigma_gas_temp, disk_params->DD, sim_opts, disk_params);
    interpol(disk_params->dpressvec, disk_params->rvec, current_r, &dpress_temp, disk_params->DD, sim_opts, disk_params);
    interpol(disk_params->ugvec, disk_params->rvec, current_r, &ugas_temp, disk_params->DD, sim_opts, disk_params);
    
    if (sigma_gas_temp < EPSILON_SIZE) {
        fprintf(stderr, "WARNING [int_step]: Interpolated sigma_gas is <= 0 for r=%.10lg. Skipping step.\n", current_r);
        *new_prad_ptr = current_r;
        *new_size_ptr = psize;
        return;
    }

    eqrhs(psize, dpress_temp, sigma_gas_temp, ugas_temp, current_r, &dy1, disk_params);

    // 2. lépés
    ytemp = current_r + 0.5 * step * dy1;
    if (ytemp < disk_params->RMIN) {
        ytemp = disk_params->RMIN;
    }
    if (ytemp > disk_params->RMAX) {
        ytemp = disk_params->RMAX;
    }
    interpol(disk_params->sigmavec, disk_params->rvec, ytemp, &sigma_gas_temp, disk_params->DD, sim_opts, disk_params);
    interpol(disk_params->dpressvec, disk_params->rvec, ytemp, &dpress_temp, disk_params->DD, sim_opts, disk_params);
    interpol(disk_params->ugvec, disk_params->rvec, ytemp, &ugas_temp, disk_params->DD, sim_opts, disk_params);
    
    if (sigma_gas_temp < EPSILON_SIZE) {
        fprintf(stderr, "WARNING [int_step]: Interpolated sigma_gas is <= 0 for r=%.10lg. Skipping step.\n", ytemp);
        *new_prad_ptr = current_r;
        *new_size_ptr = psize;
        return;
    }
    
    eqrhs(psize, dpress_temp, sigma_gas_temp, ugas_temp, ytemp, &dy2, disk_params);

    // 3. lépés
    ytemp = current_r + 0.5 * step * dy2;
    if (ytemp < disk_params->RMIN) {
        ytemp = disk_params->RMIN;
    }
    if (ytemp > disk_params->RMAX) {
        ytemp = disk_params->RMAX;
    }
    interpol(disk_params->sigmavec, disk_params->rvec, ytemp, &sigma_gas_temp, disk_params->DD, sim_opts, disk_params);
    interpol(disk_params->dpressvec, disk_params->rvec, ytemp, &dpress_temp, disk_params->DD, sim_opts, disk_params);
    interpol(disk_params->ugvec, disk_params->rvec, ytemp, &ugas_temp, disk_params->DD, sim_opts, disk_params);
    
    if (sigma_gas_temp < EPSILON_SIZE) {
        fprintf(stderr, "WARNING [int_step]: Interpolated sigma_gas is <= 0 for r=%.10lg. Skipping step.\n", ytemp);
        *new_prad_ptr = current_r;
        *new_size_ptr = psize;
        return;
    }
    
    eqrhs(psize, dpress_temp, sigma_gas_temp, ugas_temp, ytemp, &dy3, disk_params);
    
    // 4. lépés
    ytemp = current_r + step * dy3;
    if (ytemp < disk_params->RMIN) {
        ytemp = disk_params->RMIN;
    }
    if (ytemp > disk_params->RMAX) {
        ytemp = disk_params->RMAX;
    }
    interpol(disk_params->sigmavec, disk_params->rvec, ytemp, &sigma_gas_temp, disk_params->DD, sim_opts, disk_params);
    interpol(disk_params->dpressvec, disk_params->rvec, ytemp, &dpress_temp, disk_params->DD, sim_opts, disk_params);
    interpol(disk_params->ugvec, disk_params->rvec, ytemp, &ugas_temp, disk_params->DD, sim_opts, disk_params);

    if (sigma_gas_temp < EPSILON_SIZE) {
        fprintf(stderr, "WARNING [int_step]: Interpolated sigma_gas is <= 0 for r=%.10lg. Skipping step.\n", ytemp);
        *new_prad_ptr = current_r;
        *new_size_ptr = psize;
        return;
    }
    
    eqrhs(psize, dpress_temp, sigma_gas_temp, ugas_temp, ytemp, &dy4, disk_params);

    *new_prad_ptr = current_r + step * (dy1 + 2.0 * dy2 + 2.0 * dy3 + dy4) / 6.0;
    if (*new_prad_ptr < disk_params->RMIN) {
        *new_prad_ptr = disk_params->RMIN;
    }
    if (*new_prad_ptr > disk_params->RMAX) {
        *new_prad_ptr = disk_params->RMAX;
    }

    if (sim_opts->growth == 1.0) {
        interpol(aggregated_sigmad_ptr, aggregated_rdvec_ptr, current_r, &sigma_dust_temp, grid_dd, sim_opts, disk_params);
        *new_size_ptr = update_particle_size(psize,
                                             disk_params->PDENSITY,
                                             sigma_gas_temp,
                                             sigma_dust_temp,
                                             current_r,
                                             disk_params->dpressvec[(int)((current_r - disk_params->RMIN)/disk_params->DD)],
                                             dpress_temp,
                                             step,
                                             disk_params);
    }
}

// Helper function to get max particle distance (distance_au)
double get_max_particle_distance(const dust_particle_t *particles_array, int num_particles) {
    if (particles_array == NULL || num_particles <= 0) return 0.0;
    double max_dist = 0.0;
    for (int k = 0; k < num_particles; ++k) {
        if (particles_array[k].distance_au > max_dist) {
            max_dist = particles_array[k].distance_au;
        }
    }
    return max_dist;
}

// Helper function to get min particle distance (distance_au)
double get_min_particle_distance(const dust_particle_t *particles_array, int num_particles, double min_dist_cutoff) {
    if (particles_array == NULL || num_particles <= 0) return HUGE_VAL; 
    double max_reciprocal = 0.0; 
    
    for (int k = 0; k < num_particles; ++k) {
        if (particles_array[k].distance_au > min_dist_cutoff && particles_array[k].distance_au > 0.0) {
            double current_reciprocal = 1.0 / particles_array[k].distance_au;
            if (current_reciprocal > max_reciprocal) {
                max_reciprocal = current_reciprocal;
            }
        }
    }
    return (max_reciprocal > 0.0) ? (1.0 / max_reciprocal) : 0.0; 
}

// Function to find minimum of three doubles (from utils.h or similar)
double find_min_three(double val1, double val2, double val3) {
    double min_val = val1;
    if (val2 < min_val) min_val = val2;
    if (val3 < min_val) min_val = val3;
    return min_val;
}

void tIntegrate(disk_t *disk_params, const simulation_options_t *sim_opts, output_files_t *output_files) {
    const double CFL_FACTOR_DRIFT = 0.05;
    const double CFL_FACTOR_KEPLER = 0.005;
    const double MIN_TIMESTEP_YEARS = 1.0e-10;
    const double MAX_TIMESTEP_ABSOLUTE_YEARS = 1.0;

    ParticleData_t p_data;
    HeaderData_t header_data_for_files;

    double initial_max_drdt = 0.0;
    double particle_dt_suggestion_in_years;

    p_data.particles_pop1 = NULL;
    p_data.particles_pop2 = NULL;
    p_data.num_particles_pop1 = PARTICLE_NUMBER;
    p_data.num_particles_pop2 = PARTICLE_NUMBER;

    double L = 0.;
    if (disk_params == NULL) {
        fprintf(stderr, "ERROR [tIntegrate]: disk_params pointer is NULL!\n");
        exit(EXIT_FAILURE);
    }
    if (sim_opts->drift == 1.) {
        PARTICLE_NUMBER = get_particle_count(sim_opts->dust_input_filename);
    } else {
        PARTICLE_NUMBER = 0;
    }
    if (PARTICLE_NUMBER > 0) {
        int num_pop2 = (sim_opts->twopop == 1.0) ? PARTICLE_NUMBER : 0;
        allocate_particle_data(&p_data, PARTICLE_NUMBER, num_pop2, (int)sim_opts->twopop);
    }
    if (sim_opts->drift == 1.) {
        if (initialize_mass_accumulation_file(output_files, sim_opts, disk_params, &header_data_for_files) != 0) {
            fprintf(stderr, "ERROR: Failed to set up initial output files. Exiting.\n");
            free_particle_data(&p_data);
            exit(EXIT_FAILURE);
        }
        load_dust_particles(&p_data, sim_opts->dust_input_filename, disk_params, sim_opts);
    }

    double max_dist_pop1 = 0.0;
    double min_dist_pop1 = 0.0;
    double max_dist_pop2 = 0.0;
    double min_dist_pop2 = 0.0;

    char dens_name[MAX_PATH_LEN] = "";
    char dust_name[MAX_PATH_LEN] = "";
    char micron_dust_name[MAX_PATH_LEN] = "";

    double t = 0.0;
    double t_integration_in_internal_units = sim_opts->TMAX * (2.0 * M_PI);
    const double delta_r_uniform_grid = (disk_params->RMAX - disk_params->RMIN) / (double)(disk_params->NGRID - 1);
    double gas_dt_suggestion_in_years = time_step(disk_params);

    if (sim_opts->drift == 1.0 && p_data.num_particles_pop1 > 0) {
        for (int i = 0; i < p_data.num_particles_pop1; ++i) {
            // A drdt értékét használjuk, feltételezve, hogy az már be van állítva
            if (p_data.particles_pop1[i].drdt != HUGE_VAL) {
                if (fabs(p_data.particles_pop1[i].drdt) > initial_max_drdt) {
                    initial_max_drdt = fabs(p_data.particles_pop1[i].drdt);
                }
            }
        }
    }

    fprintf(stderr,"init max drdt %lg\n", initial_max_drdt);
    
    if (initial_max_drdt > 1e-15) {
        particle_dt_suggestion_in_years = CFL_FACTOR_DRIFT * delta_r_uniform_grid / initial_max_drdt;
    } else {
        particle_dt_suggestion_in_years = MAX_TIMESTEP_ABSOLUTE_YEARS;
    }

    double deltat_in_years = fmin(gas_dt_suggestion_in_years, particle_dt_suggestion_in_years);
    if (sim_opts->DT > 0.0) {
        deltat_in_years = fmin(deltat_in_years, sim_opts->DT);
    }
    deltat_in_years = fmax(deltat_in_years, MIN_TIMESTEP_YEARS);
    deltat_in_years = fmin(deltat_in_years, MAX_TIMESTEP_ABSOLUTE_YEARS);
    
    double output_interval_years = sim_opts->TMAX / sim_opts->WO;
    deltat_in_years = fmin(deltat_in_years, output_interval_years);
    double deltat = deltat_in_years * (2.0 * M_PI);
    
    double overall_min_dist = HUGE_VAL;
    double overall_max_dist = 0.0;

    do {
        if (sim_opts->drift == 1.) {
            for (int i = 0; i < p_data.num_particles_pop1; i++) {
                if (p_data.particles_pop1[i].distance_au > 0.0) {
                    p_data.particles_pop1[i].distance_au_reciprocal = 1.0 / p_data.particles_pop1[i].distance_au;
                } else {
                    p_data.particles_pop1[i].distance_au_reciprocal = 0.0;
                }
            }
            max_dist_pop1 = get_max_particle_distance(p_data.particles_pop1, p_data.num_particles_pop1);
            min_dist_pop1 = get_min_particle_distance(p_data.particles_pop1, p_data.num_particles_pop1, disk_params->RMIN);

            if (sim_opts->twopop == 1) {
                for (int i = 0; i < p_data.num_particles_pop2; i++) {
                    if (p_data.particles_pop2[i].distance_au > 0.0) {
                        p_data.particles_pop2[i].distance_au_reciprocal = 1.0 / p_data.particles_pop2[i].distance_au;
                    } else {
                        p_data.particles_pop2[i].distance_au_reciprocal = 0.0;
                    }
                }
                max_dist_pop2 = get_max_particle_distance(p_data.particles_pop2, p_data.num_particles_pop2);
                min_dist_pop2 = get_min_particle_distance(p_data.particles_pop2, p_data.num_particles_pop2, disk_params->RMIN);
                overall_min_dist = fmin(min_dist_pop1, min_dist_pop2);
                overall_max_dist = fmax(max_dist_pop1, max_dist_pop2);
            } else {
                overall_min_dist = min_dist_pop1;
                overall_max_dist = max_dist_pop1;
            }

            double current_time_years = t / (2.0 * M_PI);
            
            if ((fmod(current_time_years, output_interval_years) < deltat_in_years || current_time_years == 0) &&
                (current_time_years >= L - (deltat_in_years * 0.5))) {
                fprintf(stderr,"\n--- Simulation Time: %.2e years (Internal time: %.2e, L: %.2e) ---\n", current_time_years, t, L);

                if (sim_opts->evol == 1 || current_time_years == 0) {
                    snprintf(dens_name, MAX_PATH_LEN, "%s/%s/%s_%08d.dat", sim_opts->output_dir_name, LOGS_DIR, FILE_DENS_PREFIX, (int)L);
                }
                snprintf(dust_name, MAX_PATH_LEN, "%s/%s/%s_%08d.dat", sim_opts->output_dir_name, LOGS_DIR, FILE_DUST_PREFIX,(int)L);
                snprintf(micron_dust_name, MAX_PATH_LEN, "%s/%s/micron_%s_%08d.dat", sim_opts->output_dir_name, LOGS_DIR, FILE_DUST_PREFIX, (int)L);

                if (sim_opts->evol == 1 || current_time_years == 0) {
                    output_files->surface_file = fopen(dens_name, "w");
                    if (output_files->surface_file != NULL) {
                        HeaderData_t gas_header_data = {.current_time = current_time_years, .is_initial_data = (current_time_years == 0.0)};
                        write_file_header(output_files->surface_file, FILE_TYPE_GAS_DENSITY, &gas_header_data);
                    }
                }

                output_files->dust_file = fopen(dust_name, "w");
                if (output_files->dust_file != NULL) {
                    HeaderData_t dust_header_data = {.current_time = current_time_years, .is_initial_data = (current_time_years == 0.0)};
                    write_file_header(output_files->dust_file, FILE_TYPE_DUST_EVOL, &dust_header_data);
                }

                if (sim_opts->twopop == 1.) {
                    output_files->micron_dust_file = fopen(micron_dust_name, "w");
                    if (output_files->micron_dust_file != NULL) {
                        HeaderData_t micron_dust_header_data = {.current_time = current_time_years, .is_initial_data = (current_time_years == 0.0)};
                        write_file_header(output_files->micron_dust_file, FILE_TYPE_MICRON_DUST_EVOL, &micron_dust_header_data);
                    }
                }

                if (sim_opts->growth == 1.) {
                    calculate_dust_density_grid(&p_data, disk_params, sim_opts);
                }

                if (sim_opts->evol == 1 || current_time_years == 0) {
                    write_gas_profile_to_file(disk_params, output_files);
                }
                if (sim_opts->growth == 1.) {
                    write_dust_profile_to_file((int)L, &p_data, disk_params, sim_opts, output_files);
                }
                L = L + output_interval_years;
                close_snapshot_files(output_files, dens_name, dust_name, micron_dust_name, sim_opts);
            }

            if (sim_opts->evol == 1.) {
                get_gas_surface_density_pressure_pressure_gradient(sim_opts, disk_params, deltat);
                validate_disk_state(disk_params);  // ← IDE!

            }

            if (p_data.particles_pop1 == NULL) {
                fprintf(stderr, "ERROR [tIntegrate]: particles_pop1 is NULL before update_particle_positions call!\n");
                exit(EXIT_FAILURE);
            }
            if (disk_params == NULL) {
                fprintf(stderr, "ERROR [tIntegrate]: disk_params is NULL before update_particle_positions call!\n");
                exit(EXIT_FAILURE);
            }
            update_particle_positions(p_data.particles_pop1, PARTICLE_NUMBER, deltat, t, sim_opts, disk_params);
            if (sim_opts->twopop == 1.) {
                update_particle_positions(p_data.particles_pop2, PARTICLE_NUMBER, deltat, t, sim_opts, disk_params);
            }
            double particle_dt_drift_limit_in_years = HUGE_VAL;
            if (sim_opts->drift == 1.0 && PARTICLE_NUMBER > 0) {
                double max_abs_drdt_from_current_step = 0.0;
                for (int i = 0; i < p_data.num_particles_pop1; ++i) {
                    if (p_data.particles_pop1[i].drdt != HUGE_VAL) {
                        if (fabs(p_data.particles_pop1[i].drdt) > max_abs_drdt_from_current_step) {
                            max_abs_drdt_from_current_step = fabs(p_data.particles_pop1[i].drdt);
                        }
                    }
                }
                if (sim_opts->twopop == 1) {
                    for (int i = 0; i < p_data.num_particles_pop2; ++i) {
                        if (p_data.particles_pop2[i].drdt != HUGE_VAL) {
                            if (fabs(p_data.particles_pop2[i].drdt) > max_abs_drdt_from_current_step) {
                                max_abs_drdt_from_current_step = fabs(p_data.particles_pop2[i].drdt);
                            }
                        }
                    }
                }
                if (max_abs_drdt_from_current_step > 1e-15) {
                    particle_dt_drift_limit_in_years = CFL_FACTOR_DRIFT * delta_r_uniform_grid / max_abs_drdt_from_current_step;
                } else {
                    particle_dt_drift_limit_in_years = MAX_TIMESTEP_ABSOLUTE_YEARS;
                }
            } else {
                particle_dt_drift_limit_in_years = MAX_TIMESTEP_ABSOLUTE_YEARS;
            }

            deltat_in_years = fmin(gas_dt_suggestion_in_years, particle_dt_drift_limit_in_years);
            if (sim_opts->DT > 0.0) {
                deltat_in_years = fmin(deltat_in_years, sim_opts->DT);
            }
            deltat_in_years = fmax(deltat_in_years, MIN_TIMESTEP_YEARS);
            deltat_in_years = fmin(deltat_in_years, MAX_TIMESTEP_ABSOLUTE_YEARS);
            deltat_in_years = fmin(deltat_in_years, output_interval_years);
            deltat = deltat_in_years * (2.0 * M_PI);
            t = t + deltat;

            double current_time_years_after_step = t / (2.0 * M_PI);
            if (!(overall_max_dist >= disk_params->RMIN && overall_min_dist != overall_max_dist)) {
                fprintf(stderr,"DEBUG [tIntegrate]: Simulation termination condition met (overall_max_dist < RMIN or overall_min_dist == overall_max_dist) at time: %.2e years.\n", current_time_years_after_step);
                goto cleanup;
            }
        } else {
            double current_time_years = t / (2.0 * M_PI);
            double output_interval_years = sim_opts->TMAX / sim_opts->WO;
            if ((fmod(current_time_years, output_interval_years) < deltat_in_years || current_time_years == 0) && L - current_time_years < deltat_in_years * 0.5) {
                fprintf(stderr,"\n--- Simulation Time: %.2e years (Internal time: %.2e, L: %.2e) ---\n", current_time_years, t, L);
                fprintf(stderr,"DEBUG [tIntegrate]: Outputting data for gas-only simulation at time %.2e. L=%.2e\n", current_time_years, L);
                snprintf(dens_name, MAX_PATH_LEN, "%s/%s/%s_%08d.dat", sim_opts->output_dir_name, LOGS_DIR, FILE_DENS_PREFIX, (int)L);
                fprintf(stderr, "DEBUG [tIntegrate]: Outputting %s_%08d.dat to %s.\n", dens_name, (int)L, dens_name);
                output_files->surface_file = fopen(dens_name, "w");
                if (output_files->surface_file == NULL) {
                    fprintf(stderr, "ERROR: Could not open %s for writing in gas-only branch.\n", dens_name);
                } else {
                    fprintf(stderr, "DEBUG [tIntegrate]: Opened %s for writing in gas-only branch.\n", dens_name);
                    HeaderData_t gas_header_data = {.current_time = current_time_years, .is_initial_data = (current_time_years == 0.0)};
                    write_file_header(output_files->surface_file, FILE_TYPE_GAS_DENSITY, &gas_header_data);
                }
                write_gas_profile_to_file(disk_params, output_files);
                if (output_files->surface_file != NULL) {
                    fclose(output_files->surface_file);
                    output_files->surface_file = NULL;
                    fprintf(stderr, "DEBUG [tIntegrate]: Closed %s in gas-only branch.\n", dens_name);
                }
                L = L + output_interval_years;
                fprintf(stderr,"DEBUG [tIntegrate]: Updated L to %.2e.\n", L);
            }
            fprintf(stderr,"DEBUG [tIntegrate]: Calling get_gas_surface_density_pressure_pressure_gradient for gas-only evolution.\n");
            get_gas_surface_density_pressure_pressure_gradient(sim_opts, disk_params,deltat);
            validate_disk_state(disk_params);  // ← IDE!

            t = t + deltat;
        }
    } while (t <= t_integration_in_internal_units);

    fprintf(stderr,"\n\nDEBUG [tIntegrate]: Main simulation loop finished (t > t_integration).\n");

cleanup:
    cleanup_simulation_resources(&p_data, output_files, sim_opts);
    fprintf(stderr,"DEBUG [tIntegrate]: Cleanup completed.\n");
}


