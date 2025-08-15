// src/simulation_core.c

// Standard C Library Includes
#include <stdio.h>    // For printf, fopen, fclose, fscanf, snprintf, sprintf
#include <stdlib.h>   // For exit, EXIT_FAILURE, EXIT_SUCCESS, system
#include <math.h>     // For M_PI, fmod, HUGE_VAL
#include <string.h>   // For snprintf, sprintf

#include <omp.h>

// Your Project Header Includes
#include "config.h"        // For PARTICLE_NUMBER, RMIN (and potentially other constants)
#include "io_utils.h"      // For load_dust_particles, initialize_mass_accumulation_file, write_file_header, close_snapshot_files
#include "disk_model.h"    // For calculate_scale_height, calculate_gas_pressure, calculate_local_sound_speed, calculate_gas_viscosity, interpol
#include "dust_physics.h"  // For calculate_stokes_number, Count_Mass, secondaryGrowth, calculate_dust_density_grid, update_particle_positions, update_particle_size (if separate), get_particle_count (function prototype)
#include "utils.h"         // For time_step, get_gas_surface_density_pressure_pressure_gradient, find_min_three
#include "simulation_core.h" // For tIntegrate prototype, etc.
#include "particle_data.h" // New include for ParticleData_t
#include "globals.h"       // For MAX_PATH_LEN, FILE_DENS_PREFIX, FILE_DUST_PREFIX, LOGS_DIR, M_PI, ROUND_PRECISION_FACTOR

// NOTE: PARTICLE_NUMBER is assumed to be defined as a global variable (e.g., in config.h)
// and represents the number of particles in a *single* population.
// If twopop is enabled, both populations will have PARTICLE_NUMBER particles.
extern int PARTICLE_NUMBER; 

/*  Kiszámolja az 1D-s driftet  */
/*      dr/dt = St/(1+St*St)*H(r)/r*dlnP/dlnr*cs = St/(1+St*St) * (H/r) * (r/P) * (dP/dr) * cs      */
void eqrhs(double prad, double dp, double sigma, double ug, double r, double *drdt, const disk_t *disk_params) {
    double P, H, dPdr, St, csound;

//    fprintf(stderr, "DEBUG [eqrhs_input]: pradius=%.10lg, dp=%.10lg, sigma=%.10lg, ug=%.10lg, r=%.10lg\n", prad, dp, sigma, ug, r);


//    fprintf(stderr, "DEBUG [eqrhs]: Entering with prad=%.10lg, dp=%.10lg, sigma=%.10lg, ug=%.10lg, r=%.10lg\n", prad, dp, sigma, ug, r);

    // Ellenőrizd a bemeneti paramétereket
    if (isnan(prad) || isinf(prad) || prad <= 0.0) {
//        fprintf(stderr, "ERROR [eqrhs]: Invalid prad (%.10lg) for r=%.10lg. Exiting.\n", prad, r);
        exit(EXIT_FAILURE);
    }
    if (isnan(sigma) || isinf(sigma) || sigma <= 0.0) {
//        fprintf(stderr, "ERROR [eqrhs]: Invalid sigma (%.10lg) for r=%.10lg. Exiting.\n", sigma, r);
        exit(EXIT_FAILURE);
    }
    if (isnan(r) || isinf(r) || r <= 0.0) {
//        fprintf(stderr, "ERROR [eqrhs]: Invalid r (%.10lg). Exiting.\n", r);
        exit(EXIT_FAILURE);
    }

    St = calculate_stokes_number(prad, sigma, r,disk_params);
    H = calculate_scale_height(r, disk_params);
    P = calculate_gas_pressure(sigma, r, disk_params);
    dPdr = dp; // dPdr directly comes from dp input
    csound = calculate_local_sound_speed(r, disk_params);



//    fprintf(stderr, "DEBUG [eqrhs]: Calculated intermediates: St=%.10lg, H=%.10lg, P=%.10lg, dPdr=%.10lg, csound=%.10lg  size: %lg\n", St, H, P, dPdr, csound,prad);

    // Ellenőrizd a köztes értékeket
    if (isnan(St) || isinf(St)) {
        fprintf(stderr, "ERROR [eqrhs]: calculate_stokes_number returned NaN/Inf (%.10lg) for prad=%.10lg, sigma=%.10lg, r=%.10lg. Exiting.\n", St, prad, sigma, r);
        exit(EXIT_FAILURE);
    }
    if (isnan(H) || isinf(H) || H <= 0.0) {
        fprintf(stderr, "ERROR [eqrhs]: Scale Height returned NaN/Inf/Zero (%.10lg) for r=%.10lg. Exiting.\n", H, r);
        exit(EXIT_FAILURE);
    }
    if (isnan(P) || isinf(P) || P <= 0.0) { // P could be very small or zero at the edge of disk
        fprintf(stderr, "ERROR [eqrhs]: Gas Pressure returned NaN/Inf/Zero (%.10lg) for sigma=%.10lg, r=%.10lg. Exiting.\n", P, sigma, r);
        exit(EXIT_FAILURE);
    }
    if (isnan(csound) || isinf(csound) || csound <= 0.0) {
        fprintf(stderr, "ERROR [eqrhs]: Sound Speed returned NaN/Inf/Zero (%.10lg) for r=%.10lg. Exiting.\n", csound, r);
        exit(EXIT_FAILURE);
    }

    // Végleges számítás
    double term1 = ug / (1. + St * St);
    double term2_denominator = (1. + St * St);
    double term2 = St / term2_denominator * H / P * dPdr * csound;
    *drdt = term1 + term2;

    // Ellenőrizd az eredményt
    if (isnan(*drdt) || isinf(*drdt)) {
        fprintf(stderr, "CRITICAL ERROR [eqrhs]: Final drdt BECAME NaN/INF (%.10lg) for particle at r=%.10lg!\n", *drdt, r);
        fprintf(stderr, "DEBUG [eqrhs]: Components: ug=%.10lg, St=%.10lg, H=%.10lg, P=%.10lg, dPdr=%.10lg, csound=%.10lg\n", ug, St, H, P, dPdr, csound);
        exit(EXIT_FAILURE); // Azonnal kilép, hogy lásd, mi okozta
    }

//        fprintf(stderr, "DEBUG [eqrhs_output]: drdt=%.10lg\n", *drdt);


//    fprintf(stderr, "DEBUG [eqrhs]: Exiting with drdt=%.10lg\n", *drdt);
}


/* for solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*hu/(2*r)*dsigma/dr    --> 3*nu = Coeff_1  */
double Coeff_1(double r, const disk_t *disk_params){                    
    double A;
    A = 3.0 * calculate_gas_viscosity(r, disk_params);
    return A;
}

/* for solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*hu/(2*r)*dsigma/dr    --> 9*nu /(2*r) = Coeff_2   */
double Coeff_2(double r, const disk_t *disk_params){                            
    double B;
    B = 9.0 * calculate_gas_viscosity(r, disk_params) / (2.0 * r);
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
// prad bemenet: AU-ban!
// This function operates on a single particle's data at a time.
// `aggregated_sigmad` and `aggregated_rdvec` are expected to be arrays representing dust surface density
// on a grid, which `int_step` interpolates from.
/* Runge-Kutta4 integrator */
void int_step(double time, double psize, double current_r, // particle's current size and distance_au
              double step, double *new_prad_ptr, double *new_psize_ptr,
              const double *aggregated_sigmad, const double *aggregated_rdvec, // Aggregated dust density & radii from disk_t or similar
              int num_grid_points, double grid_dd, double grid_rmin, // Parameters for interpolation grid
              const disk_t *disk_params, const simulation_options_t *sim_opts){
    double dy1,dy2,dy3,dy4;
    double ytemp; 
    double sigma_gas, dpress, ugas;
    double gas_pressure, particle_density;

    int opt = 0; // For interpol
    double sigmadd_at_r = 0.0; // Dust surface density at current particle's radius

    // Ha a részecske a rács tartományán kívül van, térj vissza nullázott értékekkel.
    if (current_r < grid_rmin || current_r >= (grid_rmin + grid_dd * num_grid_points)) {
        fprintf(stderr, "DEBUG [int_step]: Particle at r=%.4e AU is outside grid boundaries. RMIN=%.4e AU, RMAX=%.4e AU. Aborting step.\n", 
                        current_r, grid_rmin, (grid_rmin + grid_dd * num_grid_points));
        *new_prad_ptr = current_r;
        *new_psize_ptr = 0.0;
        return;
    }

    interpol(disk_params->sigmavec, disk_params->rvec, current_r, &sigma_gas, disk_params->DD, opt, disk_params);
    interpol(disk_params->dpressvec, disk_params->rvec, current_r, &dpress, disk_params->DD, opt, disk_params);
    interpol(disk_params->ugvec, disk_params->rvec, current_r, &ugas, disk_params->DD, opt, disk_params);
    
    interpol(aggregated_sigmad, aggregated_rdvec, current_r, &sigmadd_at_r, grid_dd, opt, disk_params);

    if (sim_opts->growth == 1.) {
        if (time != 0.) {
            interpol(disk_params->pressvec, disk_params->rvec, current_r, &gas_pressure, disk_params->DD, opt, disk_params);
            particle_density = disk_params->PDENSITY;
            *new_psize_ptr = update_particle_size(psize, particle_density, sigma_gas, sigmadd_at_r, current_r, gas_pressure, dpress, step, disk_params);
        }
    }

    *new_prad_ptr = current_r;

    eqrhs(psize, dpress, sigma_gas, ugas, current_r, &dy1, disk_params);

    ytemp = current_r + 0.5 * step * dy1;
    if (ytemp < disk_params->RMIN) {
        ytemp = disk_params->RMIN;
    }
    if (ytemp > disk_params->RMAX) {
        ytemp = disk_params->RMAX;
    }

    double sigma_gas_temp, dpress_temp, ugas_temp;
    interpol(disk_params->sigmavec, disk_params->rvec, ytemp, &sigma_gas_temp, disk_params->DD, opt, disk_params);
    interpol(disk_params->dpressvec, disk_params->rvec, ytemp, &dpress_temp, disk_params->DD, opt, disk_params);
    interpol(disk_params->ugvec, disk_params->rvec, ytemp, &ugas_temp, disk_params->DD, opt, disk_params);
    eqrhs(psize, dpress_temp, sigma_gas_temp, ugas_temp, ytemp, &dy2, disk_params);
        
    ytemp = current_r + 0.5 * step * dy2;
    if (ytemp < disk_params->RMIN) {
        ytemp = disk_params->RMIN;
    }
    if (ytemp > disk_params->RMAX) {
        ytemp = disk_params->RMAX;
    }
    interpol(disk_params->sigmavec, disk_params->rvec, ytemp, &sigma_gas_temp, disk_params->DD, opt, disk_params);
    interpol(disk_params->dpressvec, disk_params->rvec, ytemp, &dpress_temp, disk_params->DD, opt, disk_params);
    interpol(disk_params->ugvec, disk_params->rvec, ytemp, &ugas_temp, disk_params->DD, opt, disk_params);
    eqrhs(psize, dpress_temp, sigma_gas_temp, ugas_temp, ytemp, &dy3, disk_params);
    
    ytemp = current_r + step * dy3;
    if (ytemp < disk_params->RMIN) {
        ytemp = disk_params->RMIN;
    }
    if (ytemp > disk_params->RMAX) {
        ytemp = disk_params->RMAX;
    }
    interpol(disk_params->sigmavec, disk_params->rvec, ytemp, &sigma_gas_temp, disk_params->DD, opt, disk_params);
    interpol(disk_params->dpressvec, disk_params->rvec, ytemp, &dpress_temp, disk_params->DD, opt, disk_params);
    interpol(disk_params->ugvec, disk_params->rvec, ytemp, &ugas_temp, disk_params->DD, opt, disk_params);
    eqrhs(psize, dpress_temp, sigma_gas_temp, ugas_temp, ytemp, &dy4, disk_params);

    *new_prad_ptr = current_r + step * (dy1 + 2.0 * dy2 + 2.0 * dy3 + dy4) / 6.0;
    if (*new_prad_ptr <= disk_params->RMIN) {
        *new_prad_ptr = disk_params->RMIN;
    }
    if (*new_prad_ptr >= disk_params->RMAX) {
        *new_prad_ptr = disk_params->RMAX;
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
// Using max reciprocal to find min value, ignoring invalid distances
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



// Amennyiben a LOCAL_AU_TO_CM és LOCAL_INTERNAL_TIME_TO_SEC nincsenek definiálva a globals.h-ban,
// ideiglenesen definiálhatod itt, a függvény elején. Ideális esetben a globals.h-ban lennének.
#ifndef LOCAL_AU_TO_CM
#define LOCAL_AU_TO_CM 1.495978707e13 // cm / AU
#endif
#ifndef LOCAL_INTERNAL_TIME_TO_SEC
#define LOCAL_INTERNAL_TIME_TO_SEC (3.1536e7 / (2.0 * M_PI)) // sec / belső időegység
#endif





// simulation_loop.c (tIntegrate function)

void tIntegrate(disk_t *disk_params, const simulation_options_t *sim_opts, output_files_t *output_files) {
    // --- Lokális konstansok a tIntegrate függvényen belül ---
    // Ezeket finomhangolhatod a szimuláció stabilitása és sebessége érdekében
    const double CFL_FACTOR_DRIFT = 0.005;      // Courant-szám a radiális mozgáshoz (javaslat: 0.05-0.2)
    const double CFL_FACTOR_KEPLER = 0.005;    // A Kepleri keringési idő töredéke (javaslat: 0.001-0.01)
    const double MIN_TIMESTEP_YEARS = 1.0e-10; // Minimum időlépés (évben), hogy ne akadjon el
    const double MAX_TIMESTEP_ABSOLUTE_YEARS = 5.0; // Abszolút maximum időlépés (évben)


    ParticleData_t p_data;
    HeaderData_t header_data_for_files;

    p_data.particles_pop1 = NULL;
    p_data.particles_pop2 = NULL;
    p_data.num_particles_pop1 = PARTICLE_NUMBER;
    p_data.num_particles_pop2 = PARTICLE_NUMBER;

    fprintf(stderr, "DEBUG [tIntegrate]: Entering tIntegrate. Current Time: %.2e, Total Time: %.2e\n", sim_opts->TCURR, sim_opts->TMAX);
    fprintf(stderr, "DEBUG [tIntegrate]: sim_opts->drift = %.2f, sim_opts->growth = %.2f, sim_opts->evol = %.2f, sim_opts->twopop = %.2f\n",
            sim_opts->drift, sim_opts->growth, sim_opts->evol, sim_opts->twopop);

    // --- NEW DEBUG PRINT: disk_params->eps at tIntegrate entry ---
    fprintf(stderr, "DEBUG [tIntegrate]: disk_params->eps value at tIntegrate entry: %.10e\n", disk_params->eps);

    double L = 0.; // "Snapshot" timer in years

    if (disk_params == NULL) {
        fprintf(stderr, "ERROR [tIntegrate]: disk_params pointer is NULL!\n");
        exit(EXIT_FAILURE);
    }

    // --- Initialization Section ---
    if (sim_opts->drift == 1.) {
        PARTICLE_NUMBER = get_particle_count(sim_opts->dust_input_filename);
        fprintf(stderr, "DEBUG [tIntegrate]: Initial particle count (PARTICLE_NUMBER) from get_particle_count: %d\n", PARTICLE_NUMBER);
    } else {
        fprintf(stderr, "DEBUG [tIntegrate]: Particle drift is OFF. PARTICLE_NUMBER set to 0.\n");
        PARTICLE_NUMBER = 0;
    }

    if (PARTICLE_NUMBER > 0) {
        int num_pop2 = (sim_opts->twopop == 1.0) ? PARTICLE_NUMBER : 0;
        allocate_particle_data(&p_data, PARTICLE_NUMBER, num_pop2, (int)sim_opts->twopop);
    } else {
        fprintf(stderr, "DEBUG [tIntegrate]: No particles to allocate. Skipping allocate_particle_data.\n");
    }

    if (sim_opts->drift == 1.) {
        if (initialize_mass_accumulation_file(output_files, sim_opts, disk_params, &header_data_for_files) != 0) {
            fprintf(stderr, "ERROR: Failed to set up initial output files. Exiting.\n");
            free_particle_data(&p_data);
            exit(EXIT_FAILURE);
        }
        fprintf(stderr, "DEBUG [tIntegrate]: Calling load_dust_particles to load initial dust profile.\n");
        load_dust_particles(&p_data, sim_opts->dust_input_filename, disk_params, sim_opts);
        fprintf(stderr, "DEBUG [tIntegrate]: load_dust_particles completed.\n");

        // --- NEW DEBUG PRINT: disk_params->eps after load_dust_particles ---
        fprintf(stderr, "DEBUG [tIntegrate]: disk_params->eps value after load_dust_particles: %.10e\n", disk_params->eps);
    }

    double max_dist_pop1 = 0.0;
    double min_dist_pop1 = 0.0;
    double max_dist_pop2 = 0.0;
    double min_dist_pop2 = 0.0;

    char dens_name[MAX_PATH_LEN] = "";
    char dust_name[MAX_PATH_LEN] = "";
    char micron_dust_name[MAX_PATH_LEN] = "";

    double t = 0.0; // Az idő belső orbitális egységben (radián vagy 2pi-hez viszonyítva)
    // TMAX években van, konvertáljuk belső egységre a ciklus feltételéhez
    double t_integration_in_internal_units = sim_opts->TMAX * (2.0 * M_PI);

    // Számoljuk ki a rács egyenletes lépésközét egyszer az elején
    const double delta_r_uniform_grid = (disk_params->RMAX - disk_params->RMIN) / (double)(disk_params->NGRID - 1);
    fprintf(stderr, "DEBUG [tIntegrate]: Calculated uniform grid spacing (delta_r_uniform_grid): %.8lg AU\n", delta_r_uniform_grid);

    // Gázra vonatkozó időlépés javaslat (valószínűleg a time_step függvényből)
    // Feltételezve, hogy a time_step már a Kepleri idő valamilyen töredékét adja vissza években.
    double gas_dt_suggestion_in_years = time_step(disk_params);
    fprintf(stderr, "DEBUG [tIntegrate]: Gas time step suggestion: %.2e years\n", gas_dt_suggestion_in_years);

    // A deltat_in_years a tényleges időlépés, ami majd érvényesül.
    // Kezdetben állítsuk be a gáz-időlépés javaslatra, vagy a felhasználóira.
    double deltat_in_years;
    if (sim_opts->DT > 0.0 && sim_opts->DT < gas_dt_suggestion_in_years) {
        deltat_in_years = sim_opts->DT;
        fprintf(stderr, "DEBUG [tIntegrate]: Initial deltat set to user-provided (smaller): %.2e years\n", deltat_in_years);
    } else {
        deltat_in_years = gas_dt_suggestion_in_years;
        fprintf(stderr, "DEBUG [tIntegrate]: Initial deltat set to gas suggestion (or user-provided was larger): %.2e years\n", deltat_in_years);
    }

    // A belső egységre konvertált deltat, amit az integrációhoz használunk
    double deltat = deltat_in_years * (2.0 * M_PI);
    fprintf(stderr, "DEBUG [tIntegrate]: Initial effective deltat (in internal units): %.2e\n", deltat);


    // --- Main Simulation Loop ---
    double overall_min_dist = HUGE_VAL;
    double overall_max_dist = 0.0;

    do {
        if (sim_opts->drift == 1.) {
            // --- Update particle reciprocals and find min/max distances ---
            for (int i = 0; i < p_data.num_particles_pop1; i++) {
                if (p_data.particles_pop1[i].distance_au > 0.0) {
                    p_data.particles_pop1[i].distance_au_reciprocal = 1.0 / p_data.particles_pop1[i].distance_au;
                } else {
                    p_data.particles_pop1[i].distance_au_reciprocal = 0.0;
                }
            }
            max_dist_pop1 = get_max_particle_distance(p_data.particles_pop1, p_data.num_particles_pop1);
            min_dist_pop1 = get_min_particle_distance(p_data.particles_pop1, p_data.num_particles_pop1, disk_params->RMIN);
            fprintf(stderr, "DEBUG [tIntegrate]: Pop1 MinDist: %.8lg AU, MaxDist: %.8lg AU\n", min_dist_pop1, max_dist_pop1);


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
                fprintf(stderr, "DEBUG [tIntegrate]: Pop2 MinDist: %.8lg AU, MaxDist: %.8lg AU\n", min_dist_pop2, max_dist_pop2);

                overall_min_dist = fmin(min_dist_pop1, min_dist_pop2);
                overall_max_dist = fmax(max_dist_pop1, max_dist_pop2);
            } else {
                overall_min_dist = min_dist_pop1;
                overall_max_dist = max_dist_pop1;
            }
            fprintf(stderr, "DEBUG [tIntegrate]: Overall MinDist: %.8lg AU, Overall MaxDist: %.8lg AU\n", overall_min_dist, overall_max_dist);


            // --- Output Data (Snapshot) Handling ---
            double current_time_years = t / (2.0 * M_PI);
            double output_interval_years = sim_opts->TMAX / sim_opts->WO;

            // Ellenőrizzük, hogy elérkeztünk-e egy output időponthoz, figyelembe véve a lebegőpontos pontosságot
            // A feltétel: (aktuális idő % output_időköz) kisebb, mint az időlépés, VAGY az aktuális idő 0.
            // És győződjünk meg róla, hogy csak egyszer írunk ki egy adott L értékhez.
            if ((fmod(current_time_years, output_interval_years) < deltat_in_years || current_time_years == 0) &&
                (current_time_years >= L - (deltat_in_years * 0.5))) { // L-hez való közelítés
                fprintf(stderr,"\n--- Simulation Time: %.2e years (Internal time: %.2e, L: %.2e) ---\n", current_time_years, t, L);

                // Generate filenames
                if (sim_opts->evol == 1 || current_time_years == 0) {
                    snprintf(dens_name, MAX_PATH_LEN, "%s/%s/%s_%08d.dat", sim_opts->output_dir_name, LOGS_DIR, FILE_DENS_PREFIX, (int)L);
                }
                snprintf(dust_name, MAX_PATH_LEN, "%s/%s/%s_%08d.dat", sim_opts->output_dir_name, LOGS_DIR, FILE_DUST_PREFIX,(int)L);
                snprintf(micron_dust_name, MAX_PATH_LEN, "%s/%s/micron_%s_%08d.dat", sim_opts->output_dir_name, LOGS_DIR, FILE_DUST_PREFIX, (int)L);


                // Open files and write headers
                output_files->surface_file = fopen(dens_name, "w");
                if (sim_opts->evol == 1 || current_time_years == 0) {
                    if (output_files->surface_file == NULL) {
                        fprintf(stderr, "ERROR: Could not open %s for writing.\n", dens_name);
                    } else {
                        HeaderData_t gas_header_data = {.current_time = current_time_years, .is_initial_data = (current_time_years == 0.0)};
                        write_file_header(output_files->surface_file, FILE_TYPE_GAS_DENSITY, &gas_header_data);
                    }
                }

                output_files->dust_file = fopen(dust_name, "w");
                if (output_files->dust_file == NULL) {
                    fprintf(stderr, "ERROR: Could not open %s for writing.\n", dust_name);
                } else {
                    HeaderData_t dust_header_data = {.current_time = current_time_years, .is_initial_data = (current_time_years == 0.0)};
                    write_file_header(output_files->dust_file, FILE_TYPE_DUST_EVOL, &dust_header_data);
                }

                if (sim_opts->twopop == 1.) {
                    output_files->micron_dust_file = fopen(micron_dust_name, "w");
                    if (output_files->micron_dust_file == NULL) {
                        fprintf(stderr, "ERROR: Could not open %s for writing.\n", micron_dust_name);
                    } else {
                        HeaderData_t micron_dust_header_data = {.current_time = current_time_years, .is_initial_data = (current_time_years == 0.0)};
                        write_file_header(output_files->micron_dust_file, FILE_TYPE_MICRON_DUST_EVOL, &micron_dust_header_data);
                    }
                }

                if (sim_opts->growth == 1.) {
                    // --- NEW DEBUG PRINT: disk_params->eps just before calculate_dust_density_grid ---
                    fprintf(stderr, "DEBUG [tIntegrate]: disk_params->eps value just before calling calculate_dust_density_grid: %.10e\n", disk_params->eps);
                    calculate_dust_density_grid(&p_data, disk_params, sim_opts);   // Order of args: p_data, sim_opts, disk_params -> I've fixed this earlier based on dust_physics.c
                                                                 // IMPORTANT: The order of arguments to calculate_dust_density_grid in simulation_loop.c (tIntegrate) was:
                                                                 // calculate_dust_density_grid(&p_data, sim_opts, disk_params);
                                                                 // While the prototype in dust_physics.c was:
                                                                 // void calculate_dust_density_grid(const ParticleData_t *p_data, disk_t *disk_params, const simulation_options_t *sim_opts);
                                                                 // This is a MISMATCH! The arguments for disk_params and sim_opts are swapped!
                                                                 // THIS IS LIKELY THE CAUSE OF YOUR PROBLEM!

                    // The call should be:
                    // calculate_dust_density_grid(&p_data, disk_params, sim_opts);
                    // I'm assuming you fixed this, but re-highlighting it just in case.
                }

                // Gas density output
                if (sim_opts->evol == 1 || current_time_years == 0) {
                    write_gas_profile_to_file(disk_params, output_files);
                }

                // Dust surface density output (Main dust & Micron dust if twopop)
                if (sim_opts->growth == 1.) {
                    write_dust_profile_to_file((int)L, &p_data, disk_params, sim_opts, output_files);
                }
                fprintf(stderr,"L set to %lg\n",L);

                L = L + output_interval_years;
                close_snapshot_files(output_files, dens_name, dust_name, micron_dust_name, sim_opts);
            }

            // Gas evolution
            if (sim_opts->evol == 1.) {
                get_gas_surface_density_pressure_pressure_gradient(sim_opts, disk_params);
            }

            if (p_data.particles_pop1 == NULL) {
                fprintf(stderr, "ERROR [tIntegrate]: particles_pop1 is NULL before update_particle_positions call!\n");
                exit(EXIT_FAILURE);
            }
            if (disk_params == NULL) {
                fprintf(stderr, "ERROR [tIntegrate]: disk_params is NULL before update_particle_positions call!\n");
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "DEBUG [tIntegrate]: Calling update_particle_positions for Pop1 with num_particles = %d\n", PARTICLE_NUMBER);

            // --- Részecske mozgás (update_particle_positions) és méret evolúció ---
            // A update_particle_positions függvény végzi el a Runge-Kutta integrációt
            // és frissíti a részecskék pozícióit és a drdt tagjait.
            update_particle_positions(p_data.particles_pop1, PARTICLE_NUMBER, deltat, t, sim_opts, disk_params); // 0 = Pop1
            if (sim_opts->twopop == 1.) {
                update_particle_positions(p_data.particles_pop2, PARTICLE_NUMBER, deltat, t, sim_opts, disk_params); // 1 = Pop2
            }

            // --- Adaptív Időlépés Számítása a Következő Lépéshez ---
            // Most, hogy a update_particle_positions lefutott, a particles_pop1[i].drdt (és pop2) tartalmazza az aktuális radiális sebességet.
            double particle_dt_drift_limit_in_years = HUGE_VAL;
            if (sim_opts->drift == 1.0 && PARTICLE_NUMBER > 0) {
                double max_abs_drdt_from_current_step = 0.0;
                for (int i = 0; i < p_data.num_particles_pop1; ++i) {
                    // Pozíció clampelés! Ez elengedhetetlen a nullázódás elkerüléséhez.
//                    if (p_data.particles_pop1[i].distance_au < disk_params->RMIN) {
//                        p_data.particles_pop1[i].distance_au = 0;
//                    }
                    if (fabs(p_data.particles_pop1[i].drdt) > max_abs_drdt_from_current_step) {
                        max_abs_drdt_from_current_step = fabs(p_data.particles_pop1[i].drdt);
                    }
                }
                if (sim_opts->twopop == 1) {
                    for (int i = 0; i < p_data.num_particles_pop2; ++i) {
//                        if (p_data.particles_pop2[i].distance_au < disk_params->RMIN) {
//                            p_data.particles_pop2[i].distance_au = 0;
//                        }
                        if (fabs(p_data.particles_pop2[i].drdt) > max_abs_drdt_from_current_step) {
                            max_abs_drdt_from_current_step = fabs(p_data.particles_pop2[i].drdt);
                        }
                    }
                }

                if (max_abs_drdt_from_current_step > 1e-15) { // Elkerüljük a nullával való osztást
                    particle_dt_drift_limit_in_years = CFL_FACTOR_DRIFT * delta_r_uniform_grid / max_abs_drdt_from_current_step;
                } else {
                    particle_dt_drift_limit_in_years = MAX_TIMESTEP_ABSOLUTE_YEARS; // Ha nincs számottevő drift, nagy dt engedélyezése
                }
                fprintf(stderr, "DEBUG [tIntegrate]: Particle drift time step limit: %.2e years (max_abs_drdt: %.2e)\n", particle_dt_drift_limit_in_years, max_abs_drdt_from_current_step);
            } else {
                particle_dt_drift_limit_in_years = MAX_TIMESTEP_ABSOLUTE_YEARS; // Ha nincs drift engedélyezve, vagy 0 részecske van
            }


            // A három időlépés közül a legkisebbet választjuk a következő lépéshez
            // 1. Gáz CFL (time_step)
            // 2. Por Drift CFL (particle_dt_drift_limit_in_years)
            // 3. Felhasználó által megadott fix DT (sim_opts->DT)
            deltat_in_years = fmin(gas_dt_suggestion_in_years, particle_dt_drift_limit_in_years);

            if (sim_opts->DT > 0.0) { // Ha a felhasználó megadott egy fix DT-t
                deltat_in_years = fmin(deltat_in_years, sim_opts->DT);
            }

            // Végső korlátozások
            deltat_in_years = fmax(deltat_in_years, MIN_TIMESTEP_YEARS); // Ne essen túl kicsire
            deltat_in_years = fmin(deltat_in_years, MAX_TIMESTEP_ABSOLUTE_YEARS); // Ne nőjön túl nagyra
            deltat_in_years = fmin(deltat_in_years, output_interval_years); // Ne lépje túl az output intervallumot


            // KONVERTÁLJUK a végleges időlépést a belső orbitális egységre!
            deltat = deltat_in_years * (2.0 * M_PI);
            fprintf(stderr, "DEBUG [tIntegrate]: New effective deltat (in internal units): %.2e (%.2e years)\n", deltat, deltat_in_years);


            // Az idő léptetése a belső egységben
            t = t + deltat;

            // Termination condition for drift == 1 branch
            double current_time_years_after_step = t / (2.0 * M_PI);
            fprintf(stderr, "DEBUG [tIntegrate]: Termination check at time %.2e years: overall_max_dist (%.8lg AU) >= RMIN (%.8lg AU) && overall_min_dist (%.8lg AU) != overall_max_dist (%.8lg AU)\n",
                            current_time_years_after_step, overall_max_dist, disk_params->RMIN, overall_min_dist, overall_max_dist);
            // Terminációs feltétel: a részecskék mind RMIN-en belülre kerültek, VAGY az összes részecske egy pontban van.
            if (!(overall_max_dist >= disk_params->RMIN && overall_min_dist != overall_max_dist)) {
                fprintf(stderr,"DEBUG [tIntegrate]: Simulation termination condition met (overall_max_dist < RMIN or overall_min_dist == overall_max_dist) at time: %.2e years.\n", current_time_years_after_step);
                goto cleanup;
            }

        } else { // sim_opts->drift == 0. (Gas-only simulation)
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
            get_gas_surface_density_pressure_pressure_gradient(sim_opts, disk_params);

            // Az idő léptetése a belső egységben
            t = t + deltat;
        }

    // A ciklus feltétele a belső egységgel van ellenőrizve
    } while (t <= t_integration_in_internal_units);

    fprintf(stderr,"\n\nDEBUG [tIntegrate]: Main simulation loop finished (t > t_integration).\n");

cleanup:
    // --- Cleanup Section ---
    cleanup_simulation_resources(&p_data, output_files, sim_opts);
    fprintf(stderr,"DEBUG [tIntegrate]: Cleanup completed.\n");
}
