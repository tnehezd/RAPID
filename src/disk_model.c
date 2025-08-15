// src/disk_model.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h> // for DBL_MAX, or use a very large number

#include "disk_model.h"
#include "config.h"
#include "simulation_types.h"
#include "globals.h"
#include "dust_physics.h"
#include "io_utils.h"
#include "utils.h"
#include "simulation_core.h"

// This file contains the implementation of functions for initializing and evolving gas disk properties.

// Egy nagyon kicsi, de nem nulla érték, amivel elkerülhető
// a nullával való osztás.
#define EPSILON 1e-10

// A tanh argumentumának maximális értéke a numerikus stabilitás
// érdekében. A tanh(50) már gyakorlatilag 1.0, így ez biztonságos.
#define TANH_MAX_ARG 100.0

/* --- Disk Parameter Initialization and Setup Implementations --- */


void initialize_grid_cells(disk_t *disk_params) {
    int i;
    for(i = 0; i <= disk_params->NGRID + 1; i++) {
        // Initialize radial grid cell positions, including ghost cells at i=0 and i=NGRID+1
        disk_params->rvec[i] = disk_params->RMIN + (i - 1) * disk_params->DD;
    }
}

void initial_gas_surface_density_profile(disk_t *disk_params){
    int i;

    // Apply a power-law initial profile: sigma = SIGMA0 * r^(SIGMAP_EXP)
    for(i = 1; i <= disk_params->NGRID; i++) {
        disk_params->sigmavec[i] = disk_params->SIGMA0 * pow(disk_params->rvec[i], disk_params->SIGMAP_EXP);
    }

    calculate_boundary(disk_params->sigmavec, disk_params);
}

void initial_gas_pressure_profile(disk_t *disk_params){
    int i;

    // Calculate initial gas pressure for each grid cell
    for(i = 1; i <= disk_params->NGRID; i++) {
        disk_params->pressvec[i] = calculate_gas_pressure(disk_params->sigmavec[i], disk_params->rvec[i], disk_params);
    }
    calculate_boundary(disk_params->pressvec, disk_params);
}

void initial_gas_pressure_gradient_profile(disk_t *disk_params){
    // Calculate initial radial gas pressure gradient
    calculate_gas_pressure_gradient(disk_params);
    calculate_boundary(disk_params->dpressvec, disk_params);
}

void initial_gas_velocity_profile(disk_t *disk_params){
    // Calculate initial radial gas velocity
    calculate_gas_velocity(disk_params);
    calculate_boundary(disk_params->ugvec, disk_params);
}

// Calculating the turbulent alpha parameter --> reducing alpha with alpha_r
double calculate_turbulent_alpha(double r, const disk_t *disk_params) {
    // Kiszámítja a turbulens alfa paramétert. Ez a verzió robusztusabb a numerikus stabilitás szempontjából,
    // megakadályozva a nullával való osztást és a nagy tanh argumentumokat.

    // A Dr_dze_i és Dr_dze_o paraméterekhez EPSILON-t adunk,
    // hogy elkerüljük a nullával való osztást.
    double drdze_i_safe = disk_params->Dr_dze_i;
    double drdze_o_safe = disk_params->Dr_dze_o;

    // Az EPSILON hozzáadása, ha az érték nullához közelít.
    if (fabs(drdze_i_safe) < EPSILON) {
        drdze_i_safe = (drdze_i_safe >= 0) ? EPSILON : -EPSILON;
    }
    if (fabs(drdze_o_safe) < EPSILON) {
        drdze_o_safe = (drdze_o_safe >= 0) ? EPSILON : -EPSILON;
    }

    double tanh_arg1 = (r - disk_params->r_dze_i) / drdze_i_safe;
    double tanh_arg2 = (disk_params->r_dze_o - r) / drdze_o_safe;

    // A tanh argumentumának korlátozása a numerikus stabilitás érdekében.
    // Ezzel elkerülhető a nagyon nagy argumentumokból adódó hiba.
    if (tanh_arg1 > TANH_MAX_ARG) tanh_arg1 = TANH_MAX_ARG;
    if (tanh_arg1 < -TANH_MAX_ARG) tanh_arg1 = -TANH_MAX_ARG;
    if (tanh_arg2 > TANH_MAX_ARG) tanh_arg2 = TANH_MAX_ARG;
    if (tanh_arg2 < -TANH_MAX_ARG) tanh_arg2 = -TANH_MAX_ARG;

    double alpha_r = 1.0 - 0.5 * (1.0 - disk_params->a_mod) * (tanh(tanh_arg1) + tanh(tanh_arg2));

    // A visszaadott értéknek mindig pozitívnak kell lennie.
    // JAVÍTVA: Ahelyett, hogy 0.0-t ad vissza, ami nullás viszkozitást eredményez,
    // egy kis pozitív értéket ad vissza, ami megakadályozza a nullával való osztást.
    if (alpha_r < EPSILON) {
        return EPSILON * disk_params->alpha_visc;
    }
    
    return alpha_r * disk_params->alpha_visc;
}


/* --- Gas Disk Property Calculation Implementations --- */

double calculate_gas_viscosity(double r, const disk_t *disk_params) {
    double nu;
    double cs, H;

    H = calculate_scale_height(r, disk_params);
    cs = calculate_local_sound_speed(r, disk_params);

    // nu = alpha * cs * H
    nu = calculate_turbulent_alpha(r, disk_params) * cs * H;
    return nu;
}

double calculate_scale_height(double r, const disk_t *disk_params) {
    if (disk_params == NULL) {
        fprintf(stderr, "ERROR [calculate_scale_height]: disk_params is NULL!\n");
        return 0.0;
    }

    // H = HASP * r^(1+FLIND)
    double calculated_result = pow(r, 1. + disk_params->FLIND) * disk_params->HASP;
    return calculated_result;
}

double calculate_keplerian_velocity(double r, const disk_t *disk_params) {
    // v_K = sqrt(G * M_star / r)
    return sqrt(G_GRAV_CONST * disk_params->STAR_MASS / r);
}

double calculate_keplerian_angular_velocity(double r, const disk_t *disk_params) {
    // omega_K = sqrt(G * M_star / r^3)
    return sqrt(G_GRAV_CONST * disk_params->STAR_MASS / r / r / r);
}

double calculate_local_sound_speed(double r, const disk_t *disk_params) {
    // c_s = omega_K * H
    return calculate_keplerian_angular_velocity(r, disk_params) * calculate_scale_height(r, disk_params);
}

double calculate_midplane_gas_density(double sigma, double r, const disk_t *disk_params) {
    // rho_midplane = (1 / sqrt(2*PI)) * Sigma / H (for Gaussian vertical profile)
    double H_val = calculate_scale_height(r, disk_params);
    if (H_val < 1e-12) { // A very small number to prevent division by zero
        fprintf(stderr, "WARNING [calculate_midplane_gas_density]: H is too small at r=%.10e. Returning 0.\n", r);
        return 0.0;
    }
    return 1. / sqrt(2.0 * M_PI) * sigma / H_val;
}

double calculate_gas_pressure(double sigma, double r, const disk_t *disk_params) {
    // p = rho_gas * c_s^2
    double cs_val = calculate_local_sound_speed(r, disk_params);
    return calculate_midplane_gas_density(sigma, r, disk_params) * cs_val * cs_val;
}

void calculate_gas_pressure_gradient(disk_t *disk_params) {
    int i;
    double ptemp, pvec[disk_params->NGRID + 2];

    // Calculate radial pressure gradient using central finite difference
    for (i = 1; i <= disk_params->NGRID; i++) {
        // Safe check for rvec to prevent division by zero
        if (fabs(disk_params->rvec[i + 1] - disk_params->rvec[i - 1]) < 1e-12) {
             pvec[i] = 0.0;
             continue;
        }

        ptemp = (disk_params->pressvec[i + 1] - disk_params->pressvec[i - 1]) / (2.0 * disk_params->DD);
        pvec[i] = ptemp;
    }
    // Copy results to the disk_params structure's dpressvec array
    for (i = 1; i <= disk_params->NGRID; i++) {
        disk_params->dpressvec[i] = pvec[i];
    }
}

double coefficient_for_gas_velocity(double sigma, double r) {
    // Coefficient for the gas radial velocity equation: -3 / (Sigma * R^0.5)
    double denominator = sigma * sqrt(r);
    if (denominator == 0.0 || isnan(denominator) || isinf(denominator)) {
        return 0.0; // Handle edge cases to prevent division by zero or sqrt of negative
    }
    return -1.0 * (3.0 / denominator);
}

void calculate_gas_velocity(disk_t *disk_params) {
    double tempug;
    double ugvec_temp_calc[disk_params->NGRID + 2]; // Holds intermediate (nu * Sigma * R^0.5) terms
    double ug_derivative_terms[disk_params->NGRID + 1]; // Holds the derivative part of the equation

    int i;

    // Calculate (nu * Sigma * R^0.5) for each grid point
    #pragma omp parallel for private(i)
    for (i = 0; i <= disk_params->NGRID + 1; i++) {
        ugvec_temp_calc[i] = disk_params->sigmavec[i] * calculate_gas_viscosity(disk_params->rvec[i], disk_params) * sqrt(disk_params->rvec[i]);
    }

    // Calculate the radial derivative using central finite difference
    #pragma omp parallel for private(i, tempug)
    for (i = 1; i <= disk_params->NGRID; i++) {
        if (fabs(disk_params->rvec[i + 1] - disk_params->rvec[i - 1]) < 1e-12) {
             ug_derivative_terms[i] = 0.0;
             continue;
        }
        
        tempug = (ugvec_temp_calc[i + 1] - ugvec_temp_calc[i - 1]) / (2.0 * disk_params->DD);
        ug_derivative_terms[i] = coefficient_for_gas_velocity(disk_params->sigmavec[i], disk_params->rvec[i]) * tempug;
    }

    // Assign the calculated radial velocities to the disk_params structure
    #pragma omp parallel for private(i)
    for (i = 1; i <= disk_params->NGRID; i++) {
        disk_params->ugvec[i] = ug_derivative_terms[i];
    }
}

void get_gas_surface_density_pressure_pressure_gradient(const simulation_options_t *sim_opts, disk_t *disk_params, double dt) {
    double u, u_bi, u_fi;
    // Ideiglenes tömb a (sigma * nu) frissített értékeinek tárolására
    double uvec_new[disk_params->NGRID + 2]; 

    int i;
    
    // DEBUG: Ellenőrzés a függvény elején
//    fprintf(stderr, "INFO: get_gas_surface_density_pressure_pressure_gradient belépés. Sigma a 10. ponton: %e\n", disk_params->sigmavec[10]);

    // Uvec tömb inicializálása a jelenlegi adatok alapján
    double uvec_old[disk_params->NGRID + 2];
//    #pragma omp parallel for private(i)
    for(i = 0; i <= disk_params->NGRID + 1; i++) {
        double r_val = disk_params->rvec[i];
        double sigma_val = disk_params->sigmavec[i];
        double nu_val = calculate_gas_viscosity(r_val, disk_params);
        uvec_old[i] = sigma_val * nu_val;
//        fprintf(stderr," INFO: r: %lg UVEC %lg\n", disk_params->rvec[i], uvec_old[i]);
        // uvec_new[i] = uvec_old[i]; // Inicializálás
    }

    // A diffúziós egyenlet megoldása az uvec_new tömbbe
    for (i = 1; i <= disk_params->NGRID; i++) {
        u = uvec_old[i];
        u_bi = uvec_old[i - 1];
        u_fi = uvec_old[i + 1];

        double temp = Coeff_1(disk_params->rvec[i], disk_params) * (u_fi - 2.0 * u + u_bi) / (disk_params->DD * disk_params->DD) +
                      Coeff_2(disk_params->rvec[i], disk_params) * (u_fi - u_bi) / (2.0 * disk_params->DD);
        
        uvec_new[i] = uvec_old[i] + dt * temp;

        if (uvec_new[i] < 0.0) {
            uvec_new[i] = 0.0;
        }

        // DEBUG: Itt kiírhatod a temp értékét, hogy lásd, nem nulla-e
        // fprintf(stderr, "DEBUG: r=%e, temp=%e, DT*temp=%e\n", disk_params->rvec[i], temp, sim_opts->DT * temp);
    }
    
    // Szegélyfeltételek beállítása az uvec_new tömbre
    uvec_new[0] = uvec_old[0];
    uvec_new[disk_params->NGRID + 1] = uvec_old[disk_params->NGRID + 1];

    // A sigmavec frissítése az uvec_new tömb alapján
    #pragma omp parallel for private(i)
    for (i = 0; i <= disk_params->NGRID + 1; i++) {
        double r_val = disk_params->rvec[i];
        double nu_val = calculate_gas_viscosity(r_val, disk_params);

        if (nu_val > DBL_EPSILON) {
            disk_params->sigmavec[i] = uvec_new[i] / nu_val;
        } else {
            disk_params->sigmavec[i] = 0.0;
        }
    }
    
    // A nyomás és a nyomás gradiens frissítése
    // Ez a lépés most már a frissített sigmavec-et használja
    #pragma omp parallel for private(i)
    for (i = 0; i <= disk_params->NGRID + 1; i++) {
         disk_params->pressvec[i] = calculate_gas_pressure(disk_params->sigmavec[i], disk_params->rvec[i], disk_params);
    }
    
    calculate_gas_pressure_gradient(disk_params);

    // Szegélyfeltételek alkalmazása
    calculate_boundary(disk_params->sigmavec, disk_params);
    calculate_boundary(disk_params->pressvec, disk_params);
    calculate_boundary(disk_params->dpressvec, disk_params);
    
    // DEBUG: Ellenőrzés a függvény végén
//    fprintf(stderr, "INFO: get_gas_surface_density_pressure_pressure_gradient kilépés. Sigma a 10. ponton: %e\n", disk_params->sigmavec[10]);
}
