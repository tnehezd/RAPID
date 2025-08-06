// src/disk_model.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "disk_model.h"
#include "config.h"
#include "simulation_types.h"
#include "globals.h"
#include "dust_physics.h"
#include "io_utils.h"
#include "utils.h"
#include "simulation_core.h"

// This file contains the implementation of functions for initializing and evolving gas disk properties.

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
    double alpha_r;
    alpha_r = 1.0 - 0.5 * (1.0 - disk_params->a_mod) * (tanh((r - disk_params->r_dze_i) / disk_params->Dr_dze_i) + tanh((disk_params->r_dze_o - r) / disk_params->Dr_dze_o));
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
    return 1. / sqrt(2.0 * M_PI) * sigma / calculate_scale_height(r, disk_params);
}

double calculate_gas_pressure(double sigma, double r, const disk_t *disk_params) {
    // p = rho_gas * c_s^2
    return calculate_midplane_gas_density(sigma, r, disk_params) * calculate_local_sound_speed(r, disk_params) * calculate_local_sound_speed(r, disk_params);
}

void calculate_gas_pressure_gradient(disk_t *disk_params) {
    int i;
    double ptemp, pvec[disk_params->NGRID + 2];

    // Calculate radial pressure gradient using central finite difference
    for (i = 1; i <= disk_params->NGRID; i++) {
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
    if (sigma == 0.0 || r <= 0.0) {
        return 0.0; // Handle edge cases to prevent division by zero or sqrt of negative
    }
    return -1.0 * (3.0 / (sigma * sqrt(r)));
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
        tempug = (ugvec_temp_calc[i + 1] - ugvec_temp_calc[i - 1]) / (2.0 * disk_params->DD);
        ug_derivative_terms[i] = coefficient_for_gas_velocity(disk_params->sigmavec[i], disk_params->rvec[i]) * tempug;
    }

    // Assign the calculated radial velocities to the disk_params structure
    #pragma omp parallel for private(i)
    for (i = 1; i <= disk_params->NGRID; i++) {
        disk_params->ugvec[i] = ug_derivative_terms[i];
    }
}

void get_gas_surface_density_pressure_pressure_gradient(const simulation_options_t *sim_opts, disk_t *disk_params) {
    double u, u_bi, u_fi;
    double sigma_temp[disk_params->NGRID + 2]; // Temporary array for updated surface density
    double uvec[disk_params->NGRID + 2];       // Temporary array for (sigma * nu) or similar terms

    int i;

    // Set boundary conditions for temporary sigma array
    sigma_temp[0] = disk_params->sigmavec[0];
    sigma_temp[disk_params->NGRID + 1] = disk_params->sigmavec[disk_params->NGRID + 1];

    // Initialize uvec at boundary points
    uvec[0] = disk_params->sigmavec[0] * calculate_gas_viscosity(disk_params->rvec[0], disk_params);
    uvec[disk_params->NGRID + 1] = disk_params->sigmavec[disk_params->NGRID + 1] * calculate_gas_viscosity(disk_params->rvec[disk_params->NGRID + 1], disk_params);

    // Calculate uvec for internal grid points
    #pragma omp parallel for private(i)
    for(i = 1; i <= disk_params->NGRID; i++) {
        uvec[i] = disk_params->sigmavec[i] * calculate_gas_viscosity(disk_params->rvec[i], disk_params);
    }

    // This loop calculates the new sigma_temp using a finite difference scheme for gas evolution.
    // It processes values sequentially due to data dependencies.
    for (i = 1; i <= disk_params->NGRID; i++) {
        u = uvec[i];
        u_bi = uvec[i - 1];
        u_fi = uvec[i + 1];

        // Numerical update for sigma_temp (e.g., discretized diffusion equation)
        double temp = Coeff_1(disk_params->rvec[i], disk_params) * (u_fi - 2.0 * u + u_bi) / (disk_params->DD * disk_params->DD) +
                      Coeff_2(disk_params->rvec[i], disk_params) * (u_fi - u_bi) / (2.0 * disk_params->DD);
        
        sigma_temp[i] = uvec[i] + sim_opts->DT * temp;
    }

    // Update the actual disk parameters based on the calculated sigma_temp
    #pragma omp parallel for private(i)
    for (i = 1; i <= disk_params->NGRID; i++) {
        disk_params->sigmavec[i] = sigma_temp[i] / calculate_gas_viscosity(disk_params->rvec[i], disk_params);
        disk_params->pressvec[i] = calculate_gas_pressure(disk_params->sigmavec[i], disk_params->rvec[i], disk_params);
    }

    // Recalculate pressure gradient and apply boundary conditions to all relevant arrays.
    calculate_gas_pressure_gradient(disk_params);
    calculate_boundary(disk_params->sigmavec, disk_params);
    calculate_boundary(disk_params->pressvec, disk_params);
    calculate_boundary(disk_params->dpressvec, disk_params);
}