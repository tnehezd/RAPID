// src/integrator.h

#include <stdio.h>    
#include <stdlib.h>   
#include <math.h>     
#include <string.h>   
#include <omp.h>
#include "config.h"      
#include "ascii_output.h"     
#include "disk_model.h"   
#include "dust_physics.h" 
#include "utils.h"        
#include "simulation_core.h"
#include "particle_data.h"
#include "gas_physics.h"
#include "boundary_conditions.h"
#include "integrator.h"

void integrateParticleRungeKutta4Structured(double actual_time, double particle_radius, double actual_timestep, double particle_distance, double *particle_distance_new, double *particle_radius_new, const DiskParameters *disk_params, const SimulationOptions *simulation_options) {

    double k1_slope, k2_slope, k3_slope, k4_slope;
    double intermediate_r;
    double gas_sigma, gas_gradP, gas_vrad;
    double actual_gas_pressure = 0.0;
    double actual_dust_sigma = 0.0;

    // --- STEP 1: k1 ---
    linearInterpolation(disk_params->gas_surface_density_vector, disk_params->radial_grid, particle_distance, &gas_sigma, disk_params->delta_r, disk_params);
    linearInterpolation(disk_params->gas_pressure_gradient_vector, disk_params->radial_grid, particle_distance, &gas_gradP, disk_params->delta_r, disk_params);
    linearInterpolation(disk_params->gas_velocity_vector, disk_params->radial_grid, particle_distance, &gas_vrad, disk_params->delta_r, disk_params);
    calculate1DDustDrift(particle_radius, gas_gradP, gas_sigma, gas_vrad, particle_distance, &k1_slope, disk_params);

    // --- STEP 2: k2 ---
    intermediate_r = particle_distance + 0.5 * actual_timestep * k1_slope;
    linearInterpolation(disk_params->gas_surface_density_vector, disk_params->radial_grid, intermediate_r, &gas_sigma, disk_params->delta_r, disk_params);
    linearInterpolation(disk_params->gas_pressure_gradient_vector, disk_params->radial_grid, intermediate_r, &gas_gradP, disk_params->delta_r, disk_params);
    linearInterpolation(disk_params->gas_velocity_vector, disk_params->radial_grid, intermediate_r, &gas_vrad, disk_params->delta_r, disk_params);
    calculate1DDustDrift(particle_radius, gas_gradP, gas_sigma, gas_vrad, intermediate_r, &k2_slope, disk_params);

    // --- STEP 3: k3 ---
    intermediate_r = particle_distance + 0.5 * actual_timestep * k2_slope;
    linearInterpolation(disk_params->gas_surface_density_vector, disk_params->radial_grid, intermediate_r, &gas_sigma, disk_params->delta_r, disk_params);
    linearInterpolation(disk_params->gas_pressure_gradient_vector, disk_params->radial_grid, intermediate_r, &gas_gradP, disk_params->delta_r, disk_params);
    linearInterpolation(disk_params->gas_velocity_vector, disk_params->radial_grid, intermediate_r, &gas_vrad, disk_params->delta_r, disk_params);
    calculate1DDustDrift(particle_radius, gas_gradP, gas_sigma, gas_vrad, intermediate_r, &k3_slope, disk_params);

    // --- STEP 4: k4 ---
    intermediate_r = particle_distance + actual_timestep * k3_slope;
    linearInterpolation(disk_params->gas_surface_density_vector, disk_params->radial_grid, intermediate_r, &gas_sigma, disk_params->delta_r, disk_params);
    linearInterpolation(disk_params->gas_pressure_gradient_vector, disk_params->radial_grid, intermediate_r, &gas_gradP, disk_params->delta_r, disk_params);
    linearInterpolation(disk_params->gas_velocity_vector, disk_params->radial_grid, intermediate_r, &gas_vrad, disk_params->delta_r, disk_params);
    calculate1DDustDrift(particle_radius, gas_gradP, gas_sigma, gas_vrad, intermediate_r, &k4_slope, disk_params);

    // --- Final Position ---
    *particle_distance_new = particle_distance + (actual_timestep / 6.0) * (k1_slope + 2.0*k2_slope + 2.0*k3_slope + k4_slope);

    // --- Dust Growth ---
    if (simulation_options->option_for_dust_growth == 1.0 && actual_time != 0.0) {
        // Interpolálás Euler-rácsról
        double r_rel = (particle_distance - disk_params->r_min) / disk_params->delta_r;
        int idx = (int)floor(r_rel);
        double f = r_rel - (double)idx; // lineáris interpoláció a két cella között

        if (idx < 0) idx = 0;
        if (idx >= disk_params->grid_number - 1) idx = disk_params->grid_number - 2;

        actual_dust_sigma = (1.0 - f) * disk_params->dust_surface_density_euler[idx] + f * disk_params->dust_surface_density_euler[idx + 1];

        // Frissített részecske méret
        *particle_radius_new = calculateDustParticleSize(particle_radius, disk_params->particle_density, gas_sigma, actual_dust_sigma, particle_distance, gas_gradP, actual_gas_pressure, actual_timestep, disk_params);
    } else {
        *particle_radius_new = particle_radius;
    }
}
