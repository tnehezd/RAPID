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

void integrateParticleRungeKutta4(double particle_radius,
                                  double actual_timestep, double particle_distance, double *particle_distance_new, double *particle_radius_new, const DiskParameters 
                                  *disk_params){

    double k1_slope, k2_slope, k3_slope, k4_slope;
    double intermediate_radius;
    double gas_sigma, gas_grad, gas_vrad;
    
    // 1. HELYI GÁZPROFIL INTERPOLÁCIÓJA A KIINDULÓPONTON
    linearInterpolation(disk_params->gas_surface_density_vector, disk_params->radial_grid, particle_distance, &gas_sigma, disk_params->delta_r, disk_params);
    linearInterpolation(disk_params->gas_pressure_gradient_vector, disk_params->radial_grid, particle_distance, &gas_grad, disk_params->delta_r, disk_params);
    linearInterpolation(disk_params->gas_velocity_vector, disk_params->radial_grid, particle_distance, &gas_vrad, disk_params->delta_r, disk_params);
    *particle_radius_new = particle_radius;

    // 3. RUNGE-KUTTA 4 INTEGRÁCIÓ (Valódi térbeli újramintavételezéssel!)
    
    // K1 - Kiindulóponton
    calculate1DDustDrift(particle_radius, gas_grad, gas_sigma, gas_vrad, particle_distance, &k1_slope, disk_params);

    // K2 - Első köztes ponton (újra-interpoláljuk a változó gázteret!)
    intermediate_radius = particle_distance + 0.5 * actual_timestep * k1_slope;
    linearInterpolation(disk_params->gas_surface_density_vector, disk_params->radial_grid, intermediate_radius, &gas_sigma, disk_params->delta_r, disk_params);
    linearInterpolation(disk_params->gas_pressure_gradient_vector, disk_params->radial_grid, intermediate_radius, &gas_grad, disk_params->delta_r, disk_params);
    linearInterpolation(disk_params->gas_velocity_vector, disk_params->radial_grid, intermediate_radius, &gas_vrad, disk_params->delta_r, disk_params);
    calculate1DDustDrift(particle_radius, gas_grad, gas_sigma, gas_vrad, intermediate_radius, &k2_slope, disk_params);
        
    // K3 - Második köztes ponton (újra-interpolálunk)
    intermediate_radius = particle_distance + 0.5 * actual_timestep * k2_slope;
    linearInterpolation(disk_params->gas_surface_density_vector, disk_params->radial_grid, intermediate_radius, &gas_sigma, disk_params->delta_r, disk_params);
    linearInterpolation(disk_params->gas_pressure_gradient_vector, disk_params->radial_grid, intermediate_radius, &gas_grad, disk_params->delta_r, disk_params);
    linearInterpolation(disk_params->gas_velocity_vector, disk_params->radial_grid, intermediate_radius, &gas_vrad, disk_params->delta_r, disk_params);
    calculate1DDustDrift(particle_radius, gas_grad, gas_sigma, gas_vrad, intermediate_radius, &k3_slope, disk_params);
    
    // K4 - Végponton (újra-interpolálunk)
    intermediate_radius = particle_distance + actual_timestep * k3_slope;
    linearInterpolation(disk_params->gas_surface_density_vector, disk_params->radial_grid, intermediate_radius, &gas_sigma, disk_params->delta_r, disk_params);
    linearInterpolation(disk_params->gas_pressure_gradient_vector, disk_params->radial_grid, intermediate_radius, &gas_grad, disk_params->delta_r, disk_params);
    linearInterpolation(disk_params->gas_velocity_vector, disk_params->radial_grid, intermediate_radius, &gas_vrad, disk_params->delta_r, disk_params);
    calculate1DDustDrift(particle_radius, gas_grad, gas_sigma, gas_vrad, intermediate_radius, &k4_slope, disk_params);

    // Súlyozott RK4 átlagolás az új pozícióhoz
    *particle_distance_new = particle_distance + actual_timestep * (k1_slope + 2.0 * k2_slope + 2.0 * k3_slope + k4_slope) / 6.0;
}
