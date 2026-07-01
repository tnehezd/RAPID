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

void integrateParticleRungeKutta4(double actual_time, double particle_radius, const double *dust_surfacedensity, const double *particle_distance_grid, 
                                  double actual_timestep, double particle_distance, double *particle_distance_new, double *particle_radius_new, const DiskParameters 
                                  *disk_params, const SimulationOptions *simulation_options){

    double k1_slope,k2_slope,k3_slope,k4_slope;
    double intermediate_radius, grid_aligned_radius;
    double actual_gas_surfacedensity, actual_gas_pressure, actual_gas_pressure_gradient, actual_gas_radial_velocity; 
    double intrinsic_particle_density, updated_particle_radius;
    double actual_dust_surfacedensity = 0.0;
    
    linearInterpolation(disk_params->gas_surface_density_vector,disk_params->radial_grid,particle_distance,&actual_gas_surfacedensity,disk_params->delta_r,disk_params);
    linearInterpolation(disk_params->gas_pressure_gradient_vector,disk_params->radial_grid,particle_distance,&actual_gas_pressure_gradient,disk_params->delta_r,disk_params);
    linearInterpolation(disk_params->gas_velocity_vector,disk_params->radial_grid,particle_distance,&actual_gas_radial_velocity,disk_params->delta_r,disk_params);

    double particle_grid_spacing = (disk_params->r_max - disk_params->r_min) / (particle_number-1);
    int particle_grid_density = (int)(1./particle_grid_spacing);//
    particle_grid_density = particle_grid_density * ROUNDING_FACTOR;
    double particle_grid_density_double = (double) particle_grid_density;
    int nearest_particle_index;

    nearest_particle_index = (int)floor(particle_distance * particle_grid_density_double+0.5);
    grid_aligned_radius = (double)nearest_particle_index / particle_grid_density_double;
    
    int i;
    for(i=0;i<particle_number;i++) {
        if(grid_aligned_radius == particle_distance_grid[i]) {
            actual_dust_surfacedensity = dust_surfacedensity[i];
            break;
        }
    }

    if(simulation_options->option_for_dust_growth == 1.) {
        if(actual_time != 0.) {
            updated_particle_radius = particle_radius;
            linearInterpolation(disk_params->gas_pressure_vector,disk_params->radial_grid,particle_distance,&actual_gas_pressure,disk_params->delta_r,disk_params);
            intrinsic_particle_density = disk_params->particle_density; 
            updated_particle_radius = calculateDustParticleSize(particle_radius,intrinsic_particle_density,actual_gas_surfacedensity,actual_dust_surfacedensity,particle_distance,actual_gas_pressure_gradient,actual_gas_pressure_gradient,actual_timestep,disk_params);
            particle_radius = updated_particle_radius;
        }
    }

    *particle_radius_new = particle_radius;

    calculate1DDustDrift(particle_radius, actual_gas_pressure_gradient, actual_gas_surfacedensity, actual_gas_radial_velocity, particle_distance, &k1_slope,disk_params);

    intermediate_radius = particle_distance + 0.5 * actual_timestep * k1_slope;
    calculate1DDustDrift(particle_radius, actual_gas_pressure_gradient, actual_gas_surfacedensity, actual_gas_radial_velocity, intermediate_radius, &k2_slope,disk_params);
        
    intermediate_radius = particle_distance + 0.5 * actual_timestep * k2_slope;
    calculate1DDustDrift(particle_radius, actual_gas_pressure_gradient, actual_gas_surfacedensity, actual_gas_radial_velocity, intermediate_radius, &k3_slope,disk_params);
    
    intermediate_radius = particle_distance + actual_timestep * k3_slope;
    calculate1DDustDrift(particle_radius, actual_gas_pressure_gradient, actual_gas_surfacedensity, actual_gas_radial_velocity, intermediate_radius, &k4_slope,disk_params);

    *particle_distance_new = particle_distance + actual_timestep * (k1_slope + 2.0 * k2_slope + 2.0 * k3_slope + k4_slope) / 6.0;

}



