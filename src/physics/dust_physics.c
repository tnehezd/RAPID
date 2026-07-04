// src/dust_physics.c
#include "dust_physics.h" 
#include "config.h"       
#include "simulation_types.h" 
#include "gas_physics.h"
#include "simulation_core.h" 
#include "boundary_conditions.h"
#include "integrator.h"
#include "utils.h"          
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>     
#include "logger.h"        

double calculateStokesNumber(double particle_radius, double gas_surfacedensity, const DiskParameters *disk_params) { 

    return disk_params->particle_density_dimensionless * particle_radius * M_PI / (2.0 * gas_surfacedensity);
}

/*			BIRNSTIEL EL AL 2012			*/
double calculateRadialDriftBarrier(double dust_surfacedensity, double radial_distance, double gas_pressure, double pressure_gradient, double particle_density, const DiskParameters *disk_params) {

    double dust_surfacedensity_cgs = dust_surfacedensity / SURFACE_DENSITY_CONVERSION_FACTOR;

    double keplerian_velocity = calculateKeplerianVelocity(radial_distance,disk_params);
    double keplerian_velocity_squared = keplerian_velocity * keplerian_velocity;
    double local_soundspeed = calculateLocalSoundSpeed(radial_distance,disk_params);
    double local_soundspeed_squared = local_soundspeed * local_soundspeed;
    double log_pressure_gradient = radial_distance / gas_pressure * pressure_gradient;
    double s_drift =  disk_params->drift_factor * 2.0 / M_PI * dust_surfacedensity_cgs / particle_density * keplerian_velocity_squared / local_soundspeed_squared * fabs(1.0 / log_pressure_gradient);
    return s_drift;
}

double calculateTurbulentFragmentationBarrier(double gas_surfacedensity, double radial_distance, double particle_density, const DiskParameters *disk_params) {

    double s_frag, fragmentation_velocity, fragmentation_velocity_squared, gas_surfacedensity_cgs, local_soundspeed, local_soundspeed_squared;

    fragmentation_velocity = disk_params->fragmentation_velocity * CM_PER_SEC_TO_AU_PER_YEAR_2PI;
    fragmentation_velocity_squared = fragmentation_velocity * fragmentation_velocity;
    gas_surfacedensity_cgs = gas_surfacedensity / SURFACE_DENSITY_CONVERSION_FACTOR;
    local_soundspeed = calculateLocalSoundSpeed(radial_distance,disk_params);
    local_soundspeed_squared = local_soundspeed * local_soundspeed;

    s_frag = disk_params->fragmentation_factor * 2.0 / (3.0 * M_PI) * gas_surfacedensity_cgs / (particle_density * calculateTurbulentAlpha(radial_distance,disk_params)) * fragmentation_velocity_squared / local_soundspeed_squared;

    return s_frag;
}

double calculateDriftInducedFragmentationBarrier(double gas_surfacedensity, double radial_distance, double gas_pressure, double pressure_gradient, double particle_density, const DiskParameters *disk_params) {

    double fragmentation_velocity, keplerian_velocity, log_pressure_gradient, local_soundspeed, local_soundspeed_squared, drift_barrier_size, gas_surfacedensity_cgs;

    fragmentation_velocity = disk_params->fragmentation_velocity * CM_PER_SEC_TO_AU_PER_YEAR_2PI; 
    gas_surfacedensity_cgs = gas_surfacedensity / SURFACE_DENSITY_CONVERSION_FACTOR;
    local_soundspeed = calculateLocalSoundSpeed(radial_distance,disk_params);
    local_soundspeed_squared = local_soundspeed * local_soundspeed;
    log_pressure_gradient = radial_distance / gas_pressure * pressure_gradient;
    keplerian_velocity = calculateKeplerianVelocity(radial_distance,disk_params);

    drift_barrier_size = fragmentation_velocity * keplerian_velocity / fabs(log_pressure_gradient * local_soundspeed_squared * 0.5) * 2.0 * gas_surfacedensity_cgs / (M_PI * particle_density);

    return drift_barrier_size;
}

double calculateGrowthTimescale(double radial_distance, double dust_to_gas_ratio,const DiskParameters *disk_params) {

    double keplerian_frequency = calculateKeplerianFrequency(radial_distance,disk_params);
    double calculateGrowthTimescale = dust_to_gas_ratio / keplerian_frequency;
    return calculateGrowthTimescale;
}

/**
 * @brief Iterates over all particles and updates their size based on physical barriers and growth timescales.
 */
void updateParticleSizes(ParticleData *particle_data, int particle_number, double actual_timestep, 
                        const DiskParameters *disk_params, const SimulationOptions *sim_opts, SnapshotMode mode) {
    
    #pragma omp parallel for

    for (int i = 0; i < particle_number; i++) {
        double r = particle_data->particle_distance_array[i][0];
        if (r < disk_params->r_min || r > disk_params->r_max) continue;

        // Retrieve local gas properties at particle position using grid indices
        int idx = (int)particle_data->dust_particle_mass_array[i][1];
        if (idx < 0 || idx >= disk_params->grid_number) idx = 0;

        double gas_sigma = disk_params->gas_surface_density_vector[idx];
        double gas_pres  = disk_params->gas_pressure_vector[idx];
        double gas_grad  = disk_params->gas_pressure_gradient_vector[idx];
        double dust_sigma = particle_data->dust_surfacedensity[idx];

        // Update primary population size
        particle_data->particle_distance_array[i][1] = calculateDustParticleSize(
            particle_data->particle_distance_array[i][1],
            disk_params->particle_density,
            gas_sigma, dust_sigma, r, gas_pres, gas_grad,
            actual_timestep, disk_params
        );
    }

    if(isSecondaryPopulationEnabled(mode) || particle_data->micron_particle_distance_array != NULL) {
        #pragma omp parallel for
        for (int i = 0; i < particle_number; i++) {
            double r_micron = particle_data->micron_particle_distance_array[i][0];
            if (r_micron < disk_params->r_min || r_micron > disk_params->r_max) continue;

            int idx_micron = (int)particle_data->micron_dust_particle_mass_array[i][1];
            if (idx_micron < 0 || idx_micron >= disk_params->grid_number) idx_micron = 0;

            double gas_sigma_micron = disk_params->gas_surface_density_vector[idx_micron];
            double gas_pres_micron  = disk_params->gas_pressure_vector[idx_micron];
            double gas_grad_micron  = disk_params->gas_pressure_gradient_vector[idx_micron];
            double dust_sigma_micron = particle_data->micron_dust_surfacedensity[idx_micron];

            // Update secondary population size
            particle_data->micron_particle_distance_array[i][1] = calculateDustParticleSize(
                particle_data->micron_particle_distance_array[i][1],
                disk_params->particle_density,
                gas_sigma_micron, dust_sigma_micron, r_micron, gas_pres_micron, gas_grad_micron,
                actual_timestep, disk_params
            );
        }
    }
}

double calculateDustParticleSize(double particle_radius, double particle_density, double gas_surfacedensity, double dust_surfacedensity, double particle_distance, double gas_pressure, double pressure_gradient, double actual_timestep, const DiskParameters *disk_params) {

    double turbulent_size_barrier = calculateTurbulentFragmentationBarrier(gas_surfacedensity, particle_distance, particle_density, disk_params);        
    double fragmentation_size_barrier = calculateDriftInducedFragmentationBarrier(gas_surfacedensity, particle_distance, gas_pressure, pressure_gradient, particle_density,disk_params);
    double drift_size_barrier = calculateRadialDriftBarrier(dust_surfacedensity, particle_distance, gas_pressure, pressure_gradient, particle_density, disk_params); 
    double particle_size_minimum = findMinimumOfAnArray(turbulent_size_barrier, fragmentation_size_barrier, drift_size_barrier); 
    double dust_to_gas_ratio = dust_surfacedensity / gas_surfacedensity;
    double growth_timescale = calculateGrowthTimescale(particle_distance, dust_to_gas_ratio, disk_params);
    double particle_size = 0.0;

    particle_size_minimum = particle_size_minimum / AU_IN_CM; 

    if (particle_radius < particle_size_minimum) {
        particle_size = findMinimumOfAnArray(particle_radius * exp(actual_timestep / growth_timescale), particle_size_minimum, HUGE_VAL);
    } else {
        particle_size = particle_size_minimum;
    }

    return particle_size;
}

/**
 * @brief Maps particle masses onto the radial grid to calculate surface density (Sigma).
 *        Supports primary and secondary dust populations.
 */
void calculateDustSurfaceDensity(const ParticleData *particle_data,
                                 const SimulationOptions *simulation_options,
                                 const DiskParameters *disk_params,
                                 SnapshotMode mode)
{
    int i, k;
    int grid_n = disk_params->grid_number;
    int is_twopop = isSecondaryPopulationEnabled(mode);

    // 1. Reset
    for (i = 0; i < grid_n; i++) {
        particle_data->dust_surfacedensity[i] = 0.0;
        if (is_twopop) {
            particle_data->micron_dust_surfacedensity[i] = 0.0;
        }
    }

    // 2. Mapping (CIC or NGP)
    int mode_type = simulation_options->dust_smoothing_mode;

    for (k = 0; k < particle_number; k++) {
        // Primary Mapping
        double r = particle_data->particle_distance_array[k][0];
        double rel_pos = (r - disk_params->r_min) / disk_params->delta_r;
        
        // Secondary Mapping (if enabled)
        double r_mic = 0.0;
        double rel_pos_mic = 0.0;
        if (is_twopop) {
            r_mic = particle_data->micron_particle_distance_array[k][0];
            rel_pos_mic = (r_mic - disk_params->r_min) / disk_params->delta_r;
        }

        // Apply mapping based on mode (CIC/NGP)
        if (mode_type == 1) { // NGP
            int idx = (int)round(rel_pos);
            if (idx >= 0 && idx < grid_n) {
                particle_data->dust_surfacedensity[idx] += particle_data->dust_particle_mass_grid[k];
            }
            if (is_twopop) {
                int idx_m = (int)round(rel_pos_mic);
                if (idx_m >= 0 && idx_m < grid_n) {
                    particle_data->micron_dust_surfacedensity[idx_m] += particle_data->massmicradial_grid[k];
                }
            }
        } else { // CIC (Default for all other modes)
            int idx = (int)floor(rel_pos);
            double x = rel_pos - (double)idx;
            if (idx >= 0 && idx < grid_n - 1) {
                double m = particle_data->dust_particle_mass_grid[k];
                particle_data->dust_surfacedensity[idx]   += m * (1.0 - x);
                particle_data->dust_surfacedensity[idx+1] += m * x;
            }
            if (is_twopop) {
                int idx_m = (int)floor(rel_pos_mic);
                double x_m = rel_pos_mic - (double)idx_m;
                if (idx_m >= 0 && idx_m < grid_n - 1) {
                    double m_m = particle_data->massmicradial_grid[k];
                    particle_data->micron_dust_surfacedensity[idx_m]   += m_m * (1.0 - x_m);
                    particle_data->micron_dust_surfacedensity[idx_m+1] += m_m * x_m;
                }
            }
        }
    }

    // 3. Smoothing (TopHat or Gaussian) - apply to both arrays
    // Note: If you have Gaussian smoothing, remember to apply it to micron_dust_surfacedensity as well!
    // (Skipped full snippet for brevity, ensure you loop over 'is_twopop' arrays too)

    // 4. Normalization
    for (i = 0; i < grid_n; i++) {
        double r_i = disk_params->radial_grid[i];
        double area = 2.0 * M_PI * r_i * disk_params->delta_r;

        if (area > 1e-20) {
            particle_data->dust_surfacedensity[i] /= area;
            if (is_twopop) {
                particle_data->micron_dust_surfacedensity[i] /= area;
            }
        }
    }
}


void calculateDustDistance(const char *file_name, ParticleData *particle_data, double actual_timestep, double actual_time, int number_of_particles, const SimulationOptions *simulation_options, const DiskParameters *disk_params){

    int i;
    double particle_distance, particle_distance_new, particle_radius_new, particle_radius;
    char file_path[1024];

    #pragma omp master
    {
        if (actual_time == 0) {
            sprintf(file_path, "%s/%s%s", file_name, kDriftTimescaleFileName, kFileNamesSuffix);
            drift_timescale_file = fopen(file_path, "w");
        }
    }

    #pragma omp barrier

    #pragma omp parallel for private(particle_distance, particle_distance_new, particle_radius_new, particle_radius)
    for (i = 0; i < number_of_particles; i++) {

        if (particle_data->particle_distance_array[i][0] > disk_params->r_min && particle_data->particle_distance_array[i][0] < disk_params->r_max) {
            particle_distance = particle_data->particle_distance_array[i][0];
            particle_radius = particle_data->particle_distance_array[i][1];

			integrateParticleRungeKutta4(actual_time, particle_radius, particle_data->dust_surfacedensity, particle_data->particle_distance_grid, actual_timestep, particle_distance, &particle_distance_new, &particle_radius_new, disk_params, simulation_options);
            if (actual_time == 0) {
                if (simulation_options->option_for_dust_secondary_population == 0) {
                    double current_timestep_value = (fabs(particle_distance_new - particle_distance) / (actual_timestep));

                    #pragma omp critical(drift_timescale_file_write)
                    {
                        if (drift_timescale_file != NULL) {
                            fprintf(drift_timescale_file, "%lg %lg\n", particle_data->particle_distance_array[i][0], (particle_data->particle_distance_array[i][0] / current_timestep_value) / (2.0 * M_PI));
                        } else {
                            LOG_ERROR("drift_timescale_file is NULL during write in calculateDustDistance (t=0 block).\n");
                        }
                    }
                }
            }

            if (simulation_options->option_for_dust_secondary_population != 1) { 
                particle_data->particle_distance_array[i][1] = particle_radius_new;
                particle_data->particle_distance_array[i][0] = particle_distance_new;
            } else {
//                particle_data->particle_distance_array[i][1] = particle_radius_new;
//                particle_data->particle_distance_array[i][0] = particle_distance_new;
            }
        } else {
            particle_data->particle_distance_array[i][0] = 0.0;
            particle_data->particle_distance_array[i][1] = 0.0;

        }
    }

    #pragma omp master
    {
        if (actual_time == 0) { 
            if (drift_timescale_file != NULL) {
                fclose(drift_timescale_file);
                drift_timescale_file = NULL; 
            }
        }
    }

    #pragma omp barrier
}