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


void calculateDustSurfaceDensity(const ParticleData *particle_data, const SimulationOptions *simulation_options, const DiskParameters *disk_params) {

        
    double grid_step = (disk_params->r_max - disk_params->r_min) / (particle_number - 1);
    int i;
    double temporary_dust_surfacedensity[particle_number][3];
    double temporary_micron_dust_surfacedensity[particle_number][3];

    for(i=0; i<particle_number; i++){
        temporary_dust_surfacedensity[i][0] = 0.0; temporary_dust_surfacedensity[i][1] = 0.0; temporary_dust_surfacedensity[i][2] = 0.0;
        temporary_micron_dust_surfacedensity[i][0] = 0.0; temporary_micron_dust_surfacedensity[i][1] = 0.0; temporary_micron_dust_surfacedensity[i][2] = 0.0;
        particle_data->particle_distance_grid[i] = 0.0;
        particle_data->micron_particle_distance_grid[i] = 0.0;
        particle_data->dust_surfacedensity[i] = 0.0;
        particle_data->micron_dust_surfacedensity[i] = 0.0;
    }

    calculateDustSurfaceDensityFromRepresentativeMass(particle_data->particle_distance_array, particle_data->dust_particle_mass_grid, temporary_dust_surfacedensity, particle_number,disk_params);
    if (simulation_options->option_for_dust_secondary_population == 1.0) { 
        calculateDustSurfaceDensityFromRepresentativeMass(particle_data->micron_particle_distance_array, particle_data->massmicradial_grid, temporary_micron_dust_surfacedensity, particle_number,disk_params);
    }

    mergeParticlesByRadius(temporary_dust_surfacedensity, grid_step, particle_number,disk_params);
    if (simulation_options->option_for_dust_secondary_population == 1.0) { 
        mergeParticlesByRadius(temporary_micron_dust_surfacedensity, grid_step, particle_number,disk_params);
    }

    #pragma omp parallel for private(i)
    for (i = 0; i < particle_number; i++) {
        particle_data->particle_distance_grid[i] = temporary_dust_surfacedensity[i][1];
        particle_data->dust_surfacedensity[i] = temporary_dust_surfacedensity[i][0];

        if (simulation_options->option_for_dust_secondary_population == 1.0) { 
            particle_data->micron_particle_distance_grid[i] = temporary_micron_dust_surfacedensity[i][1];
            particle_data->micron_dust_surfacedensity[i] = temporary_micron_dust_surfacedensity[i][0];
        }
    }
}



void calculateDustDistanceStructured(const char *file_name, StructuredParticleData *data, double actual_timestep, double actual_time, const SimulationOptions *simulation_options, const DiskParameters *disk_params) {

    char file_path[1024];
    size_t n_r = data->n_r;
    size_t n_z = data->n_z;

    // --- 1. Timescale fájl megnyitása (Csak t=0-nál és csak a master szálon) ---
    #pragma omp master
    {
        if (actual_time == 0) {
            sprintf(file_path, "%s/%s%s", file_name, kDriftTimescaleFileName, kFileNamesSuffix);
            drift_timescale_file = fopen(file_path, "w");
            if (drift_timescale_file == NULL) {
                fprintf(stderr, "ERROR: Could not open drift timescale file: %s\n", file_path);
            }
        }
    }
    // Megvárjuk, amíg a master végez a fájlművelettel
    #pragma omp barrier

    // --- 2. Párhuzamosított Lagrange-i léptetés ---
    #pragma omp parallel
    {
        double r_new, radius_new;

        #pragma omp for collapse(2) schedule(static)
        for (size_t i = 0; i < n_r; i++) {
            for (size_t j = 0; j < n_z; j++) {
                DustParticle *p = &data->particles[i][j];

                // Csak a korong határain belüli részecskéket számoljuk
                if (p->r_au > disk_params->r_min && p->r_au < disk_params->r_max) {

                    // RK4 integráció a fix Euler-rács (disk_params->radial_grid) használatával
                    integrateParticleRungeKutta4(
                        actual_time,
                        p->radius,
                        disk_params->dust_surface_density_euler,
                        disk_params->radial_grid, // Fix hivatkozási rács
                        actual_timestep,
                        p->r_au,
                        &r_new,
                        &radius_new,
                        disk_params,
                        simulation_options
                    );

                    // --- Timescale írás: Csak t=0-nál, a midplane rétegből (j=n_z/2) ---
                    if (actual_time == 0 && j == n_z / 2 && simulation_options->option_for_dust_secondary_population == 0) {
                        double v_drift = fabs(r_new - p->r_au) / actual_timestep;
                        if (v_drift > 0) {
                            #pragma omp critical(drift_timescale_file_write)
                            {
                                if (drift_timescale_file != NULL) {
                                    // Drift időskála években kifejezve
                                    fprintf(drift_timescale_file, "%lg %lg\n", p->r_au, (p->r_au / v_drift) / (2.0 * M_PI));
                                }
                            }
                        }
                    }

                    // Fizikai állapot frissítése (Lagrange-i követés)
                    p->r_au = r_new;
                    p->radius = radius_new;

                } else {
                    // Ha elhagyta a tartományt, deaktiváljuk
                    p->r_au = 0.0;
                    p->radius = 0.0;
                    p->mass_g = 0.0;
                }
            }
        }
    } // Parallel régió vége

    // --- 3. Fájl lezárása (Csak t=0-nál a master szálon) ---
    #pragma omp master
    {
        if (actual_time == 0 && drift_timescale_file != NULL) {
            fclose(drift_timescale_file);
            drift_timescale_file = NULL;
        }
    }
    // Biztosítjuk, hogy a fájl lezáruljon, mielőtt visszatérünk a fő ciklusba
    #pragma omp barrier
}