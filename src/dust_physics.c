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

/**
 * @brief Kiszámítja a porszemcsék új pozícióját és méretét strukturált rácson.
 * Ez a függvény összefogja a sűrűség-leképezést, a driftet és a növekedést.
 */
void calculateDustDistanceStructured(const char *file_name, StructuredParticleData *data, 
                                     double actual_timestep, double actual_time, 
                                     const SimulationOptions *simulation_options, 
                                     DiskParameters *disk_params) {




    
    char file_path[1024];
    size_t n_r = data->n_r;
    size_t n_z = data->n_z;

    // --- 1. LEKÉPEZÉS (MAPPING): Részecskék tömegének vetítése az Euler-rácsra ---
    // Kiválasztjuk, melyik algoritmus frissítse a disk_params->dust_surface_density_euler tömböt.
    // Ezt soros módban futtatjuk, mert a belső függvényeid már tartalmazzák a logikát.
    switch(simulation_options->dust_mapping_mode) {
        case 0: 
            updateDustSurfaceDensityEulerian(data, disk_params); 
            break;
        case 1: 
            updateDustSurfaceDensityEulerianCIC(data, disk_params); 
            break;
        case 2: 
            updateDustSurfaceDensityEulerianTSC(data, disk_params); 
            break;
        case 3: 
            updateDustSurfaceDensitySmart(data, disk_params); // TSC + Lyukkitöltés
            break;
        default:
            updateDustSurfaceDensitySmart(data, disk_params);
    }

    // --- 2. Timescale fájl megnyitása (Csak t=0-nál és csak a master szálon) ---
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
    #pragma omp barrier

    // --- 3. MOZGATÁS ÉS NÖVEKEDÉS (Párhuzamosított Lagrange-i léptetés) ---
    #pragma omp parallel
    {
        double r_new, radius_new;

        #pragma omp for collapse(2) schedule(static)
        for (size_t i = 0; i < n_r; i++) {
            for (size_t j = 0; j < n_z; j++) {
                DustParticle *p = &data->particles[i][j];

                // Csak a korong határain belüli aktív részecskéket számoljuk
                if (p->r_au > disk_params->r_min && p->r_au < disk_params->r_max) {

                    // RK4 integráció: Az előbb frissített porsűrűség rácsot használjuk
                    integrateParticleRungeKutta4Structured(
                        actual_time, 
                        p->radius, 
                        actual_timestep, 
                        p->r_au, 
                        &r_new, 
                        &radius_new, 
                        disk_params, 
                        simulation_options
                    );


                    // --- Timescale írás (t=0, midplane, primer populáció) ---
                    if (actual_time == 0 && j == n_z / 2 && simulation_options->option_for_dust_secondary_population == 0) {
                        double v_drift = fabs(r_new - p->r_au) / actual_timestep;
                        if (v_drift > 0) {
                            #pragma omp critical(drift_timescale_file_write)
                            {
                                if (drift_timescale_file != NULL) {
                                    // Sugár [AU] és Drift időskála [év]
                                    fprintf(drift_timescale_file, "%lg %lg\n", p->r_au, (p->r_au / v_drift) / (2.0 * M_PI));
                                }
                            }
                        }
                    }

                    // Lagrange-i frissítés
                    p->r_au = r_new;
                    p->radius = radius_new;

                } else {
                    // Határon kívül: deaktiválás és tömegvesztés
                    p->r_au = 0.0;
                    p->radius = 0.0;
                    p->mass_g = 0.0;
                }
            }
        }
    } // Parallel régió vége

    // --- 4. Fájl lezárása (Csak t=0-nál a master szálon) ---
    #pragma omp master
    {
        if (actual_time == 0 && drift_timescale_file != NULL) {
            fclose(drift_timescale_file);
            drift_timescale_file = NULL;
        }
    }
    #pragma omp barrier
}