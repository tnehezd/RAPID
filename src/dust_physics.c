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

/*  Stokes number in the Epstein drag regime  */
double calculateStokesNumber(double particle_radius, double gas_surfacedensity, const DiskParameters *disk_params) { 

    return disk_params->particle_density_dimensionless * particle_radius * M_PI / (2.0 * gas_surfacedensity);
}

void calculateParticleMass(int number_of_particles, double (*dust_particle_mass_array)[5], int indii, int indio, int indoi, int indoo, double *massiout, double *massoout, const SimulationOptions *simulation_options) {

    double massitemp = 0.0;
    double massotemp = 0.0;
    int i;

    // simulation_options->flag_for_deadzone (ez helyettesíti az optdze-t): 1.0 = dinamikus DZE (flag-elt), 0.0 = fix DZE (nem flag-elt)
    if(simulation_options->flag_for_deadzone == 1.0) { // Dinamikus DZE: flag-ek használatával (dust_particle_mass_array[i][3] és [4])
        #pragma omp parallel for private(i) reduction(+:massitemp, massotemp)
        for (i = 0; i < number_of_particles; i++) {
            // A részecske aktuális grid indexe (dust_particle_mass_array[i][1]-ből)
            int current_r_index = (int)dust_particle_mass_array[i][1]; 

            // --- Belső DZE ---
            // Ellenőrizzük, hogy a részecske grid indexe a belső DZE tartományában van-e
            if ((current_r_index >= indii) && (current_r_index <= indio)) {
                if (dust_particle_mass_array[i][3] == 0.0) { // Belső DZE flag ellenőrzése
                    #pragma omp critical(inner_dze_update)
                    {
                        dust_particle_mass_array[i][3] = 1.0;
                        massitemp = massitemp + dust_particle_mass_array[i][0]; // Tömeg hozzáadása a [0] indexről
                    }
                }
            }

            // --- Külső DZE ---
            // Ellenőrizzük, hogy a részecske grid indexe a külső DZE tartományában van-e
            if ((current_r_index >= indoi) && (current_r_index <= indoo)) {
                if (dust_particle_mass_array[i][4] == 0.0) { // Külső DZE flag ellenőrzése (Flag[4])
                    #pragma omp critical(outer_dze_update)
                    {
                        dust_particle_mass_array[i][4] = 1.0;
                        massotemp = massotemp + dust_particle_mass_array[i][0]; // Tömeg hozzáadása a [0] indexről
                    }
                }
            }
        }
    } else { // Fix DZE (simulation_options->flag_for_deadzone == 0.0): Nincsenek flag-ek a tömeg felhalmozáshoz
        #pragma omp parallel for private(i) reduction(+:massitemp, massotemp)
        for (i = 0; i < number_of_particles; i++) {
            int current_r_index = (int)dust_particle_mass_array[i][1]; // A részecske grid indexe

            // --- Belső DZE (Fix) ---
            // Az eredeti kódod else ága nem tartalmazta a belső DZE gyűjtését. 
            // Ha szeretnéd, akkor ide kell tenni a logikát:
            if ((current_r_index >= indii) && (current_r_index <= indio)) {
                #pragma omp critical(inner_dze_update_fixed)
                {
                    massitemp = massitemp + dust_particle_mass_array[i][0]; 
                }
            }

            // --- Külső DZE (Fix) ---
            if ((current_r_index >= indoi) && (current_r_index <= indoo)) {
                #pragma omp critical(outer_dze_update_fixed)
                {
                    massotemp = massotemp + dust_particle_mass_array[i][0];
                }
            }
        }
    }

    *massiout = massitemp;
    *massoout = massotemp;
}


/*			BIRNSTIEL EL AL 2012			*/

//reprezentativ reszecske kezdeti meretenek meghatarozasa
// 1. radialis drift altal meghatarozott maximalis meret			--> kimenet cm-ben!
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

// 2. a kis skalaju turbulencia altal okozott fragmentacio szerinti maximalis meret	--> kimenet cm-ben!
double calculateTurbulentFragmentationBarrier(double gas_surfacedensity, double radial_distance, double particle_density, const DiskParameters *disk_params) {

    double s_frag, fragmentation_velocity, fragmentation_velocity_squared, gas_surfacedensity_cgs, local_soundspeed, local_soundspeed_squared;

    fragmentation_velocity = disk_params->fragmentation_velocity * CM_PER_SEC_TO_AU_PER_YEAR_2PI; /*	cm/sec --> AU / (yr/2pi)	*/
    fragmentation_velocity_squared = fragmentation_velocity * fragmentation_velocity;
    gas_surfacedensity_cgs = gas_surfacedensity / SURFACE_DENSITY_CONVERSION_FACTOR;
    local_soundspeed = calculateLocalSoundSpeed(radial_distance,disk_params); // / CM_PER_SEC_TO_AU_PER_YEAR_2PI; // Komment ki, ha a calculateLocalSoundSpeed már megfelelő mértékegységben van
    local_soundspeed_squared = local_soundspeed * local_soundspeed;

    s_frag = disk_params->fragmentation_factor * 2.0 / (3.0 * M_PI) * gas_surfacedensity_cgs / (particle_density * calculateTurbulentAlpha(radial_distance,disk_params)) * fragmentation_velocity_squared / local_soundspeed_squared;

    return s_frag;
}

// 3. radialis drift altal okozott fragmentacio szerinti maximalis meret		--> kimenet cm-ben!
double calculateDriftInducedFragmentationBarrier(double gas_surfacedensity, double radial_distance, double gas_pressure, double pressure_gradient, double particle_density, const DiskParameters *disk_params) {

    double fragmentation_velocity, keplerian_velocity, log_pressure_gradient, local_soundspeed, local_soundspeed_squared, drift_barrier_size, gas_surfacedensity_cgs;

    fragmentation_velocity = disk_params->fragmentation_velocity * CM_PER_SEC_TO_AU_PER_YEAR_2PI; /*	cm/sec --> AU / (yr/2pi)	*/
    gas_surfacedensity_cgs = gas_surfacedensity / SURFACE_DENSITY_CONVERSION_FACTOR;
    local_soundspeed = calculateLocalSoundSpeed(radial_distance,disk_params);
    local_soundspeed_squared = local_soundspeed * local_soundspeed;
    log_pressure_gradient = radial_distance / gas_pressure * pressure_gradient;
    keplerian_velocity = calculateKeplerianVelocity(radial_distance,disk_params);

    drift_barrier_size = fragmentation_velocity * keplerian_velocity / fabs(log_pressure_gradient * local_soundspeed_squared * 0.5) * 2.0 * gas_surfacedensity_cgs / (M_PI * particle_density);

    return drift_barrier_size;
}

/*	a reszecskek novekedesenek idoskalaja	*/
double calculateGrowthTimescale(double radial_distance, double dust_to_gas_ratio,const DiskParameters *disk_params) {

    double keplerian_frequency = calculateKeplerianFrequency(radial_distance,disk_params);
    double calculateGrowthTimescale = dust_to_gas_ratio / keplerian_frequency;
    return calculateGrowthTimescale;
}

/*	kiszamolja az adott helyen a reszecske meretet --> BIRNSTIEL CIKK	*/
double calculateDustParticleSize(double particle_radius, double particle_density, double gas_surfacedensity, double dust_surfacedensity, double particle_distance, double gas_pressure, double pressure_gradient, double actual_timestep, const DiskParameters *disk_params) {

    double turbulent_size_barrier = calculateTurbulentFragmentationBarrier(gas_surfacedensity, particle_distance, particle_density, disk_params);           // cm-ben
    double fragmentation_size_barrier = calculateDriftInducedFragmentationBarrier(gas_surfacedensity, particle_distance, gas_pressure, pressure_gradient, particle_density,disk_params); // cm-ben
    double drift_size_barrier = calculateRadialDriftBarrier(dust_surfacedensity, particle_distance, gas_pressure, pressure_gradient, particle_density, disk_params); // cm-ben
    double particle_size_minimum = findMinimumOfAnArray(turbulent_size_barrier, fragmentation_size_barrier, drift_size_barrier);         // cm-ben -- megadja, hogy a fenti ket reszecske korlatbol melyik ad kisebb meretet (az a reszecskenovekedes felso korlatja
    //	double dust_to_gas_ratio = gas_surfacedensity / 100.;
    double dust_to_gas_ratio = dust_surfacedensity / gas_surfacedensity; // A korábbi kódban fordítva volt, feltételezem, hogy dust_to_gas_ratio = (por sűrűség) / (gáz sűrűség)
    double growth_timescale = calculateGrowthTimescale(particle_distance, dust_to_gas_ratio, disk_params);
    double particle_size = 0.0;

    particle_size_minimum = particle_size_minimum / AU_IN_CM; // AU-ban

    /*	kiszamolja, hogy a fenti particle_size_minimum, vagy a novekedesi idoskalabol szarmazo meret korlatozza a reszecske meretet	*/
    if (particle_radius < particle_size_minimum) {
        particle_size = findMinimumOfAnArray(particle_radius * exp(actual_timestep / growth_timescale), particle_size_minimum, HUGE_VAL);
    } else { // prad >= particle_size_minimum
        particle_size = particle_size_minimum;
    }

    return particle_size;
}


void calculateDustSurfaceDensity(double max_param, double min_param, const ParticleData *particle_data, const SimulationOptions *simulation_options, const DiskParameters *disk_params) {

    // Suppress unused parameter warnings
    (void)max_param;
    (void)min_param;
    
    
    double grid_step = (disk_params->r_max - disk_params->r_min) / (particle_number - 1);
    int i;

    // A temp tömbök deklarálását érdemes a scope tetejére tenni
    double temporary_dust_surfacedensity[particle_number][3];
    double temporary_micron_dust_surfacedensity[particle_number][3];

    // Inicializálás, ha szükséges (bár a calculateDustSurfaceDensity valószínűleg felülírja)
    for(i=0; i<particle_number; i++){
        temporary_dust_surfacedensity[i][0] = 0.0; temporary_dust_surfacedensity[i][1] = 0.0; temporary_dust_surfacedensity[i][2] = 0.0;
        temporary_micron_dust_surfacedensity[i][0] = 0.0; temporary_micron_dust_surfacedensity[i][1] = 0.0; temporary_micron_dust_surfacedensity[i][2] = 0.0;
        particle_data->particle_distance_grid[i] = 0.0;
        particle_data->micron_particle_distance_grid[i] = 0.0;
        particle_data->dust_surfacedensity[i] = 0.0;
        particle_data->micron_dust_surfacedensity[i] = 0.0;
    }


    // calculateDustSurfaceDensity és mergeParticlesByRadius függvények hívásai:
    // Ezek valószínűleg szekvenciálisak, hacsak a függvények belsejében nincs OpenMP.
    // Ha ezek a függvények valamilyen globális állapotot módosítanak, akkor kritikusak.
    // Feltételezve, hogy a 'temporary_dust_surfacedensity' és 'temporary_micron_dust_surfacedensity' kizárólagosan a hívásaikban vannak feldolgozva,
    // és nem ütköznek más szálakkal globális adatokon keresztül.
    calculateDustSurfaceDensityFromRepresentativeMass(particle_data->particle_distance_array, particle_data->dust_particle_mass_grid, temporary_dust_surfacedensity, particle_number,disk_params);
    if (simulation_options->option_for_dust_secondary_population == 1.0) { // Használjunk double összehasonlítást
        calculateDustSurfaceDensityFromRepresentativeMass(particle_data->micron_particle_distance_array, particle_data->massmicradial_grid, temporary_micron_dust_surfacedensity, particle_number,disk_params);
    }

    mergeParticlesByRadius(temporary_dust_surfacedensity, grid_step, particle_number,disk_params);
    if (simulation_options->option_for_dust_secondary_population == 1.0) { // Használjunk double összehasonlítást
        mergeParticlesByRadius(temporary_micron_dust_surfacedensity, grid_step, particle_number,disk_params);
    }

    // Utolsó másoló ciklus: Ez is jól párhuzamosítható.
    #pragma omp parallel for private(i)
    for (i = 0; i < particle_number; i++) {
        particle_data->particle_distance_grid[i] = temporary_dust_surfacedensity[i][1];
        particle_data->dust_surfacedensity[i] = temporary_dust_surfacedensity[i][0];

        if (simulation_options->option_for_dust_secondary_population == 1.0) { // double összehasonlítás
            particle_data->micron_particle_distance_grid[i] = temporary_micron_dust_surfacedensity[i][1];
            particle_data->micron_dust_surfacedensity[i] = temporary_micron_dust_surfacedensity[i][0];
        }
    }
}


/*	Fuggveny a porszemcsek uj tavolsaganak elraktarozasara		*/
void calculateDustDistance(const char *file_name, ParticleData *particle_data, double actual_timestep, double actual_time, int number_of_particles, const SimulationOptions *simulation_options, const DiskParameters *disk_params){

    int i;
    double particle_distance, particle_distance_new, particle_radius_new, particle_radius;
    char file_path[1024];


    // Fájlkezelés t==0 esetén: ez valószínűleg egyszer történik meg a szimuláció elején,
    // még mielőtt az igazi párhuzamosítás elkezdődne a fő ciklusban.
    // Biztosítjuk, hogy csak egy szál nyissa meg/zárja be a fájlt.
    #pragma omp master
    {
        if (actual_time == 0) {
            sprintf(file_path, "%s/%s%s", file_name, kDriftTimescaleFileName, kFileNamesSuffix);
            drift_timescale_file = fopen(file_path, "w");
        }
    }
    // Szinkronizálás a master szál befejezéséig.
    #pragma omp barrier

    // **ITT A PÁRHUZAMOSÍTÁS!**
    // Az `i` ciklus független iterációkkal rendelkezik, minden szál a saját `radius[i]` elemen dolgozik.
    #pragma omp parallel for private(particle_distance, particle_distance_new, particle_radius_new, particle_radius)
    for (i = 0; i < number_of_particles; i++) {
        // Csak a r_min és r_max közötti részecskékkel foglalkozunk.
        // A 0.0-ra állítás kívül esik a párhuzamos részen, ha az if feltétel nem teljesül.
        if (particle_data->particle_distance_array[i][0] > disk_params->r_min && particle_data->particle_distance_array[i][0] < disk_params->r_max) {
            particle_distance = particle_data->particle_distance_array[i][0];
            particle_radius = particle_data->particle_distance_array[i][1];

			integrateParticleRungeKutta4(actual_time, particle_radius, particle_data->dust_surfacedensity, particle_data->particle_distance_grid, actual_timestep, particle_distance, &particle_distance_new, &particle_radius_new, disk_params, simulation_options);
            if (actual_time == 0) {
                if (simulation_options->option_for_dust_secondary_population == 0) {
                    double current_timestep_value = (fabs(particle_distance_new - particle_distance) / (actual_timestep));
                    // Azért kell a critical szekció, mert az drift_timescale_file fájlba írunk.
                    // Ez a critical szekció biztosítja, hogy egyszerre csak egy szál írjon a fájlba.

                    #pragma omp critical(drift_timescale_file_write)
                    {
                        // Ellenőrizzük, hogy a fájlmutató nem NULL
                        if (drift_timescale_file != NULL) {
                            fprintf(drift_timescale_file, "%lg %lg\n", particle_data->particle_distance_array[i][0], (particle_data->particle_distance_array[i][0] / current_timestep_value) / (2.0 * M_PI));
                        } else {
                            fprintf(stderr, "ERROR: drift_timescale_file is NULL during write in calculateDustDistance (t=0 block).\n");
                        }
                    }
                }
            }

            if (simulation_options->option_for_dust_secondary_population != 1) { // Ha növekedés engedélyezett vagy valami más mód
                particle_data->particle_distance_array[i][1] = particle_radius_new;
                particle_data->particle_distance_array[i][0] = particle_distance_new;
            } else { // simulation_options->option_for_dust_secondary_population == 1, csak drift
                particle_data->particle_distance_array[i][0] = particle_distance_new;
            }
        } else {
            // Ha a részecske r_min vagy r_max kívülre kerül, 0.0-ra állítjuk a pozícióját.
            // Ez a hozzárendelés iterációnként független, így párhuzamosítható.
            particle_data->particle_distance_array[i][0] = 0.0;
        }
    }

    // A fájl bezárása ismételten egy szál által kell, hogy történjen.
    #pragma omp master
    {
        if (actual_time == 0) { // Csak akkor zárjuk be, ha meg is nyitottuk a t=0 blokkban
            if (drift_timescale_file != NULL) {
                fclose(drift_timescale_file);
                drift_timescale_file = NULL; // Fontos, hogy NULL-ra állítsuk, miután bezártuk.
            }
        }
    }
    // Szinkronizálás a master szál befejezéséig.
    #pragma omp barrier
}