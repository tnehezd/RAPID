
// src/simulation_core.c

// Standard C Library Includes
#include <stdio.h>    // For printf, fopen, fclose, fscanf, snprintf, sprintf
#include <stdlib.h>   // For exit, EXIT_FAILURE, EXIT_SUCCESS, system
#include <math.h>     // For M_PI, fmod, HUGE_VAL (and pow if used by other functions)
#include <string.h>   // For snprintf, sprintf

#include <omp.h>
#include <time.h>


// Your Project Header Includes
#include "config.h"       
#include "io_utils.h"     
#include "disk_model.h"   
#include "dust_physics.h" 
#include "utils.h"        
#include "simulation_core.h"
#include "particle_data.h" // Új include
#include "gas_physics.h"
#include "boundary_conditions.h"
#include "integrator.h"


/*	Kiszamolja az 1D-s driftet	*/
/*  	dr/dt = St/(1+St*St)*H(r)/r*dlnP/dlnr*cs = St/(1+St*St) * (H/r) * (r/P) * (dP/dr) * cs		*/
void calculate1DDustDrift(double particle_radius, double pressure_gradient, double gas_surface_density, double gas_velocity, double radial_distance, double *drift_velocity, const DiskParameters *disk_params) {

    double local_pressure, local_pressure_scaleheight, local_pressure_gradient, stokes_number, local_soundspeed;
      
    stokes_number = calculateStokesNumber(particle_radius,gas_surface_density,disk_params);
    local_pressure_scaleheight = calculatePressureScaleHeight(radial_distance,disk_params);   
    local_pressure = calculateGasPressure(gas_surface_density,radial_distance,disk_params);
    local_pressure_gradient = pressure_gradient;
    local_soundspeed = calculateLocalSoundSpeed(radial_distance,disk_params); 

    *drift_velocity = gas_velocity / (1. + stokes_number * stokes_number) + stokes_number / (1. + stokes_number * stokes_number) * local_pressure_scaleheight / local_pressure * local_pressure_gradient * local_soundspeed;	/* bearamlas sebessege: Birnstiel PHD	*/
}



double calculateTimeStep(const DiskParameters *disk_params) { // Add const here too
    double max_diffusion_coefficient, time_step;
    int i;

// IT WOULD BE BETTER TO CALCULATE FOR BOTH PARTS OF THE FUNCTION!!!!

    max_diffusion_coefficient = -10000.0;
    
    for(i = 0; i < disk_params->grid_number; i++) {
        if(ftcsSecondDerivativeCoefficient(disk_params->radial_grid[i], disk_params) > max_diffusion_coefficient) {
            max_diffusion_coefficient = ftcsSecondDerivativeCoefficient(disk_params->radial_grid[i], disk_params);
        }
    }
    time_step = disk_params->delta_r * disk_params->delta_r / (2.0 * max_diffusion_coefficient);
    fprintf(stderr," Actual calculateTimeStep: delta_r = %.2e, time_step = %.2e\n", disk_params->delta_r, time_step);

    return time_step;
}


static void handleSnapshot(
    double t,
    double current_time_years,
    double *snapshot,
    double deltat,
    ParticleData *particle_data,
    DiskParameters *disk_params,
    const SimulationOptions *sim_opts,
    OutputFiles *output_files,
    double *masstempiin, double *massmtempiin,
    double *masstempoin, double *massmtempoin,
    double *masstempiout, double *massmtempiout,
    double *masstempoout, double *massmtempoout,
    double *tavin, double *tavout,
    double min_radius, double max_radius,
    char *dens_name,
    char *dust_name,
    char *dust_name2,
    char *size_name
){

    fprintf(stderr,
        "\n--- Simulation Time: %.2e years (Internal time: %.2e, snapshot: %.2e) ---\n",
        current_time_years, t, *snapshot);

    buildSnapshotFilenames(dens_name, dust_name, dust_name2, size_name, sim_opts, (int)(*snapshot)); 


    // Fájlok megnyitása és fejlécek írása
    output_files->surface_file = openSnapshotFile(dens_name, FILE_TYPE_GAS_DENSITY, current_time_years);
    output_files->dust_file = openSnapshotFile(dust_name, FILE_TYPE_DUST_MOTION, current_time_years);

    if (sim_opts->option_for_dust_secondary_population == 1.) {
        output_files->micron_dust_file = openSnapshotFile(dust_name2, FILE_TYPE_MICRON_MOTION, current_time_years);
    }


}

static void snapshotInitAtT0(double t, double current_time_years, ParticleData *particle_data, DiskParameters *disk_params, const SimulationOptions *sim_opts, int particle_number, double min_radius, double max_radius) {
    // Eredeti t==0 logika
    if (current_time_years == 0) {
        updateParticleGridIndices(particle_data->radius,  particle_data->partmassind,  particle_data->massvec, t, particle_number, disk_params);
        if (sim_opts->option_for_dust_secondary_population == 1) updateParticleGridIndices( particle_data->radiusmicr,  particle_data->partmassmicrind,  particle_data->massmicradial_grid, t, particle_number, disk_params);

            if (sim_opts->option_for_dust_growth == 1.) {
                calculateDustSurfaceDensity(max_radius, min_radius,  particle_data->radius,  particle_data->radiusmicr,  particle_data->sigmad,  particle_data->sigmadm,  particle_data->massvec,  particle_data->massmicradial_grid,  particle_data->rdvec,  particle_data->rmicvec, sim_opts, disk_params);
            }
        }
}

static void snapshotPrintGas(double current_time_years, DiskParameters *disk_params, OutputFiles *output_files,const SimulationOptions *sim_opts, double snapshot) {

    if (sim_opts->option_for_evolution == 1 || current_time_years == 0) {
        if(snapshot != 0) printGasSurfaceDensityPressurePressureDerivateFile(disk_params, output_files);
    }

    
}

static void snapshotPrintDust(int snapshot, ParticleData *particle_data, DiskParameters *disk_params, const SimulationOptions *sim_opts, OutputFiles *output_files, char *size_name) {
    if (sim_opts->option_for_dust_drift == 1) {
        printDustParticleSizeFile(size_name, snapshot, particle_data->radius, particle_data->radiusmicr, disk_params, sim_opts, output_files);
    }
}

static void snapshotResetMasses(ParticleData *particle_data, int particle_number, const SimulationOptions *sim_opts, double *masstempiout, double *massmtempiout, double *masstempoout, double *massmtempoout) {
    // 1) mass accumulation nullázása
    *masstempiout = 0.0;
    *massmtempiout = 0.0;
    *masstempoout = 0.0;
    *massmtempoout = 0.0;

    // 2) deadzone reset
    if (sim_opts->flag_for_deadzone == 1.0 && particle_number > 0) {
        for (int k = 0; k < particle_number; k++) {
            particle_data->partmassind[k][3] = 0.0;
            particle_data->partmassind[k][4] = 0.0;
        }
    }
}

static void snapshotMassGrowthAndSigma(double t, double snapshot, ParticleData *particle_data, DiskParameters *disk_params, const SimulationOptions *sim_opts, OutputFiles *output_files,
                                       double *masstempiin, double *massmtempiin, double *masstempoin, double *massmtempoin, double *masstempiout, double *massmtempiout,
                                       double *masstempoout, double *massmtempoout, double *tavin, double *tavout, double min_radius, double max_radius) {
    // 1) mass growth print
    printMassGrowthAtDZEFile(snapshot,
                             particle_data->partmassind,
                             particle_data->partmassmicrind,
                             t,
                             *masstempiin, *masstempoin,
                             *massmtempiin, *massmtempoin,
                             masstempiout, masstempoout,
                             massmtempiout, massmtempoout,
                             tavin, tavout,
                             disk_params, sim_opts, output_files);

    // 2) input mass update
    *masstempiin  = *masstempiout;
    *massmtempiin = *massmtempiout;
    *masstempoin  = *masstempoout;
    *massmtempoin = *massmtempoout;

    // 3) dust surface density print
    if (sim_opts->option_for_dust_growth == 1.) {
        printDustSurfaceDensityPressurePressureDerivateFile(particle_data->rdvec, particle_data->rmicvec, particle_data->sigmad, particle_data->sigmadm, disk_params, sim_opts, output_files, (int)snapshot);
    }
}

static void snapshotAdvance(double *snapshot, const SimulationOptions *sim_opts) {
    *snapshot += (double)(sim_opts->maximum_simulation_time / sim_opts->output_frequency);
}

static int isSnapshotDue(double current_time_years, double snapshot, double deltat,const SimulationOptions *sim_opts) {
    double interval = sim_opts->maximum_simulation_time / sim_opts->output_frequency;

    int periodic_snapshot = (fmod(current_time_years, interval) < deltat);
    int initial_snapshot  = (current_time_years == 0.0);
    int snapshot_sync = ((snapshot - current_time_years) < deltat);

    return (periodic_snapshot || initial_snapshot) && snapshot_sync;
}



void timeIntegrationForTheSystem(DiskParameters *disk_params, const SimulationOptions *sim_opts, OutputFiles *output_files) {
    ParticleData particle_data;
    HeaderData header_data_for_files; // Később inicializáljuk a setupInitialOutputFiles-ban


    double snapshot = 0.; // Években mért "pillanatfelvétel" időzítő

    if (disk_params == NULL) {
        fprintf(stderr, "ERROR [timeIntegrationForTheSystem]: disk_params_ptr is NULL!\n");
        exit(1); // Program leállítása, ha kritikus hiba van
    }

    // --- Inicializálási szakasz ---
    if (sim_opts->option_for_dust_drift == 1.) {
        particle_number = calculateNumbersOfParticles(sim_opts->dust_input_filename);
    } else {
        fprintf(stderr, "ERROR [timeIntegrationForTheSystem]: Particle drift is OFF. particle_number set to 0.\n");
        particle_number = 0;
    }

    if (particle_number > 0 && allocateParticleData(&particle_data, particle_number, (int)sim_opts->option_for_dust_secondary_population) != 0) {
        fprintf(stderr, "ERROR: Failed to allocate particle data. Exiting.\n");
        exit(EXIT_FAILURE);
    }

    // Fájl inicializálás a meglévő io_utils függvény hívásával
    if (sim_opts->option_for_dust_drift == 1.) {
        if (setupInitialOutputFiles(output_files, sim_opts, disk_params, &header_data_for_files) != 0) {
            fprintf(stderr, "ERROR: Failed to set up initial output files. Exiting.\n");
            exit(EXIT_FAILURE);
        }
        // loadDustParticlesFromFile hívása a részecskeadatok beolvasására
        loadDustParticlesFromFile(particle_data.radius, particle_data.radiusmicr, particle_data.massvec, particle_data.massmicradial_grid, sim_opts->dust_input_filename);
    }

    int i; // Hagyjuk meg ezt a ciklusváltozót a C89 kompatibilitás kedvéért, ha szükséges
    
    // Ideiglenes puffer a fájlneveknek a ciklusban
    char dens_name[MAX_PATH_LEN] = "";
    char dust_name[MAX_PATH_LEN] = "";
    char dust_name2[MAX_PATH_LEN] = "";
    char size_name[MAX_PATH_LEN] = "";

    double t = 0.0;
    double t_integration_in_internal_units = sim_opts->maximum_simulation_time * 2.0 * M_PI;
    double deltat = calculateTimeStep(disk_params) / 5.0;

    // user_defined_time_step felülbírálása, ha a felhasználó megadott kisebb értéket
    if (sim_opts->user_defined_time_step > 0.0 && sim_opts->user_defined_time_step < deltat) {
        ((SimulationOptions *)sim_opts)->user_defined_time_step = deltat; // Az eredeti kód hibásan a deltat-t vette át, ha sim_opts->user_defined_time_step nagyobb volt
    } else {
        ((SimulationOptions *)sim_opts)->user_defined_time_step = deltat; // Az eredeti kód szerint, ha nem kisebb, akkor deltat a user_defined_time_step
    }

    // Mass accumulation változók
    double masstempiin = 0, massmtempiin = 0, masstempoin = 0, massmtempoin = 0;
    double masstempiout = 0, massmtempiout = 0, masstempoout = 0, massmtempoout = 0;
    double tavin = 0, tavout = 0; // Távolságok a printMassGrowthAtDZEFile-hoz


    if (sim_opts->option_for_dust_secondary_population == 0 && particle_number > 0) {
        for (i = 0; i < particle_number; i++) {
            particle_data.radiusmicr[i][0] = 0;
            particle_data.radiusmicr[i][1] = 0;
            particle_data.partmassmicrind[i][0] = 0;
            particle_data.partmassmicrind[i][1] = 0;
            particle_data.massmicradial_grid[i] = 0;
        }
    }

    // --- Fő szimulációs ciklus ---
    do {
        if (sim_opts->option_for_dust_drift == 1.) {

            double min_radius, max_radius;

            computeParticleRadiusRange(&particle_data,particle_number,sim_opts->option_for_dust_secondary_population,&min_radius,&max_radius);

            // --- Kimeneti adatok (pillanatfelvétel) kezelése ---
            double current_time_years = t / (2.0 * M_PI);
            if (isSnapshotDue(current_time_years, snapshot, deltat, sim_opts)) {

                handleSnapshot(t, current_time_years, &snapshot, deltat, &particle_data, disk_params, sim_opts, output_files, &masstempiin, &massmtempiin, &masstempoin, &massmtempoin, &masstempiout, &massmtempiout, &masstempoout, &massmtempoout,&tavin, &tavout, min_radius, max_radius, dens_name, dust_name, dust_name2, size_name);
                snapshotInitAtT0(t, current_time_years, &particle_data, disk_params, sim_opts, particle_number, min_radius, max_radius);
                snapshotPrintGas(current_time_years, disk_params, output_files, sim_opts, snapshot);
                snapshotPrintDust((int)snapshot, &particle_data, disk_params, sim_opts, output_files, size_name);
                snapshotResetMasses(&particle_data, particle_number, sim_opts, &masstempiout, &massmtempiout, &masstempoout, &massmtempoout);
                snapshotMassGrowthAndSigma(t, snapshot,&particle_data, disk_params, sim_opts,output_files,&masstempiin, &massmtempiin,&masstempoin, &massmtempoin,&masstempiout, &massmtempiout,&masstempoout, &massmtempoout,&tavin, &tavout,min_radius, max_radius);

                fprintf(stderr,"snapshot set to %lg\n",snapshot);
                snapshotAdvance(&snapshot, sim_opts);
                closeSnapshotFiles(output_files, dens_name, dust_name, dust_name2, sim_opts);
            }

            // Gas option_for_evolutionution
            if (sim_opts->option_for_evolution == 1.) {
                refreshGasSurfaceDensityPressurePressureGradient(sim_opts, disk_params);
            }

            // Count masses and get sigma_d for the next step (always done)
            updateParticleGridIndices(particle_data.radius, particle_data.partmassind, particle_data.massvec, t, particle_number, disk_params);
            if (sim_opts->option_for_dust_secondary_population == 1) updateParticleGridIndices(particle_data.radiusmicr, particle_data.partmassmicrind, particle_data.massmicradial_grid, t, particle_number, disk_params);

            if (sim_opts->option_for_dust_growth == 1.) {
                calculateDustSurfaceDensity(max_radius, min_radius, particle_data.radius, particle_data.radiusmicr, particle_data.sigmad, particle_data.sigmadm, particle_data.massvec, particle_data.massmicradial_grid, particle_data.rdvec, particle_data.rmicvec, sim_opts, disk_params);
            }

            // Get radii for next step
            int optsize = 0;
            calculateDustDistance(sim_opts->output_dir_name, optsize, particle_data.radius, particle_data.sigmad, particle_data.rdvec, deltat, t, particle_number, sim_opts, disk_params);

            if (sim_opts->option_for_dust_secondary_population == 1.) {
                optsize = 1;
                calculateDustDistance(sim_opts->output_dir_name, optsize, particle_data.radiusmicr, particle_data.sigmadm, particle_data.rmicvec, deltat, t, particle_number, sim_opts, disk_params);

            }

            t = t + deltat;

            // Kilépési feltétel a drift == 1 ágon
            // Fontos: absolute_maximum_radius és absolute_minimum_radius frissül az előző szakaszban, azt használjuk itt.
//            if (!(absolute_maximum_radius >= disk_params->r_min && absolute_minimum_radius != absolute_maximum_radius)) {
//                fprintf(stderr,"DEBUG [timeIntegrationForTheSystem]: Simulation termination condition met (absolute_maximum_radius < r_min or absolute_minimum_radius == absolute_maximum_radius) at time: %lg.\n",snapshot);

//                goto cleanup; // Ugrás a tisztításra
//            }

        } else { // sim_opts->option_for_dust_drift == 0. (Gas-only simulation)
            double current_time_years = t / (2.0 * M_PI);

            if ((fmod(current_time_years, (sim_opts->maximum_simulation_time / sim_opts->output_frequency)) < deltat || current_time_years == 0) && snapshot - current_time_years < deltat) {
                fprintf(stderr,"\n--- Simulation Time: %.2e years (Internal time: %.2e, snapshot: %.2e) ---\n", current_time_years, t, snapshot);

                fprintf(stderr,"DEBUG [timeIntegrationForTheSystem]: Outputting data for gas-only simulation at time %.2e. snapshot=%.2e\n", current_time_years, snapshot);
                snprintf(dens_name, MAX_PATH_LEN, "%s/%s/%s_%08d.dat", sim_opts->output_dir_name, kLogFilesDirectory, kGasDensityProfileFilePrefix, (int)snapshot);
                fprintf(stderr, "DEBUG [timeIntegrationForTheSystem]: Outputting %s_%08d.dat to %s.\n", kGasDensityProfileFilePrefix, (int)snapshot, dens_name);

                output_files->surface_file = fopen(dens_name, "w");
                if (output_files->surface_file == NULL) {
                    fprintf(stderr, "ERROR: Could not open %s for writing in gas-only branch.\n", dens_name);
                } else {
                    fprintf(stderr, "DEBUG [timeIntegrationForTheSystem]: Opened %s for writing in gas-only branch.\n", dens_name);
                    HeaderData gas_header_data = {.current_time = current_time_years, .is_initial_data = (current_time_years == 0.0)};
                    printFileHeader(output_files->surface_file, FILE_TYPE_GAS_DENSITY, &gas_header_data);
                }

                printGasSurfaceDensityPressurePressureDerivateFile(disk_params, output_files);


                if (output_files->surface_file != NULL) {
                    fclose(output_files->surface_file);
                    output_files->surface_file = NULL;
                    fprintf(stderr, "DEBUG [timeIntegrationForTheSystem]: Closed %s in gas-only branch.\n", dens_name);
                }

                snapshot = snapshot + (double)(sim_opts->maximum_simulation_time / sim_opts->output_frequency);
                fprintf(stderr,"DEBUG [timeIntegrationForTheSystem]: Updated snapshot to %.2e.\n", snapshot);
            }

            fprintf(stderr,"DEBUG [timeIntegrationForTheSystem]: Calling refreshGasSurfaceDensityPressurePressureGradient for gas-only evolution.\n");
            refreshGasSurfaceDensityPressurePressureGradient(sim_opts, disk_params);
            
            t = t + deltat;

        }


    } while (t <= t_integration_in_internal_units);

    fprintf(stderr,"\n\nDEBUG [timeIntegrationForTheSystem]: Main simulation loop finished (t > t_integration_in_internal_units).\n");

cleanup:
    // --- Tisztítási szakasz ---
    cleanupSimulationResources(&particle_data, output_files, sim_opts);
    fprintf(stderr,"DEBUG [timeIntegrationForTheSystem]: Cleanup completed.\n");
}

