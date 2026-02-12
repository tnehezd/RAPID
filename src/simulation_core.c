// src/simulation_core.c
#include <stdio.h>    
#include <stdlib.h>   
#include <stdint.h>
#include <math.h>     
#include <string.h> 
#include <hdf5.h>  
#include <omp.h>
#include <time.h>
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
#include "hdf5_output.h"


void calculate1DDustDrift(double particle_radius, double pressure_gradient, double gas_surface_density, double gas_velocity, double radial_distance, double *drift_velocity, const DiskParameters *disk_params) {

    double local_pressure, local_pressure_scaleheight, local_pressure_gradient, stokes_number, local_soundspeed;
      
    stokes_number = calculateStokesNumber(particle_radius,gas_surface_density,disk_params);
    local_pressure_scaleheight = calculatePressureScaleHeight(radial_distance,disk_params);   
    local_pressure = calculateGasPressure(gas_surface_density,radial_distance,disk_params);
    local_pressure_gradient = pressure_gradient;
    local_soundspeed = calculateLocalSoundSpeed(radial_distance,disk_params); 

    *drift_velocity = gas_velocity / (1. + stokes_number * stokes_number) + stokes_number / (1. + stokes_number * stokes_number) * local_pressure_scaleheight / local_pressure * local_pressure_gradient * local_soundspeed;
}

double calculateTimeStep(const DiskParameters *disk_params) {

    double max_diffusion_coefficient, time_step;
    int i;
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


static void handleSnapshot(double actual_time, double current_time_years, double *snapshot, const SimulationOptions *sim_opts, OutputFiles *output_files, 
                           char *dens_name, char *dust_name, char *dust_name2, char *size_name){

    fprintf(stderr, "\n--- Simulation Time: %.2e years (Internal time: %.2e, snapshot: %.2e) ---\n", current_time_years, actual_time, *snapshot);
    buildSnapshotFilenames(dens_name, dust_name, dust_name2, size_name, sim_opts, (int)(*snapshot)); 
    output_files->surface_file = openSnapshotFile(dens_name, FILE_TYPE_GAS_DENSITY, current_time_years);
    output_files->dust_file = openSnapshotFile(dust_name, FILE_TYPE_DUST_DENSITY, current_time_years);
    if (sim_opts->option_for_dust_secondary_population == 1.0) {
        output_files->micron_dust_file = openSnapshotFile(dust_name2, FILE_TYPE_DUST_MICRON_DENSITY, current_time_years);
    }
    output_files->size_file = openSnapshotFile(size_name, FILE_TYPE_PARTICLE_SIZE, current_time_years);
}

static void snapshotInitAtT0(double t, double current_time_years, ParticleData *particle_data, DiskParameters *disk_params, const SimulationOptions *sim_opts, int particle_number) {

    if (current_time_years == 0) {
        updateParticleGridIndices(particle_data, t, particle_number, disk_params);
        if (sim_opts->option_for_dust_secondary_population == 1) updateParticleGridIndices(particle_data, t, particle_number, disk_params);
            if (sim_opts->option_for_dust_growth == 1.) {
                calculateDustSurfaceDensity(particle_data, sim_opts, disk_params);
            }
        }
}

static void snapshotPrintGas(DiskParameters *disk_params, OutputFiles *output_files,const SimulationOptions *sim_opts) {

    if (sim_opts->option_for_evolution == 1 ) {
        printGasSurfaceDensityPressurePressureDerivateFile(disk_params, output_files);
    }

    
}

static void snapshotPrintDust(int snapshot, ParticleData *particle_data, DiskParameters *disk_params, const SimulationOptions *sim_opts, OutputFiles *output_files, char *size_name) {

    if (sim_opts->option_for_dust_drift == 1) {
        printDustParticleSizeFile(size_name, snapshot, particle_data->particle_distance_array, particle_data->micron_particle_distance_array, disk_params, sim_opts, output_files);
    }
}

static void snapshotResetMasses(ParticleData *particle_data, int particle_number, const SimulationOptions *sim_opts) {

    if (sim_opts->flag_for_deadzone == 1.0 && particle_number > 0) {
        for (int k = 0; k < particle_number; k++) {
            particle_data->dust_particle_mass_array[k][3] = 0.0;
            particle_data->dust_particle_mass_array[k][4] = 0.0;
        }
    }
}

static void snapshotDustSurfacedensity(double snapshot, ParticleData *particle_data, DiskParameters *disk_params, const SimulationOptions *sim_opts, OutputFiles *output_files) {

    if (sim_opts->option_for_dust_growth == 1.) {
        printDustSurfaceDensityPressurePressureDerivateFile(particle_data->particle_distance_grid, particle_data->micron_particle_distance_grid, particle_data->dust_surfacedensity, particle_data->micron_dust_surfacedensity, disk_params, sim_opts, output_files, (int)snapshot);
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

static void handleSnapshotASCII(double t, double current_time_years, double *snapshot, ParticleData *particle_data, int particle_number, DiskParameters *disk_params, const SimulationOptions *sim_opts,
                                OutputFiles *output_files, char *dens_name, char *dust_name, char *dust_name2, char *size_name) {

    handleSnapshot(t, current_time_years, snapshot, sim_opts, output_files,dens_name, dust_name, dust_name2, size_name);
    snapshotInitAtT0(t, current_time_years,particle_data, disk_params, sim_opts,particle_number);
    snapshotPrintGas(disk_params, output_files,sim_opts);
    snapshotPrintDust((int)(*snapshot), particle_data, disk_params,sim_opts, output_files, size_name);
    snapshotResetMasses(particle_data, particle_number, sim_opts);
    snapshotDustSurfacedensity(*snapshot, particle_data, disk_params, sim_opts,output_files);
    fprintf(stderr, "snapshot set to %lg\n", *snapshot);

    PressureTrap current_traps[3];
    int num_found = identifyPressureTraps(disk_params, current_traps, 3);

    for (int k = 0; k < num_found; k++) {
        if (current_traps[k].radial_position > 0) {
            double local_H = calculatePressureScaleHeight(current_traps[k].radial_position, disk_params);
            current_traps[k].inner_boundary = current_traps[k].radial_position - 1.0 * local_H;
            current_traps[k].outer_boundary = current_traps[k].radial_position + 1.0 * local_H;
            calculateMassInSpecificTrap(&current_traps[k], particle_data, particle_number, sim_opts);
        }
    }

    printTrapMassEvolution(*snapshot, num_found, current_traps, output_files);
    snapshotAdvance(snapshot, sim_opts);
    closeSnapshotFiles(output_files, sim_opts);

}

static void handleSnapshotHDF5(double *t, double current_time_years, double *snapshot,
                               const SimulationOptions *sim_opts, OutputFiles *output_files,
                               DiskParameters *disk_params, ParticleData *particle_data)
{
    char *filename = NULL;
    asprintf(&filename, "%s/%s/snapshot_%08d.h5",
             sim_opts->output_dir_name, kLogFilesDirectory, (int)(*snapshot));

    // Init file
    if (initHDF5File(filename, output_files) != 0) {
        fprintf(stderr, "ERROR: Could not initialize HDF5 file %s\n", filename);
        return;
    }

    // Write datasets to the already opened file
    if (particle_data != NULL) {
        hid_t file_id = (hid_t)(intptr_t)output_files->hdf5_file;
        writeHDF5SnapshotToFile(file_id, sim_opts, disk_params, particle_data);
    }

    // Close file
    closeHDF5File(output_files);

    (*snapshot) += (double)(sim_opts->maximum_simulation_time / sim_opts->output_frequency);
    fprintf(stderr, "HDF5 snapshot %d done at time %.2e years\n", (int)(*snapshot), current_time_years);
}



static void simulateDustDriftStep(double *t, double deltat, double *snapshot, ParticleData *particle_data, int particle_number, DiskParameters *disk_params, const SimulationOptions *sim_opts,
                                  OutputFiles *output_files, char *dens_name, char *dust_name, char *dust_name2, char *size_name) {

    double min_radius, max_radius;
    double current_time_years = *t / (2.0 * M_PI);

    computeParticleRadiusRange(particle_data, particle_number, sim_opts->option_for_dust_secondary_population, &min_radius, &max_radius);


    if (isSnapshotDue(current_time_years, *snapshot, deltat, sim_opts)) {

        if (sim_opts->output_format == OUTPUT_ASCII) {
            handleSnapshotASCII(*t, current_time_years, snapshot,particle_data,particle_number,disk_params,sim_opts,
                                output_files,dens_name,dust_name,dust_name2,size_name);
        } else {
            handleSnapshotHDF5(t, current_time_years, snapshot, sim_opts, output_files, disk_params, particle_data);
        }

    }

    if (sim_opts->option_for_evolution == 1.) {
        refreshGasSurfaceDensityPressurePressureGradient(sim_opts, disk_params);
    }

    updateParticleGridIndices(particle_data,*t, particle_number, disk_params);

    if (sim_opts->option_for_dust_secondary_population == 1) {
        updateParticleGridIndices(particle_data,*t, particle_number, disk_params);
    }

    if (sim_opts->option_for_dust_growth == 1.) {
        calculateDustSurfaceDensity(particle_data,sim_opts, disk_params);
    }

    calculateDustDistance(sim_opts->output_dir_name,particle_data, deltat, *t,particle_number, sim_opts, disk_params);

    if (sim_opts->option_for_dust_secondary_population == 1.) {
        calculateDustDistance(sim_opts->output_dir_name, particle_data,deltat, *t,particle_number, sim_opts, disk_params);
    }

    *t += deltat;
}



static void simulateGasOnlyStep(double *t,double deltat,double *snapshot,DiskParameters *disk_params,const SimulationOptions *sim_opts,OutputFiles *output_files,char *dens_name){

    double current_time_years = *t / (2.0 * M_PI);

    if (isSnapshotDue(current_time_years, *snapshot, deltat, sim_opts)) {

        if (sim_opts->output_format == OUTPUT_ASCII) {

            asprintf(&dens_name, "%s/%s/%s_%08d%s",sim_opts->output_dir_name,kLogFilesDirectory,kGasDensityProfileFilePrefix,(int)(*snapshot),kFileNamesSuffix);
            fprintf(stderr,
                    "\n--- Simulation Time: %.2e years (Internal time: %.2e, snapshot: %.2e) ---\n",
                    current_time_years, *t, *snapshot);

            output_files->surface_file = fopen(dens_name, "w");

            if (!output_files->surface_file) {
                fprintf(stderr, "ERROR: Could not open %s for writing.\n", dens_name);
            } else {
                HeaderData gas_header_data = {.current_time = current_time_years,.is_initial_data = (current_time_years == 0.0)};
                printFileHeader(output_files->surface_file,FILE_TYPE_GAS_DENSITY,&gas_header_data);
            }

            printGasSurfaceDensityPressurePressureDerivateFile(disk_params, output_files);

            if (output_files->surface_file) {
                fclose(output_files->surface_file);
                output_files->surface_file = NULL;
            } else {
            
                handleSnapshotHDF5(t, current_time_years, snapshot, sim_opts, output_files, disk_params, NULL);
            
            }
        }
        snapshotAdvance(snapshot, sim_opts);
        fprintf(stderr, "snapshot updated to %lg\n", *snapshot);
    }

    refreshGasSurfaceDensityPressurePressureGradient(sim_opts, disk_params);

    *t += deltat;
}


void timeIntegrationForTheSystem(SnapshotMode mode, DiskParameters *disk_params, const SimulationOptions *sim_opts, OutputFiles *output_files) {

    ParticleData particle_data;
    HeaderData header_data_for_files; 

    double snapshot = 0.; 

    if (disk_params == NULL) {
        fprintf(stderr, "ERROR [timeIntegrationForTheSystem]: disk_params_ptr is NULL!\n");
        exit(1);
    }

    memset(&particle_data, 0, sizeof(ParticleData));

    if (mode > 2) {
        particle_number = calculateNumbersOfParticles(sim_opts->dust_input_filename);
    } else {
        fprintf(stderr, "ERROR [timeIntegrationForTheSystem]: Particle drift is OFF. particle_number set to 0.\n");
        particle_number = 0;
    }

    if (particle_number > 0 && allocateParticleData(&particle_data, particle_number, (int)sim_opts->option_for_dust_secondary_population) != 0) {
        fprintf(stderr, "ERROR: Failed to allocate particle data. Exiting.\n");
        exit(EXIT_FAILURE);
    }

    if (mode>2) {
        if (setupInitialOutputFiles(output_files, sim_opts, disk_params, &header_data_for_files) != 0) {
            fprintf(stderr, "ERROR: Failed to set up initial output files. Exiting.\n");
            exit(EXIT_FAILURE);
        }
        loadDustParticlesFromFile(&particle_data, sim_opts->dust_input_filename);
    }

    int i;
    char dens_name[MAX_PATH_LEN] = "";
    char dust_name[MAX_PATH_LEN] = "";
    char dust_name2[MAX_PATH_LEN] = "";
    char size_name[MAX_PATH_LEN] = "";
    double t = 0.0;
    double t_integration_in_internal_units = sim_opts->maximum_simulation_time * 2.0 * M_PI;
    double deltat = calculateTimeStep(disk_params) / 5.0;

    if (sim_opts->user_defined_time_step > 0.0 && sim_opts->user_defined_time_step < deltat) {
        ((SimulationOptions *)sim_opts)->user_defined_time_step = deltat;
    } else {
        ((SimulationOptions *)sim_opts)->user_defined_time_step = deltat; 
    }

    if (sim_opts->option_for_dust_secondary_population == 0 && particle_number > 0) {
        for (i = 0; i < particle_number; i++) {
            particle_data.micron_particle_distance_array[i][0] = 0;
            particle_data.micron_particle_distance_array[i][1] = 0;
            particle_data.micron_dust_particle_mass_array[i][0] = 0;
            particle_data.micron_dust_particle_mass_array[i][1] = 0;
            particle_data.massmicradial_grid[i] = 0;
        }
    }

    do {
        if (mode > 1) {
            simulateDustDriftStep(&t, deltat, &snapshot, &particle_data, particle_number, disk_params, sim_opts, output_files, dens_name, dust_name, dust_name2, size_name );
        } else { 
            simulateGasOnlyStep(&t, deltat, &snapshot,disk_params, sim_opts, output_files,dens_name);
        }    
    } while (t <= t_integration_in_internal_units);

    fprintf(stderr,"\n\nDEBUG [timeIntegrationForTheSystem]: Main simulation loop finished (t > t_integration_in_internal_units).\n");
    cleanupSimulationResources(&particle_data, output_files);
    fprintf(stderr,"DEBUG [timeIntegrationForTheSystem]: Cleanup completed.\n");
}


