// src/simulation_core.c
#include <stdio.h>    
#include <stdlib.h>   
#include <stdint.h>
#include <math.h>     
#include <string.h> 
#include <hdf5.h>  
#include <omp.h>
#include <time.h>
#include "simulation_core.h"
#include "config.h"       
#include "ascii_output.h"     
#include "disk_model.h"   
#include "dust_physics.h" 
#include "utils.h"        
#include "particle_data.h" 
#include "gas_physics.h"
#include "boundary_conditions.h"
#include "integrator.h"
#include "hdf5_output.h"
#include "print_terminal.h"
#include "logger.h"


void calculate1DDustDrift(double particle_radius, double pressure_gradient, double gas_surface_density, double gas_velocity, double radial_distance, double *drift_velocity, const DiskParameters *disk_params) {
    double local_pressure, local_pressure_scaleheight, local_pressure_gradient, stokes_number, local_soundspeed;
      
    stokes_number = calculateStokesNumber(particle_radius,gas_surface_density,disk_params);
    local_pressure_scaleheight = calculatePressureScaleHeight(radial_distance,disk_params);   
    local_pressure = calculateGasPressure(gas_surface_density,radial_distance,disk_params);
    local_pressure_gradient = pressure_gradient;
    local_soundspeed = calculateLocalSoundSpeed(radial_distance,disk_params); 

    *drift_velocity = gas_velocity / (1. + stokes_number * stokes_number) + stokes_number / (1. + stokes_number * stokes_number) * local_pressure_scaleheight / local_pressure * local_pressure_gradient * local_soundspeed;
}

double getMaximumDriftVelocity(const ParticleData *particle_data, int particle_number, const DiskParameters *disk_params, SnapshotMode mode) {
    double max_v = 0.0;
    double v_drift = 0.0;

    // 1. Iterate over the primary (cm-sized) dust population
    for (int i = 0; i < particle_number; i++) {
        double r = particle_data->particle_distance_array[i][0];
        if (r < disk_params->r_min || r > disk_params->r_max) continue;

        int idx = (int)particle_data->dust_particle_mass_array[i][1]; // Grid index
        // Boundary safety check to prevent segmentation faults
        if (idx < 0) {
            idx = 0;
        } else if (idx >= disk_params->grid_number) {
            idx = disk_params->grid_number - 1;
        }
        
        calculate1DDustDrift(
            particle_data->particle_distance_array[i][1], // Particle radius
            disk_params->gas_pressure_gradient_vector[idx],
            disk_params->gas_surface_density_vector[idx],
            disk_params->gas_velocity_vector[idx],
            r,
            &v_drift,
            disk_params
        );

        if (fabs(v_drift) > max_v) max_v = fabs(v_drift);
    }

    // 2. Iterate over the secondary (micron-sized) dust population if enabled
    if (isSecondaryPopulationEnabled(mode)) {
        for (int i = 0; i < particle_number; i++) {
            double r_mic = particle_data->micron_particle_distance_array[i][0];
            if (r_mic < disk_params->r_min || r_mic > disk_params->r_max) continue;

            int idx_mic = (int)particle_data->micron_dust_particle_mass_array[i][1]; // Micron grid index
            
            calculate1DDustDrift(
                particle_data->micron_particle_distance_array[i][1], // Micron radius
                disk_params->gas_pressure_gradient_vector[idx_mic],
                disk_params->gas_surface_density_vector[idx_mic],
                disk_params->gas_velocity_vector[idx_mic],
                r_mic,
                &v_drift,
                disk_params
            );

            if (fabs(v_drift) > max_v) max_v = fabs(v_drift);
        }
    }

    return max_v;
}

double calculateTimeStep(const DiskParameters *disk_params, double max_drift_v) {
    double max_viscosity = -1e10;
    double min_photo_dt = 1e10; 
    int i;
    int N = disk_params->grid_number;
    double WIND_CRIT = 1e-20;
    
    // 1. Gas/Viscous and Photoevaporation loops (already existing)
    for(i = 1; i <= N; i++) {
        double current_nu = calculateKinematicViscosity(disk_params->radial_grid[i], disk_params);
        if(current_nu > max_viscosity) {
            max_viscosity = current_nu;
        }

        if (disk_params->enable_photoevaporation && disk_params->sigma_dot_photoevap != NULL) {
            double sigma = disk_params->gas_surface_density_vector[i];
            double sigma_dot = disk_params->sigma_dot_photoevap[i];
            
            if (sigma_dot > 1e-20 && sigma > WIND_CRIT) {
                double cell_photo_dt = 0.15 * (sigma / sigma_dot);
                if (cell_photo_dt < min_photo_dt) {
                    min_photo_dt = cell_photo_dt;
                }
            }
        }
    }

    double safety_factor = 0.35;
    double viscous_dt = safety_factor * (disk_params->delta_r * disk_params->delta_r) / (2.0 * max_viscosity);
    
    // 2. Dust Drift CFL condition (The "Brake" you wanted)
    double dust_dt = 1e10; // Default: effectively infinity
    if (max_drift_v > 1e-15) { // Avoid division by zero
        double dust_safety_factor = 0.4;
        dust_dt = dust_safety_factor * disk_params->delta_r / max_drift_v;
    }
    
    // 3. Select the most restrictive timestep
    double time_step = viscous_dt;
    if (disk_params->enable_photoevaporation && min_photo_dt < time_step) {
        time_step = min_photo_dt;
    }
    if (dust_dt < time_step) {
        time_step = dust_dt;
    }

    return time_step;
}

static void handleSnapshot(double current_time_years, double *output_time, const SimulationOptions *sim_opts, OutputFiles *output_files, 
                           char *dens_name, char *dust_name, char *dust_name2, char *size_name, char *size_name2, SnapshotMode mode) {
    buildSnapshotFilenames(dens_name, dust_name, dust_name2, size_name, size_name2, sim_opts, (int)(*output_time),mode); 
    if(isGasEvolutionEnabled(mode)) {
        output_files->surface_file = openSnapshotFile(dens_name, FILE_TYPE_GAS_DENSITY, current_time_years);
    }
    output_files->dust_file = openSnapshotFile(dust_name, FILE_TYPE_DUST_DENSITY, current_time_years);
    if (isSecondaryPopulationEnabled(mode)) {
        output_files->micron_dust_file = openSnapshotFile(dust_name2, FILE_TYPE_DUST_MICRON_DENSITY, current_time_years);
    }
    output_files->size_file = openSnapshotFile(size_name, FILE_TYPE_PARTICLE_SIZE, current_time_years);
    if (isSecondaryPopulationEnabled(mode)) {
        output_files->size_micron_file = openSnapshotFile(size_name2, FILE_TYPE_MICRON_PARTICLE_SIZE, current_time_years);
    }
}

static void snapshotInitAtT0(double t, double current_time_years, ParticleData *particle_data, DiskParameters *disk_params, const SimulationOptions *sim_opts, int particle_number, SnapshotMode mode) {
    if (current_time_years == 0) {
        updateParticleGridIndices(particle_data, t, particle_number, disk_params, isSecondaryPopulationEnabled(mode));

        if (isDustEnabled(mode)) {
            calculateDustSurfaceDensity(particle_data, sim_opts, disk_params, mode);
        }
    }
}

static void snapshotPrintGas(DiskParameters *disk_params, OutputFiles *output_files, int is_gas_evolution_enabled) {
    
    if(is_gas_evolution_enabled) {
        printGasSurfaceDensityPressurePressureDerivateFile(disk_params, output_files);
    }
}

static void snapshotPrintDust(int output_time, ParticleData *particle_data, DiskParameters *disk_params, const SimulationOptions *sim_opts, OutputFiles *output_files, char *size_name, char *size_name2, SnapshotMode mode) {
    if (isDustEnabled(mode)) {
        printDustParticleSizeFile(size_name, size_name2, output_time, particle_data, disk_params, output_files, mode);
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

static void snapshotDustSurfacedensity(double output_time, ParticleData *particle_data, DiskParameters *disk_params, const SimulationOptions *sim_opts, OutputFiles *output_files, SnapshotMode mode) {
    calculateDustSurfaceDensity(particle_data, sim_opts, disk_params, mode);

    printDustSurfaceDensityPressurePressureDerivateFile(
        disk_params->radial_grid,
        disk_params->radial_grid,
        particle_data->dust_surfacedensity,
        particle_data->micron_dust_surfacedensity,       // Második populáció (micron)
        disk_params, 
        sim_opts, 
        output_files, 
        (int)output_time
    );

}

static void snapshotAdvance(double *output_time, const SimulationOptions *sim_opts) {
    *output_time += (double)(sim_opts->maximum_simulation_time / sim_opts->output_frequency);
}

static int isSnapshotDue(double current_time_years, double output_time, double deltat,const SimulationOptions *sim_opts) {
    double interval = sim_opts->maximum_simulation_time / sim_opts->output_frequency;
    int periodic_output_time = (fmod(current_time_years, interval) < deltat);
    int initial_output_time  = (current_time_years == 0.0);
    int output_time_sync = ((output_time - current_time_years) < deltat);

    return (periodic_output_time || initial_output_time) && output_time_sync;
}

static void handleSnapshotASCII(double t, double current_time_years, double *output_time, ParticleData *particle_data, int particle_number, DiskParameters *disk_params, const SimulationOptions *sim_opts,
                                OutputFiles *output_files, char *dens_name, char *dust_name, char *dust_name2, char *size_name, char *size_name2, SnapshotMode mode) {

    handleSnapshot(current_time_years, output_time, sim_opts, output_files, dens_name, dust_name, dust_name2, size_name, size_name2, mode);
    snapshotInitAtT0(t, current_time_years, particle_data, disk_params, sim_opts, particle_number, mode);
    snapshotPrintGas(disk_params, output_files, isGasEvolutionEnabled(mode));
    snapshotPrintDust((int)(*output_time), particle_data, disk_params, sim_opts, output_files, size_name, size_name2, mode);
    snapshotResetMasses(particle_data, particle_number, sim_opts);
    snapshotDustSurfacedensity(*output_time, particle_data, disk_params, sim_opts, output_files, mode);

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

    printTrapMassEvolution(*output_time, num_found, current_traps, output_files);
    snapshotAdvance(output_time, sim_opts);
    closeSnapshotFiles(output_files, mode);
}

static void handleSnapshotHDF5(double output_time, const SimulationOptions *sim_opts, OutputFiles *output_files, DiskParameters *disk_params, ParticleData *particle_data) {
    char *filename = NULL;
    asprintf(&filename, "%s/%s/%s_%08d%s",
             sim_opts->output_dir_name, kLogFilesDirectory, kSnapshotOutputFileNamePrefix, (int)(output_time), kFileNamesHDF5Suffix);

    if (initHDF5File(filename, output_files) != 0) {
        LOG_ERROR("Could not initialize HDF5 file %s", filename);
        return;
    }

    hid_t file_id = (hid_t)(intptr_t)output_files->hdf5_file;

    // --- FIX: If in gas-only mode: create a dummy particle data structure ---
    if (particle_data == NULL) {
        ParticleData dummy_particle_data;
        memset(&dummy_particle_data, 0, sizeof(ParticleData));

        writeHDF5SnapshotToFile(output_time, file_id, sim_opts, disk_params, &dummy_particle_data);
    } else {
        writeHDF5SnapshotToFile(output_time, file_id, sim_opts, disk_params, particle_data);
    }

    closeHDF5File(output_files);
    free(filename);
}

static void simulateDustDriftStep(double *t, double deltat, double *output_time, ParticleData *particle_data, 
                                  int particle_number, DiskParameters *disk_params, SimulationOptions *sim_opts,
                                  OutputFiles *output_files, char *dens_name, char *dust_name, char *dust_name2, 
                                  char *size_name, char *size_name2, SnapshotMode mode) {

    double min_radius, max_radius;
    double current_time_years = *t / (2.0 * M_PI);


    computeParticleRadiusRange(particle_data, particle_number, sim_opts->option_for_dust_secondary_population, &min_radius, &max_radius);

    if (isSnapshotDue(current_time_years, *output_time, deltat, sim_opts)) {
        if (sim_opts->output_format == OUTPUT_ASCII) {
            handleSnapshotASCII(*t, current_time_years, output_time, particle_data, particle_number, disk_params, sim_opts,
                               output_files, dens_name, dust_name, dust_name2, size_name, size_name2, mode);
        } else {

            handleSnapshotHDF5(*output_time, sim_opts, output_files, disk_params, particle_data);

            PressureTrap current_traps[MAX_TRAPS];
            int num_found = identifyPressureTraps(disk_params, current_traps, MAX_TRAPS);
            double trap_pos[MAX_TRAPS] = {0.0};  

            for (int i = 0; i < num_found && i < MAX_TRAPS; i++) {
                trap_pos[i] = current_traps[i].radial_position;
            }

            double primary_mass[MAX_TRAPS] = {0};
            double secondary_mass[MAX_TRAPS] = {0};
            double total_mass[MAX_TRAPS] = {0};

            for (int k = 0; k < num_found; k++) {
                if (current_traps[k].radial_position > 0.0) {
                    double local_H = calculatePressureScaleHeight(current_traps[k].radial_position, disk_params);
                    current_traps[k].inner_boundary = current_traps[k].radial_position - local_H;
                    current_traps[k].outer_boundary = current_traps[k].radial_position + local_H;

                    calculateMassInSpecificTrap(&current_traps[k], particle_data, particle_number, sim_opts);

                    primary_mass[k]   = current_traps[k].primary_dust_mass;
                    secondary_mass[k] = current_traps[k].secondary_dust_mass;
                    total_mass[k]     = primary_mass[k] + secondary_mass[k];
                }
            }

            appendMassTimeSeries(*output_time, primary_mass, secondary_mass, total_mass, trap_pos);
            snapshotAdvance(output_time, sim_opts);
        }
    }

    if (isGasEvolutionEnabled(mode)) {
        refreshGasSurfaceDensityPressurePressureGradient(sim_opts, disk_params);
    }


    // 1. Update grid indices based on current positions
    updateParticleGridIndices(particle_data, *t, particle_number, disk_params, isSecondaryPopulationEnabled(mode));

    // 2. Refresh dust density distribution
    if (isDustEnabled(mode)) {
        calculateDustSurfaceDensity(particle_data, sim_opts, disk_params, mode);
    }

    // 3. Update particle sizes if growth is enabled
    if (isDustGrowthEnabled(mode)) {
        updateParticleSizes(particle_data, particle_number, deltat, disk_params, sim_opts, mode);
    }

    // 4. Evolve particle positions
    calculateDustDistance(sim_opts->output_dir_name, particle_data, deltat, *t, particle_number, sim_opts, disk_params);
    if(isSecondaryPopulationEnabled(mode)) {
        applyCoagulationMassTransfer(particle_data, disk_params, sim_opts, deltat);
    }

    *t += deltat;
}

static void simulateGasOnlyStep(double *t, double deltat, double *output_time, DiskParameters *disk_params,
                                SimulationOptions *sim_opts, OutputFiles *output_files, char *dens_name)
{
    double current_time_years = *t / (2.0 * M_PI);

    // 1. ELŐSZÖR alkalmazzuk a peremfeltételt, hogy a szellemmezők tiszta adatok legyenek!
    sim_opts->current_bc_target = 0; // Sigma
    applyBoundaryConditions(disk_params->gas_surface_density_vector, disk_params, sim_opts);

    // 2. UTÁNA jöhet a gázállapot és a gradiensek frissítése
    refreshGasSurfaceDensityPressurePressureGradient(sim_opts, disk_params);

    applyBoundaryConditions(disk_params->gas_surface_density_vector, disk_params, sim_opts);

    // --- Ellenőrzés CSAK UTÁNA ---
    for (int i = 0; i <= disk_params->grid_number + 1; i++) {
        if (isnan(disk_params->gas_surface_density_vector[i]) || 
            isnan(disk_params->gas_pressure_vector[i]) || 
            isnan(disk_params->gas_pressure_gradient_vector[i])) {
            LOG_ERROR("NAN DETECTED at index %d! (Sigma: %e, P: %e, dP/dr: %e)", 
                      i, disk_params->gas_surface_density_vector[i], 
                         disk_params->gas_pressure_vector[i], 
                         disk_params->gas_pressure_gradient_vector[i]);
            exit(1);
        }
    }

    if (isSnapshotDue(current_time_years, *output_time, deltat, sim_opts)) {
        if (sim_opts->output_format == OUTPUT_ASCII) {
            asprintf(&dens_name, "%s/%s/%s_%08d%s", sim_opts->output_dir_name, kLogFilesDirectory, kGasDensityProfileFilePrefix, (int)(*output_time), kFileNamesSuffix);
            output_files->surface_file = fopen(dens_name, "w");

            if (!output_files->surface_file) {
                fprintf(stderr, "\nERROR: Could not open %s for writing.\n", dens_name);
            } else {
                HeaderData gas_header_data = { .current_time = current_time_years, .is_initial_data = (current_time_years == 0.0) };
                printFileHeader(output_files->surface_file, FILE_TYPE_GAS_DENSITY, &gas_header_data);
                printGasSurfaceDensityPressurePressureDerivateFile(disk_params, output_files);
                fclose(output_files->surface_file);
                output_files->surface_file = NULL;
            }
        } else {    
            handleSnapshotHDF5(*output_time, sim_opts, output_files, disk_params, NULL);
        }
        snapshotAdvance(output_time, sim_opts);
    }

    *t += deltat;
}

void timeIntegrationForTheSystem(SnapshotMode mode, DiskParameters *disk_params, SimulationOptions *sim_opts, OutputFiles *output_files) {
    ParticleData particle_data;
    HeaderData header_data_for_files;     

    double output_time = 0.; 
    double initial_disk_mass = 0.0;
    for (int i = 1; i <= disk_params->grid_number; i++) {
        initial_disk_mass += 2.0 * M_PI * disk_params->radial_grid[i] * disk_params->delta_r * disk_params->gas_surface_density_vector[i];
    }
    LOG_INFO("Numerically integrated initial disk mass is: %.5e M_Sun\n", initial_disk_mass);
    LOG_INFO("YAML provided initial disk mass is: %.5e M_Sun\n", disk_params->total_disk_mass);
    LOG_INFO("Initial calculated disk mass is %.5e M_Sun\n", initial_disk_mass);

    if (disk_params == NULL) {
        LOG_ERROR("disk_params_ptr is NULL!\n");
        exit(1);
    }

    memset(&particle_data, 0, sizeof(ParticleData));

    if (mode >= SnapshotDrift) {
        particle_number = calculateNumbersOfParticles(sim_opts->dust_input_filename);
    } else {
        LOG_INFO("Gas-only mode active. Particle count set to 0.\n");
        particle_number = 0;
    }

    if (particle_number > 0 && allocateParticleData(&particle_data, particle_number, (int)sim_opts->option_for_dust_secondary_population) != 0) {
        LOG_ERROR("Failed to allocate particle data. Exiting.\n");
        exit(EXIT_FAILURE);
    }

    if (mode >= SnapshotDrift) {
        if (setupInitialOutputFiles(output_files, sim_opts, disk_params, &header_data_for_files) != 0) {
            LOG_ERROR("Failed to set up initial output files. Exiting.\n");
            exit(EXIT_FAILURE);
        }
        loadDustParticlesFromFile(&particle_data, sim_opts->dust_input_filename);
    }

    char dens_name[MAX_PATH_LEN] = "";
    char dust_name[MAX_PATH_LEN] = "";
    char dust_name2[MAX_PATH_LEN] = "";
    char size_name[MAX_PATH_LEN] = "";
    char size_name2[MAX_PATH_LEN] = "";
    double t = 0.0;
    double t_integration_in_internal_units = sim_opts->maximum_simulation_time * 2.0 * M_PI;

    if (sim_opts->output_format == OUTPUT_HDF5) {
        char *ts_filename = NULL;
        asprintf(&ts_filename, "%s/%s/%s%s", sim_opts->output_dir_name, kLogFilesDirectory, kTimeSeriesForMassAccumulatinFileName, kFileNamesHDF5Suffix);
        initializeMassTimeSeries(ts_filename);
        free(ts_filename);
    }


// --- MAIN TEMPORAL INTEGRATION LOOP ---
    int step_counter = 0;
    double current_disk_mass = initial_disk_mass;
    double target_termination_mass = initial_disk_mass * 0.001;

    static double last_snapshot_time = 0.0; // static, hogy megjegyezze két lépés között
    double snapshot_interval = sim_opts->maximum_simulation_time / sim_opts->output_frequency;

    do {
        static double dt_old = 0.0;

        double max_drift_velocity = getMaximumDriftVelocity(&particle_data, particle_number, disk_params, mode);        
        double dt_new = calculateTimeStep(disk_params, max_drift_velocity) / 5.0;

        if (dt_old == 0.0) dt_old = dt_new;

        /* smoothing */
        double deltat = 0.7 * dt_old + 0.3 * dt_new;

        dt_old = deltat;

        ((SimulationOptions *)sim_opts)->user_defined_time_step = deltat;

        double current_time_years = t / (2.0 * M_PI);
        int snapshot_done = isSnapshotDue(current_time_years, output_time, deltat, sim_opts);

        if (isDustEnabled(mode)) {
            simulateDustDriftStep(&t, deltat, &output_time, &particle_data, particle_number, disk_params, sim_opts, output_files, dens_name, dust_name, dust_name2, size_name, size_name2, mode);
        } else { 
            simulateGasOnlyStep(&t, deltat, &output_time, disk_params, sim_opts, output_files, dens_name);
        }    
        step_counter++;

        // --- Calcualte mass in each step for termination ---
        current_disk_mass = 0.0;
        int has_nan = 0;
        int fully_evaporated_cells = 0;

        for (int i = 1; i <= disk_params->grid_number; i++) {
            // Don't allow negative surface densities, set them to zero if they occur
            if (disk_params->gas_surface_density_vector[i] < 0.0) {
                disk_params->gas_surface_density_vector[i] = 0.0;
            }

            if (isnan(disk_params->gas_surface_density_vector[i])) {
                has_nan = 1;
            }

            if (disk_params->gas_surface_density_vector[i] <= disk_params->density_floor) {
                fully_evaporated_cells++;
            }

            double cell_mass = 2.0 * M_PI * disk_params->radial_grid[i] * disk_params->delta_r * disk_params->gas_surface_density_vector[i];
            current_disk_mass += cell_mass;
        }

        // --- Diagnostic output ---
        if (step_counter % 10 == 0 || snapshot_done) {
            const char *mode_str = snapshotModeToString(mode);
            printStatus(step_counter, deltat, current_time_years, t, output_time, mode_str, snapshot_done,
                        current_disk_mass, target_termination_mass, initial_disk_mass, last_snapshot_time, snapshot_interval, sim_opts);
        }

        // --- Automatic termination ---
        // If the mass is lower than the threshold, or if NaN is detected, or the has_nan fag is triggered
        if (current_disk_mass < target_termination_mass || isnan(current_disk_mass) || has_nan || (fully_evaporated_cells >= (int)(disk_params->grid_number * 0.95))) { 
            LOG_INFO("Early exit triggered!\n");
            if (isnan(current_disk_mass) || has_nan) {
                LOG_ERROR("Reason: Numerical instability detected (NaN). Disk is effectively evaporated.\n");
            } else if (fully_evaporated_cells >= (int)(disk_params->grid_number * 0.95)) {
                LOG_ERROR("Reason: 95%% of the radial grid cells are fully evaporated (below %.1e M_Sun).\n", disk_params->density_floor);
            } else {
                LOG_ERROR("Reason: Disk mass reached the 0.1%% threshold: %.4e M_Sun (Target: %.4e M_Sun)\n", 
                        current_disk_mass, target_termination_mass);
            }
            break; 
        }

        if (snapshot_done) {
            last_snapshot_time = current_time_years;
        }


    } while (t <= t_integration_in_internal_units);

    if (sim_opts->output_format == OUTPUT_HDF5) {
        closeMassTimeSeries();
    }

    LOG_INFO("Main simulation loop finished.\n");
    cleanupSimulationResources(&particle_data, output_files);
    LOG_INFO("Cleanup completed.\n");
}