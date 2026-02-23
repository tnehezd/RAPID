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
    #include "vertical_settling.h"



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


    static void handleSnapshot(double actual_time, double current_time_years, double *output_time, const SimulationOptions *sim_opts, OutputFiles *output_files, 
                               char *dens_name, char *dust_name, char *dust_name2, char *size_name){

        fprintf(stderr, "\n--- Simulation Time: %.2e years (Internal time: %.2e, output_time: %.2e in ASCII mode) ---\n", current_time_years, actual_time, *output_time);
        buildSnapshotFilenames(dens_name, dust_name, dust_name2, size_name, sim_opts, (int)(*output_time)); 
        output_files->surface_file = openSnapshotFile(dens_name, FILE_TYPE_GAS_DENSITY, current_time_years);
        output_files->dust_file = openSnapshotFile(dust_name, FILE_TYPE_DUST_DENSITY, current_time_years);
        if (sim_opts->option_for_dust_secondary_population == 1.0) {
            output_files->micron_dust_file = openSnapshotFile(dust_name2, FILE_TYPE_DUST_MICRON_DENSITY, current_time_years);
        }
        output_files->size_file = openSnapshotFile(size_name, FILE_TYPE_PARTICLE_SIZE, current_time_years);
    }

    static void snapshotInitAtT0(double current_time_years, StructuredParticleData *structured_particle_data, DiskParameters *disk_params, const SimulationOptions *sim_opts) {

        if (current_time_years == 0) {
            updateStructuredParticleGridIndices(structured_particle_data, disk_params);
            if (sim_opts->option_for_dust_secondary_population == 1) {
                updateStructuredParticleGridIndices(structured_particle_data, disk_params);
                if (sim_opts->option_for_dust_growth == 1.) {
                    updateDustSurfaceDensitySmart(structured_particle_data, disk_params);
                }
            }
        }
    }

    static void snapshotPrintGas(DiskParameters *disk_params, OutputFiles *output_files,const SimulationOptions *sim_opts) {

        if (sim_opts->option_for_evolution == 1 ) {
            printGasSurfaceDensityPressurePressureDerivateFile(disk_params, output_files);
        }

        
    }

    static void snapshotPrintDust(int output_time, StructuredParticleData *structured_particle_data, DiskParameters *disk_params, const SimulationOptions *sim_opts, OutputFiles *output_files, char *size_name) {

        if (sim_opts->option_for_dust_drift == 1) {
            printDustParticleSizeFileStructured(size_name, output_time, structured_particle_data, disk_params, sim_opts, output_files);
        }
    }

    static void snapshotDustSurfacedensity(double output_time, StructuredParticleData *data, DiskParameters *disk_params, const SimulationOptions *sim_opts, OutputFiles *output_files) {

        if (sim_opts->option_for_dust_growth == 1.0) {
            
            // --- 1. LEKÉPEZÉS (MAPPING) ---
            // Frissítjük a disk_params->dust_surface_density_euler-t az új struktúrából
            switch(sim_opts->dust_mapping_mode) {
                case 0:  updateDustSurfaceDensityEulerian(data, disk_params);    break;
                case 1:  updateDustSurfaceDensityEulerianCIC(data, disk_params); break;
                case 2:  updateDustSurfaceDensityEulerianTSC(data, disk_params); break;
                case 3:  updateDustSurfaceDensitySmart(data, disk_params);       break;
                default: updateDustSurfaceDensitySmart(data, disk_params);       break;
            }

            // --- 2. KIÍRÁS ---
            // A leképezett (már integrált és normalizált) adatokat küldjük a fájlba
            printDustSurfaceDensityPressurePressureDerivateFile(
                disk_params->radial_grid,               // X-tengely rácsa
                disk_params->radial_grid,               // Koordináták
                disk_params->dust_surface_density_euler, // Gáz/Por sűrűség (opcionális)
                disk_params->dust_surface_density_euler, // Az imént kiszámolt por felületi sűrűség
                disk_params, 
                sim_opts, 
                output_files, 
                (int)output_time
            );
        }
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

    static void handleSnapshotASCII(double t, double current_time_years, double *output_time, StructuredParticleData *structured_particle_data, DiskParameters *disk_params, const SimulationOptions *sim_opts,
                                    OutputFiles *output_files, char *dens_name, char *dust_name, char *dust_name2, char *size_name) {


        handleSnapshot(t, current_time_years, output_time, sim_opts, output_files,dens_name, dust_name, dust_name2, size_name);
        snapshotInitAtT0(current_time_years,structured_particle_data,disk_params, sim_opts);
        snapshotPrintGas(disk_params, output_files,sim_opts);
        snapshotPrintDust((int)(*output_time), structured_particle_data, disk_params,sim_opts, output_files, size_name);
        snapshotDustSurfacedensity(*output_time, structured_particle_data, disk_params, sim_opts, output_files);


        if(sim_opts->dimension == 2) {
            writeDustField2D(structured_particle_data, sim_opts->output_dir_name, *output_time, NULL);
        }

        PressureTrap current_traps[3];
        int num_found = identifyPressureTraps(disk_params, current_traps, 3);


        for (int i=0;i<disk_params->grid_number;i++) {
            if(disk_params->gas_pressure_gradient_vector[i]>0)
            fprintf(stderr,"DEBUG: gradient[%d] = %.6e\n", i, disk_params->gas_pressure_gradient_vector[i]);
        }



        for (int k = 0; k < num_found; k++) {
            if (current_traps[k].radial_position > 0) {
                double local_H = calculatePressureScaleHeight(current_traps[k].radial_position, disk_params);
                current_traps[k].inner_boundary = current_traps[k].radial_position - 1.0 * local_H;
                current_traps[k].outer_boundary = current_traps[k].radial_position + 1.0 * local_H;
                calculateMassInSpecificTrapStructured(&current_traps[k], structured_particle_data, sim_opts);
            }
        }

        printTrapMassEvolution(*output_time, num_found, current_traps, output_files);
        snapshotAdvance(output_time, sim_opts);
        closeSnapshotFiles(output_files, sim_opts);

    }

    static void handleSnapshotHDF5(double output_time, const SimulationOptions *sim_opts, OutputFiles *output_files, DiskParameters *disk_params, StructuredParticleData *structured_particle_data) {
        char *filename = NULL;
        asprintf(&filename, "%s/%s/%s_%08d%s",
                 sim_opts->output_dir_name, kLogFilesDirectory, kSnapshotOutputFileNamePrefix, (int)(output_time), kFileNamesHDF5Suffix);

        // Init file
        if (initHDF5File(filename, output_files) != 0) {
            fprintf(stderr, "ERROR: Could not initialize HDF5 file %s\n", filename);
            free(filename); // Ne felejtsd el felszabadítani az asprintf által foglalt memóriát!
            return;
        }

        // Write datasets to the already opened file
        if (structured_particle_data != NULL) {
            hid_t file_id = (hid_t)(intptr_t)output_files->hdf5_file;
            
            // Itt hívjuk meg az új verziót!
            writeHDF5SnapshotToFile(output_time, file_id, sim_opts, disk_params, structured_particle_data);
        }

        // Close file
        closeHDF5File(output_files);
        free(filename);
    }

    static void simulateDustDriftStep(double *t, double deltat, double *output_time, DiskParameters *disk_params, const SimulationOptions *sim_opts,
                                      OutputFiles *output_files, char *dens_name, char *dust_name, char *dust_name2, char *size_name, StructuredParticleData *structured_particle_data) {

        double current_time_years = *t / (2.0 * M_PI);


        if (isSnapshotDue(current_time_years, *output_time, deltat, sim_opts)) {


            if (sim_opts->output_format == OUTPUT_ASCII) {
                handleSnapshotASCII(*t, current_time_years, output_time, structured_particle_data, disk_params, sim_opts,
                                    output_files, dens_name, dust_name, dust_name2, size_name);
            } else {
                handleSnapshotHDF5(*output_time, sim_opts, output_files, disk_params, structured_particle_data);
                
                PressureTrap current_traps[MAX_TRAPS];
                int num_found = identifyPressureTraps(disk_params, current_traps, MAX_TRAPS);
                double trap_pos[MAX_TRAPS] = {0.0};
                for (int i = 0; i < num_found && i < MAX_TRAPS; i++) trap_pos[i] = current_traps[i].radial_position;

                double primary_mass[MAX_TRAPS] = {0}, secondary_mass[MAX_TRAPS] = {0}, total_mass[MAX_TRAPS] = {0};

                for (int k = 0; k < num_found; k++) {
                    if (current_traps[k].radial_position > 0.0) {
                        double local_H = calculatePressureScaleHeight(current_traps[k].radial_position, disk_params);
                        current_traps[k].inner_boundary = current_traps[k].radial_position - local_H;
                        current_traps[k].outer_boundary = current_traps[k].radial_position + local_H;
                        calculateMassInSpecificTrapStructured(&current_traps[k], structured_particle_data, sim_opts);
                        primary_mass[k]   = current_traps[k].primary_dust_mass;
                        secondary_mass[k] = current_traps[k].secondary_dust_mass;
                        total_mass[k]     = primary_mass[k] + secondary_mass[k];
                    }
                }
                appendMassTimeSeries(*output_time, primary_mass, secondary_mass, total_mass, trap_pos);
                snapshotAdvance(output_time, sim_opts);
            }
        }

        // --- FIZIKAI SZÁMÍTÁSOK (STRUKTURÁLT MOTOR) ---

        if (sim_opts->option_for_evolution == 1.) {
            refreshGasSurfaceDensityPressurePressureGradient(deltat, disk_params);
        }

        updateDustSurfaceDensitySmart(structured_particle_data, disk_params);

        // 2. Vertikális ülepedés
        applyVerticalSettlingDeterministic(structured_particle_data, disk_params, deltat); 

        // 3. Radiális drift és növekedés (A timescale fájlt ez kezeli t=0-nál!)
        calculateDustDistanceStructured(sim_opts->output_dir_name, structured_particle_data, deltat, *t, sim_opts, disk_params);


        // 5. Régi indexek frissítése a snapshotokhoz (a szinkronizált adatokból)
        updateStructuredParticleGridIndices(structured_particle_data, disk_params);

        *t += deltat;
    }


    static void simulateGasOnlyStep(double *t,double deltat,double *output_time,DiskParameters *disk_params,const SimulationOptions *sim_opts,OutputFiles *output_files,char *dens_name){

        double current_time_years = *t / (2.0 * M_PI);

        if (isSnapshotDue(current_time_years, *output_time, deltat, sim_opts)) {

            if (sim_opts->output_format == OUTPUT_ASCII) {

                asprintf(&dens_name, "%s/%s/%s_%08d%s",sim_opts->output_dir_name,kLogFilesDirectory,kGasDensityProfileFilePrefix,(int)(*output_time),kFileNamesSuffix);

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
                
                    handleSnapshotHDF5(*output_time, sim_opts, output_files, disk_params, NULL);
                
                }
            }
            snapshotAdvance(output_time, sim_opts);
        }

        refreshGasSurfaceDensityPressurePressureGradient(deltat, disk_params);

        *t += deltat;
    }

    void timeIntegrationForTheSystem(SnapshotMode mode, DiskParameters *disk_params, const SimulationOptions *sim_opts, OutputFiles *output_files, StructuredParticleData *structured_particle_data) {

        HeaderData header_data_for_files; 
        double output_time = 0.0;

        if (disk_params == NULL) {
            fprintf(stderr, "ERROR [timeIntegrationForTheSystem]: disk_params_ptr is NULL!\n");
            exit(1);
        }

        // --- Particle number calculation ---
        if (mode > 2) {
            particle_number = sim_opts->number_of_dust_particles;
        } else {
            fprintf(stderr, "DEBUG [timeIntegrationForTheSystem]: Particle drift is OFF. particle_number set to 0.\n");
            particle_number = 0;
        }

        // --- Output files setup ---
        if (mode > 2) {
            if (setupInitialOutputFiles(output_files, sim_opts, disk_params, &header_data_for_files) != 0) {
                fprintf(stderr, "ERROR: Failed to set up initial output files. Exiting.\n");
                exit(EXIT_FAILURE);
            }
        }

        char dens_name[MAX_PATH_LEN] = "";
        char dust_name[MAX_PATH_LEN] = "";
        char dust_name2[MAX_PATH_LEN] = "";
        char size_name[MAX_PATH_LEN] = "";
        double t = 0.0;
        double t_integration_in_internal_units = sim_opts->maximum_simulation_time * 2.0 * M_PI;




        double dt_cfl = calculateTimeStep(disk_params) / 5.0;
        double dt_user = sim_opts->user_defined_time_step;

        double deltat;

        /* user nem adott meg dt-t → használjuk a stabilat */
        if (dt_user <= 0.0) {
            deltat = dt_cfl;
        }
        /* user adott meg → limitáljuk */
        else {
            deltat = fmin(dt_cfl, dt_user);
        }

        /* opcionálisan logoljuk */
        fprintf(stderr,
        "TIMESTEP: dt_cfl=%.3e  dt_user=%.3e  dt_used=%.3e\n",
        dt_cfl, dt_user, deltat);




        // --- StructuredParticleData: át kell venni az init_tool-ból, nem malloc-olni újra ---
        if (structured_particle_data == NULL || structured_particle_data->particles == NULL) {
            fprintf(stderr, "ERROR [timeIntegrationForTheSystem]: StructuredParticleData not initialized! Exiting.\n");
            exit(EXIT_FAILURE);
        }

        // ---- Time series init (HDF5) ----
        if (sim_opts->output_format == OUTPUT_HDF5) {
            char *ts_filename = NULL;

            if (asprintf(&ts_filename,"%s/%s/%s%s", sim_opts->output_dir_name, kLogFilesDirectory, kTimeSeriesForMassAccumulatinFileName, kFileNamesHDF5Suffix) == -1) {
                fprintf(stderr, "ERROR: asprintf failed for time_series filename\n");
                exit(EXIT_FAILURE);
            }

            initializeMassTimeSeries(ts_filename);
            free(ts_filename);
        }

        // ---- Main integration loop ----
        do {
            if (mode > 1) {
                simulateDustDriftStep(&t, deltat, &output_time, disk_params, sim_opts, output_files, dens_name, dust_name, dust_name2, size_name, structured_particle_data);
            } else { 
                simulateGasOnlyStep(&t, deltat, &output_time, disk_params, sim_opts, output_files, dens_name);
            }    
        } while (t <= t_integration_in_internal_units);

        if (sim_opts->output_format == OUTPUT_HDF5) {
            closeMassTimeSeries();
        }



        // ---- Cleanup ----
        cleanupSimulationResources(structured_particle_data, output_files);

        fprintf(stderr, "\n\nDEBUG [timeIntegrationForTheSystem]: Main simulation loop finished (t > t_integration_in_internal_units).\n");
        fprintf(stderr, "DEBUG [timeIntegrationForTheSystem]: Cleanup completed.\n");
    }
