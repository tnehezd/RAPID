
// src/simulation_core.c

// Standard C Library Includes
#include <stdio.h>    // For printf, fopen, fclose, fscanf, snprintf, sprintf
#include <stdlib.h>   // For exit, EXIT_FAILURE, EXIT_SUCCESS, system
#include <math.h>     // For M_PI, fmod, HUGE_VAL (and pow if used by other functions)
#include <string.h>   // For snprintf, sprintf

#include <omp.h>

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
void calculate1DDustDrift(double pradius, double dp, double sigma, double ug, double r, double *drdt, const disk_t *disk_params) {

    double P, H, dPdr, St, csound;
      
    St = calculateStokesNumber(pradius,sigma,disk_params);
    H = calculatePressureScaleHeight(r,disk_params);   
    P = calculateGasPressure(sigma,r,disk_params);
    dPdr = dp;
    csound = calculateLocalSoundSpeed(r,disk_params); 

    *drdt = ug / (1. + St * St) + St / (1. + St * St) * H / P * dPdr * csound;	/* bearamlas sebessege: Birnstiel PHD	*/
}


/* for solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*hu/(2*r)*dsigma/dr	--> 3*nu  	*/
double ftcsSecondDerivativeCoefficient(double r, const disk_t *disk_params){					
    double A;
    A = 3.0 * calculateKinematicViscosity(r, disk_params);
    return A;
}

/* for solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*hu/(2*r)*dsigma/dr	--> 9*nu /(2*r)	*/
double ftcsFirstDerivativeCoefficient(double r, const disk_t *disk_params){							
    double B;
    B = 9.0 * calculateKinematicViscosity(r,disk_params) / (2.0 * r);
    return B;
}

double calculateTimeStep(const disk_t *disk_params) { // Add const here too
    double A_max, stepping;
    int i;

// IT WOULD BE BETTER TO CALCULATE FOR BOTH PARTS OF THE FUNCTION!!!!

    A_max = -10000.0;
    
    for(i = 0; i < disk_params->grid_number; i++) {
        if(ftcsSecondDerivativeCoefficient(disk_params->rvec[i], disk_params) > A_max) {
            A_max = ftcsSecondDerivativeCoefficient(disk_params->rvec[i], disk_params);
        }
    }
    stepping = disk_params->DD * disk_params->DD / (2.0 * A_max);
    fprintf(stderr," Actual calculateTimeStep: DD = %.2e, stepping = %.2e\n", disk_params->DD, stepping);

    return stepping;
}


void timeIntegrationForTheSystem(disk_t *disk_params, const simulation_options_t *sim_opts, output_files_t *output_files) {
    ParticleData_t p_data;
    HeaderData_t header_data_for_files; // Később inicializáljuk a setupInitialOutputFiles-ban


    double L = 0.; // Években mért "pillanatfelvétel" időzítő

    if (disk_params == NULL) {
        fprintf(stderr, "ERROR [timeIntegrationForTheSystem]: disk_params_ptr is NULL!\n");
        exit(1); // Program leállítása, ha kritikus hiba van
    }

    // --- Inicializálási szakasz ---
    if (sim_opts->drift == 1.) {
        particle_number = calculateNumbersOfParticles(sim_opts->dust_input_filename);
    } else {
        fprintf(stderr, "ERROR [timeIntegrationForTheSystem]: Particle drift is OFF. particle_number set to 0.\n");
        particle_number = 0;
    }

    if (particle_number > 0 && allocateParticleData(&p_data, particle_number, (int)sim_opts->twopop) != 0) {
        fprintf(stderr, "ERROR: Failed to allocate particle data. Exiting.\n");
        exit(EXIT_FAILURE);
    }

    // Fájl inicializálás a meglévő io_utils függvény hívásával
    if (sim_opts->drift == 1.) {
        if (setupInitialOutputFiles(output_files, sim_opts, disk_params, &header_data_for_files) != 0) {
            fprintf(stderr, "ERROR: Failed to set up initial output files. Exiting.\n");
            exit(EXIT_FAILURE);
        }
        // loadDustParticlesFromFile hívása a részecskeadatok beolvasására
        loadDustParticlesFromFile(p_data.radius, p_data.radiusmicr, p_data.massvec, p_data.massmicrvec, sim_opts->dust_input_filename);
    }

    // További inicializálások
    double max = 0.0, min = 0.0, max2 = 0.0, min2 = 0.0;
    int i; // Hagyjuk meg ezt a ciklusváltozót a C89 kompatibilitás kedvéért, ha szükséges
    
    // Ideiglenes puffer a fájlneveknek a ciklusban
    char dens_name[MAX_PATH_LEN] = "";
    char dust_name[MAX_PATH_LEN] = "";
    char dust_name2[MAX_PATH_LEN] = "";
    char size_name[MAX_PATH_LEN] = "";

    double t = 0.0;
    double t_integration_in_internal_units = sim_opts->TMAX * 2.0 * M_PI;
    double deltat = calculateTimeStep(disk_params) / 5.0;

    // DT felülbírálása, ha a felhasználó megadott kisebb értéket
    if (sim_opts->DT > 0.0 && sim_opts->DT < deltat) {
        ((simulation_options_t *)sim_opts)->DT = deltat; // Az eredeti kód hibásan a deltat-t vette át, ha sim_opts->DT nagyobb volt
    } else {
        ((simulation_options_t *)sim_opts)->DT = deltat; // Az eredeti kód szerint, ha nem kisebb, akkor deltat a DT
    }

    // Mass accumulation változók
    double masstempiin = 0, massmtempiin = 0, masstempoin = 0, massmtempoin = 0;
    double masstempiout = 0, massmtempiout = 0, masstempoout = 0, massmtempoout = 0;
    double tavin = 0, tavout = 0; // Távolságok a printMassGrowthAtDZEFile-hoz


    if (sim_opts->twopop == 0 && particle_number > 0) {
        for (i = 0; i < particle_number; i++) {
            p_data.radiusmicr[i][0] = 0;
            p_data.radiusmicr[i][1] = 0;
            p_data.partmassmicrind[i][0] = 0;
            p_data.partmassmicrind[i][1] = 0;
            p_data.massmicrvec[i] = 0;
        }
    }

    // --- Fő szimulációs ciklus ---
    do {
        if (sim_opts->drift == 1.) {

            // Radius reciprok számítása a min/max kereséshez
            for (i = 0; i < particle_number; i++) {
                if (p_data.radius[i][0] > 0. && p_data.radius[i][0] > disk_params->r_min) {
                    p_data.radius_rec[i][0] = 1. / p_data.radius[i][0];
                } else {
                    p_data.radius_rec[i][0] = 0.; // Vagy valamilyen "érvénytelen" érték, ami kizárja a min/max keresésből
                }
            }

            max = findMaximumOfAnArray(p_data.radius, particle_number);
            min = findMaximumOfAnArray(p_data.radius_rec, particle_number); // Megkeresi a távolság reciprokának maximumát
            min = 1. / min; // Ebből lesz a távolság minimuma

            double mint, maxt;

            if (sim_opts->twopop == 1) {
                // Micron részecskék radius reciprok számítása
                for (i = 0; i < particle_number; i++) {
                    if (p_data.radiusmicr[i][0] > 0. && p_data.radiusmicr[i][0] > disk_params->r_min) {
                        p_data.radius_rec[i][0] = 1. / p_data.radiusmicr[i][0];
                    } else {
                        p_data.radius_rec[i][0] = 0.;
                    }
                }

                max2 = findMaximumOfAnArray(p_data.radiusmicr, particle_number);
                min2 = findMaximumOfAnArray(p_data.radius_rec, particle_number);
                min2 = 1. / min2;

                mint = findMinimumOfAnArray(min, min2, HUGE_VAL); // Itt a te findMinimumOfAnArray(s1, s2, s3) függvényedet használjuk
                // Ahogy korábban beszéltük, ha findMaximumOfAnArray(s1,s2,s3) lenne, az jobb lenne,
                // de a meglévő findMinimumOfAnArray-t használva a reciprok trükkkel:
                maxt = findMinimumOfAnArray(1. / max, 1. / max2, HUGE_VAL);
                maxt = 1. / maxt;
            } else {
                mint = min;
                maxt = max;
            }

            // --- Kimeneti adatok (pillanatfelvétel) kezelése ---
            double current_time_years = t / (2.0 * M_PI);
            if ((fmod(current_time_years, (sim_opts->TMAX / sim_opts->WO)) < deltat || current_time_years == 0) && L - current_time_years < deltat) {
                fprintf(stderr,"\n--- Simulation Time: %.2e years (Internal time: %.2e, L: %.2e) ---\n", current_time_years, t, L);


                if (current_time_years != 0) {
                    if (sim_opts->evol == 1) {
                        snprintf(dens_name, MAX_PATH_LEN, "%s/%s/%s_%08d%s", sim_opts->output_dir_name, kLogFilesDirectory, kGasDensityProfileFilePrefix, (int)L,kFileNamesSuffix);
                    }
                }

                snprintf(dust_name, MAX_PATH_LEN, "%s/%s/%s_%08d%s", sim_opts->output_dir_name, kLogFilesDirectory,kDustDensityProfileFilePrefix, (int)L,kFileNamesSuffix);
                snprintf(dust_name2, MAX_PATH_LEN, "%s/%s/%s_%08d%s", sim_opts->output_dir_name, kLogFilesDirectory,kDustDensityProfileFilePrefix,(int)L,kFileNamesSuffix);
                snprintf(size_name, MAX_PATH_LEN, "%s/%s/%s_%08d%s", sim_opts->output_dir_name, kLogFilesDirectory,kDustParticleSizeFileName, (int)L,kFileNamesSuffix);

                // Fájlok megnyitása és fejlécek írása
                output_files->surface_file = fopen(dens_name, "w");
                if(L != 0) {
                    if (output_files->surface_file == NULL) {
                        fprintf(stderr, "ERROR: Could not open %s for writing.\n", dens_name);
                    } else {
                        HeaderData_t gas_header_data = {.current_time = current_time_years, .is_initial_data = (current_time_years == 0.0)};
                        printFileHeader(output_files->surface_file, FILE_TYPE_GAS_DENSITY, &gas_header_data);
                    }
                }

                output_files->dust_file = fopen(dust_name, "w");
                if (output_files->dust_file == NULL) {
                    fprintf(stderr, "ERROR: Could not open %s for writing.\n", dust_name);
                } else {
                    HeaderData_t dust_header_data = {.current_time = current_time_years, .is_initial_data = (current_time_years == 0.0)};
                    printFileHeader(output_files->dust_file, FILE_TYPE_DUST_MOTION, &dust_header_data);
                }

                if (sim_opts->twopop == 1.) {
                    output_files->micron_dust_file = fopen(dust_name2, "w");
                    if (output_files->micron_dust_file == NULL) {
                        fprintf(stderr, "ERROR: Could not open %s for writing.\n", dust_name2);
                    } else {
                        HeaderData_t micron_dust_header_data = {.current_time = current_time_years, .is_initial_data = (current_time_years == 0.0)};
                        printFileHeader(output_files->micron_dust_file, FILE_TYPE_MICRON_MOTION, &micron_dust_header_data);
                    }
                }

                // Eredeti t==0 logika
                if (current_time_years == 0) {
                    updateParticleGridIndices(p_data.radius, p_data.partmassind, p_data.massvec, t, particle_number, disk_params);
                    if (sim_opts->twopop == 1) updateParticleGridIndices(p_data.radiusmicr, p_data.partmassmicrind, p_data.massmicrvec, t, particle_number, disk_params);

                    if (sim_opts->growth == 1.) {
                        calculateDustSurfaceDensity(max, min, p_data.radius, p_data.radiusmicr, p_data.sigmad, p_data.sigmadm, p_data.massvec, p_data.massmicrvec, p_data.rdvec, p_data.rmicvec, sim_opts, disk_params);
                    }
                }

                // Gas density output
                if (sim_opts->evol == 1 || current_time_years == 0) {
                    if(L != 0) printGasSurfaceDensityPressurePressureDerivateFile(disk_params, output_files);
                }

                // Particle position and size output
                if (sim_opts->drift == 1) {
                    printDustParticleSizeFile(size_name, (int)L, p_data.radius, p_data.radiusmicr, disk_params, sim_opts, output_files);
                }

                // Reset mass accumulation variables for next interval
                masstempiout = 0;
                massmtempiout = 0;
                masstempoout = 0;
                massmtempoout = 0;


                // Resetting partmassind[k][3] and [k][4]
                if (sim_opts->dzone == 1.0 && particle_number > 0) { // Ellenőrzés particle_number-re
                    for (int k = 0; k < particle_number; k++) {
                        p_data.partmassind[k][3] = 0.0;
                        p_data.partmassind[k][4] = 0.0;
                    }
                    
                }

                printMassGrowthAtDZEFile(L, p_data.partmassind, p_data.partmassmicrind, t, masstempiin, masstempoin, massmtempiin, massmtempoin, &masstempiout, &masstempoout, &massmtempiout, &massmtempoout, &tavin, &tavout, disk_params, sim_opts, output_files);
                // Update input mass for next printMassGrowthAtDZEFile call
                masstempiin = masstempiout;
                massmtempiin = massmtempiout;
                masstempoin = masstempoout;
                massmtempoin = massmtempoout;

                if (sim_opts->growth == 1.) {
                    printDustSurfaceDensityPressurePressureDerivateFile(p_data.rdvec, p_data.rmicvec, p_data.sigmad, p_data.sigmadm, disk_params, sim_opts, output_files, (int)L);
                }
                fprintf(stderr,"L set to %lg\n",L);

                L = L + (double)(sim_opts->TMAX / sim_opts->WO);
                // Fájlok bezárása, amelyek csak ezen időintervallumban voltak nyitva
                closeSnapshotFiles(output_files, dens_name, dust_name, dust_name2, sim_opts);
            }

            // Gas evolution
            if (sim_opts->evol == 1.) {
                refreshGasSurfaceDensityPressurePressureGradient(sim_opts, disk_params);
            }

            // Count masses and get sigma_d for the next step (always done)
            updateParticleGridIndices(p_data.radius, p_data.partmassind, p_data.massvec, t, particle_number, disk_params);
            if (sim_opts->twopop == 1) updateParticleGridIndices(p_data.radiusmicr, p_data.partmassmicrind, p_data.massmicrvec, t, particle_number, disk_params);

            if (sim_opts->growth == 1.) {
                calculateDustSurfaceDensity(max, min, p_data.radius, p_data.radiusmicr, p_data.sigmad, p_data.sigmadm, p_data.massvec, p_data.massmicrvec, p_data.rdvec, p_data.rmicvec, sim_opts, disk_params);
            }

            // Get radii for next step
            int optsize = 0;
            calculateDustDistance(sim_opts->output_dir_name, optsize, p_data.radius, p_data.sigmad, p_data.rdvec, deltat, t, particle_number, sim_opts, disk_params);

            if (sim_opts->twopop == 1.) {
                optsize = 1;
                calculateDustDistance(sim_opts->output_dir_name, optsize, p_data.radiusmicr, p_data.sigmadm, p_data.rmicvec, deltat, t, particle_number, sim_opts, disk_params);

            }

            t = t + deltat;

            // Kilépési feltétel a drift == 1 ágon
            // Fontos: maxt és mint frissül az előző szakaszban, azt használjuk itt.
//            if (!(maxt >= disk_params->r_min && mint != maxt)) {
//                fprintf(stderr,"DEBUG [timeIntegrationForTheSystem]: Simulation termination condition met (maxt < r_min or mint == maxt) at time: %lg.\n",L);

//                goto cleanup; // Ugrás a tisztításra
//            }

        } else { // sim_opts->drift == 0. (Gas-only simulation)
            double current_time_years = t / (2.0 * M_PI);

            if ((fmod(current_time_years, (sim_opts->TMAX / sim_opts->WO)) < deltat || current_time_years == 0) && L - current_time_years < deltat) {
                fprintf(stderr,"\n--- Simulation Time: %.2e years (Internal time: %.2e, L: %.2e) ---\n", current_time_years, t, L);

                fprintf(stderr,"DEBUG [timeIntegrationForTheSystem]: Outputting data for gas-only simulation at time %.2e. L=%.2e\n", current_time_years, L);
                snprintf(dens_name, MAX_PATH_LEN, "%s/%s/%s_%08d.dat", sim_opts->output_dir_name, kLogFilesDirectory, kGasDensityProfileFilePrefix, (int)L);
                fprintf(stderr, "DEBUG [timeIntegrationForTheSystem]: Outputting %s_%08d.dat to %s.\n", kGasDensityProfileFilePrefix, (int)L, dens_name);

                output_files->surface_file = fopen(dens_name, "w");
                if (output_files->surface_file == NULL) {
                    fprintf(stderr, "ERROR: Could not open %s for writing in gas-only branch.\n", dens_name);
                } else {
                    fprintf(stderr, "DEBUG [timeIntegrationForTheSystem]: Opened %s for writing in gas-only branch.\n", dens_name);
                    HeaderData_t gas_header_data = {.current_time = current_time_years, .is_initial_data = (current_time_years == 0.0)};
                    printFileHeader(output_files->surface_file, FILE_TYPE_GAS_DENSITY, &gas_header_data);
                }

                printGasSurfaceDensityPressurePressureDerivateFile(disk_params, output_files);


                if (output_files->surface_file != NULL) {
                    fclose(output_files->surface_file);
                    output_files->surface_file = NULL;
                    fprintf(stderr, "DEBUG [timeIntegrationForTheSystem]: Closed %s in gas-only branch.\n", dens_name);
                }

                L = L + (double)(sim_opts->TMAX / sim_opts->WO);
                fprintf(stderr,"DEBUG [timeIntegrationForTheSystem]: Updated L to %.2e.\n", L);
            }

            fprintf(stderr,"DEBUG [timeIntegrationForTheSystem]: Calling refreshGasSurfaceDensityPressurePressureGradient for gas-only evolution.\n");
            refreshGasSurfaceDensityPressurePressureGradient(sim_opts, disk_params);
            
            t = t + deltat;
        }

    } while (t <= t_integration_in_internal_units);

    fprintf(stderr,"\n\nDEBUG [timeIntegrationForTheSystem]: Main simulation loop finished (t > t_integration_in_internal_units).\n");

cleanup:
    // --- Tisztítási szakasz ---
    cleanupSimulationResources(&p_data, output_files, sim_opts);
    fprintf(stderr,"DEBUG [timeIntegrationForTheSystem]: Cleanup completed.\n");
}

