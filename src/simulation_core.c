
// src/simulation_core.c

// Standard C Library Includes
#include <stdio.h>    // For printf, fopen, fclose, fscanf, snprintf, sprintf
#include <stdlib.h>   // For exit, EXIT_FAILURE, EXIT_SUCCESS, system
#include <math.h>     // For M_PI, fmod, HUGE_VAL (and pow if used by other functions)
#include <string.h>   // For snprintf, sprintf

#include <omp.h>

// Your Project Header Includes
#include "config.h"       // For PARTICLE_NUMBER, TMAX, WO, RMIN, DT, optdr, sim_opts->twopop, sim_opts->growth, optev, r_dze_i, r_dze_o
#include "io_utils.h"     // For timePar (though not called in tIntegrate, it's io-related), reszecskek_szama, por_be, Print_Sigma, Print_Pormozg_Size, Print_Mass, Print_Sigmad. Also for globals: filenev1, filenev3, fout, foutmicr, massfil
#include "disk_model.h"   // If any disk_model functions are called (e.g., Perem indirectly if sigma/press depend on it) - Though not directly visible in tIntegrate, often needed for global disk parameters. Add if you hit implicit declaration for disk_model functions.
#include "dust_physics.h" // For Count_Mass, secondaryGrowth, find_max, find_min, Get_Sigmad, Get_Radius
#include "utils.h"        // For time_step, Get_Sigma_P_dP, and potentially other utility functions
#include "simulation_core.h"
#include "particle_data.h" // Új include


/*	Kiszamolja az 1D-s driftet	*/
/*  	dr/dt = St/(1+St*St)*H(r)/r*dlnP/dlnr*cs = St/(1+St*St) * (H/r) * (r/P) * (dP/dr) * cs		*/
void eqrhs(double pradius, double dp, double sigma, double ug, double r, double *drdt, const disk_t *disk_params) {
    // DEBUG [eqrhs]: Entry point
    // printf("DEBUG [eqrhs]: r=%.2e, pradius=%.2e, dp=%.2e, sigma=%.2e, ug=%.2e\n", r, pradius, dp, sigma, ug);

    double P, H, dPdr, St, csound;
      
    St = Stokes_Number(pradius,sigma,disk_params);
    H = scale_height(r,disk_params);   
    P = press(sigma,r,disk_params);
    dPdr = dp;
    csound = c_sound(r,disk_params); 

    *drdt = ug / (1. + St * St) + St / (1. + St * St) * H / P * dPdr * csound;	/* bearamlas sebessege: Birnstiel PHD	*/
//    printf("DEBUG [eqrhs]: St=%.2e, H=%.2e, P=%.2e, dPdr=%.2e, csound=%.2e, drdt=%.2e\n", St, H, P, dPdr, csound, *drdt);
}


/* for solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*hu/(2*r)*dsigma/dr	--> 3*nu = Coeff_1 	*/
double Coeff_1(double r, const disk_t *disk_params){					
    // printf("DEBUG [Coeff_1]: r=%.2e\n", r); // This might be called very frequently, uncomment only if needed
    double A;
    A = 3.0 * visc(r, disk_params);
    return A;
}

/* for solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*hu/(2*r)*dsigma/dr	--> 9*nu /(2*r) = Coeff_2 	*/
double Coeff_2(double r, const disk_t *disk_params){							
    // printf("DEBUG [Coeff_2]: r=%.2e\n", r); // This might be called very frequently, uncomment only if needed
    double B;
    B = 9.0 * visc(r,disk_params) / (2.0 * r);
    return B;
}

double time_step(const disk_t *disk_params) { // Add const here too
    printf("DEBUG [time_step]: Calculating minimum time step.\n");
    fprintf(stderr, "DEBUG [time_step]: Entry. disk_params address=%p, FLIND=%.2f, HASP=%.2f\n",
            (void*)disk_params, disk_params->FLIND, disk_params->HASP);

    double A_max, stepping;
    int i;

    A_max = -10000.0;
    
    for(i = 0; i < disk_params->NGRID; i++) {
        if(Coeff_1(disk_params->rvec[i], disk_params) > A_max) {
            A_max = Coeff_1(disk_params->rvec[i], disk_params);
        }
    }
    printf("DEBUG [time_step]: A_max = %.2e\n", A_max);

    stepping = disk_params->DD * disk_params->DD / (2.0 * A_max);
    printf("DEBUG [time_step]: DD = %.2e, stepping = %.2e\n", disk_params->DD, stepping);

    return stepping;
}


/*	Runge-Kutta4 integrator	*/
// prad bemenet: AU-ban!
void int_step(double time, double prad, const double *sigmad, const double *rdvec, double step, double y, double *ynew, double *pradnew, const disk_t *disk_params, const simulation_options_t *sim_opts){
    // printf("DEBUG [int_step]: Called for time=%.2e, y=%.2e, prad=%.2e, step=%.2e\n", time, y, prad, step);
    double dy1,dy2,dy3,dy4;
    double ytemp, ytemp2;
    double sigma, dpress, ugas; 
    double pdens, p;
    double pradtemp;
    int opt = 0;
    double sigmadd = 0.0;
    
/*	Mivel a kulongozo parametereket csak a megadott gridcella pontokban ismerjuk, de ez nem feltetlen egyezik meg a reszecskek pozicijaval, ezert minden fontos parametert interpolalunk a reszecskek tavolsagara	*/
    interpol(disk_params->sigmavec,disk_params->rvec,y,&sigma,disk_params->DD,opt,disk_params);
    interpol(disk_params->dpressvec,disk_params->rvec,y,&dpress,disk_params->DD,opt,disk_params);
    interpol(disk_params->ugvec,disk_params->rvec,y,&ugas,disk_params->DD,opt,disk_params);
    // printf("DEBUG [int_step]: Interpolated: sigma=%.2e, dpress=%.2e, ugas=%.2e\n", sigma, dpress, ugas);

    double dd = (disk_params->RMAX - disk_params->RMIN) / (PARTICLE_NUMBER-1);
    int dker = (int)(1./dd);//
    dker = dker * KEREK;
    double ddker = (double) dker;
    int temp;

    temp = (int)floor(y * ddker+0.5);
    ytemp2 = (double)temp / ddker;
    
    int i;
    for(i=0;i<PARTICLE_NUMBER;i++) {
        if(ytemp2 == rdvec[i]) {
            sigmadd = sigmad[i];
            break;
        }
    }
    // if (i == PARTICLE_NUMBER) {
    //     printf("DEBUG [int_step]: Warning: ytemp2 (%.2e) not found in rdvec.\n", ytemp2);
    // } else {
    //     printf("DEBUG [int_step]: sigmadd found: %.2e at index %d\n", sigmadd, i);
    // }

    if(sim_opts->growth == 1.) {		// ha van reszecskenovekedes
        // printf("DEBUG [int_step]: Particle growth is ON.\n");
        if(time != 0.) {	// ha nem t0 idopontban vagyunk
            // printf("DEBUG [int_step]: Time is not zero (%.2e). Calculating particle growth.\n", time);
            pradtemp = prad;
            interpol(disk_params->pressvec,disk_params->rvec,y,&p,disk_params->DD,opt,disk_params);
            pdens = disk_params->PDENSITY; 
            pradtemp = getSize(prad,pdens,sigma,sigmadd,y,p,dpress,step,disk_params);	// itt szamolja a reszecskenovekedest
            prad = pradtemp;
            // printf("DEBUG [int_step]: New prad (after growth): %.2e\n", prad);
        }
    }

    *pradnew = prad;

/*	Itt szamolja a reszecske poziciojat	*/
    eqrhs(prad, dpress, sigma, ugas, y, &dy1,disk_params);
    // printf("DEBUG [int_step]: RK4 dy1=%.2e\n", dy1);

    ytemp = y + 0.5 * step * dy1;
    eqrhs(prad, dpress, sigma, ugas, ytemp, &dy2,disk_params);
    // printf("DEBUG [int_step]: RK4 dy2=%.2e\n", dy2);
        
    ytemp = y + 0.5 * step * dy2;
    eqrhs(prad, dpress, sigma, ugas, ytemp, &dy3,disk_params);
    // printf("DEBUG [int_step]: RK4 dy3=%.2e\n", dy3);
    
    ytemp = y + step * dy3;
    eqrhs(prad, dpress, sigma, ugas, ytemp, &dy4,disk_params);
    // printf("DEBUG [int_step]: RK4 dy4=%.2e\n", dy4);

    *ynew = y + step * (dy1 + 2.0 * dy2 + 2.0 * dy3 + dy4) / 6.0;
    // printf("DEBUG [int_step]: New ynew=%.2e\n", *ynew);
}



void tIntegrate(disk_t *disk_params, const simulation_options_t *sim_opts, output_files_t *output_files) {
    ParticleData_t p_data;
    HeaderData_t header_data_for_files; // Később inicializáljuk a setup_initial_output_files-ban

    if (disk_params == NULL) {
        fprintf(stderr, "ERROR [tIntegrate]: disk_params_ptr is NULL!\n");
        exit(1); // Program leállítása, ha kritikus hiba van
    }
    fprintf(stderr, "DEBUG_TINTEGRATE_DISK_PARAMS_PTR: Address of disk_params_ptr: %p\n", (void*)disk_params);
    fprintf(stderr, "DEBUG_TINTEGRATE_DISK_PARAMS_H_F: HASP=%.10lg, FLIND=%.10lg\n",
            disk_params->HASP, disk_params->FLIND);

    // --- Inicializálási szakasz ---
    if (sim_opts->drift == 1.) {
        PARTICLE_NUMBER = reszecskek_szama(sim_opts->dust_input_filename);
    } else {
        fprintf(stderr, "ERROR [tIntegrate]: Particle drift is OFF. PARTICLE_NUMBER set to 0.\n");
        PARTICLE_NUMBER = 0;
    }

    if (PARTICLE_NUMBER > 0 && allocate_particle_data(&p_data, PARTICLE_NUMBER, (int)sim_opts->twopop) != 0) {
        fprintf(stderr, "ERROR: Failed to allocate particle data. Exiting.\n");
        exit(EXIT_FAILURE);
    }

    // Fájl inicializálás a meglévő io_utils függvény hívásával
    if (sim_opts->drift == 1.) {
        if (setup_initial_output_files(output_files, sim_opts, disk_params, &header_data_for_files) != 0) {
            fprintf(stderr, "ERROR: Failed to set up initial output files. Exiting.\n");
            exit(EXIT_FAILURE);
        }
        // por_be hívása a részecskeadatok beolvasására
        por_be(p_data.radius, p_data.radiusmicr, p_data.massvec, p_data.massmicrvec, sim_opts->dust_input_filename);
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
    double t_integration = sim_opts->TMAX * 2.0 * M_PI;
    double deltat = time_step(disk_params) / 5.0;

    // DT felülbírálása, ha a felhasználó megadott kisebb értéket
    if (sim_opts->DT > 0.0 && sim_opts->DT < deltat) {
        ((simulation_options_t *)sim_opts)->DT = deltat; // Az eredeti kód hibásan a deltat-t vette át, ha sim_opts->DT nagyobb volt
    } else {
        ((simulation_options_t *)sim_opts)->DT = deltat; // Az eredeti kód szerint, ha nem kisebb, akkor deltat a DT
    }

    // Mass accumulation változók
    double L = 0.; // Években mért "pillanatfelvétel" időzítő
    double masstempiin = 0, massmtempiin = 0, masstempoin = 0, massmtempoin = 0;
    double masstempiout = 0, massmtempiout = 0, masstempoout = 0, massmtempoout = 0;
    double tavin = 0, tavout = 0; // Távolságok a Print_Mass-hoz


    if (sim_opts->twopop == 0 && PARTICLE_NUMBER > 0) {
        for (i = 0; i < PARTICLE_NUMBER; i++) {
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
            for (i = 0; i < PARTICLE_NUMBER; i++) {
                if (p_data.radius[i][0] > 0. && p_data.radius[i][0] > disk_params->RMIN) {
                    p_data.radius_rec[i][0] = 1. / p_data.radius[i][0];
                } else {
                    p_data.radius_rec[i][0] = 0.; // Vagy valamilyen "érvénytelen" érték, ami kizárja a min/max keresésből
                }
            }

            max = find_max(p_data.radius, PARTICLE_NUMBER);
            min = find_max(p_data.radius_rec, PARTICLE_NUMBER); // Megkeresi a távolság reciprokának maximumát
            min = 1. / min; // Ebből lesz a távolság minimuma

            double mint, maxt;

            if (sim_opts->twopop == 1) {
                // Micron részecskék radius reciprok számítása
                for (i = 0; i < PARTICLE_NUMBER; i++) {
                    if (p_data.radiusmicr[i][0] > 0. && p_data.radiusmicr[i][0] > disk_params->RMIN) {
                        p_data.radius_rec[i][0] = 1. / p_data.radiusmicr[i][0];
                    } else {
                        p_data.radius_rec[i][0] = 0.;
                    }
                }

                max2 = find_max(p_data.radiusmicr, PARTICLE_NUMBER);
                min2 = find_max(p_data.radius_rec, PARTICLE_NUMBER);
                min2 = 1. / min2;

                mint = find_min(min, min2, HUGE_VAL); // Itt a te find_min(s1, s2, s3) függvényedet használjuk
                // Ahogy korábban beszéltük, ha find_max(s1,s2,s3) lenne, az jobb lenne,
                // de a meglévő find_min-t használva a reciprok trükkkel:
                maxt = find_min(1. / max, 1. / max2, HUGE_VAL);
                maxt = 1. / maxt;
            } else {
                mint = min;
                maxt = max;
            }

            // --- Kimeneti adatok (pillanatfelvétel) kezelése ---
            double current_time_years = t / (2.0 * M_PI);
            if ((fmod(current_time_years, (sim_opts->TMAX / sim_opts->WO)) < deltat || current_time_years == 0) && L - current_time_years < deltat) {
                printf("\n--- Simulation Time: %.2e years (Internal time: %.2e, L: %.2e) ---\n", current_time_years, t, L);

                if (current_time_years != 0) {
                    if (sim_opts->evol == 1) {
                        snprintf(dens_name, MAX_PATH_LEN, "%s/%s/%s_%08d.dat", sim_opts->output_dir_name, LOGS_DIR, FILE_DENS_PREFIX, (int)L);
                    }
                }

                snprintf(dust_name, MAX_PATH_LEN, "%s/%s/dust.%i.dat", sim_opts->output_dir_name, LOGS_DIR, (int)L);
                snprintf(dust_name2, MAX_PATH_LEN, "%s/%s/dustmic.%i.dat", sim_opts->output_dir_name, LOGS_DIR, (int)L);
                snprintf(size_name, MAX_PATH_LEN, "%s/%s/size.%d.dat", sim_opts->output_dir_name, LOGS_DIR, (int)L);

                fprintf(stderr, "DEBUG [tIntegrate]: Output file names set: %s, %s, %s.\n", dust_name, dust_name2, size_name);

                // Fájlok megnyitása és fejlécek írása
                output_files->surface_file = fopen(dens_name, "w");
                if (output_files->surface_file == NULL) {
                    fprintf(stderr, "ERROR: Could not open %s for writing.\n", dens_name);
                } else {
                    fprintf(stderr, "DEBUG [tIntegrate]: Opened %s for writing.\n", dens_name);
                    HeaderData_t gas_header_data = {.current_time = current_time_years, .is_initial_data = (current_time_years == 0.0)};
                    print_file_header(output_files->surface_file, FILE_TYPE_GAS_DENSITY, &gas_header_data);
                }

                output_files->dust_file = fopen(dust_name, "w");
                if (output_files->dust_file == NULL) {
                    fprintf(stderr, "ERROR: Could not open %s for writing.\n", dust_name);
                } else {
                    fprintf(stderr, "DEBUG [tIntegrate]: Opened %s for writing.\n", dust_name);
                    HeaderData_t dust_header_data = {.current_time = current_time_years, .is_initial_data = (current_time_years == 0.0)};
                    print_file_header(output_files->dust_file, FILE_TYPE_DUST_MOTION, &dust_header_data);
                }

                if (sim_opts->twopop == 1.) {
                    output_files->micron_dust_file = fopen(dust_name2, "w");
                    if (output_files->micron_dust_file == NULL) {
                        fprintf(stderr, "ERROR: Could not open %s for writing.\n", dust_name2);
                    } else {
                        fprintf(stderr, "DEBUG [tIntegrate]: Opened %s for writing.\n", dust_name2);
                        HeaderData_t micron_dust_header_data = {.current_time = current_time_years, .is_initial_data = (current_time_years == 0.0)};
                        print_file_header(output_files->micron_dust_file, FILE_TYPE_MICRON_MOTION, &micron_dust_header_data);
                    }
                }

                // Eredeti t==0 logika
                if (current_time_years == 0) {
                    Count_Mass(p_data.radius, p_data.partmassind, p_data.massvec, t, PARTICLE_NUMBER, disk_params);
                    if (sim_opts->twopop == 1) Count_Mass(p_data.radiusmicr, p_data.partmassmicrind, p_data.massmicrvec, t, PARTICLE_NUMBER, disk_params);

                    if (sim_opts->growth == 1.) {
                        Get_Sigmad(max, min, p_data.radius, p_data.radiusmicr, p_data.sigmad, p_data.sigmadm, p_data.massvec, p_data.massmicrvec, p_data.rdvec, p_data.rmicvec, sim_opts, disk_params);
                    }
                }

                // Gas density output
                if (sim_opts->evol == 1 || current_time_years == 0) {
                    Print_Sigma(disk_params, output_files);
                }

                // Particle position and size output
                if (sim_opts->drift == 1) {
                    fprintf(stderr, "DEBUG [tIntegrate]: sim_opts->drift is ON. Calling Print_Pormozg_Size.\n");
                    Print_Pormozg_Size(size_name, (int)L, p_data.radius, p_data.radiusmicr, disk_params, sim_opts, output_files);
                    fprintf(stderr, "DEBUG [tIntegrate]: Print_Pormozg_Size completed.\n");
                }

                // Reset mass accumulation variables for next interval
                masstempiout = 0;
                massmtempiout = 0;
                masstempoout = 0;
                massmtempoout = 0;


                // Resetting partmassind[k][3] and [k][4]
                if (sim_opts->dzone == 1.0 && PARTICLE_NUMBER > 0) { // Ellenőrzés PARTICLE_NUMBER-re
                    for (int k = 0; k < PARTICLE_NUMBER; k++) {
                        p_data.partmassind[k][3] = 0.0;
                        p_data.partmassind[k][4] = 0.0;
                    }
                    
                }

                fprintf(stderr, "DEBUG [tIntegrate]: Calling Print_Mass.\n");

                Print_Mass(L, p_data.partmassind, p_data.partmassmicrind, t, masstempiin, masstempoin, massmtempiin, massmtempoin, &masstempiout, &masstempoout, &massmtempiout, &massmtempoout, &tavin, &tavout, disk_params, sim_opts, output_files);

                fprintf(stderr, "DEBUG [tIntegrate]: Print_Mass completed. Outputs: masstempiout=%.2e, massmtempiout=%.2e\n", masstempiout, massmtempiout);

                printf("masstempiin %lg mtempiin %lg oin %lg  moin %lg ue o-ra: %lg  %lg  %lg  %lg tavin %lg out %lg\n", masstempiin, massmtempiin, masstempoin, massmtempoin,masstempiout, massmtempiout, masstempoout, massmtempoout,tavin, tavout);

                // Update input mass for next Print_Mass call
                masstempiin = masstempiout;
                massmtempiin = massmtempiout;
                masstempoin = masstempoout;
                massmtempoin = massmtempoout;

                if (sim_opts->growth == 1.) {
                    fprintf(stderr, "DEBUG [tIntegrate]: sim_opts->growth is ON. Calling Print_Sigmad (refreshing).\n");
                    Print_Sigmad(p_data.rdvec, p_data.rmicvec, p_data.sigmad, p_data.sigmadm, disk_params, sim_opts, output_files);
                    fprintf(stderr, "DEBUG [tIntegrate]: Print_Sigmad (refreshing) completed.\n");
                }

                L = L + (double)(sim_opts->TMAX / sim_opts->WO);
                printf("DEBUG [tIntegrate]: Updated L to %.2e.\n", L);

                // Fájlok bezárása, amelyek csak ezen időintervallumban voltak nyitva
                close_snapshot_files(output_files, dens_name, dust_name, dust_name2, sim_opts);
            }

            // Gas evolution
            if (sim_opts->evol == 1.) {
                Get_Sigma_P_dP(sim_opts, disk_params);
            }

            // Count masses and get sigma_d for the next step (always done)
            Count_Mass(p_data.radius, p_data.partmassind, p_data.massvec, t, PARTICLE_NUMBER, disk_params);
            if (sim_opts->twopop == 1) Count_Mass(p_data.radiusmicr, p_data.partmassmicrind, p_data.massmicrvec, t, PARTICLE_NUMBER, disk_params);

            if (sim_opts->growth == 1.) {
                Get_Sigmad(max, min, p_data.radius, p_data.radiusmicr, p_data.sigmad, p_data.sigmadm, p_data.massvec, p_data.massmicrvec, p_data.rdvec, p_data.rmicvec, sim_opts, disk_params);
            }

            // Get radii for next step
            int optsize = 0;
            Get_Radius(sim_opts->output_dir_name, optsize, p_data.radius, p_data.sigmad, p_data.rdvec, deltat, t, PARTICLE_NUMBER, sim_opts, disk_params);

            if (sim_opts->twopop == 1.) {
                optsize = 1;
                Get_Radius(sim_opts->output_dir_name, optsize, p_data.radiusmicr, p_data.sigmadm, p_data.rmicvec, deltat, t, PARTICLE_NUMBER, sim_opts, disk_params);

            }

            t = t + deltat;

            // Kilépési feltétel a drift == 1 ágon
            // Fontos: maxt és mint frissül az előző szakaszban, azt használjuk itt.
            if (!(maxt >= disk_params->RMIN && mint != maxt)) {
                printf("DEBUG [tIntegrate]: Simulation termination condition met (maxt < RMIN or mint == maxt).\n");
                printf("A program sikeresen lefutott az integralasi ido vege elott (t: %lg). \n\nNyomj ENTER-t a kilepeshez!\n", L);
                goto cleanup; // Ugrás a tisztításra
            }

        } else { // sim_opts->drift == 0. (Gas-only simulation)
            double current_time_years = t / (2.0 * M_PI);

            if ((fmod(current_time_years, (sim_opts->TMAX / sim_opts->WO)) < deltat || current_time_years == 0) && L - current_time_years < deltat) {
                printf("\n--- Simulation Time: %.2e years (Internal time: %.2e, L: %.2e) ---\n", current_time_years, t, L);

                printf("DEBUG [tIntegrate]: Outputting data for gas-only simulation at time %.2e. L=%.2e\n", current_time_years, L);
                snprintf(dens_name, MAX_PATH_LEN, "%s/%s/%s_%08d.dat", sim_opts->output_dir_name, LOGS_DIR, FILE_DENS_PREFIX, (int)L);
                fprintf(stderr, "DEBUG [tIntegrate]: Outputting %s_%08d.dat to %s.\n", FILE_DENS_PREFIX, (int)L, dens_name);

                output_files->surface_file = fopen(dens_name, "w");
                if (output_files->surface_file == NULL) {
                    fprintf(stderr, "ERROR: Could not open %s for writing in gas-only branch.\n", dens_name);
                } else {
                    fprintf(stderr, "DEBUG [tIntegrate]: Opened %s for writing in gas-only branch.\n", dens_name);
                    HeaderData_t gas_header_data = {.current_time = current_time_years, .is_initial_data = (current_time_years == 0.0)};
                    print_file_header(output_files->surface_file, FILE_TYPE_GAS_DENSITY, &gas_header_data);
                }

                Print_Sigma(disk_params, output_files);
                printf("DEBUG [tIntegrate]: Print_Sigma completed.\n");

                if (output_files->surface_file != NULL) {
                    fclose(output_files->surface_file);
                    output_files->surface_file = NULL;
                    fprintf(stderr, "DEBUG [tIntegrate]: Closed %s in gas-only branch.\n", dens_name);
                }

                L = L + (double)(sim_opts->TMAX / sim_opts->WO);
                printf("DEBUG [tIntegrate]: Updated L to %.2e.\n", L);
            }

            printf("DEBUG [tIntegrate]: Calling Get_Sigma_P_dP for gas-only evolution.\n");
            Get_Sigma_P_dP(sim_opts, disk_params);
            printf("DEBUG [tIntegrate]: Get_Sigma_P_dP completed.\n");
            t = t + deltat;
        }

    } while (t <= t_integration);

    printf("\n\nDEBUG [tIntegrate]: Main simulation loop finished (t > t_integration).\n");
    printf("A program sikeresen lefutott, azonban elkepzelheto, hogy az integralasi ido nem volt elegendo. A legtavolabbi reszecske %lg CsE tavolsagra van a kozponti csillagtol. \n\nNyomj ENTER-t a kilepeshez!\n", max);

cleanup:
    // --- Tisztítási szakasz ---
    cleanup_simulation_resources(&p_data, output_files, sim_opts);
    printf("DEBUG [tIntegrate]: Cleanup completed.\n");
}

