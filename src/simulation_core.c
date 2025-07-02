
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
    // printf("DEBUG [eqrhs]: St=%.2e, H=%.2e, P=%.2e, dPdr=%.2e, csound=%.2e, drdt=%.2e\n", St, H, P, dPdr, csound, *drdt);
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

    fprintf(stderr, "DEBUG [tIntegrate]: Function entry. disk_params address=%p, sim_opts address=%p\n",
            (void*)disk_params, (void*)sim_opts);
    fprintf(stderr, "DEBUG [tIntegrate]:   sim_opts->output_dir_name = '%s'\n", sim_opts->output_dir_name);
    fprintf(stderr, "DEBUG [tIntegrate]:   disk_params->FLIND = %.2f, disk_params->HASP = %.2f\n",
            disk_params->FLIND, disk_params->HASP);

    // A globális PARTICLE_NUMBER változó beállítása a szimulációs opciók alapján
    if(sim_opts->drift == 1.) {
        fprintf(stderr, "DEBUG [tIntegrate]: Particle drift is ON. Counting particles from '%s'.\n", sim_opts->dust_input_filename);
        PARTICLE_NUMBER = reszecskek_szama(sim_opts->dust_input_filename);
        fprintf(stderr, "DEBUG [tIntegrate]: PARTICLE_NUMBER set to %d.\n", PARTICLE_NUMBER);
    } else {
        fprintf(stderr, "DEBUG [tIntegrate]: Particle drift is OFF. PARTICLE_NUMBER set to 0.\n");
        PARTICLE_NUMBER = 0;
    }

    // --- DINAMIKUS MEMÓRIAFOGLALÁS A RÉSZECSKE ADATOKNAK ---
    double (*radius)[2] = NULL;
    double (*radiusmicr)[2] = NULL;
    double (*radiussec)[2] = NULL;
    double (*radius_rec)[2] = NULL; // Temp array for inverse radii calculations
    double *massvec = NULL;
    double *massmicrvec = NULL;
    double *masssecvec = NULL;
    double (*partmassind)[4] = NULL;
    double (*partmassmicrind)[4] = NULL;
    double (*partmasssecind)[4] = NULL;
    double *sigmad = NULL;
    double *sigmadm = NULL;
    double *sigmads = NULL; // Surface density for secondary particles
    double *rdvec = NULL;
    double *rmicvec = NULL;
    double *rsvec = NULL; // Radii for secondary particle surface density

    // Csak akkor foglalunk memóriát, ha van részecskeszám
    if (PARTICLE_NUMBER > 0) {
        radius = malloc(PARTICLE_NUMBER * sizeof(*radius));
        radiusmicr = malloc(PARTICLE_NUMBER * sizeof(*radiusmicr));
        radius_rec = malloc(PARTICLE_NUMBER * sizeof(*radius_rec));
        massvec = malloc(PARTICLE_NUMBER * sizeof(double));
        massmicrvec = malloc(PARTICLE_NUMBER * sizeof(double));
        partmassind = malloc(PARTICLE_NUMBER * sizeof(*partmassind));
        partmassmicrind = malloc(PARTICLE_NUMBER * sizeof(*partmassmicrind));
        sigmad = malloc(PARTICLE_NUMBER * sizeof(double));
        sigmadm = malloc(PARTICLE_NUMBER * sizeof(double));
        rdvec = malloc(PARTICLE_NUMBER * sizeof(double));
        rmicvec = malloc(PARTICLE_NUMBER * sizeof(double));

        // A secondary particles tömbök mérete 4*PARTICLE_NUMBER (ha two-pop van és növekedés)
        // A malloc hívások mérete a num_particles paramétertől függ a Get_Radius-ban
        radiussec = malloc(4 * PARTICLE_NUMBER * sizeof(*radiussec));
        masssecvec = malloc(4 * PARTICLE_NUMBER * sizeof(double));
        partmasssecind = malloc(4 * PARTICLE_NUMBER * sizeof(*partmasssecind));
        sigmads = malloc(4 * PARTICLE_NUMBER * sizeof(double));
        rsvec = malloc(4 * PARTICLE_NUMBER * sizeof(double));


        // Hibaellenőrzés minden malloc hívás után
        if (!radius || !radiusmicr || !radius_rec || !massvec || !massmicrvec ||
            !partmassind || !partmassmicrind || !sigmad || !sigmadm || !rdvec || !rmicvec ||
            !radiussec || !masssecvec || !partmasssecind || !sigmads || !rsvec) {
            fprintf(stderr, "ERROR: Memory allocation failed in tIntegrate! Exiting.\n");
            exit(EXIT_FAILURE);
        }
        fprintf(stderr, "DEBUG [tIntegrate]: Particle arrays dynamically allocated. Size based on PARTICLE_NUMBER: %d\n", PARTICLE_NUMBER);
    } else {
        fprintf(stderr, "DEBUG [tIntegrate]: PARTICLE_NUMBER is 0. No particle arrays allocated.\n");
    }

    double max = 0.0, min = 0.0, max2 = 0.0, min2 = 0.0;
    int i; // Az i változó deklarációja

    // Karaktertömbök inicializálása üres sztringgel a biztonság kedvéért
    char porout[MAX_PATH_LEN] = "";
    char poroutmicr[MAX_PATH_LEN] = ""; // Mikronos por mozgás fájlneve
    char massout[MAX_PATH_LEN] = "";
    char dens_name[MAX_PATH_LEN] = ""; // Gáz felületi sűrűség fájlnév
    char dust_name[MAX_PATH_LEN] = ""; // Fő por felületi sűrűség fájlnév
    char dust_name2[MAX_PATH_LEN] = ""; // Mikronos por felületi sűrűség fájlnév
    char size_name[MAX_PATH_LEN] = ""; // Részecskeméret fájlnév

/* A részecskék adatainak beolvasása fájlból, ha a drift opció értéke 1, egyébként nem számol a program driftet */
    if(sim_opts->drift == 1.) {
        fprintf(stderr, "DEBUG [tIntegrate]: Reading particle data from file '%s'.\n", sim_opts->dust_input_filename);
        por_be(radius, radiusmicr, massvec, massmicrvec, sim_opts->dust_input_filename);    /* porrészecskék adatainak beolvasása */

        // --- ELTÁVOLÍTVA: A 'cp' parancs, ami duplikálta az initial_dust_profile.dat fájlt ---
        // Ez a másolás szükségtelen, mivel az init_tool már a megfelelő helyre generálja a fájlt.
        // A main.c-ben lévő cp parancsok az 'initial' könyvtárba való archiválást végzik.
        // --- ELTÁVOLÍTÁS VÉGE ---


/* az aktuális mappában a dust_particle_evolution.dat fájl létrehozása: ebbe kerül be a porrészecske távolsága és if(sim_opts->drift == 1.)dexe, valamint az adott időlépés */
        snprintf(porout,MAX_PATH_LEN,"%s/%s/dust_particle_evolution.dat",sim_opts->output_dir_name,LOGS_DIR);


/* ha 2 populációs a futás, akkor a mikronos pornak is létrehoz egy pormozgás fájlt, ebbe kerül be a távolság, index és az idő */
        if(sim_opts->twopop == 1.) {
            // --- JAVÍTVA: A FILE* pointer helyett a konstans fájlnevet használjuk ---
            snprintf(poroutmicr,MAX_PATH_LEN,"%s/%s/micron_particle_evolution.dat",sim_opts->output_dir_name, LOGS_DIR);
            // --- JAVÍTÁS VÉGE ---
        }
/* tömegnövekedési fájl létrehozása az aktuális mappába - ez lehet, hogy egy külön opció lesz a kimeneti adatok méretének csökkentésére */
        
        // --- MÓDOSÍTÁS: mass_accumulation_dze_edge.dat a LOGS mappába ---
        snprintf(massout,MAX_PATH_LEN,"%s/%s/mass_accumulation_dze_edge.dat",sim_opts->output_dir_name,LOGS_DIR);
        // --- MÓDOSÍTÁS VÉGE ---

        fprintf(stderr, "DEBUG [tIntegrate]: Opening output files: %s, %s (if 2pop), %s\n", porout, poroutmicr, massout);

        output_files->por_motion_file = fopen(porout,"w");
        if (output_files->por_motion_file == NULL) { fprintf(stderr, "ERROR: Could not open %s\n", porout); /*exit(EXIT_FAILURE);*/ }
        if (output_files->por_motion_file != NULL) {
            fprintf(output_files->por_motion_file, "# Main Dust Particle Motion and Properties\n");
            fprintf(output_files->por_motion_file, "# Columns: 1. Time Step, 2. Particle Index, 3. Radius [AU]\n");
            fprintf(output_files->por_motion_file, "# Data generated by Dust Drift Simulation\n");
            fflush(output_files->por_motion_file);
        }

        if(sim_opts->twopop == 1.) {
            output_files->micron_motion_file = fopen(poroutmicr,"w");
            if (output_files->micron_motion_file == NULL) { fprintf(stderr, "ERROR: Could not open %s\n", poroutmicr); /*exit(EXIT_FAILURE);*/ }
            if (output_files->micron_motion_file != NULL) {
                fprintf(output_files->micron_motion_file, "# Micron Dust Particle Motion and Properties\n");
                fprintf(output_files->micron_motion_file, "# Columns: 1. Time Step, 2. Particle Index, 3. Radius [AU]\n");
                fprintf(output_files->micron_motion_file, "# Data generated by Dust Drift Simulation (Two-Population Model)\n");
                fflush(output_files->micron_motion_file);
            }
        }

        output_files->mass_file = fopen(massout,"w");
        if (output_files->mass_file == NULL) { fprintf(stderr, "ERROR: Could not open %s\n", massout); /*exit(EXIT_FAILURE);*/ }

        // Fejléc írása a mass_accumulation_dze_edge.dat fájlba
        if (output_files->mass_file != NULL) {
            fprintf(output_files->mass_file, "# This file contains the time evolution of dust mass within specified disk regions.\n");
            fprintf(output_files->mass_file, "#\n");
            fprintf(output_files->mass_file, "# Columns:\n");
            fprintf(output_files->mass_file, "# 1. Time Step (arbitrary unit, iteration count)\n");
            fprintf(output_files->mass_file, "# 2. Inner Boundary Radius (R_in) of the inner region [AU]\n");
            fprintf(output_files->mass_file, "# 3. Total Dust Mass (non-micron + micron) within R < R_in [M_Jupiter]\n"); // JAVÍTVA: jelfut helyett output_files->mass_file
            fprintf(output_files->mass_file, "# 4. Outer Boundary Radius (R_out) of the outer region [AU]\n");
            fprintf(output_files->mass_file, "# 5. Total Dust Mass (non-micron + micron) within R > R_out [M_Jupiter]\n");
            fprintf(output_files->mass_file, "#\n");
            fprintf(output_files->mass_file, "# All masses include both 'cm-sized' (larger) and 'micron-sized' dust populations.\n");
            fprintf(output_files->mass_file, "# Data is updated periodically during the simulation.\n");
            fflush(output_files->mass_file);
        }
    }
    double t = 0.0;
    double t_integration = sim_opts->TMAX * 2.0 * M_PI;  // TMAX a sim_opts-ból
    double deltat = time_step(disk_params) / 5.0; // rvec a disk_params-ból.

    // Döntés az időlépésről:
    // Ha a fixed_time_step (sim_opts->DT) pozitív ÉS kisebb, mint a számított deltat,
    // AKKOR a sim_opts->DT-t használjuk.
    if(sim_opts->DT > 0.0 && sim_opts->DT < deltat) {
        // sim_opts->DT már eleve a helyes, kisebb, pozitív értéket tartalmazza a configból.
        // Nincs szükség felülírásra, a deltat változó itt nem releváns a továbbiakban.
    } else {
        // Ha a fixed_time_step 0, vagy nagyobb/egyenlő a számított deltat-tal,
        // AKKOR a számított deltat-ot használjuk, és BEÁLLÍTJUK a sim_opts->DT-t erre az értékre.
        ((simulation_options_t *)sim_opts)->DT = deltat; // Itt frissítjük a sim_opts->DT-t!
    }


    double L = 0.;
    double masstempiin = 0, massmtempiin = 0, masstempoin = 0, massmtempoin = 0; // a külső és belső dze-n felgyülemlett por mennyisége -- bemeneti adat (Print_Mass függvénybe)
    double masstempiout = 0, masstempoout = 0, massmtempiout = 0, massmtempoout = 0; // a külső és belső dze-n felgyülemlett por mennyisége -- kimeneti adat (Print_Mass függvényből)

    double tavin, tavout;


    // A secondary particles inicializálása
    if (PARTICLE_NUMBER > 0) {
        fprintf(stderr, "DEBUG [tIntegrate]: Initializing secondary particle arrays (size 4*PARTICLE_NUMBER = %d).\n", 4 * PARTICLE_NUMBER);
        for(i = 0; i < 4*PARTICLE_NUMBER; i++) {
            radiussec[i][0] = 0;
            radiussec[i][1] = 0;
            partmasssecind[i][0] = 0;
            partmasssecind[i][1] = 0;
            masssecvec[i] = 0;
            sigmads[i] = 0;
            rsvec[i] = 0;
        }

        // Ha nincs két populáció, a mikronos részecsketömbök nullázása
        if(sim_opts->twopop == 0) {
            fprintf(stderr, "DEBUG [tIntegrate]: Two-population simulation is OFF. Initializing micron particle arrays to zero.\n");
            for(i = 0; i < PARTICLE_NUMBER; i++) {
                radiusmicr[i][0] = 0;
                radiusmicr[i][1] = 0;
                partmassmicrind[i][0] = 0;
                partmassmicrind[i][1] = 0;
                massmicrvec[i] = 0;
            }
        }
    }


    printf("DEBUG [tIntegrate]: Entering main simulation time loop (t=%.2e, t_integration=%.2e).\n", t, t_integration);
    do {
/* Ha van drift: */
        if(sim_opts->drift == 1.) {
            if(sim_opts->twopop == 1) {
                fprintf(stderr, "DEBUG [tIntegrate]: sim_opts->twopop is ON. Calling secondaryGrowth.\n");
                secondaryGrowth(radius,radiusmicr,radiussec,partmassmicrind,partmasssecind,massmicrvec,masssecvec,disk_params);
                fprintf(stderr, "DEBUG [tIntegrate]: secondaryGrowth completed.\n");
            }
/* A minimum kereséshez létrehozza a cm-es részecskék távolságának reciprokát */
            for (i=0; i < PARTICLE_NUMBER; i++) {
                if (radius[i][0] > 0. && radius[i][0] > disk_params->RMIN) { // RMIN a disk_params-ból
                    radius_rec[i][0] = 1. / radius[i][0];
                } else {
                    radius_rec[i][0] = 0.;
                }
            }

            max = find_max(radius,PARTICLE_NUMBER);     /* Megkeresi, hogy melyik a legtávolabbi cm-es részecske a központi csillagtól */
            min = find_max(radius_rec,PARTICLE_NUMBER);     /* Megkeresi a távolság reciprokának maximumát, azaz a legkisebb távolságra lévő cm-es részecskét */
            min = 1. / min;

            double mint, maxt;

/* ha 2 populációs a szimuláció, a fentihez hasonlóan megkeresi a legnagyobb és a legkisebb távolságra lévő mikronos részecskét */
            if(sim_opts->twopop == 1) {
                for (i=0; i < PARTICLE_NUMBER; i++) {
                    if (radiusmicr[i][0] > 0. && radiusmicr[i][0] > disk_params->RMIN) { // RMIN a disk_params-ból
                        radius_rec[i][0] = 1. / radiusmicr[i][0];
                    } else {
                        radius_rec[i][0] = 0.;
                    }
                }

                max2 = find_max(radiusmicr,PARTICLE_NUMBER);    /* Megkeresi, hogy melyik a legtávolabbi részecske a központi csillagtól */
                min2 = find_max(radius_rec,PARTICLE_NUMBER);
                min2 = 1. / min2;

/* megnézi, hogy mely részecske van a legközelebb, illetve legtávolabb (mikronos, vagy cm-es) */
                mint = find_min(min,min2,HUGE_VAL);
                maxt = find_min(1. / max, 1./max2,HUGE_VAL);
                maxt = 1./ maxt;
            } else {
                mint = min;
                maxt = max;
            }

/* Ha a legtávolabbi részecske távolsága nagyobb, mint RMIN (azaz még a szimuláció tartományában van), és a min és a max nem egyenlő (azaz nem gyűlt pl. össze az összes porrészecske 1 helyen), akkor a program tovább számol, egyébként a futás leáll -- ez persze 1 dze esetén működik, meg kell oldani, hogy 2 dze esetén is leálljon akkor, ha az összes por összegyűlt a nyomási maximumokban - meg kell persze csinálni, hogy ez is opcionális legyen */
            if(maxt >= disk_params->RMIN && mint != maxt) { // RMIN a disk_params-ból
                double time = t / 2.0 / M_PI;

                if((fmod(time, (sim_opts->TMAX/sim_opts->WO)) < deltat || time == 0) && L-time < deltat){ // TMAX és WO a sim_opts-ból
                    printf("\n--- Simulation Time: %.2e years (Internal time: %.2e, L: %.2e) ---\n", time, t, L);

                    fprintf(stderr, "DEBUG [tIntegrate]: Outputting data at time %.2e. L=%.2e\n", time, L);

/* Az adatok kiírásához szükséges fájlok neveinek elmentése */
                    if (t==0) {
                        // --- MÓDOSÍTÁS: dens_name az initial surface density fájlra mutat CONFIG_DIR-ban ---
//                        snprintf(dens_name,MAX_PATH_LEN,"%s/%s/%s",sim_opts->output_dir_name,CONFIG_DIR,INITIAL_SURFACE_DENSITY_FILE);
                        // --- MÓDOSÍTÁS VÉGE ---

//                        fprintf(stderr, "DEBUG [tIntegrate]: Outputting initial gas surface density for t=0 to %s.\n", dens_name);
                    } else {
                        if(sim_opts->evol == 1) {
                            // --- MÓDOSÍTÁS: dens.X.dat a LOGS mappába (8 vezető nullával) ---
                            snprintf(dens_name,MAX_PATH_LEN,"%s/%s/density_profile_%08d.dat",sim_opts->output_dir_name,LOGS_DIR,(int)L);
                            // --- MÓDOSÍTÁS VÉGE ---
                            fprintf(stderr, "DEBUG [tIntegrate]: Outputting density_profile_%08d.dat to %s.\n", (int)L, dens_name);
                        }
                    }

                    // --- MÓDOSÍTÁS: dust.X.dat és dustmic.X.dat a LOGS mappába ---
                    snprintf(dust_name,MAX_PATH_LEN,"%s/%s/dust.%i.dat",sim_opts->output_dir_name,LOGS_DIR,(int)L);
                    snprintf(dust_name2,MAX_PATH_LEN,"%s/%s/dustmic.%i.dat",sim_opts->output_dir_name,LOGS_DIR,(int)L);
                    // --- MÓDOSÍTÁS VÉGE ---

                    // --- MÓDOSÍTÁS: size.X.dat a LOGS mappába ---
                    snprintf(size_name,MAX_PATH_LEN,"%s/%s/size.%d.dat",sim_opts->output_dir_name,LOGS_DIR,(int)L);
                    // --- MÓDOSÍTÁS VÉGE ---

                    fprintf(stderr, "DEBUG [tIntegrate]: Output file names set: %s, %s, %s.\n", dust_name, dust_name2, size_name);

                    // *******************************************************************
                    // FÁJLOK MEGNYITÁSA AZ ADATKIÍRÁSHOZ
                    // *******************************************************************
                    output_files->surface_file = fopen(dens_name, "w");
                    if (output_files->surface_file == NULL) {
                        fprintf(stderr, "ERROR: Could not open %s for writing.\n", dens_name);
                        // exit(EXIT_FAILURE); // Fontos lehet!
                    } else {
                        fprintf(stderr, "DEBUG [tIntegrate]: Opened %s for writing.\n", dens_name);
                        // FEJLÉC Hozzáadása ide
                        fprintf(output_files->surface_file, "# Gas Surface Density Data at t = %lg years\n",time);
                        fprintf(output_files->surface_file, "# Generated by init_tool_module (Date: %s %s)\n", __DATE__, __TIME__);
                        fprintf(output_files->surface_file, "#--------------------------------------------------------------------------\n");
                        fprintf(output_files->surface_file, "# %-15s %-15s %-15s %-15s\n",
             "Radius_AU", "GasSurfDensity", "GasPressure", "GasPressureDeriv");
                        fprintf(output_files->surface_file, "#--------------------------------------------------------------------------\n");
                        fflush(output_files->surface_file); // Frissítsd a fájlt azonnal
                    }

                    output_files->dust_file = fopen(dust_name, "w");
                    if (output_files->dust_file == NULL) {
                        fprintf(stderr, "ERROR: Could not open %s for writing.\n", dust_name);
                        // exit(EXIT_FAILURE); // Fontos lehet!
                    } else {
                        fprintf(stderr, "DEBUG [tIntegrate]: Opened %s for writing.\n", dust_name);
                    }

                    if(sim_opts->twopop == 1.) {
                        output_files->micron_dust_file = fopen(dust_name2, "w");
                        if (output_files->micron_dust_file == NULL) {
                            fprintf(stderr, "ERROR: Could not open %s for writing.\n", dust_name2);
                            // exit(EXIT_FAILURE); // Fontos lehet!
                        } else {
                            fprintf(stderr, "DEBUG [tIntegrate]: Opened %s for writing.\n", dust_name2);
                        }
                    }

                    // --- ELTÁVOLÍTVA: output_files->size_file megnyitása itt, mivel Print_Pormozg_Size nyitja meg ---
                    fprintf(stderr, "DEBUG [tIntegrate]: Removed redundant fopen for size file. Print_Pormozg_Size handles it.\n");
                    // *******************************************************************
                    // FÁJL MEGNYITÁS VÉGE
                    // *******************************************************************


                    if(t==0) {
//                        fprintf(stderr, "DEBUG [tIntegrate]: Initializing mass for t=0.\n");
/* A részecskék tömegét tartalmazó tömb inicializálása */
                        Count_Mass(radius,partmassind,massvec,t,PARTICLE_NUMBER,disk_params);
                        if(sim_opts->twopop==1) Count_Mass(radiusmicr,partmassmicrind,massmicrvec,t,PARTICLE_NUMBER,disk_params);
                        if(sim_opts->twopop==1) Count_Mass(radiussec,partmasssecind,masssecvec,t,4*PARTICLE_NUMBER,disk_params);
//                        fprintf(stderr, "DEBUG [tIntegrate]: Count_Mass completed for t=0.\n");

/* Ha van tömegnövekedés, akkor a por felületsűrűségének kiszámolása itt történik */
                        if(sim_opts->growth == 1.) {
//                            fprintf(stderr, "DEBUG [tIntegrate]: sim_opts->growth is ON. Calling Get_Sigmad for t=0.\n");
                            // JAVÍTVA: A Get_Sigmad függvény paraméterezése
                            Get_Sigmad(max, min, radius, radiusmicr, radiussec, sigmad, sigmadm, sigmads, massvec, massmicrvec, masssecvec, rdvec, rmicvec, rsvec, sim_opts, disk_params);
//                             fprintf(stderr, "DEBUG [tIntegrate]: Get_Sigmad completed for t=0.\n");
                        }
                    }

/* A sigma, p, dp kiírása egy fájlba */
                    if(sim_opts->evol == 1 || time == 0) {
                        Print_Sigma(disk_params, output_files);
                    } 

/* Ha számol a futás driftet, itt írja ki a részecskék távolságát és méretét */
                    if(sim_opts->drift == 1) {
                        fprintf(stderr, "DEBUG [tIntegrate]: sim_opts->drift is ON. Calling Print_Pormozg_Size.\n");
                        Print_Pormozg_Size(size_name, (int)L, radius, radiusmicr, disk_params, sim_opts, output_files);
                        fprintf(stderr, "DEBUG [tIntegrate]: Print_Pormozg_Size completed.\n");
                    }

/* A tömegnövekedési fájlba az adatok kiírása */
                    masstempiout = 0;
                    massmtempiout = 0;
                    masstempoout = 0;
                    massmtempoout = 0;
                    fprintf(stderr, "DEBUG [tIntegrate]: Calling Print_Mass.\n");
                    Print_Mass(L, disk_params->rvec, partmassind, partmassmicrind, partmasssecind, disk_params->dpressvec, masstempiin, masstempoin,massmtempiin, massmtempoin, &masstempiout, &masstempoout, &massmtempiout, &massmtempoout,&tavin, &tavout, disk_params, sim_opts, output_files);
                    fprintf(stderr, "DEBUG [tIntegrate]: Print_Mass completed. Outputs: masstempiout=%.2e, massmtempiout=%.2e\n", masstempiout, massmtempiout);

                    if(disk_params->r_dze_i != tavin) { // r_dze_i a disk_params-ból
                        masstempiin = masstempiout;
                        massmtempiin = massmtempiout;
                        fprintf(stderr, "DEBUG [tIntegrate]: r_dze_i (%.2f) != tavin (%.2f). Updating masstempiin/massmtempiin.\n", disk_params->r_dze_i, tavin);
                    }
                    if(disk_params->r_dze_o != tavout) { // r_dze_o a disk_params-ból
                        masstempoin = masstempoout;
                        massmtempoin = massmtempoout;
                        fprintf(stderr, "DEBUG [tIntegrate]: r_dze_o (%.2f) != tavout (%.2f). Updating masstempoin/massmtempoin.\n", disk_params->r_dze_o, tavout);
                    }

/* Ha van pornövekedés, kiírja a por felületsűrűségét egy fájlba --> a pornövekedéshez szükséges egyáltalán ezt kiszámolni! */
                    if(sim_opts->growth == 1.) {
                        fprintf(stderr, "DEBUG [tIntegrate]: sim_opts->growth is ON. Calling Print_Sigmad (refreshing).\n");
                        // JAVÍTVA: A Print_Sigmad függvény paraméterezése
                        Print_Sigmad(rdvec, rmicvec, sigmad, sigmadm, disk_params, sim_opts, output_files);
                        fprintf(stderr, "DEBUG [tIntegrate]: Print_Sigmad (refreshing) completed.\n");
                    }

                    L = L+(double)(sim_opts->TMAX/sim_opts->WO); // TMAX és WO a sim_opts-ból
                    printf("DEBUG [tIntegrate]: Updated L to %.2e.\n", L);

                    // *******************************************************************
                    // FÁJLOK BEZÁRÁSA AZ ADATKIÍRÁS UTÁN
                    // *******************************************************************
                    if (output_files->surface_file != NULL) { fclose(output_files->surface_file); output_files->surface_file = NULL; fprintf(stderr, "DEBUG [tIntegrate]: Closed %s.\n", dens_name); }
                    if (output_files->dust_file != NULL) { fclose(output_files->dust_file); output_files->dust_file = NULL; fprintf(stderr, "DEBUG [tIntegrate]: Closed %s.\n", dust_name); }
                    if (sim_opts->twopop == 1 && output_files->micron_dust_file != NULL) { fclose(output_files->micron_dust_file); output_files->micron_dust_file = NULL; fprintf(stderr, "DEBUG [tIntegrate]: Closed %s.\n", dust_name2); }
                    // --- ELTÁVOLÍTVA: output_files->size_file bezárása itt, mivel Print_Pormozg_Size kezeli ---
                    fprintf(stderr, "DEBUG [tIntegrate]: Removed redundant fclose for size file. Print_Pormozg_Size handles it.\n");
                    // *******************************************************************
                    // FÁJL BEZÁRÁS VÉGE
                    // *******************************************************************

                }

/* Ha az evol opció értéke 1, akkor megoldja minden lépésben a sigmára vonatkozó diffúziós egyenletet */
                if(sim_opts->evol == 1.) {
                    Get_Sigma_P_dP(sim_opts, disk_params); // Átadja a sim_opts és a disk_params címét
                }

/* A részecskék tömegét tartalmazó tömb adatainak frissítése */
                Count_Mass(radius,partmassind,massvec,t,PARTICLE_NUMBER,disk_params);
                if(sim_opts->twopop==1) Count_Mass(radiusmicr,partmassmicrind,massmicrvec,t,PARTICLE_NUMBER,disk_params);
                if(sim_opts->twopop==1) Count_Mass(radiussec,partmasssecind,masssecvec,t,4*PARTICLE_NUMBER,disk_params);

/* Ha van részecskenövekedés, akkor kiszámolja a por felületsűrűségét */
                if(sim_opts->growth == 1.) {
                    // JAVÍTVA: A Get_Sigmad függvény paraméterezése
                    Get_Sigmad(max, min, radius, radiusmicr, radiussec, sigmad, sigmadm, sigmads, massvec, massmicrvec, masssecvec, rdvec, rmicvec, rsvec, sim_opts, disk_params);
                }

                int optsize = 0;
/* A cm-es részecskék esetén az optsize értéke 0 */
/* Itt számolja ki a program a cm-es részecskék új távolságát (és méretét, ha kell) */
                Get_Radius(sim_opts->output_dir_name,optsize,radius,sigmad,rdvec,deltat,t,PARTICLE_NUMBER,sim_opts,disk_params); // output_dir_name a sim_opts-ból

/* Ha a futás 2 populációs, akkor az optsize értéke 1 */
/* Itt számolja ki a program a mikronos részecskék új távolságát */
                if(sim_opts->twopop == 1.) {
                    optsize = 1;
                    printf("DEBUG [tIntegrate]: sim_opts->twopop is ON. Calling Get_Radius for Micron particles (optsize=%d, PARTICLE_NUMBER=%d).\n", optsize, PARTICLE_NUMBER);
                    // JAVÍTVA: radiusmicr és sigmadm használata
                    Get_Radius(sim_opts->output_dir_name,optsize,radiusmicr,sigmadm,rmicvec,deltat,t,PARTICLE_NUMBER,sim_opts,disk_params);
                    printf("DEBUG [tIntegrate]: Get_Radius for Micron particles completed.\n");
                    optsize = 2;
                    printf("DEBUG [tIntegrate]: Calling Get_Radius for Secondary particles (optsize=%d, PARTICLE_NUMBER=%d).\n", optsize, 4 * PARTICLE_NUMBER);
                    // JAVÍTVA: radiussec és sigmads használata
                    Get_Radius(sim_opts->output_dir_name,optsize,radiussec,sigmads,rsvec,deltat,t,4 * PARTICLE_NUMBER,sim_opts,disk_params);
                    printf("DEBUG [tIntegrate]: Get_Radius for Secondary particles completed.\n");
                }

                t = t + deltat;     /* Időléptetés */
                printf("DEBUG [tIntegrate]: Time advanced to t=%.2e.\n", t);

            } else {    /* Ha a legmesszebbi részecske távolsága már nem nagyobb, vagy egyenlő, mint RMIN, vagy a legkisebb távolságra lévő részecske távolsága, akkor a program "figyelmeztető szöveg" mellett sikeresen kilép, nem fut "feleslegesen" tovább. */
                printf("DEBUG [tIntegrate]: Simulation termination condition met (maxt < RMIN or mint == maxt).\n");
                printf("A program sikeresen lefutott az integralasi ido vege elott (t: %lg). \n\nNyomj ENTER-t a kilepeshez!\n",L);
                exit(EXIT_SUCCESS);
            }

        } else { /* Ez az az eset, ha a program nem számol driftet, azaz csak a gáz felületsűrűségének fejlődésére vagyunk kíváncsiak */
            double time = t / 2.0 / M_PI;

            if((fmod(time, (sim_opts->TMAX/sim_opts->WO)) < deltat || time == 0) && L-time < deltat){ // TMAX és WO a sim_opts-ból
                printf("\n--- Simulation Time: %.2e years (Internal time: %.2e, L: %.2e) ---\n", time, t, L);

                printf("DEBUG [tIntegrate]: Outputting data for gas-only simulation at time %.2e. L=%.2e\n", time, L);
                if (t==0) {
                    // --- MÓDOSÍTÁS: dens_name az initial surface density fájlra mutat CONFIG_DIR-ban (gáz-only ág) ---
//                    snprintf(dens_name,MAX_PATH_LEN,"%s/%s/%s",sim_opts->output_dir_name,CONFIG_DIR,INITIAL_SURFACE_DENSITY_FILE);
                    // --- MÓDOSÍTÁS VÉGE ---
//                    fprintf(stderr, "DEBUG [tIntegrate]: Outputting initial gas surface density for t=0 to %s.\n", dens_name);
                } else {
                    // --- MÓDOSÍTÁS: dens.X.dat a LOGS mappába (gáz-only ág, 8 vezető nullával) ---
                    snprintf(dens_name,MAX_PATH_LEN,"%s/%s/density_profile_%08d.dat",sim_opts->output_dir_name,LOGS_DIR,(int)L);
                    // --- MÓDOSÍTÁS VÉGE ---
                    fprintf(stderr, "DEBUG [tIntegrate]: Outputting density_profile_%08d.dat to %s.\n", (int)L, dens_name);
                }

                output_files->surface_file = fopen(dens_name, "w"); // Fájl megnyitása a Print_Sigma hívása előtt
                if (output_files->surface_file == NULL) {
                    fprintf(stderr, "ERROR: Could not open %s for writing in gas-only branch.\n", dens_name);
                } else {
                    fprintf(stderr, "DEBUG [tIntegrate]: Opened %s for writing in gas-only branch.\n", dens_name);
                    // FEJLÉC Hozzáadása ide
                    fprintf(output_files->surface_file, "# Gas Surface Density Data at t = %lg years\n",time);
                    fprintf(output_files->surface_file, "# Generated by init_tool_module (Date: %s %s)\n", __DATE__, __TIME__);
                    fprintf(output_files->surface_file, "#--------------------------------------------------------------------------\n");
                    fprintf(output_files->surface_file, "# %-15s %-15s %-15s %-15s\n",
             "Radius_AU", "GasSurfDensity", "GasPressure", "GasPressureDeriv");
                    fprintf(output_files->surface_file, "#--------------------------------------------------------------------------\n");
                    fflush(output_files->surface_file); // Frissítsd a fájlt azonnal
                }

                Print_Sigma(disk_params, output_files);
                printf("DEBUG [tIntegrate]: Print_Sigma completed.\n");
                
                if (output_files->surface_file != NULL) { // Fájl bezárása a Print_Sigma hívása után
                    fclose(output_files->surface_file);
                    output_files->surface_file = NULL;
                    fprintf(stderr, "DEBUG [tIntegrate]: Closed %s in gas-only branch.\n", dens_name);
                }

                L = L+(double)(sim_opts->TMAX/sim_opts->WO); // TMAX és WO a sim_opts-ból
                printf("DEBUG [tIntegrate]: Updated L to %.2e.\n", L);
            }

            printf("DEBUG [tIntegrate]: Calling Get_Sigma_P_dP for gas-only evolution.\n");
            Get_Sigma_P_dP(sim_opts, disk_params);
            printf("DEBUG [tIntegrate]: Get_Sigma_P_dP completed.\n");
            t = t + deltat;     /* Időléptetés */
            printf("DEBUG [tIntegrate]: Time advanced to t=%.2e.\n", t);
        }

    } while (t <= t_integration);


/* Az időléptetés leteltével a program sikeresen kilép */
    // Memória felszabadítása
    if (PARTICLE_NUMBER > 0) {
        free(radius);
        free(radiusmicr);
        free(radius_rec);
        free(massvec);
        free(massmicrvec);
        free(partmassind);
        free(partmassmicrind);
        free(sigmad);
        free(sigmadm);
        free(rdvec);
        free(rmicvec);

        free(radiussec);
        free(masssecvec);
        free(partmasssecind);
        free(sigmads);
        free(rsvec);

        fprintf(stderr, "DEBUG [tIntegrate]: All dynamically allocated particle arrays freed.\n");
    }

    // Fájlok bezárása, amelyek az integrációs ciklus elején nyíltak meg (egyszer).
    if (output_files->por_motion_file != NULL) {
        fclose(output_files->por_motion_file);
        fprintf(stderr, "DEBUG [tIntegrate]: Closed dust_particle_evolution.dat\n");
    }
    if (sim_opts->twopop == 1 && output_files->micron_motion_file != NULL) {
        fclose(output_files->micron_motion_file);
        fprintf(stderr, "DEBUG [tIntegrate]: Closed micron_particle_evolution.dat\n");
    }
    if (output_files->mass_file != NULL) {
        fclose(output_files->mass_file);
        fprintf(stderr, "DEBUG [tIntegrate]: Closed mass_accumulation_dze_edge.dat\n");
    }
    // Az itt bezárandó fájlok listájából KIVETTÜK a surface/dens, dust, size fájlokat,
    // mivel azok már be lettek zárva minden iteráció végén a ciklusban, amikor az adatok kiírásra kerültek.

    printf("\n\nDEBUG [tIntegrate]: Main simulation loop finished (t > t_integration).\n");
    printf("A program sikeresen lefutott, azonban elkepzelheto, hogy az integralasi ido nem volt elegendo. A legtavolabbi reszecske %lg CsE tavolsagra van a kozponti csillagtol. \n\nNyomj ENTER-t a kilepeshez!\n", max); // Visszaállítottam a 'max' változó használatát
}

void secondaryGrowth(double rad[][2], double radmicr[][2], double radsec[][2], double partmicind[][4], double partsecind[][4], double *massmicvec, double *masssecvec, const disk_t *disk_params) {
//    printf("DEBUG [secondaryGrowth]: Entering secondaryGrowth function.\n");
//    fprintf(stderr, "DEBUG [secondaryGrowth]: Function entry. disk_params address=%p\n", (void*)disk_params);
//    fprintf(stderr, "DEBUG [secondaryGrowth]:   disk_params->FLIND = %.2f, disk_params->HASP = %.2f\n",
//            disk_params->FLIND, disk_params->HASP);

    int i, j; //, k=0;
    int step = 0;
    double temp[PARTICLE_NUMBER], temp2[4*PARTICLE_NUMBER], tempmic[PARTICLE_NUMBER], temp2ch1[PARTICLE_NUMBER], temp2ch2[PARTICLE_NUMBER], temp2ch3[PARTICLE_NUMBER], temp2ch4[PARTICLE_NUMBER];
    int hist[PARTICLE_NUMBER], histmic[PARTICLE_NUMBER];
//    printf("DEBUG [secondaryGrowth]: Declared local arrays (hist, temp, etc.) based on PARTICLE_NUMBER: %d.\n", PARTICLE_NUMBER);

/*	A porszemcse generalo program gridfelbontasa -- nem feltetlen egyezik meg a jelenlegi gridfelbontassal!	*/
    double dd = (disk_params->RMAX - disk_params->RMIN) / (PARTICLE_NUMBER-1), rtempvec[PARTICLE_NUMBER];
    double temptemp[4*PARTICLE_NUMBER], masstemp[4*PARTICLE_NUMBER][2], sizetemp[4*PARTICLE_NUMBER],mtp[4*PARTICLE_NUMBER], indtemp[4*PARTICLE_NUMBER][2];
//    printf("DEBUG [secondaryGrowth]: Calculated dd=%.2e based on RMAX=%.2f, RMIN=%.2f, PARTICLE_NUMBER=%d.\n", dd, disk_params->RMAX, disk_params->RMIN, PARTICLE_NUMBER);

//    printf("DEBUG [secondaryGrowth]: Initializing radsec, partsecind, masssecvec, sizetemp, mtp, indtemp arrays.\n");
    for(i=0; i < 4*PARTICLE_NUMBER; i++) {	
        temp2[i] = radsec[i][0];			// masodlagosan novesztett porszemcsek tavolsaga
        sizetemp[i] = radsec[i][1];
        temptemp[i] = 0;				// atmeneti vektor, amely elmeni, hogy hol noveszt uj reszecsket a program
        masstemp[i][0] = 0;				// atemeneti tomb, amely a novesztett reszecske tomeget tarolja el
        masstemp[i][1] = 0;				// atmeneti tomb, amely a novesztett reszecske tavolsagat tarolja el
        mtp[i] = masssecvec[i];
        if(temp2[i] != 0) step++;			// megszamolja, hogy hany novesztett reszecske van a fuggvenybe valo belepeskor
        indtemp[i][0] = partsecind[i][3];
    }
//    printf("DEBUG [secondaryGrowth]: Initial step (count of existing secondary particles) = %d.\n", step);

//    printf("DEBUG [secondaryGrowth]: Initializing temp, tempmic, rtempvec, hist, histmic, temp2chX arrays.\n");

    #pragma omp parallel for

    for(i=0; i < PARTICLE_NUMBER; i++) {
        temp[i] = rad[i][0];				// cm-es porszemcsek tavolsaga
        tempmic[i] = radmicr[i][0];			// mikoronos porszemcsek tavolsaga
        rtempvec[i] = disk_params->RMIN + i * dd;			// atmeneti r vektor
        rtempvec[i] = rtempvec[i] + dd / 2.0;		// ez adja a reszecskek eredeti helyet (a generalo programban ide helyeztuk el a reszecskeket	
        hist[i] = 0;					// cm-es porszemcsek hisztogramja
        histmic[i] = 0;					// mikoronos porszemcsek hisztogramja

        temp2ch1[i] = temp2[i];
        temp2ch2[i] = temp2[i+PARTICLE_NUMBER];
        temp2ch3[i] = temp2[i+2*PARTICLE_NUMBER];
        temp2ch4[i] = temp2[i+3*PARTICLE_NUMBER];
    }
//    printf("DEBUG [secondaryGrowth]: Array initializations complete.\n");


/*	Megvizsgalja, hogy az adott tavolsagon hany reszecske talalhato	*/
/*	Valamiert csak ugy ad jo kepet, ha az if-es felteleket beirom.	*/
//    printf("DEBUG [secondaryGrowth]: Calling histogram for various particle types.\n");
    for(i=0;i<PARTICLE_NUMBER;i++) {
        if(temp[i] != 0) {
            // printf("DEBUG [secondaryGrowth]: Calling histogram for temp[%d]=%.2e\n", i, temp[i]);
            histogram(temp[i],hist,dd,disk_params);
        }
        if(tempmic[i] != 0) {
            // printf("DEBUG [secondaryGrowth]: Calling histogram for tempmic[%d]=%.2e\n", i, tempmic[i]);
            histogram(tempmic[i],histmic,dd,disk_params);
        }
        if(temp2ch1[i] != 0) {
            // printf("DEBUG [secondaryGrowth]: Calling histogram for temp2ch1[%d]=%.2e\n", i, temp2ch1[i]);
            histogram(temp2ch1[i],hist,dd,disk_params);
        }
        if(temp2ch2[i] != 0) {
            // printf("DEBUG [secondaryGrowth]: Calling histogram for temp2ch2[%d]=%.2e\n", i, temp2ch2[i]);
            histogram(temp2ch2[i],hist,dd,disk_params);
        }
        if(temp2ch3[i] != 0) {
            // printf("DEBUG [secondaryGrowth]: Calling histogram for temp2ch3[%d]=%.2e\n", i, temp2ch3[i]);
            histogram(temp2ch3[i],hist,dd,disk_params);
        }
        if(temp2ch4[i] != 0) {
            // printf("DEBUG [secondaryGrowth]: Calling histogram for temp2ch4[%d]=%.2e\n", i, temp2ch4[i]);
            histogram(temp2ch4[i],hist,dd,disk_params);
        }
    }
//    printf("DEBUG [secondaryGrowth]: All histogram calls completed.\n");

    j=0;
//    printf("DEBUG [secondaryGrowth]: Checking for new particle growth conditions.\n");
    for(i=0;i<PARTICLE_NUMBER;i++) {
/*	Megvizsgalja, hogy kifogyott-e az adott tavolsagrol a cm-es por. Ha igen, viszont az adott tavolsagon van mikronos, akkor novel egy ujabb reszecsket, ha a mikronos por altal kepviselt tomeg az eredeti helyen levo tomeg szazadanal nagyobb	*/
        if(hist[i] == 0 && histmic[i] != 0) {
            // printf("DEBUG [secondaryGrowth]: Condition met for i=%d (hist[i]=%d, histmic[i]=%d).\n", i, hist[i], histmic[i]);
            if(rtempvec[i] != 0) {
                double sigcurr, siginit;
                sigcurr = partmicind[i][0]/(2.0 * M_PI * dd * (temp[i] - dd/2.));
                siginit = partmicind[i][2]/(2.0 * M_PI * dd * (temp[i] - dd/2.));
                // printf("DEBUG [secondaryGrowth]: sigcurr=%.2e, siginit=%.2e\n", sigcurr, siginit);
                if(sigcurr >= siginit / 100.) {
//                    printf("DEBUG [secondaryGrowth]: New particle grown at rtempvec[%d]=%.2e. j=%d\n", i, rtempvec[i], j);
/*	Uj reszecske novelese --> az adott helyen levo mikronos reszecske tomegenek 75%-at kapja az uj reszecske, ezzel egyutt a mikronos reszecske altal kepviselt reszecske tomege 75%-kal csokken	*/
                    massmicvec[i] = partmicind[i][0] * (1.-.75);	// a mikronos reszecske tomegenek csokkentese
                    masstemp[j][1] = partmicind[i][0] * .75;	// novesztett reszecske altal kepviselt tomeg
                    masstemp[j][0] = rtempvec[i];			// novesztett reszecske altal kepviselt tavolsag
                    temptemp[j] = -1*rtempvec[i];			// atmeneti vektor, amely elmenti, hogy mely helyen noveszt uj reszecsket a program
                    j++;
                }
            }
        }
    }
//    printf("DEBUG [secondaryGrowth]: Particle growth checks completed. Total new particles tracked by j=%d.\n", j);

/*	A novesztett reszecske tavolsagat tartalmazo atmeneti vektor sorbarendezese	*/
//    printf("DEBUG [secondaryGrowth]: Calling sort for temptemp (4*PARTICLE_NUMBER = %d).\n", 4 * PARTICLE_NUMBER);
    sort(temptemp,4*PARTICLE_NUMBER);	
//    printf("DEBUG [secondaryGrowth]: Sort completed.\n");

//    printf("DEBUG [secondaryGrowth]: Processing sorted temptemp to radsec.\n");

    #pragma omp parallel for
    for(i=0; i < 4*PARTICLE_NUMBER-step; i++) {
        if((temptemp[i] != 0)) {
            radsec[i][0] = -1. * temptemp[i];
            // printf("DEBUG [secondaryGrowth]: Moving temptemp[%d]=%.2e to radsec[%d][0]=%.2e\n", i, temptemp[i], i, radsec[i][0]);
            temptemp[i] = 0;
        } else {
            radsec[i][0] = 0;
            radsec[i][1] = 0;
        }
    }
//    printf("DEBUG [secondaryGrowth]: First pass of radsec population completed.\n");

//    printf("DEBUG [secondaryGrowth]: Copying radsec to temptemp and re-sorting.\n");

    #pragma omp parallel for

    for(i=0; i < 4*PARTICLE_NUMBER; i++) {
        temptemp[i] = radsec[i][0];
        radsec[i][0] = 0;
        radsec[i][1] = 0;
    }

    sort(temptemp,4*PARTICLE_NUMBER);
//    printf("DEBUG [secondaryGrowth]: Re-sort of temptemp completed.\n");
//	double index = find_max(indtemp,4*PARTICLE_NUMBER);

//    printf("DEBUG [secondaryGrowth]: Final population of radsec and masssecvec.\n");
    for(i=0; i < 4*PARTICLE_NUMBER; i++) {
        if(temptemp[i] >= disk_params->RMIN) {
            radsec[i][0] = temptemp[i];
            for(j = 0; j < 4*PARTICLE_NUMBER; j++) {
                if(masstemp[j][0] == radsec[i][0]) {
                    masssecvec[i] = masstemp[j][1];
                    partsecind[i][0] = masssecvec[i];
                }
                if(temp2[j] == radsec[i][0]) {
                    masssecvec[i] = mtp[j];
                    radsec[i][1] = sizetemp[j];
                    partsecind[i][3] = indtemp[j][0];
                    partsecind[i][0] = masssecvec[i];
                } else {
                    radsec[i][1] = 1e-6/AU2CM; // This will overwrite if the previous condition wasn't met.
                }
            }
        }
    }
//    printf("DEBUG [secondaryGrowth]: Exiting secondaryGrowth function.\n");
}