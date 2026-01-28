// src/dust_physics.c
#include "dust_physics.h" 
#include "config.h"       
#include "simulation_types.h" 
#include "gas_physics.h"
#include "simulation_core.h" 
#include "boundary_conditions.h"
#include "utils.h"          
#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <omp.h>             // OpenMP támogatáshoz

// Globális változó deklarációk, ha nem lennének meg máshol (pl. config.h)
// Fontos: ezeknek a típusoknak egyezniük kell a config.h-ban deklaráltakkal!
// Ha már szerepelnek a config.h-ban, akkor ezeket innen törölni kell,
// vagy csak az extern kulcsszót meghagyni!



/*	Calculates the Stokes number for each particle	*/
/*	St = rho_particle * radius_particle * PI / (2 * sigma)	*/
double calculateStokesNumber(double pradius, double sigma, disk_t *disk_params) { /*	in the Epstein drag regime	*/
    return disk_params->PDENSITYDIMLESS * pradius * M_PI / (2.0 * sigma);
}

void calculateParticleMass(int n, double (*partmassind)[5], int indii, int indio, int indoi, int indoo, double *massiout, double *massoout, const simulation_options_t *sim_opts) {

    // Debug üzenet frissítve az indexekre

    double massitemp = 0.0;
    double massotemp = 0.0;
    int i;

    // sim_opts->dzone (ez helyettesíti az optdze-t): 1.0 = dinamikus DZE (flag-elt), 0.0 = fix DZE (nem flag-elt)
    if(sim_opts->dzone == 1.0) { // Dinamikus DZE: flag-ek használatával (partmassind[i][3] és [4])
        #pragma omp parallel for private(i) reduction(+:massitemp, massotemp)
        for (i = 0; i < n; i++) {
            // A részecske aktuális grid indexe (partmassind[i][1]-ből)
            int current_r_index = (int)partmassind[i][1]; 

            // --- Belső DZE ---
            // Ellenőrizzük, hogy a részecske grid indexe a belső DZE tartományában van-e
            if ((current_r_index >= indii) && (current_r_index <= indio)) {
                if (partmassind[i][3] == 0.0) { // Belső DZE flag ellenőrzése
                    #pragma omp critical(inner_dze_update)
                    {
                        partmassind[i][3] = 1.0;
                        massitemp = massitemp + partmassind[i][0]; // Tömeg hozzáadása a [0] indexről
                    }
                }
            }

            // --- Külső DZE ---
            // Ellenőrizzük, hogy a részecske grid indexe a külső DZE tartományában van-e
            if ((current_r_index >= indoi) && (current_r_index <= indoo)) {
                if (partmassind[i][4] == 0.0) { // Külső DZE flag ellenőrzése (Flag[4])
                    #pragma omp critical(outer_dze_update)
                    {
                        partmassind[i][4] = 1.0;
                        massotemp = massotemp + partmassind[i][0]; // Tömeg hozzáadása a [0] indexről
                    }
                }
            }
        }
    } else { // Fix DZE (sim_opts->dzone == 0.0): Nincsenek flag-ek a tömeg felhalmozáshoz
        #pragma omp parallel for private(i) reduction(+:massitemp, massotemp)
        for (i = 0; i < n; i++) {
            int current_r_index = (int)partmassind[i][1]; // A részecske grid indexe

            // --- Belső DZE (Fix) ---
            // Az eredeti kódod else ága nem tartalmazta a belső DZE gyűjtését. 
            // Ha szeretnéd, akkor ide kell tenni a logikát:
            if ((current_r_index >= indii) && (current_r_index <= indio)) {
                #pragma omp critical(inner_dze_update_fixed)
                {
                    massitemp = massitemp + partmassind[i][0]; 
                }
            }

            // --- Külső DZE (Fix) ---
            if ((current_r_index >= indoi) && (current_r_index <= indoo)) {
                #pragma omp critical(outer_dze_update_fixed)
                {
                    massotemp = massotemp + partmassind[i][0];
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
double calculateRadialDriftBarrier(double sigmad, double r, double p, double dp, double rho_p, const disk_t *disk_params) {

    double Sigmad_cgs = sigmad / SURFACE_DENSITY_CONVERSION_FACTOR;

    double vkep = calculateKeplerianVelocity(r,disk_params);
    double vkep2 = vkep * vkep;
    double c_s = calculateLocalSoundSpeed(r,disk_params);
    double c_s2 = c_s * c_s;
    double dlnPdlnr = r / p * dp;
    double s_drift =  disk_params->fDrift * 2.0 / M_PI * Sigmad_cgs / rho_p * vkep2 / c_s2 * fabs(1.0 / dlnPdlnr);
    return s_drift;
}

// 2. a kis skalaju turbulencia altal okozott fragmentacio szerinti maximalis meret	--> kimenet cm-ben!
double calculateTurbulentFragmentationBarrier(double sigma, double r, double rho_p, const disk_t *disk_params) {

    double s_frag, u_frag, u_frag2, Sigma_cgs, c_s, c_s2;

    u_frag = disk_params->uFrag * CM_PER_SEC_TO_AU_PER_YEAR_2PI; /*	cm/sec --> AU / (yr/2pi)	*/
    u_frag2 = u_frag * u_frag;
    Sigma_cgs = sigma / SURFACE_DENSITY_CONVERSION_FACTOR;
    c_s = calculateLocalSoundSpeed(r,disk_params); // / CM_PER_SEC_TO_AU_PER_YEAR_2PI; // Komment ki, ha a calculateLocalSoundSpeed már megfelelő mértékegységben van
    c_s2 = c_s * c_s;

    s_frag = disk_params->fFrag * 2.0 / (3.0 * M_PI) * Sigma_cgs / (rho_p * calculateTurbulentAlpha(r,disk_params)) * u_frag2 / c_s2;

    return s_frag;
}

// 3. radialis drift altal okozott fragmentacio szerinti maximalis meret		--> kimenet cm-ben!
double calculateDriftInducedFragmentationBarrier(double sigma, double r, double p, double dp, double rho_p, const disk_t *disk_params) {

    double u_frag, vkep, dlnPdlnr, c_s, c_s2, s_df, Sigma_cgs;

    u_frag = disk_params->uFrag * CM_PER_SEC_TO_AU_PER_YEAR_2PI; /*	cm/sec --> AU / (yr/2pi)	*/
    Sigma_cgs = sigma / SURFACE_DENSITY_CONVERSION_FACTOR;
    c_s = calculateLocalSoundSpeed(r,disk_params);
    c_s2 = c_s * c_s;
    dlnPdlnr = r / p * dp;
    vkep = calculateKeplerianVelocity(r,disk_params);

    s_df = u_frag * vkep / fabs(dlnPdlnr * c_s2 * 0.5) * 2.0 * Sigma_cgs / (M_PI * rho_p);

    return s_df;
}

/*	a reszecskek novekedesenek idoskalaja	*/
double calculateGrowthTimescale(double r, double eps,const disk_t *disk_params) {
    double omega = calculateKeplerianFrequency(r,disk_params);
    double calculateGrowthTimescale = eps / omega;
    return calculateGrowthTimescale;
}

/*	kiszamolja az adott helyen a reszecske meretet --> BIRNSTIEL CIKK	*/
double calculateDustParticleSize(double prad, double pdens, double sigma, double sigmad, double y, double p, double dpress_val, double dt, const disk_t *disk_params) {

    double sturb = calculateTurbulentFragmentationBarrier(sigma, y, pdens, disk_params);           // cm-ben
    double sdf = calculateDriftInducedFragmentationBarrier(sigma, y, p, dpress_val, pdens,disk_params); // cm-ben
    double srdf = calculateRadialDriftBarrier(sigmad, y, p, dpress_val, pdens, disk_params); // cm-ben
    double smin = findMinimumOfAnArray(sturb, sdf, srdf);         // cm-ben -- megadja, hogy a fenti ket reszecske korlatbol melyik ad kisebb meretet (az a reszecskenovekedes felso korlatja
    //	double eps = sigma / 100.;
    double eps = sigmad / sigma; // A korábbi kódban fordítva volt, feltételezem, hogy eps = (por sűrűség) / (gáz sűrűség)
    double tau_gr = calculateGrowthTimescale(y, eps, disk_params);
    double rt = 0.0;

    smin = smin / AU_IN_CM; // AU-ban

    /*	kiszamolja, hogy a fenti smin, vagy a novekedesi idoskalabol szarmazo meret korlatozza a reszecske meretet	*/
    if (prad < smin) {
        rt = findMinimumOfAnArray(prad * exp(dt / tau_gr), smin, HUGE_VAL);
    } else { // prad >= smin
        rt = smin;
    }

    return rt;
}


void calculateDustSurfaceDensity(double max_param, double min_param, double rad[][2], double radmicr[][2], 
                double *sigma_d, double *sigma_dm,  double *massvec, double *massmicrvec,  
                double *rd, double *rmic, const simulation_options_t *sim_opts, const disk_t *disk_params) {

    // Suppress unused parameter warnings
    (void)max_param;
    (void)min_param;

    double dd = (disk_params->RMAX - disk_params->RMIN) / (particle_number - 1);
    int i;

    // A temp tömbök deklarálását érdemes a scope tetejére tenni
    double sigdtemp[particle_number][3];
    double sigdmicrtemp[particle_number][3];

    // Inicializálás, ha szükséges (bár a calculateDustSurfaceDensity valószínűleg felülírja)
    for(i=0; i<particle_number; i++){
        sigdtemp[i][0] = 0.0; sigdtemp[i][1] = 0.0; sigdtemp[i][2] = 0.0;
        sigdmicrtemp[i][0] = 0.0; sigdmicrtemp[i][1] = 0.0; sigdmicrtemp[i][2] = 0.0;
        rd[i] = 0.0;
        rmic[i] = 0.0;
        sigma_d[i] = 0.0;
        sigma_dm[i] = 0.0;
    }


    // calculateDustSurfaceDensity és mergeParticlesByRadius függvények hívásai:
    // Ezek valószínűleg szekvenciálisak, hacsak a függvények belsejében nincs OpenMP.
    // Ha ezek a függvények valamilyen globális állapotot módosítanak, akkor kritikusak.
    // Feltételezve, hogy a 'sigdtemp' és 'sigdmicrtemp' kizárólagosan a hívásaikban vannak feldolgozva,
    // és nem ütköznek más szálakkal globális adatokon keresztül.
    calculateInitialDustSurfaceDensity(rad, massvec, sigdtemp, particle_number,disk_params);
    if (sim_opts->twopop == 1.0) { // Használjunk double összehasonlítást
        calculateInitialDustSurfaceDensity(radmicr, massmicrvec, sigdmicrtemp, particle_number,disk_params);
    }

    mergeParticlesByRadius(sigdtemp, dd, particle_number,disk_params);
    if (sim_opts->twopop == 1.0) { // Használjunk double összehasonlítást
        mergeParticlesByRadius(sigdmicrtemp, dd, particle_number,disk_params);
    }

    // Utolsó másoló ciklus: Ez is jól párhuzamosítható.
    #pragma omp parallel for private(i)
    for (i = 0; i < particle_number; i++) {
        rd[i] = sigdtemp[i][1];
        sigma_d[i] = sigdtemp[i][0];

        if (sim_opts->twopop == 1.0) { // double összehasonlítás
            rmic[i] = sigdmicrtemp[i][1];
            sigma_dm[i] = sigdmicrtemp[i][0];
        }
    }
}


/*	Fuggveny a porszemcsek uj tavolsaganak elraktarozasara		*/
void calculateDustDistance(const char *nev, int opt, double radius[][2], const double *sigmad, const double *rdvec,
                double deltat, double t, int n, const simulation_options_t *sim_opts, const disk_t *disk_params){

    int i;
    double y, y_out, prad_new, particle_radius;
    char scout[1024];

    // Fájlkezelés t==0 esetén: ez valószínűleg egyszer történik meg a szimuláció elején,
    // még mielőtt az igazi párhuzamosítás elkezdődne a fő ciklusban.
    // Biztosítjuk, hogy csak egy szál nyissa meg/zárja be a fájlt.
    #pragma omp master
    {
        if (t == 0) {
            sprintf(scout, "%s/%s%s", nev, kDriftTimescaleFileName, kFileNamesSuffix);
            drift_timescale_file = fopen(scout, "w");
        }
    }
    // Szinkronizálás a master szál befejezéséig.
    #pragma omp barrier

    // **ITT A PÁRHUZAMOSÍTÁS!**
    // Az `i` ciklus független iterációkkal rendelkezik, minden szál a saját `radius[i]` elemen dolgozik.
    #pragma omp parallel for private(y, y_out, prad_new, particle_radius)
    for (i = 0; i < n; i++) {
        // Csak a RMIN és RMAX közötti részecskékkel foglalkozunk.
        // A 0.0-ra állítás kívül esik a párhuzamos részen, ha az if feltétel nem teljesül.
        if (radius[i][0] > disk_params->RMIN && radius[i][0] < disk_params->RMAX) {
            y = radius[i][0];
            particle_radius = radius[i][1];

			integrateParticleRungeKutta4(t, particle_radius, sigmad, rdvec, deltat, y, &y_out, &prad_new, disk_params, sim_opts);
            if (t == 0) {
                if (sim_opts->twopop == 0) {
                    double current_drdt_val = (fabs(y_out - y) / (deltat));
                    // Azért kell a critical szekció, mert az drift_timescale_file fájlba írunk.
                    // Ez a critical szekció biztosítja, hogy egyszerre csak egy szál írjon a fájlba.

                    #pragma omp critical(drift_timescale_file_write)
                    {
                        // Ellenőrizzük, hogy a fájlmutató nem NULL
                        if (drift_timescale_file != NULL) {
                            fprintf(drift_timescale_file, "%lg %lg\n", radius[i][0], (radius[i][0] / current_drdt_val) / (2.0 * M_PI));
                        } else {
                            fprintf(stderr, "ERROR: drift_timescale_file is NULL during write in calculateDustDistance (t=0 block).\n");
                        }
                    }
                }
            }

            if (sim_opts->twopop != 1) { // Ha növekedés engedélyezett vagy valami más mód
                radius[i][1] = prad_new;
                radius[i][0] = y_out;
            } else { // sim_opts->twopop == 1, csak drift
                radius[i][0] = y_out;
            }
        } else {
            // Ha a részecske RMIN vagy RMAX kívülre kerül, 0.0-ra állítjuk a pozícióját.
            // Ez a hozzárendelés iterációnként független, így párhuzamosítható.
            radius[i][0] = 0.0;
        }
    }

    // A fájl bezárása ismételten egy szál által kell, hogy történjen.
    #pragma omp master
    {
        if (t == 0) { // Csak akkor zárjuk be, ha meg is nyitottuk a t=0 blokkban
            if (drift_timescale_file != NULL) {
                fclose(drift_timescale_file);
                drift_timescale_file = NULL; // Fontos, hogy NULL-ra állítsuk, miután bezártuk.
            }
        }
    }
    // Szinkronizálás a master szál befejezéséig.
    #pragma omp barrier
}