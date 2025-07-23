// src/dust_physics.c
#include "dust_physics.h" // A saját headerét mindig includolni kell
#include "config.h"       // Szükséges lehet a globális konstansokhoz
#include "simulation_types.h" // Például output_files_t, disk_t struktúrákhoz

#include "simulation_core.h" // int_step, Perem, find_num_zero, find_zero, find_r_annulus függvényekhez
#include "utils.h"           // find_min függvényhez
#include <stdio.h>
#include <stdlib.h>

#include <math.h>     // For math functions like powl, sqrtl, fabsl, floorl, tanhl
#include <omp.h>      // OpenMP támogatáshoz

// FONTOS: Győződj meg róla, hogy az összes hívott függvény deklarációja
// és definíciója (pl. Perem, int_step, find_min, Coeff_1, Coeff_2, press, dpress, u_gas, stb.)
// szintén long double paramétereket és visszatérési értékeket használ,
// amennyiben lebegőpontos számokat érintenek!
// Ugyanígy a disk_t struktúra minden releváns tagja is long double kell, hogy legyen.

/*	alpha turbulens paraméter kiszámolása --> alfa csökkentése alpha_r-rel	*/
long double calculate_turbulent_alpha(long double r, const disk_t *disk_params) {
    long double alpha_r;
    alpha_r = 1.0L - 0.5L * (1.0L - disk_params->a_mod) * (tanhl((r - disk_params->r_dze_i) / disk_params->Dr_dze_i) + tanhl((disk_params->r_dze_o - r) / disk_params->Dr_dze_o));
    return alpha_r * disk_params->alpha_visc;
}

/*	kiszamolja az adott reszecskehez tartozo Stokes szamot	*/
/*	St = rho_particle * radius_particle * PI / (2 * sigma)	*/
long double Stokes_Number(long double pradius, long double sigma, disk_t *disk_params) { /*	in the Epstein drag regime	*/
    return disk_params->PDENSITYDIMLESS * pradius * M_PIl / (2.0L * sigma);
}

/*	Lokalis viszkozitas erteke	*/
long double visc(long double r, const disk_t *disk_params) {
    long double nu;
    long double cs, H;

    H = scale_height(r, disk_params);
    cs = c_sound(r, disk_params);

    nu = calculate_turbulent_alpha(r, disk_params) * cs * H;
    return nu;
}

/*	local scale height	*/
long double scale_height(long double r, const disk_t *disk_params) {

    if (disk_params == NULL) {
        fprintf(stderr, "ERROR [scale_height]: disk_params is NULL!\n");
        return 0.0L; // Vagy valamilyen hibakód/NaN
    }

    long double calculated_result = powl(r, 1.0L + disk_params->FLIND) * disk_params->HASP;
    return calculated_result;
}

/*	lokális kepleri sebesség	*/
long double v_kep(long double r, const disk_t *disk_params) {
    return sqrtl(G_GRAV_CONST * disk_params->STAR_MASS / r);
}

/*	lokalis kepleri korfrekvencia	*/
long double kep_freq(long double r, const disk_t *disk_params) {
    return sqrtl(G_GRAV_CONST * disk_params->STAR_MASS / r / r / r);
}

/*	local sound speed		*/
long double c_sound(long double r, const disk_t *disk_params) {
    return kep_freq(r, disk_params) * scale_height(r, disk_params);
}

/*	Suruseg a midplane-ben	*/
long double rho_mp(long double sigma, long double r, const disk_t *disk_params) {
    return 1.0L / sqrtl(2.0L * M_PIl) * sigma / scale_height(r, disk_params);
}

/* local pressure of the gas p = rho_gas * cs * cs kepletbol!!	*/
long double press(long double sigma, long double r, const disk_t *disk_params) {
    return rho_mp(sigma, r, disk_params) * c_sound(r, disk_params) * c_sound(r, disk_params);
}

/*	a nyomas derivaltja	*/
void dpress(disk_t *disk_params) {
    int i;
    long double ptemp; // Changed to long double
    long double pvec[disk_params->NGRID + 2]; // Changed to long double array

    for (i = 1; i <= disk_params->NGRID; i++) {
        ptemp = (disk_params->pressvec[i + 1] - disk_params->pressvec[i - 1]) / (2.0L * disk_params->DD);
        pvec[i] = ptemp;
    }
    for (i = 1; i <= disk_params->NGRID; i++) {
        disk_params->dpressvec[i] = pvec[i];
    }
}

/*	u_gas kiszamolasahoz eltarolt koefficiens	*/
long double Coeff_3(long double sigma, long double r) {
    return -1.0L * (3.0L / (sigma * sqrtl(r)));
}

/*	u_gas = -3/(Sigma*R^0.5)*(d/dR)(nu*Sigma*R^0.5) kiszamolasa	*/
void u_gas(disk_t *disk_params) {

    long double tempug; // Changed to long double
    // Lokális tömbök, méret NGRID-hez igazítva disk_params-ból
    long double ugvec[disk_params->NGRID + 2]; // Changed to long double
    long double ugvectemp[disk_params->NGRID + 1]; // Changed to long double

    int i;

    // Első ciklus: feltölti a lokális ugvec tömböt
    #pragma omp parallel for private(i)
    for (i = 0; i <= disk_params->NGRID + 1; i++) {
        ugvec[i] = disk_params->sigmavec[i] * visc(disk_params->rvec[i], disk_params) * sqrtl(disk_params->rvec[i]);
    }

    // Második ciklus: feltölti a lokális ugvectemp tömböt
    #pragma omp parallel for private(i, tempug)
    for (i = 1; i <= disk_params->NGRID; i++) {
        tempug = (ugvec[i + 1] - ugvec[i - 1]) / (2.0L * disk_params->DD);
        ugvectemp[i] = Coeff_3(disk_params->sigmavec[i], disk_params->rvec[i]) * tempug;
    }

    // Harmadik ciklus: Az eredményt bemásolja a disk_params->ugvec-be
    for (i = 1; i <= disk_params->NGRID; i++) {
        disk_params->ugvec[i] = ugvectemp[i];
    }
}

// GetMass függvény, minden double-t long double-re alakít
void GetMass(int n, long double (*partmassind)[5], int indii, int indio, int indoi, int indoo, int indoi_buff, int indoo_buff, long double *massiout, long double *massoout, long double *massoout_buff, const simulation_options_t *sim_opts) {

    long double massitemp = 0.0L;
    long double massotemp = 0.0L;
    long double massobufftemp = 0.0L;
    int i;

    #pragma omp parallel for private(i)
    for (i = 0; i < n; i++) {
        partmassind[i][3] = 0.0L; // Reseteli a 'már besorolva' flag-et (long double)
    }

    if(sim_opts->dzone == 1.0L) { // Dinamikus DZE
        #pragma omp parallel for private(i) reduction(+:massitemp, massotemp, massobufftemp)
        for (i = 0; i < n; i++) {
            int current_r_index = (int)partmassind[i][1];

            if ((current_r_index >= indii) && (current_r_index <= indio)) {
                partmassind[i][3] = 1.0L;
                massitemp += partmassind[i][0];
            }
            else if ((current_r_index >= indoi) && (current_r_index <= indoo)) {
                partmassind[i][3] = 1.0L;
                massotemp += partmassind[i][0];
            }
            else if ((current_r_index >= indoi_buff) && (current_r_index <= indoo_buff)) {
                partmassind[i][3] = 1.0L;
                massobufftemp += partmassind[i][0];
            }
        }
    } else { // Fix DZE
        #pragma omp parallel for private(i) reduction(+:massitemp, massotemp, massobufftemp)
        for (i = 0; i < n; i++) {
            int current_r_index = (int)partmassind[i][1];

            if ((current_r_index >= indii) && (current_r_index <= indio)) {
                massitemp += partmassind[i][0];
            }
            if ((current_r_index >= indoi) && (current_r_index <= indoo)) {
                massotemp += partmassind[i][0];
            }
            if ((current_r_index >= indoi_buff) && (current_r_index <= indoo_buff)) {
                massobufftemp += partmassind[i][0];
            }
        }
    }

    *massiout = massitemp;
    *massoout = massotemp;
    *massoout_buff = massobufftemp;
}

/*			BIRNSTIEL EL AL 2012			*/

//reprezentativ reszecske kezdeti meretenek meghatarozasa
// 1. radialis drift altal meghatarozott maximalis meret			--> kimenet cm-ben!
long double a_drift(long double sigmad, long double r, long double p, long double dp, long double rho_p, const disk_t *disk_params) {

    long double Sigmad_cgs = sigmad / SDCONV;

    long double vkep = v_kep(r, disk_params);
    long double vkep2 = vkep * vkep;
    long double c_s = c_sound(r, disk_params);
    long double c_s2 = c_s * c_s;
    long double dlnPdlnr = r / p * dp;
    long double s_drift = disk_params->fDrift * 2.0L / M_PIl * Sigmad_cgs / rho_p * vkep2 / c_s2 * fabsl(1.0L / dlnPdlnr);
    return s_drift;
}

// 2. a kis skalaju turbulencia altal okozott fragmentacio szerinti maximalis meret	--> kimenet cm-ben!
long double a_turb(long double sigma, long double r, long double rho_p, const disk_t *disk_params) {

    long double s_frag, u_frag, u_frag2, Sigma_cgs, c_s, c_s2;

    u_frag = disk_params->uFrag * CMPSECTOAUPYRP2PI; /*	cm/sec --> AU / (yr/2pi)	*/
    u_frag2 = u_frag * u_frag;
    Sigma_cgs = sigma / SDCONV;
    c_s = c_sound(r, disk_params);
    c_s2 = c_s * c_s;

    s_frag = disk_params->fFrag * 2.0L / (3.0L * M_PIl) * Sigma_cgs / (rho_p * calculate_turbulent_alpha(r, disk_params)) * u_frag2 / c_s2;

    return s_frag;
}

// 3. radialis drift altal okozott fragmentacio szerinti maximalis meret		--> kimenet cm-ben!
long double a_df(long double sigma, long double r, long double p, long double dp, long double rho_p, const disk_t *disk_params) {

    long double u_frag, vkep, dlnPdlnr, c_s, c_s2, s_df, Sigma_cgs;

    u_frag = disk_params->uFrag * CMPSECTOAUPYRP2PI; /*	cm/sec --> AU / (yr/2pi)	*/
    Sigma_cgs = sigma / SDCONV;
    c_s = c_sound(r, disk_params);
    c_s2 = c_s * c_s;
    dlnPdlnr = r / p * dp;
    vkep = v_kep(r, disk_params);

    s_df = u_frag * vkep / fabsl(dlnPdlnr * c_s2 * 0.5L) * 2.0L * Sigma_cgs / (M_PIl * rho_p);

    return s_df;
}

/*	a reszecskek novekedesenek idoskalaja	*/
long double tauGr(long double r, long double eps, const disk_t *disk_params) {
    long double omega = kep_freq(r, disk_params);
    long double taugr = eps / omega;
    return taugr;
}

/*	kiszamolja az adott helyen a reszecske meretet --> BIRNSTIEL CIKK	*/
long double getSize(long double prad, long double pdens, long double sigma, long double sigmad, long double y, long double p, long double dpress_val, long double dt, const disk_t *disk_params) {

    long double sturb = a_turb(sigma, y, pdens, disk_params);    // cm-ben
    long double sdf = a_df(sigma, y, p, dpress_val, pdens, disk_params); // cm-ben
    long double srdf = a_drift(sigmad, y, p, dpress_val, pdens, disk_params); // cm-ben
    long double smin = find_min(sturb, sdf, srdf); // cm-ben -- megadja, hogy a fenti ket reszecske korlatbol melyik ad kisebb meretet (az a reszecskenovekedes felso korlatja
    long double eps = sigmad / sigma; // A korábbi kódban fordítva volt, feltételezem, hogy eps = (por sűrűség) / (gáz sűrűség)
    long double tau_gr = tauGr(y, eps, disk_params);
    long double rt = 0.0L;

    smin = smin / AU2CM; // AU-ban

    /*	kiszamolja, hogy a fenti smin, vagy a novekedesi idoskalabol szarmazo meret korlatozza a reszecske meretet	*/
    if (prad < smin) {
        rt = find_min(prad * expl(dt / tau_gr), smin, HUGE_VALL); // expl for long double exponential, HUGE_VALL for long double infinity
    } else { // prad >= smin
        rt = smin;
    }

    return rt;
}

void Get_Sigmad(long double max_param, long double min_param,
                long double rad[][2], long double radmicr[][2], // Módosítva long double-re
                // THESE PARAMETERS MUST MATCH THE HEADER EXACTLY IN ORDER AND TYPE
                long double *massvec,      // Now 5th parameter, matches header
                long double *massmicrvec,  // Now 6th parameter, matches header
                long double *sigma_d,      // Reordered and Módosítva long double-re
                long double *sigma_dm,     // Reordered and Módosítva long double-re
                long double *rd, long double *rmic, // Módosítva long double-re
                const simulation_options_t *sim_opts, const disk_t *disk_params) {

    // Suppress unused parameter warnings
    (void)max_param;
    (void)min_param;

    long double dd = (disk_params->RMAX - disk_params->RMIN) / (long double)(PARTICLE_NUMBER - 1); // PARTICLE_NUMBER int, castoljuk
    int i;

    long double sigdtemp[PARTICLE_NUMBER][3];
    long double sigdmicrtemp[PARTICLE_NUMBER][3];

    #pragma omp parallel for private(i)
    for(i = 0; i < PARTICLE_NUMBER; i++){
        sigdtemp[i][0] = 0.0L; sigdtemp[i][1] = 0.0L; sigdtemp[i][2] = 0.0L;
        sigdmicrtemp[i][0] = 0.0L; sigdmicrtemp[i][1] = 0.0L; sigdmicrtemp[i][2] = 0.0L;
        rd[i] = 0.0L; // Módosítva long double literálra
        rmic[i] = 0.0L; // Módosítva long double literálra
        sigma_d[i] = 0.0L; // Módosítva long double literálra
        sigma_dm[i] = 0.0L; // Módosítva long double literálra
    }

    loadSigDust(rad, massvec, sigdtemp, PARTICLE_NUMBER, disk_params);
    if (sim_opts->twopop == 1.0L) { // Hasonlítsunk long double-hez
        loadSigDust(radmicr, massmicrvec, sigdmicrtemp, PARTICLE_NUMBER, disk_params);
    }

    contract(sigdtemp, dd, PARTICLE_NUMBER, disk_params);
    if (sim_opts->twopop == 1.0L) { // Hasonlítsunk long double-hez
        contract(sigdmicrtemp, dd, PARTICLE_NUMBER, disk_params);
    }

    #pragma omp parallel for private(i)
    for (i = 0; i < PARTICLE_NUMBER; i++) {
        rd[i] = sigdtemp[i][1]; // Már long double
        sigma_d[i] = sigdtemp[i][0]; // Már long double

        if (sim_opts->twopop == 1.0L) { // Hasonlítsunk long double-hez
            rmic[i] = sigdmicrtemp[i][1]; // Már long double
            sigma_dm[i] = sigdmicrtemp[i][0]; // Már long double
        }
    }
}


/*	Fuggveny a sigma, p, dp kiszamolasara	*/
void Get_Sigma_P_dP(const simulation_options_t *sim_opts, disk_t *disk_params) { // Added sim_opts

    long double u, u_bi, u_fi; // Módosítva long double-re
    long double sigma_temp[disk_params->NGRID + 2]; // Módosítva long double-re
    long double uvec[disk_params->NGRID + 2];    // Módosítva long double-re

    int i;

    // Boundary conditions - access via disk_params
    sigma_temp[0] = disk_params->sigmavec[0];
    sigma_temp[disk_params->NGRID + 1] = disk_params->sigmavec[disk_params->NGRID + 1];

    // uvec temporary array initialization
    uvec[0] = disk_params->sigmavec[0] * visc(disk_params->rvec[0], disk_params);
    uvec[disk_params->NGRID + 1] = disk_params->sigmavec[disk_params->NGRID + 1] * visc(disk_params->rvec[disk_params->NGRID + 1], disk_params);

    #pragma omp parallel for
    for(i = 1; i <= disk_params->NGRID; i++) {
        uvec[i] = disk_params->sigmavec[i] * visc(disk_params->rvec[i], disk_params);
    }

    // This loop is critical due to data dependencies. Keep it sequential for correctness
    for (i = 1; i <= disk_params->NGRID; i++) {
        u = uvec[i];
        u_bi = uvec[i - 1];
        u_fi = uvec[i + 1];

        long double temp = Coeff_1(disk_params->rvec[i], disk_params) * (u_fi - 2.0L * u + u_bi) / (disk_params->DD * disk_params->DD) +
                           Coeff_2(disk_params->rvec[i], disk_params) * (u_fi - u_bi) / (2.0L * disk_params->DD);
        
        sigma_temp[i] = uvec[i] + sim_opts->DT * temp;
    }

    // This loop is parallelizable
    #pragma omp parallel for
    for (i = 1; i <= disk_params->NGRID; i++) {
        disk_params->sigmavec[i] = sigma_temp[i] / visc(disk_params->rvec[i], disk_params);
        disk_params->pressvec[i] = press(disk_params->sigmavec[i], disk_params->rvec[i], disk_params);
    }

    dpress(disk_params);
    Perem(disk_params->sigmavec, disk_params);
    Perem(disk_params->pressvec, disk_params);
    Perem(disk_params->dpressvec, disk_params);
}

/*	Fuggveny a porszemcsek uj tavolsaganak elraktarozasara		*/
void Get_Radius(const char *nev, int opt, long double radius[][2], const long double *sigmad, const long double *rdvec, // Módosítva long double-re
                long double deltat, long double t, int n, const simulation_options_t *sim_opts, const disk_t *disk_params){ // Módosítva long double-re

    int i;
    long double y, y_out, prad_new, particle_radius; // Módosítva long double-re
    char scout[1024];

    // Fájlkezelés t==0 esetén: ez valószínűleg egyszer történik meg a szimuláció elején,
    // még mielőtt az igazi párhuzamosítás elkezdődne a fő ciklusban.
    // Biztosítjuk, hogy csak egy szál nyissa meg/zárja be a fájlt.
    #pragma omp master
    {
        if (t == 0.0L) { // long double összehasonlítás
            sprintf(scout, "%s/timescale.dat", nev);
            fout2 = fopen(scout, "w");
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

            int_step(t, particle_radius, sigmad, rdvec, deltat, y, &y_out, &prad_new, disk_params, sim_opts);
            if (t == 0.0L) { // long double összehasonlítás
                if (sim_opts->twopop == 0.0L) { // long double összehasonlítás
                    long double current_drdt_val = (fabsl(y_out - y) / (deltat)); // fabsl for long double
                    // Azért kell a critical szekció, mert az fout2 fájlba írunk.
                    // Ez a critical szekció biztosítja, hogy egyszerre csak egy szál írjon a fájlba.

                    #pragma omp critical(fout2_write)
                    {
                        // Ellenőrizzük, hogy a fájlmutató nem NULL
                        if (fout2 != NULL) {
                            // %Lg a long double kiírására
                            fprintf(fout2, "%Lg %Lg\n", radius[i][0], (radius[i][0] / current_drdt_val) / (2.0L * M_PIl));
                        } else {
                            fprintf(stderr, "ERROR: fout2 is NULL during write in Get_Radius (t=0 block).\n");
                        }
                    }
                }
            }

            if (sim_opts->twopop != 1.0L) { // Ha növekedés engedélyezett vagy valami más mód (long double összehasonlítás)
                radius[i][1] = prad_new;
                radius[i][0] = y_out;
            } else { // sim_opts->twopop == 1.0L, csak drift (long double összehasonlítás)
                radius[i][0] = y_out;
            }
        } else {
            // Ha a részecske RMIN vagy RMAX kívülre kerül, 0.0-ra állítjuk a pozícióját.
            // Ez a hozzárendelés iterációnként független, így párhuzamosítható.
            radius[i][0] = 0.0L; // Módosítva long double literálra
        }
    }

    // A fájl bezárása ismételten egy szál által kell, hogy történjen.
    #pragma omp master
    {
        if (t == 0.0L) { // Csak akkor zárjuk be, ha meg is nyitottuk a t=0 blokkban
            if (fout2 != NULL) {
                fclose(fout2);
                fout2 = NULL; // Fontos, hogy NULL-ra állítsuk, miután bezártuk.
            }
        }
    }
    // Szinkronizálás a master szál befejezéséig.
    #pragma omp barrier
}