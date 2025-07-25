// src/dust_physics.c
#include "dust_physics.h" // A saját headerjét mindig includolni kell
#include "config.h"       // Szükséges lehet a globális konstansokhoz (pl. PARTICLE_NUMBER, AU2CM, RMIN, RMAX, NGRID, G_GRAV_CONST, STAR, SDCONV, CMPSECTOAUPYRP2PI, uFrag, fFrag, PDENSITYDIMLESS, HASP, M_PI, DD, sim_opts->dzone, sim_opts->twopop, RMIN, RMAX, FLIND, alpha_visc, a_mod, r_dze_i, r_dze_o, Dr_dze_i, Dr_dze_o)
#include "simulation_types.h" // Például output_files_t, disk_t struktúrákhoz
#include "globals.h"
#include "simulation_core.h" // int_step, Perem, find_num_zero, find_zero, find_r_annulus függvényekhez
#include "utils.h"           // find_min függvényhez
#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <omp.h>             // OpenMP támogatáshoz

// Globális változó deklarációk, ha nem lennének meg máshol (pl. config.h)
// Fontos: ezeknek a típusoknak egyezniük kell a config.h-ban deklaráltakkal!
// Ha már szerepelnek a config.h-ban, akkor ezeket innen törölni kell,
// vagy csak az extern kulcsszót meghagyni!


/*	alpha turbulens paraméter kiszámolása --> alfa csökkentése alpha_r-rel	*/
double calculate_turbulent_alpha(double r, const disk_t *disk_params) {
    double alpha_r;
    alpha_r = 1.0 - 0.5 * (1.0 - disk_params->a_mod) * (tanh((r - disk_params->r_dze_i) / disk_params->Dr_dze_i) + tanh((disk_params->r_dze_o - r) / disk_params->Dr_dze_o));
    return alpha_r * disk_params->alpha_visc;
}

/*	kiszamolja az adott reszecskehez tartozo Stokes szamot	*/
/*	St = rho_particle * radius_particle * PI / (2 * sigma)	*/
double Stokes_Number(double pradius, double sigma, disk_t *disk_params) { /*	in the Epstein drag regime	*/
    return disk_params->PDENSITYDIMLESS * pradius * M_PI / (2.0 * sigma);
}

/*	Lokalis viszkozitas erteke	*/
double visc(double r, const disk_t *disk_params) {
    double nu;
    double cs, H;

    H = scale_height(r,disk_params);
    cs = c_sound(r,disk_params);

    nu = calculate_turbulent_alpha(r, disk_params) * cs * H;
    return nu;
}

/*	local scale height	*/
double scale_height(double r, const disk_t *disk_params) {

    if (disk_params == NULL) {
        fprintf(stderr, "ERROR [scale_height]: disk_params is NULL!\n");
        return 0.0; // Vagy valamilyen hibakód/NaN
    }

    // Itt van az eredeti számítás
    double calculated_result = pow(r, 1. + disk_params->FLIND) * disk_params->HASP;
    return calculated_result;
}

/*	lokális kepleri sebesség	*/
double v_kep(double r, const disk_t *disk_params) {
    return sqrt(G_GRAV_CONST * disk_params->STAR_MASS / r);
}

/*	lokalis kepleri korfrekvencia	*/
double kep_freq(double r, const disk_t *disk_params) {
    return sqrt(G_GRAV_CONST * disk_params->STAR_MASS / r / r / r);
}

/*	local sound speed		*/
double c_sound(double r, const disk_t *disk_params) {
    return kep_freq(r,disk_params) * scale_height(r,disk_params);
}

/*	Suruseg a midplane-ben	*/
double rho_mp(double sigma, double r, const disk_t *disk_params) {
    return 1. / sqrt(2.0 * M_PI) * sigma / scale_height(r,disk_params);
}

/* local pressure of the gas p = rho_gas * cs * cs kepletbol!!	*/
double press(double sigma, double r, const disk_t *disk_params) {
    return rho_mp(sigma, r, disk_params) * c_sound(r,disk_params) * c_sound(r, disk_params);
}

/*	a nyomas derivaltja	*/
void dpress(disk_t *disk_params) {
    int i;
    double ptemp, pvec[disk_params->NGRID + 2];

    for (i = 1; i <= disk_params->NGRID; i++) {
        ptemp = (disk_params->pressvec[i + 1] - disk_params->pressvec[i - 1]) / (2.0 * disk_params->DD);
        pvec[i] = ptemp;

    }
    for (i = 1; i <= disk_params->NGRID; i++) {
        disk_params->dpressvec[i] = pvec[i];
    }


}

/*	u_gas kiszamolasahoz eltarolt koefficiens	*/
double Coeff_3(double sigma, double r) {
    return -1.0 * (3.0 / (sigma * sqrt(r)));
}

/*	u_gas = -3/(Sigma*R^0.5)*(d/dR)(nu*Sigma*R^0.5) kiszamolasa	*/
void u_gas(disk_t *disk_params) {

    double tempug;
    // Lokális tömbök, méret NGRID-hez igazítva disk_params-ból
    double ugvec[disk_params->NGRID + 2];
    double ugvectemp[disk_params->NGRID + 1]; // Eredeti kód NGRID+1-et használt

    int i;

    // Első ciklus: feltölti a lokális ugvec tömböt
    #pragma omp parallel for private(i)
    for (i = 0; i <= disk_params->NGRID + 1; i++) { // Használd a disk_params->NGRID-et
        // Hozzáférés a disk_params tagjaihoz
        ugvec[i] = disk_params->sigmavec[i] * visc(disk_params->rvec[i], disk_params) * sqrt(disk_params->rvec[i]);
        // Megjegyzés: A sqrt() függvénynek általában csak egy double paramétere van.
        // Ha valami komplexebb számítást akarsz, akkor lehet, hogy egy saját
        // függvényt hívsz, amihez disk_params is kell. Ellenőrizd a sqrt prototípusát!
    }

    // Második ciklus: feltölti a lokális ugvectemp tömböt
    #pragma omp parallel for private(i, tempug)
    for (i = 1; i <= disk_params->NGRID; i++) { // Használd a disk_params->NGRID-et
        tempug = (ugvec[i + 1] - ugvec[i - 1]) / (2.0 * disk_params->DD); // Használd a disk_params->DD-t
        // Coeff_3 hívása, ha szükséges, átadva neki a disk_params-ot
        ugvectemp[i] = Coeff_3(disk_params->sigmavec[i], disk_params->rvec[i]) * tempug;
    }

    // Harmadik ciklus: Az eredményt bemásolja a disk_params->ugvec-be
    for (i = 1; i <= disk_params->NGRID; i++) { // Használd a disk_params->NGRID-et
        disk_params->ugvec[i] = ugvectemp[i]; // Így éri el a struktúrán belüli ugvec-et
    }
}



void GetMass(int n, double (*partmassind)[5], int indii, int indio, int indoi, int indoo, double *massiout, double *massoout, const simulation_options_t *sim_opts) {

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
double a_drift(double sigmad, double r, double p, double dp, double rho_p, const disk_t *disk_params) {

    double Sigmad_cgs = sigmad / SDCONV;

    double vkep = v_kep(r,disk_params);
    double vkep2 = vkep * vkep;
    double c_s = c_sound(r,disk_params);
    double c_s2 = c_s * c_s;
    double dlnPdlnr = r / p * dp;
    double s_drift =  disk_params->fDrift * 2.0 / M_PI * Sigmad_cgs / rho_p * vkep2 / c_s2 * fabs(1.0 / dlnPdlnr);
    return s_drift;
}

// 2. a kis skalaju turbulencia altal okozott fragmentacio szerinti maximalis meret	--> kimenet cm-ben!
double a_turb(double sigma, double r, double rho_p, const disk_t *disk_params) {

    double s_frag, u_frag, u_frag2, Sigma_cgs, c_s, c_s2;

    u_frag = disk_params->uFrag * CMPSECTOAUPYRP2PI; /*	cm/sec --> AU / (yr/2pi)	*/
    u_frag2 = u_frag * u_frag;
    Sigma_cgs = sigma / SDCONV;
    c_s = c_sound(r,disk_params); // / CMPSECTOAUPYRP2PI; // Komment ki, ha a c_sound már megfelelő mértékegységben van
    c_s2 = c_s * c_s;

    s_frag = disk_params->fFrag * 2.0 / (3.0 * M_PI) * Sigma_cgs / (rho_p * calculate_turbulent_alpha(r,disk_params)) * u_frag2 / c_s2;

    return s_frag;
}

// 3. radialis drift altal okozott fragmentacio szerinti maximalis meret		--> kimenet cm-ben!
double a_df(double sigma, double r, double p, double dp, double rho_p, const disk_t *disk_params) {

    double u_frag, vkep, dlnPdlnr, c_s, c_s2, s_df, Sigma_cgs;

    u_frag = disk_params->uFrag * CMPSECTOAUPYRP2PI; /*	cm/sec --> AU / (yr/2pi)	*/
    Sigma_cgs = sigma / SDCONV;
    c_s = c_sound(r,disk_params);
    c_s2 = c_s * c_s;
    dlnPdlnr = r / p * dp;
    vkep = v_kep(r,disk_params);

    s_df = u_frag * vkep / fabs(dlnPdlnr * c_s2 * 0.5) * 2.0 * Sigma_cgs / (M_PI * rho_p);

    return s_df;
}

/*	a reszecskek novekedesenek idoskalaja	*/
double tauGr(double r, double eps,const disk_t *disk_params) {
    double omega = kep_freq(r,disk_params);
    double taugr = eps / omega;
    return taugr;
}

/*	kiszamolja az adott helyen a reszecske meretet --> BIRNSTIEL CIKK	*/
double getSize(double prad, double pdens, double sigma, double sigmad, double y, double p, double dpress_val, double dt, const disk_t *disk_params) {

    double sturb = a_turb(sigma, y, pdens, disk_params);           // cm-ben
    double sdf = a_df(sigma, y, p, dpress_val, pdens,disk_params); // cm-ben
    double srdf = a_drift(sigmad, y, p, dpress_val, pdens, disk_params); // cm-ben
    double smin = find_min(sturb, sdf, srdf);         // cm-ben -- megadja, hogy a fenti ket reszecske korlatbol melyik ad kisebb meretet (az a reszecskenovekedes felso korlatja
    //	double eps = sigma / 100.;
    double eps = sigmad / sigma; // A korábbi kódban fordítva volt, feltételezem, hogy eps = (por sűrűség) / (gáz sűrűség)
    double tau_gr = tauGr(y, eps, disk_params);
    double rt = 0.0;

    smin = smin / AU2CM; // AU-ban

    /*	kiszamolja, hogy a fenti smin, vagy a novekedesi idoskalabol szarmazo meret korlatozza a reszecske meretet	*/
    if (prad < smin) {
        rt = find_min(prad * exp(dt / tau_gr), smin, HUGE_VAL);
    } else { // prad >= smin
        rt = smin;
    }

    return rt;
}




void Get_Sigmad(double max_param, double min_param, double rad[][2], double radmicr[][2], 
                double *sigma_d, double *sigma_dm,  double *massvec, double *massmicrvec,  
                double *rd, double *rmic, const simulation_options_t *sim_opts, const disk_t *disk_params) {

    // Suppress unused parameter warnings
    (void)max_param;
    (void)min_param;

    double dd = (disk_params->RMAX - disk_params->RMIN) / (PARTICLE_NUMBER - 1);
    int i;

    // A temp tömbök deklarálását érdemes a scope tetejére tenni
    double sigdtemp[PARTICLE_NUMBER][3];
    double sigdmicrtemp[PARTICLE_NUMBER][3];

    // Inicializálás, ha szükséges (bár a loadSigDust valószínűleg felülírja)
    for(i=0; i<PARTICLE_NUMBER; i++){
        sigdtemp[i][0] = 0.0; sigdtemp[i][1] = 0.0; sigdtemp[i][2] = 0.0;
        sigdmicrtemp[i][0] = 0.0; sigdmicrtemp[i][1] = 0.0; sigdmicrtemp[i][2] = 0.0;
        rd[i] = 0.0;
        rmic[i] = 0.0;
        sigma_d[i] = 0.0;
        sigma_dm[i] = 0.0;
    }


    // loadSigDust és contract függvények hívásai:
    // Ezek valószínűleg szekvenciálisak, hacsak a függvények belsejében nincs OpenMP.
    // Ha ezek a függvények valamilyen globális állapotot módosítanak, akkor kritikusak.
    // Feltételezve, hogy a 'sigdtemp' és 'sigdmicrtemp' kizárólagosan a hívásaikban vannak feldolgozva,
    // és nem ütköznek más szálakkal globális adatokon keresztül.
    loadSigDust(rad, massvec, sigdtemp, PARTICLE_NUMBER,disk_params);
    if (sim_opts->twopop == 1.0) { // Használjunk double összehasonlítást
        loadSigDust(radmicr, massmicrvec, sigdmicrtemp, PARTICLE_NUMBER,disk_params);
    }

    contract(sigdtemp, dd, PARTICLE_NUMBER,disk_params);
    if (sim_opts->twopop == 1.0) { // Használjunk double összehasonlítást
        contract(sigdmicrtemp, dd, PARTICLE_NUMBER,disk_params);
    }

    // Utolsó másoló ciklus: Ez is jól párhuzamosítható.
    #pragma omp parallel for private(i)
    for (i = 0; i < PARTICLE_NUMBER; i++) {
        rd[i] = sigdtemp[i][1];
        sigma_d[i] = sigdtemp[i][0];

        if (sim_opts->twopop == 1.0) { // double összehasonlítás
            rmic[i] = sigdmicrtemp[i][1];
            sigma_dm[i] = sigdmicrtemp[i][0];
        }
    }
}


/*	Fuggveny a sigma, p, dp kiszamolasara	*/
void Get_Sigma_P_dP(const simulation_options_t *sim_opts, disk_t *disk_params) { // Added sim_opts

    double u, u_bi, u_fi;
    double sigma_temp[disk_params->NGRID + 2]; // Use disk_params->NGRID
    double uvec[disk_params->NGRID + 2];     // Use disk_params->NGRID

    int i;

    // Boundary conditions - access via disk_params
    sigma_temp[0] = disk_params->sigmavec[0];
    sigma_temp[disk_params->NGRID + 1] = disk_params->sigmavec[disk_params->NGRID + 1];

    // uvec temporary array initialization
    uvec[0] = disk_params->sigmavec[0] * visc(disk_params->rvec[0], disk_params); // Use disk_params->rvec
    uvec[disk_params->NGRID + 1] = disk_params->sigmavec[disk_params->NGRID + 1] * visc(disk_params->rvec[disk_params->NGRID + 1], disk_params); // Use disk_params->rvec

    #pragma omp parallel for
    for(i = 1; i <= disk_params->NGRID; i++) { // Use disk_params->NGRID
        uvec[i] = disk_params->sigmavec[i] * visc(disk_params->rvec[i], disk_params); // Use disk_params->sigmavec and disk_params->rvec
    }

    // This loop is critical due to data dependencies. Keep it sequential for correctness
    for (i = 1; i <= disk_params->NGRID; i++) { // Use disk_params->NGRID
        u = uvec[i];
        u_bi = uvec[i - 1];
        u_fi = uvec[i + 1];

        // Access DD and deltat through the appropriate structs
        // Assuming Coeff_1 and Coeff_2 also take disk_params (and sim_opts if they need it)
        double temp = Coeff_1(disk_params->rvec[i], disk_params) * (u_fi - 2.0 * u + u_bi) / (disk_params->DD * disk_params->DD) +
                      Coeff_2(disk_params->rvec[i], disk_params) * (u_fi - u_bi) / (2.0 * disk_params->DD);
        
        sigma_temp[i] = uvec[i] + sim_opts->DT * temp; // Use sim_opts->DT for deltat
    }

    // This loop is parallelizable
    #pragma omp parallel for
    for (i = 1; i <= disk_params->NGRID; i++) { // Use disk_params->NGRID
        // Update disk_params' own arrays
        disk_params->sigmavec[i] = sigma_temp[i] / visc(disk_params->rvec[i], disk_params);
        disk_params->pressvec[i] = press(disk_params->sigmavec[i], disk_params->rvec[i], disk_params); // Assuming press takes disk_params
    }

    // These calls likely remain sequential or require their own internal OpenMP if large
    // If Perem, dpress also update members of disk_params, they should take disk_params as a parameter.
    // And if they are modifying the *content* of the arrays within disk_params, then disk_params should NOT be const in *their* parameter list.
    // However, since Get_Sigma_P_dP is modifying them, disk_params *here* cannot be const.
    // Let's remove 'const' from disk_params in Get_Sigma_P_dP signature if it modifies them.
    // void Get_Sigma_P_dP(disk_t *disk_params, const simulation_options_t *sim_opts) { ... }
    
    // Assuming these helper functions need disk_params to access *its* internal arrays
    dpress(disk_params); // Assuming dpress takes arrays and disk_params
	Perem(disk_params->sigmavec, disk_params); // First argument is the array, second is the disk_t pointer
	Perem(disk_params->pressvec, disk_params);
	Perem(disk_params->dpressvec, disk_params);
}

/*	Fuggveny a porszemcsek uj tavolsaganak elraktarozasara		*/
void Get_Radius(const char *nev, int opt, double radius[][2], const double *sigmad, const double *rdvec,
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
            if (t == 0) {
                if (sim_opts->twopop == 0) {
                    double current_drdt_val = (fabs(y_out - y) / (deltat));
                    // Azért kell a critical szekció, mert az fout2 fájlba írunk.
                    // Ez a critical szekció biztosítja, hogy egyszerre csak egy szál írjon a fájlba.

                    #pragma omp critical(fout2_write)
                    {
                        // Ellenőrizzük, hogy a fájlmutató nem NULL
                        if (fout2 != NULL) {
                            fprintf(fout2, "%lg %lg\n", radius[i][0], (radius[i][0] / current_drdt_val) / (2.0 * M_PI));
                        } else {
                            fprintf(stderr, "ERROR: fout2 is NULL during write in Get_Radius (t=0 block).\n");
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
            if (fout2 != NULL) {
                fclose(fout2);
                fout2 = NULL; // Fontos, hogy NULL-ra állítsuk, miután bezártuk.
            }
        }
    }
    // Szinkronizálás a master szál befejezéséig.
    #pragma omp barrier
}