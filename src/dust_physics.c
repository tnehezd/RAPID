// src/dust_physics.c
#include <stdio.h>
#include <stdlib.h> // For calloc, free
#include <string.h>
#include <errno.h>
#include <math.h>   // For M_PI, floor, fmod
#include <omp.h>    // For OpenMP directives

#include "disk_model.h"
#include "dust_physics.h"
#include "config.h"
#include "simulation_types.h"
#include "globals.h"
#include "io_utils.h"
#include "simulation_core.h"
#include "utils.h"

// Por felületi sűrűség számítása (Lagrange-ról Euler-re)
void calculate_dust_surface_density_profile(
    double *output_sigma_d_grid,       // Kimenet: A kiszámított por felületi sűrűség (gridenként)
    double *output_r_grid_centers,     // Kimenet: A grid cellák középpontjainak radiális pozíciói
    const double radin[][2],           // Bemenet: Részecskék radiális pozíciója (const)
    const double *massin,              // Bemenet: Részecskék tömege (const)
    int n_particles,                   // Bemenet: A részecskék száma
    int n_grid_cells,                  // Bemenet: A grid cellák száma, amire a sűrűséget számoljuk
    const disk_t *disk_params          // Bemenet: A korong paraméterei (const)
) {
    // Memória allokáció és inicializálás

    // Ideiglenes tömb a binekben összegyűjtött tömeg tárolására.
    // Fontos: a CIC miatt szükség van a NGRID + 2 méretre (egy-egy plusz bin a széleken)
    double *total_mass_in_bins = (double *)calloc(n_grid_cells + 2, sizeof(double));
    if (total_mass_in_bins == NULL) {
        fprintf(stderr, "ERROR: Memory allocation failed for total_mass_in_bins in calculate_dust_surface_density_profile.\n");
        exit(EXIT_FAILURE);
    }

    // Tömb a bin-ek pontos területének tárolására.
    double *bin_areas = (double *)calloc(n_grid_cells + 2, sizeof(double));
    if (bin_areas == NULL) {
        fprintf(stderr, "ERROR: Memory allocation failed for bin_areas in calculate_dust_surface_density_profile.\n");
        free(total_mass_in_bins); // Felszabadítjuk az első allokált memóriát is hiba esetén
        exit(EXIT_FAILURE);
    }

    // A rács celláinak szélessége (delta r).
    double dr_cell_width = (disk_params->RMAX - disk_params->RMIN) / (double)n_grid_cells;

    // 1. lépés: A grid cellák területének kiszámítása és a középpontok feltöltése
    for (int j = 0; j < n_grid_cells; j++) {
        double r_inner = disk_params->RMIN + j * dr_cell_width;
        double r_outer = disk_params->RMIN + (j + 1) * dr_cell_width;

        // A gyűrű területe: pi * (R_outer^2 - R_inner^2)
        bin_areas[j] = M_PI * (r_outer * r_outer - r_inner * r_inner);

        // Kitöltjük a kimeneti rácspontok sugaraival (cella középpontok)
        output_r_grid_centers[j] = r_inner + 0.5 * dr_cell_width;

        // Inicializáljuk a kimeneti sűrűség tömböt 0-val.
        output_sigma_d_grid[j] = 0.0;
    }

    // 2. lépés: Részecskék tömegének szétosztása a binek között (Cloud-in-Cell - CIC)
    // Végigmegyünk minden egyes porrészecskén. Párhuzamosítható.
    #pragma omp parallel for
    for (int i = 0; i < n_particles; i++) {
        double current_r = radin[i][0]; // A részecske aktuális radiális pozíciója
        double current_mass = massin[i]; // A részecske tömege

        // Csak azokat a részecskéket vesszük figyelembe, amelyek a szimulációs tartományon belül vannak.
        // Ha egy részecske RMIN-en kívülre került, az kiesett, és nem járul hozzá a sűrűséghez.
        if (current_r < disk_params->RMIN || current_r >= disk_params->RMAX) {
            continue; // Ugrás a következő részecskére
        }

        // Kiszámítjuk a részecske "normált" radiális pozícióját a rácshoz képest.
        double r_normalized = (current_r - disk_params->RMIN) / dr_cell_width;

        // Megkeressük az alsó bin indexét, amiben a részecske elhelyezkedik.
        int lower_bin_idx = (int)floor(r_normalized);

        // Kiszámítjuk, hogy a részecske mennyire van "belül" az alsó binben (0 és 1 közötti érték).
        double fraction_in_lower_bin = r_normalized - lower_bin_idx;

        // Hozzáadjuk a tömeg egy részét az alsó binhez.
        // Critical szekció szükséges, mert több szál is írhat ugyanabba a total_mass_in_bins elembe.
        if (lower_bin_idx >= 0 && lower_bin_idx < n_grid_cells) {
            #pragma omp atomic
            total_mass_in_bins[lower_bin_idx] += current_mass * (1.0 - fraction_in_lower_bin);
        }

        // Hozzáadjuk a tömeg másik részét a felső binhez (ha van ilyen).
        // Fontos: a total_mass_in_bins tömb mérete n_grid_cells + 2, így az indexek
        // 0-tól n_grid_cells + 1-ig érvényesek. lower_bin_idx + 1 maximálisan n_grid_cells lehet.
        // Tehát a feltételnek azt kell biztosítania, hogy lower_bin_idx + 1 ne lépje túl a tömb végét.
        if ((lower_bin_idx + 1) < (n_grid_cells + 2)) { // Javított felső határ az tömb eléréshez
            #pragma omp atomic
            total_mass_in_bins[lower_bin_idx + 1] += current_mass * fraction_in_lower_bin;
        }
    }

    // 3. lépés: Sűrűség számítása a gyűjtött tömegekből és a kimeneti tömb feltöltése
    // Végigmegyünk minden grid cellán. Párhuzamosítható.
    #pragma omp parallel for
    for (int j = 0; j < n_grid_cells; j++) {
        // Ellenőrizzük, hogy a bin_areas[j] ne legyen nulla, elkerülve a lebegőpontos hibát.
        if (bin_areas[j] > 1e-18) { // Nagyon kicsi pozitív értékkel hasonlítjuk össze
            output_sigma_d_grid[j] = total_mass_in_bins[j] / bin_areas[j];
        } else {
            output_sigma_d_grid[j] = 0.0; // Ha a terület nulla, a sűrűség is nulla
        }
    }

    // Memória felszabadítás
    free(total_mass_in_bins);
    free(bin_areas);
}

// Gáz turbulencia és por interakciók
/* alpha turbulens paraméter kiszámolása --> alfa csökkentése alpha_r-rel    */
double calculate_turbulent_alpha(double r, const disk_t *disk_params) {
    double alpha_r;
    alpha_r = 1.0 - 0.5 * (1.0 - disk_params->a_mod) * (tanh((r - disk_params->r_dze_i) / disk_params->Dr_dze_i) + tanh((disk_params->r_dze_o - r) / disk_params->Dr_dze_o));
    return alpha_r * disk_params->alpha_visc;
}

/* kiszamolja az adott reszecskehez tartozo Stokes szamot    */
/* St = rho_particle * radius_particle * PI / (2 * sigma)    */
double Stokes_Number(double pradius, double sigma, disk_t *disk_params) { /* in the Epstein drag regime    */
    return disk_params->PDENSITYDIMLESS * pradius * M_PI / (2.0 * sigma);
}

// Tömeggyűjtés a holt zónákban
void GetMass(int n, double (*partmassind)[5], int indii, int indio, int indoi, int indoo, double *massiout, double *massoout, const simulation_options_t *sim_opts) {

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
                    // Critical szekció, mert módosítjuk a partmassind[i][3] értéket és a massitemp-et
                    #pragma omp critical(inner_dze_update)
                    {
                        partmassind[i][3] = 1.0;
                        massitemp += partmassind[i][0]; // Tömeg hozzáadása a [0] indexről
                    }
                }
            }

            // --- Külső DZE ---
            // Ellenőrizzük, hogy a részecske grid indexe a külső DZE tartományában van-e
            if ((current_r_index >= indoi) && (current_r_index <= indoo)) {
                if (partmassind[i][4] == 0.0) { // Külső DZE flag ellenőrzése (Flag[4])
                    // Critical szekció, mert módosítjuk a partmassind[i][4] értéket és a massotemp-et
                    #pragma omp critical(outer_dze_update)
                    {
                        partmassind[i][4] = 1.0;
                        massotemp += partmassind[i][0]; // Tömeg hozzáadása a [0] indexről
                    }
                }
            }
        }
    } else { // Fix DZE (sim_opts->dzone == 0.0): Nincsenek flag-ek a tömeg felhalmozáshoz
        #pragma omp parallel for private(i) reduction(+:massitemp, massotemp)
        for (i = 0; i < n; i++) {
            int current_r_index = (int)partmassind[i][1]; // A részecske grid indexe

            // --- Belső DZE (Fix) ---
            if ((current_r_index >= indii) && (current_r_index <= indio)) {
                massitemp += partmassind[i][0];
            }

            // --- Külső DZE (Fix) ---
            if ((current_r_index >= indoi) && (current_r_index <= indoo)) {
                massotemp += partmassind[i][0];
            }
        }
    }

    *massiout = massitemp;
    *massoout = massotemp;
}

// Részecskenövekedés (Birnstiel et al. 2012)
//reprezentativ reszecske kezdeti meretenek meghatarozasa
// 1. radialis drift altal meghatarozott maximalis meret            --> kimenet cm-ben!
double a_drift(double sigmad, double r, double p, double dp, double rho_p, const disk_t *disk_params) {

    double Sigmad_cgs = sigmad / GAS_SD_CONV_RATE;

    double vkep = calculate_keplerian_velocity(r,disk_params);
    double vkep2 = vkep * vkep;
    double c_s = calculate_local_sound_speed(r,disk_params);
    double c_s2 = c_s * c_s;
    double dlnPdlnr = r / p * dp;
    double s_drift =  disk_params->fDrift * 2.0 / M_PI * Sigmad_cgs / rho_p * vkep2 / c_s2 * fabs(1.0 / dlnPdlnr);
    return s_drift;
}

// 2. a kis skalaju turbulencia altal okozott fragmentacio szerinti maximalis meret    --> kimenet cm-ben!
double a_turb(double sigma, double r, double rho_p, const disk_t *disk_params) {

    double s_frag, u_frag, u_frag2, Sigma_cgs, c_s, c_s2;

    u_frag = disk_params->uFrag * CM_PER_SEC_TO_AU_PER_YEAR_OVER_2PI; /* cm/sec --> AU / (yr/2pi)    */
    u_frag2 = u_frag * u_frag;
    Sigma_cgs = sigma / GAS_SD_CONV_RATE;
    c_s = calculate_local_sound_speed(r,disk_params);
    c_s2 = c_s * c_s;

    s_frag = disk_params->fFrag * 2.0 / (3.0 * M_PI) * Sigma_cgs / (rho_p * calculate_turbulent_alpha(r,disk_params)) * u_frag2 / c_s2;

    return s_frag;
}

// 3. radialis drift altal okozott fragmentacio szerinti maximalis meret        --> kimenet cm-ben!
double a_df(double sigma, double r, double p, double dp, double rho_p, const disk_t *disk_params) {

    double u_frag, vkep, dlnPdlnr, c_s, c_s2, s_df, Sigma_cgs;

    u_frag = disk_params->uFrag * CM_PER_SEC_TO_AU_PER_YEAR_OVER_2PI; /* cm/sec --> AU / (yr/2pi)    */
    Sigma_cgs = sigma / GAS_SD_CONV_RATE;
    c_s = calculate_local_sound_speed(r,disk_params);
    c_s2 = c_s * c_s;
    dlnPdlnr = r / p * dp;
    vkep = calculate_keplerian_velocity(r,disk_params);

    s_df = u_frag * vkep / fabs(dlnPdlnr * c_s2 * 0.5) * 2.0 * Sigma_cgs / (M_PI * rho_p);

    return s_df;
}

/* a reszecskek novekedesenek idoskalaja    */
double tauGr(double r, double eps,const disk_t *disk_params) {
    double omega = calculate_keplerian_angular_velocity(r,disk_params);
    double taugr = eps / omega;
    return taugr;
}

/* kiszamolja az adott helyen a reszecske meretet --> BIRNSTIEL CIKK    */
double getSize(double prad, double pdens, double sigma, double sigmad, double y, double p, double dpress_val, double dt, const disk_t *disk_params) {

    double sturb = a_turb(sigma, y, pdens, disk_params);     // cm-ben
    double sdf = a_df(sigma, y, p, dpress_val, pdens,disk_params); // cm-ben
    double srdf = a_drift(sigmad, y, p, dpress_val, pdens, disk_params); // cm-ben
    double smin = find_min(sturb, sdf, srdf);         // cm-ben -- megadja, hogy a fenti ket reszecske korlatbol melyik ad kisebb meretet (az a reszecskenovekedes felso korlatja
    //    double eps = sigma / 100.;
    double eps = sigmad / sigma; // A korábbi kódban fordítva volt, feltételezem, hogy eps = (por sűrűség) / (gáz sűrűség)
    double tau_gr = tauGr(y, eps, disk_params);
    double rt = 0.0;

    smin = smin / AU_TO_CM; // AU-ban

    /* kiszamolja, hogy a fenti smin, vagy a novekedesi idoskalabol szarmazo meret korlatozza a reszecske meretet    */
    if (prad < smin) {
        rt = find_min(prad * exp(dt / tau_gr), smin, HUGE_VAL);
    } else { // prad >= smin
        rt = smin;
    }

    return rt;
}

// Por felületi sűrűség lekérése és frissítése
void Get_Sigmad(double max_param, double min_param, double rad[][2], double radmicr[][2],
                double *sigma_d, double *sigma_dm,  double *massvec, double *massmicrvec,
                double *rd, double *rmic, const simulation_options_t *sim_opts, disk_t *disk_params) { // disk_params már nem const!

    // Suppress unused parameter warnings
    (void)max_param;
    (void)min_param;

    // A `dd` már nem szükséges itt, mivel a calculate_dust_surface_density_profile számolja
    // a dr_cell_width-et a disk_params->RMAX, RMIN, NGRID alapján.

    // Calculate dust surface density profile for the main dust population
    calculate_dust_surface_density_profile(sigma_d, rd, rad, massvec, PARTICLE_NUMBER, disk_params->NGRID, disk_params);

    // If two-population simulation is enabled, calculate for micron-sized dust as well
    if (sim_opts->twopop == 1.0) {
        calculate_dust_surface_density_profile(sigma_dm, rmic, radmicr, massmicrvec, PARTICLE_NUMBER, disk_params->NGRID, disk_params);
    }

    // A `contract` függvényre már nincs szükség, mivel a `calculate_dust_surface_density_profile`
    // már a CIC módszerrel számolja a gridenkénti sűrűséget.

    // FONTOS LÉPÉS: Másoljuk át a kiszámolt sűrűséget a disk_params->sigmadustvec-be.
    // A disk_params->rvec már a grid pontokat tartalmazza, azt nem kell másolni ide.
    for (int i = 0; i < disk_params->NGRID; i++) {
        disk_params->sigmadustvec[i] = sigma_d[i];
        // Ha van külön tömb a mikronos pornak a disk_t-ben, azt is itt töltenénk fel:
        // pl: disk_params->sigmadustmicrvec[i] = sigma_dm[i];
        // Mivel nincs ilyen a disk_t-ben, a sigma_dm külön marad, mint Get_Sigmad kimenet.
    }
}

// Részecske radiális mozgásának frissítése
/* Fuggveny a porszemcsek uj tavolsaganak elraktarozasara        */
void Get_Radius(const char *nev, int opt, double radius[][2], const double *sigmad, const double *rdvec,
                double deltat, double t, int n, const simulation_options_t *sim_opts, const disk_t *disk_params){

    int i;
    double y, y_out, prad_new, particle_radius;
    char timescale_out[MAX_PATH_LEN]; // Használjuk a MAX_PATH_LEN definíciót
    HeaderData_t header_data_for_files; // Később inicializáljuk a setup_initial_output_files-ban

    // Fájlkezelés t==0 esetén: ez valószínűleg egyszer történik meg a szimuláció elején,
    // még mielőtt az igazi párhuzamosítás elkezdődne a fő ciklusban.
    // Biztosítjuk, hogy csak egy szál nyissa meg/zárja be a fájlt.
    #pragma omp master
    {
        if (t == 0) {
            HeaderData_t timescale_header_data = {.current_time = 0};
            snprintf(timescale_out, sizeof(timescale_out), "%s/%s", nev, FILE_TIMESCALE);
            timescale_output_file = fopen(timescale_out, "w");
            print_file_header(timescale_output_file, FILE_TYPE_TIMESCALE, &timescale_header_data);

            if (timescale_output_file == NULL) {
                fprintf(stderr, "ERROR: Could not open timescale output file '%s': %s\n", timescale_out, strerror(errno));
                // Itt érdemes valamilyen hibakezelést végezni, pl. kilépni
                exit(EXIT_FAILURE);
            }
        }
    }
    // Szinkronizálás a master szál befejezéséig.
    #pragma omp barrier

    // **ITT A PÁRHUZAMOSÍTÁS!**
    // Az `i` ciklus független iterációkkal rendelkezik, minden szál a saját `radius[i]` elemen dolgozik.
    #pragma omp parallel for private(y, y_out, prad_new, particle_radius)
    for (i = 0; i < n; i++) {
        // Ha a részecske RMIN vagy RMAX kívülre kerül, 0.0-ra állítjuk a pozícióját, és kihagyjuk a további számításokat.
        if (radius[i][0] <= disk_params->RMIN || radius[i][0] >= disk_params->RMAX) {
            radius[i][0] = 0.0;
            // Ha a részecske "kiesett", a méretét is érdemes nullázni,
            // hogy ne vegyen részt a növekedési számításokban.
            radius[i][1] = 0.0;
            continue; // Ugrás a következő részecskére
        }

        y = radius[i][0];
        particle_radius = radius[i][1];

        // Az `int_step` függvénynek szüksége van a `sigmad` és `rdvec` (gridenkénti) tömbre,
        // hogy interpolálja a lokális por sűrűséget.
        int_step(t, particle_radius, sigmad, rdvec, deltat, y, &y_out, &prad_new, disk_params, sim_opts);

        if (t == 0) {
            if (sim_opts->twopop == 0.0) { // Használjunk double összehasonlítást
                double current_drdt_val = (fabs(y_out - y) / (deltat));
                // Azért kell a critical szekció, mert az timescale_output_file fájlba írunk.
                #pragma omp critical(timescale_output_file_write)
                {
                    if (timescale_output_file != NULL) {
                        fprintf(timescale_output_file, " %-15.6lg  %-15.6lg\n", radius[i][0], (radius[i][0] / current_drdt_val) / (2.0 * M_PI));
                    } else {
                        fprintf(stderr, "ERROR: timescale_output_file is NULL during write in Get_Radius (t=0 block).\n");
                    }
                }
            }
        }

        // Frissítjük a részecske pozícióját és méretét
        // sim_opts->growth == 1.0 (feltételezve, hogy a growth flag kontrollálja a növekedést)
        if (sim_opts->growth == 1.0) { // Ha növekedés engedélyezett
            radius[i][1] = prad_new;
            radius[i][0] = y_out;
        } else { // Ha növekedés nincs engedélyezve, csak drift
            radius[i][0] = y_out;
        }
    }

    // A fájl bezárása ismételten egy szál által kell, hogy történjen.
    #pragma omp master
    {
        if (t == 0) { // Csak akkor zárjuk be, ha meg is nyitottuk a t=0 blokkban
            if (timescale_output_file != NULL) {
                fclose(timescale_output_file);
                timescale_output_file = NULL; // Fontos, hogy NULL-ra állítsuk, miután bezártuk.
            }
        }
    }
    // Szinkronizálás a master szál befejezéséig.
    #pragma omp barrier
}