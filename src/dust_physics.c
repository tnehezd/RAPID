// src/dust_physics.c

#include <stdio.h>
#include <stdlib.h> // For calloc, free
#include <string.h>
#include <errno.h>
#include <math.h>   // For M_PI, floor, fmod
#include <omp.h>    // For OpenMP directives

#include "disk_model.h"
#include "dust_physics.h" // Ez az új Get_Radius prototípust fogja tartalmazni
#include "config.h"
#include "simulation_types.h"
#include "globals.h"
#include "io_utils.h"
#include "simulation_core.h"
#include "utils.h"
#include "particle_data.h" // Szükséges a dust_particle_t és ParticleData_t miatt

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


// AU -> cm
#define LOCAL_AU_TO_CM 1.495978707e13 // cm / AU

// Naptömeg -> gramm
#define LOCAL_SUN_MASS_TO_GRAMS 1.989e33 // g / M_Sun

// "Belső időegység" (1 év / 2pi) -> másodperc
// Az 1 AU-nál az 1 Naptömeg körüli keringési idő 1 év.
// Omega_K(1AU, 1M_Sun) = sqrt(G*M_Sun/AU^3).
// Ha G=1 és M_Sun=1, AU=1, akkor Omega_K = 1.
// A periódus T = 2pi / Omega_K. Ha Omega_K = 1, akkor T = 2pi (belső időegység).
// 1 belső időegység = 1 év / (2pi)
// 1 év = 3.1536e7 másodperc
// Tehát 1 belső időegység = 3.1536e7 / (2 * M_PI) másodperc
#define LOCAL_INTERNAL_TIME_TO_SEC (3.1536e7 / (2.0 * M_PI)) // sec / belső időegység


// A Stokes_Number függvény definíciója
// r_au: A sugár AU-ban érkezik (ahogy az rvec-ben is tárolva van)
// sigma_msun_au2: A gáz felületi sűrűsége M_Sun/AU^2-ben érkezik (ahogy a sigmavec-ben is tárolva van)
// pradius_cm: A részecske sugara CM-ben érkezik (ahogy a disk_params->PDENSITY is g/cm^3)
// disk_params: A diszk paraméterei, amiben a G_GRAV_CONST, STAR_MASS, PDENSITY vannak.
double Stokes_Number(double pradius_cm, double sigma_msun_au2, double r_au, const disk_t *disk_params) {

    // --- Konverziók CGS egységrendszerbe a számításhoz ---

    // Konvertáljuk a sugárt AU-ból cm-re
    double r_cm = r_au * LOCAL_AU_TO_CM;

    // Konvertáljuk a gáz felületi sűrűségét M_Sun/AU^2-ből g/cm^2-be
    double sigma_g_cm2 = sigma_msun_au2 * LOCAL_SUN_MASS_TO_GRAMS / (LOCAL_AU_TO_CM * LOCAL_AU_TO_CM);

    // A részecske belső sűrűsége már g/cm^3-ben van (disk_params->PDENSITY), tehát ok.
    double rho_s_g_cm3 = disk_params->PDENSITY;

    // --- Számoljuk ki a szükséges gáz paramétereket CGS-ben ---
    // Ehhez a calculate_keplerian_angular_velocity és calculate_local_sound_speed
    // függvényeknek is CGS-ben kellene eredményt adniuk,
    // VAGY belül kellene konvertálnunk őket.

    // HA a calculate_keplerian_angular_velocity (és a többi) G=1, M_star=1 egységben számol:
    // Akkor az eredménye rad / belső időegység lesz.
    // Ezt kell átkonvertálni rad/s-be.
    double omega_k_internal = calculate_keplerian_angular_velocity(r_au, disk_params); // rad / belső időegység
    double omega_k_rad_s = omega_k_internal / LOCAL_INTERNAL_TIME_TO_SEC; // rad / s

    // Hasonlóan a hangsebességgel. Ha AU/belső időegységben van, konvertálni cm/s-be:
    double csound_internal = calculate_local_sound_speed(r_au, disk_params); // AU / belső időegység
    double csound_cm_s = csound_internal * LOCAL_AU_TO_CM / LOCAL_INTERNAL_TIME_TO_SEC; // cm / s

    // A skála magasság is AU-ban van (calculate_scale_height). Konvertálni cm-be:
    double H_scale_au = calculate_scale_height(r_au, disk_params); // AU
    double H_scale_cm = H_scale_au * LOCAL_AU_TO_CM; // cm

    // --- Ellenőrzés nullával való osztás előtt ---
    if (sigma_g_cm2 <= 0.0 || csound_cm_s <= 0.0 || H_scale_cm <= 0.0 || omega_k_rad_s <= 0.0) {
        fprintf(stderr, "ERROR [Stokes_Number]: Invalid gas/disk parameters detected in CGS (sigma_g_cm2: %.2e, csound_cm_s: %.2e, H_scale_cm: %.2e, omega_k_rad_s: %.2e, r_au: %.2e).\n",
                sigma_g_cm2, csound_cm_s, H_scale_cm, omega_k_rad_s, r_au);
        return 0.0; // Hibakezelés: visszaad 0.0-t vagy egy nagyon nagy értéket
    }

    // --- A Stokes-szám Képlete (Epstein-ellenállás) CGS egységekben ---
    // St = (rho_s * a / Sigma_g) * (Omega_K * r / c_s)
    // Megjegyzés: a (Omega_K * r / c_s) tag kb. (v_Kepler / c_s).
    // Mivel c_s = Omega_K * H, ezért (Omega_K * r / c_s) = (r / H) is.
    // Használhatjuk az (r / H) formát is, ha az egyszerűbb a konverziók után.
    // Most az Omega_K * r / c_s formát használom, mivel az is számolt paraméter.

    double Stokes = (rho_s_g_cm3 * pradius_cm / sigma_g_cm2) * (omega_k_rad_s * r_cm / csound_cm_s);

    return Stokes;
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

// A Get_Sigmad függvény módosított implementációja
void Get_Sigmad(const ParticleData_t *p_data, disk_t *disk_params, const simulation_options_t *sim_opts) {
    // Ideiglenes tömbök a calculate_dust_surface_density_profile számára.
    // Ezeket a disk_params->sigmadustvec fogja felvenni a végén.
    // A disk_params->rvec már adott.
    double *sigma_d_temp = (double *)calloc(disk_params->NGRID, sizeof(double));
    double *sigma_dm_temp = (double *)calloc(disk_params->NGRID, sizeof(double));

    if (sigma_d_temp == NULL || sigma_dm_temp == NULL) {
        fprintf(stderr, "ERROR: Memory allocation failed for sigma_d_temp or sigma_dm_temp in Get_Sigmad.\n");
        free(sigma_d_temp);
        free(sigma_dm_temp);
        exit(EXIT_FAILURE);
    }

    // A radin[][2] és massin[] tömbök előkészítése a calculate_dust_surface_density_profile számára.
    // Ennek a függvénynek a bemenete a régi formátumban van (double[][2]),
    // ezért átmenetileg konvertálnunk kell a dust_particle_t adatokból.

    // Csak akkor allokálunk és töltünk fel, ha valóban vannak részecskék.
    if (p_data->num_particles_pop1 > 0) {
        double (*rad_pop1_2d)[2] = (double (*)[2])calloc(p_data->num_particles_pop1, sizeof(double[2]));
        double *mass_pop1_1d = (double *)calloc(p_data->num_particles_pop1, sizeof(double));

        if (rad_pop1_2d == NULL || mass_pop1_1d == NULL) {
            fprintf(stderr, "ERROR: Memory allocation failed for Pop1 temp arrays in Get_Sigmad.\n");
            free(sigma_d_temp); free(sigma_dm_temp);
            free(rad_pop1_2d); free(mass_pop1_1d);
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < p_data->num_particles_pop1; i++) {
            rad_pop1_2d[i][0] = p_data->particles_pop1[i].distance_au;
            rad_pop1_2d[i][1] = p_data->particles_pop1[i].current_size_cm; // Vagy valamilyen dummy, ha nem használnánk a [][1] oszlopot
            mass_pop1_1d[i] = p_data->particles_pop1[i].initial_mass_msun;
        }

        calculate_dust_surface_density_profile(sigma_d_temp, disk_params->rvec,
                                               rad_pop1_2d, mass_pop1_1d,
                                               p_data->num_particles_pop1, disk_params->NGRID, disk_params);

        free(rad_pop1_2d);
        free(mass_pop1_1d);
    }

    // Ha two-pop szimuláció van, a mikronos porra is kiszámoljuk.
    if (sim_opts->twopop == 1.0 && p_data->num_particles_pop2 > 0) {
        double (*rad_pop2_2d)[2] = (double (*)[2])calloc(p_data->num_particles_pop2, sizeof(double[2]));
        double *mass_pop2_1d = (double *)calloc(p_data->num_particles_pop2, sizeof(double));

        if (rad_pop2_2d == NULL || mass_pop2_1d == NULL) {
            fprintf(stderr, "ERROR: Memory allocation failed for Pop2 temp arrays in Get_Sigmad.\n");
            free(sigma_d_temp); free(sigma_dm_temp);
            free(rad_pop2_2d); free(mass_pop2_1d);
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < p_data->num_particles_pop2; i++) {
            rad_pop2_2d[i][0] = p_data->particles_pop2[i].distance_au;
            rad_pop2_2d[i][1] = p_data->particles_pop2[i].current_size_cm;
            mass_pop2_1d[i] = p_data->particles_pop2[i].initial_mass_msun;
        }

        calculate_dust_surface_density_profile(sigma_dm_temp, disk_params->rvec, // disk_params->rvec is good for this too
                                               rad_pop2_2d, mass_pop2_1d,
                                               p_data->num_particles_pop2, disk_params->NGRID, disk_params);

        free(rad_pop2_2d);
        free(mass_pop2_1d);
    }

    // Másoljuk át a kiszámolt sűrűségeket a disk_params megfelelő tömbjébe.
    // Itt feltételezem, hogy a disk_params tartalmazza a sigmadustvec tömböt.
    // Ha van külön tömb a mikronos pornak, azt is itt kell feltölteni (pl. sigmadustmicrvec).
    for (int i = 0; i < disk_params->NGRID; i++) {
        disk_params->sigmadustvec[i] = sigma_d_temp[i];
        // Ha van külön tömb a mikronos pornak a disk_t-ben, azt is itt töltenénk fel:
        // disk_params->sigmadustmicrvec[i] = sigma_dm_temp[i];
    }

    free(sigma_d_temp);
    free(sigma_dm_temp);
}


// Módosított Get_Radius függvény implementációja
void Get_Radius(dust_particle_t *particles_array, int num_particles, double deltat, double t,
                const simulation_options_t *sim_opts, const disk_t *disk_params) {

    // DEBUG PRINTS - Keep these for now to ensure correct pointer values
    fprintf(stderr, "DEBUG [Get_Radius]: Entering function.\n");
    fprintf(stderr, "DEBUG [Get_Radius]: particles_array address: %p, num_particles: %d\n",
            (void*)particles_array, num_particles);
    fprintf(stderr, "DEBUG [Get_Radius]: sim_opts address: %p\n", (void*)sim_opts);
    fprintf(stderr, "DEBUG [Get_Radius]: disk_params address: %p\n", (void*)disk_params);

    if (disk_params != NULL) {
        fprintf(stderr, "DEBUG [Get_Radius]: disk_params->RMIN: %lg, disk_params->RMAX: %lg\n",
                disk_params->RMIN, disk_params->RMAX);
        fprintf(stderr, "DEBUG [Get_Radius]: disk_params->sigmadustvec: %p\n", (void*)disk_params->sigmadustvec);
    } else {
        fprintf(stderr, "ERROR [Get_Radius]: disk_params pointer IS NULL on entry!\n");
        exit(EXIT_FAILURE); // Exit if disk_params is NULL
    }
    // End DEBUG PRINTS

    if (particles_array == NULL || num_particles <= 0) {
        // fprintf(stderr, "WARNING [Get_Radius]: No particles or null array passed. Skipping calculations.\n");
        return;
    }

    const double *aggregated_sigmad_ptr = disk_params->sigmadustvec;
    const double *aggregated_rdvec_ptr = disk_params->rvec; // A grid sugártömb
    int num_grid_points = disk_params->NGRID;
    double grid_dd = disk_params->DD; // disk_params->DD a rács lépésköze
    double grid_rmin = disk_params->RMIN;

    // **ITT A PÁRHUZAMOSÍTÁS!**
    #pragma omp parallel for
    for (int i = 0; i < num_particles; ++i) {
        // If a particle falls outside RMIN or RMAX, or its size/distance is invalid (0.0 or less),
        // we set its properties to 0.0 and SKIP further calculations for THIS particle in THIS timestep.
        if (particles_array[i].distance_au <= disk_params->RMIN ||
            particles_array[i].distance_au >= disk_params->RMAX ||
            particles_array[i].current_size_cm <= 0.0 || // Added check for zero or negative size
            particles_array[i].distance_au <= 0.0) {     // Added check for zero or negative distance

            particles_array[i].distance_au = 0.0;
            particles_array[i].current_size_cm = 0.0;
            particles_array[i].size_reciprocal = 0.0;
            // No need to set initial_mass_msun to 0 here, as it's typically set once at init.
            // If a particle is "dead" its mass contribution to density calculation should be zero.
            // That's handled in Get_Sigmad when preparing rad_pop_2d and mass_pop_1d arrays.
            continue; // IMPORTANT: Skip to the next particle in the loop
        }

        double current_r = particles_array[i].distance_au;
        double current_particle_size = particles_array[i].current_size_cm;

        double new_r = 0.0; // int_step will update this
        double new_size = current_particle_size; // int_step will update this (if growth=1)

        // Calling int_step with current data and freshly calculated sigmad values
        int_step(t,                                  // double time
                 current_particle_size,              // double psize
                 current_r,                          // double current_r
                 deltat,                             // double step
                 &new_r,                             // double *new_prad_ptr (pointer to new radial position)
                 &new_size,                          // double *new_psize_ptr (pointer to new particle size)
                 aggregated_sigmad_ptr,              // const double *aggregated_sigmad
                 aggregated_rdvec_ptr,               // const double *aggregated_rdvec
                 num_grid_points,                    // int num_grid_points
                 grid_dd,                            // double grid_dd
                 grid_rmin,                          // double grid_rmin
                 disk_params,                        // const disk_t *disk_params
                 sim_opts);                          // const simulation_options_t *sim_opts


        // Update particle's position and size
        particles_array[i].distance_au = new_r;

        // Only update size if growth is enabled
        if (sim_opts->growth == 1.0) {
            particles_array[i].current_size_cm = new_size;
            // Update reciprocal value too
            if (particles_array[i].current_size_cm > 0.0) {
                particles_array[i].size_reciprocal = 1.0 / particles_array[i].current_size_cm;
            } else {
                particles_array[i].size_reciprocal = 0.0;
            }
        }
        // If growth is not enabled, particle size remains unchanged.
    }
}