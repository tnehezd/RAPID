// src/dust_physics.c

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <omp.h>
#include <float.h> // for DBL_MAX, or use a very large number

#include "disk_model.h"
#include "dust_physics.h"
#include "config.h"
#include "simulation_types.h"
#include "globals.h"
#include "io_utils.h"
#include "simulation_core.h"
#include "utils.h"
#include "particle_data.h"

// A pici érték (epsilon) definíciója, a nullával való osztás
// és a nem fizikai, negatív értékek elkerülésére.
// A nagyobb stabilitás érdekében megnöveltük az értékét.
#define EPSILON_SIZE 1e-20

// A függvények implementációja

void calculate_dust_surface_density_profile(double *output_sigma_d_grid, double *output_r_grid_centers, const double *particle_radius_au, const double *particle_mass, int n_particles, int n_grid_cells, const disk_t *disk_params) {
    // Memória allokáció és inicializálás a por tömegek számára a cellákban
    double *total_mass_in_bins = (double *)calloc(n_grid_cells + 1, sizeof(double));
    if (total_mass_in_bins == NULL) {
        fprintf(stderr, "HIBA: A memóriafoglalás sikertelen a total_mass_in_bins számára.\n");
        exit(EXIT_FAILURE);
    }

    double dr_cell_width = (disk_params->RMAX - disk_params->RMIN) / (double)n_grid_cells;

    // Initialize output arrays and calculate grid centers and bin areas
    for (int j = 0; j < n_grid_cells; j++) {
        double r_inner = disk_params->RMIN + j * dr_cell_width;
        double r_outer = disk_params->RMIN + (j + 1) * dr_cell_width;
        output_r_grid_centers[j] = r_inner + 0.5 * dr_cell_width;
        output_sigma_d_grid[j] = M_PI * (r_outer * r_outer - r_inner * r_inner); // Store bin areas temporarily here
    }

    // CIC módszer: részecskék tömegének interpolálása a rácsra
    #pragma omp parallel for
    for (int i = 0; i < n_particles; i++) {
        double current_r = particle_radius_au[i];
        
        // Helyes ellenőrzés: A részecske a lemez határain belül van-e.
        if (current_r >= disk_params->RMIN && current_r < disk_params->RMAX) {
            
            // A részecske rácsindexének és a relatív pozíciójának kiszámítása
            double normalized_r = (current_r - disk_params->RMIN) / dr_cell_width;
            int bin_idx = (int)floor(normalized_r);
            double fraction = normalized_r - bin_idx;

            // Az interpoláció csak akkor történik, ha a bin indexek érvényesek
            if (bin_idx >= 0 && bin_idx < n_grid_cells) {
                // Tömeg hozzárendelése az alsó binhez
                double mass_to_lower_bin = particle_mass[i] * (1.0 - fraction);
                #pragma omp atomic
                total_mass_in_bins[bin_idx] += mass_to_lower_bin;
                
                // Tömeg hozzárendelése a felső binhez (ha létezik)
                if (bin_idx + 1 < n_grid_cells) {
                    double mass_to_upper_bin = particle_mass[i] * fraction;
                    #pragma omp atomic
                    total_mass_in_bins[bin_idx + 1] += mass_to_upper_bin;
                }
            }
        }
    }

    // Végleges sűrűség kiszámítása a bin tömegével és területével
    #pragma omp parallel for
    for (int j = 0; j < n_grid_cells; j++) {
        double bin_area = output_sigma_d_grid[j]; // bin areas were stored here
        
        if (bin_area > EPSILON_SIZE) {
            output_sigma_d_grid[j] = total_mass_in_bins[j] / bin_area;
        } else {
            // Ez az eset csak akkor fordulhat elő, ha a rács degenerált, de jó biztonsági ellenőrzés.
            output_sigma_d_grid[j] = 0.0;
        }
        // JAVÍTÁS: A negatív értékek elkerülése, ami a korábbi kimenetben is felbukkant
        if (output_sigma_d_grid[j] < 0.0) {
            output_sigma_d_grid[j] = 0.0;
        }
    }

    free(total_mass_in_bins);
}


double calculate_stokes_number(double pradius_au, double sigma_msun_au2, double r_au, const disk_t *disk_params) {
    // --- Konverzió CGS egységrendszerre a számításhoz ---
    double pradius_cm = pradius_au * AU_TO_CM;
    double r_cm = r_au * AU_TO_CM;
    double sigma_g_cm2 = sigma_msun_au2 * SUN_MASS_TO_GRAMS / (AU_TO_CM * AU_TO_CM);
    double rho_s_g_cm3 = disk_params->PDENSITY;
    double omega_k_internal_units = calculate_keplerian_angular_velocity(r_au, disk_params);
    double omega_k_rad_s = omega_k_internal_units / INTERNAL_TIME_TO_SEC;
    double csound_au_per_timeunit = calculate_local_sound_speed(r_au, disk_params);
    double csound_cm_s = csound_au_per_timeunit * AU_TO_CM / INTERNAL_TIME_TO_SEC;

    // ELLENŐRZÉS: Hozzáadott robusztusabb ellenőrzés az érvénytelen paraméterekre.
    if (sigma_g_cm2 <= EPSILON_SIZE || csound_cm_s <= EPSILON_SIZE || omega_k_rad_s <= EPSILON_SIZE) {
        fprintf(stderr, "HIBA [calculate_stokes_number]: Érvénytelen CGS paraméterek a számításhoz (sigma_g_cm2: %.2e, csound_cm_s: %.2e, omega_k_rad_s: %.2e, r_au: %.2e). Visszatérés 0.0.\n",
                sigma_g_cm2, csound_cm_s, omega_k_rad_s, r_au);
        return 0.0;
    }

    double Stokes = (rho_s_g_cm3 * pradius_cm / sigma_g_cm2) * (omega_k_rad_s * r_cm / csound_cm_s);
    
    // ELLENŐRZÉS: Visszatérés előtt ellenőrizzük, hogy az eredmény NaN-e.
    if (isnan(Stokes) || isinf(Stokes) || Stokes < 0.0) {
        fprintf(stderr, "HIBA [calculate_stokes_number]: NaN, Inf vagy negatív érték keletkezett (%.2e). Bemenet: pradius_au=%.2e, sigma_msun_au2=%.2e, r_au=%.2e.\n", Stokes, pradius_au, sigma_msun_au2, r_au);
        return 0.0;
    }

    return Stokes;
}


double calculate_max_size_from_drift(double sigmad, double r, double p, double dp, double rho_p, const disk_t *disk_params) {
    double Sigmad_cgs = sigmad / GAS_SD_CONV_RATE;
    double vkep = calculate_keplerian_velocity(r,disk_params);
    double vkep2 = vkep * vkep;
    double c_s = calculate_local_sound_speed(r,disk_params);
    double c_s2 = c_s * c_s;
    double dlnPdlnr = r / p * dp;

    // ELLENŐRZÉS: Elkerüljük a 0-val való osztást.
    if (c_s2 <= EPSILON_SIZE || fabs(dlnPdlnr) <= EPSILON_SIZE) {
        fprintf(stderr, "HIBA [calculate_max_size_from_drift]: A nevező túl kicsi vagy nulla. c_s2: %.2e, dlnPdlnr: %.2e.\n", c_s2, dlnPdlnr);
        return 0.0;
    }

    double s_drift = disk_params->fDrift * 2.0 / M_PI * Sigmad_cgs / rho_p * vkep2 / c_s2 * fabs(1.0 / dlnPdlnr);

    if (isnan(s_drift) || isinf(s_drift) || s_drift < 0.0) {
        fprintf(stderr, "HIBA [calculate_max_size_from_drift]: NaN, Inf vagy negatív érték keletkezett (%.2e). Bemenet: sigmad=%.2e, r=%.2e.\n", s_drift, sigmad, r);
        return 0.0;
    }

    return s_drift;
}


double calculate_max_size_from_turbulence(double sigma, double r, double rho_p, const disk_t *disk_params) {
    double s_frag, u_frag, u_frag2, Sigma_cgs, c_s, c_s2, alpha_turb;
    u_frag = disk_params->uFrag * CM_PER_SEC_TO_AU_PER_YEAR_OVER_2PI;
    u_frag2 = u_frag * u_frag;
    Sigma_cgs = sigma / GAS_SD_CONV_RATE;
    c_s = calculate_local_sound_speed(r,disk_params);
    c_s2 = c_s * c_s;

    alpha_turb = calculate_turbulent_alpha(r,disk_params);

    // ELLENŐRZÉS: Elkerüljük a 0-val való osztást.
    if (alpha_turb <= EPSILON_SIZE || c_s2 <= EPSILON_SIZE) {
        fprintf(stderr, "HIBA [calculate_max_size_from_turbulence]: A nevező túl kicsi vagy nulla. alpha_turb: %.2e, c_s2: %.2e.\n", alpha_turb, c_s2);
        return 0.0;
    }
    
    s_frag = disk_params->fFrag * 2.0 / (3.0 * M_PI) * Sigma_cgs / (rho_p * alpha_turb) * u_frag2 / c_s2;

    if (isnan(s_frag) || isinf(s_frag) || s_frag < 0.0) {
        fprintf(stderr, "HIBA [calculate_max_size_from_turbulence]: NaN, Inf vagy negatív érték keletkezett (%.2e). Bemenet: sigma=%.2e, r=%.2e.\n", s_frag, sigma, r);
        return 0.0;
    }

    return s_frag;
}


double calculate_max_size_from_drift_fragmentation(double sigma, double r, double p, double dp, double rho_p, const disk_t *disk_params) {

    double u_frag, vkep, dlnPdlnr, c_s, c_s2, s_df, Sigma_cgs;
    u_frag = disk_params->uFrag * CM_PER_SEC_TO_AU_PER_YEAR_OVER_2PI;
    Sigma_cgs = sigma / GAS_SD_CONV_RATE;
    c_s = calculate_local_sound_speed(r,disk_params);
    c_s2 = c_s * c_s;
    
    // JAVÍTÁS: Robusztusabb ellenőrzés a nullával való osztás elkerülésére a dlnPdlnr számításában.
    if (p <= EPSILON_SIZE || dp == 0.0) {
        fprintf(stderr, "HIBA [calculate_max_size_from_drift_fragmentation]: Érvénytelen nyomás (p=%.2e) vagy nyomásgradiens (dp=%.2e).\n", p, dp);
        return 0.0;
    }
    dlnPdlnr = r / p * dp;
    
    vkep = calculate_keplerian_velocity(r,disk_params);

    // ELLENŐRZÉS: Elkerüljük a 0-val való osztást.
    if (fabs(dlnPdlnr * c_s2 * 0.5) <= EPSILON_SIZE || rho_p <= EPSILON_SIZE) {
        fprintf(stderr, "HIBA [calculate_max_size_from_drift_fragmentation]: A nevező túl kicsi vagy nulla. dlnPdlnr: %.2e, c_s2: %.2e, rho_p: %.2e.\n", dlnPdlnr, c_s2, rho_p);
        return 0.0;
    }

    s_df = u_frag * vkep / fabs(dlnPdlnr * c_s2 * 0.5) * 2.0 * Sigma_cgs / (M_PI * rho_p);

    if (isnan(s_df) || isinf(s_df) || s_df < 0.0) {
        fprintf(stderr, "HIBA [calculate_max_size_from_drift_fragmentation]: NaN, Inf vagy negatív érték keletkezett (%.2e). Bemenet: sigma=%.2e, r=%.2e.\n", s_df, sigma, r);
        return 0.0;
    }

    return s_df;
}


double calculate_growth_timescale(double r, double eps, const disk_t *disk_params) {
    // Ha nincs por (vagy az eps rendkívül kicsi), nincs növekedés.
    // A növekedési időt végtelennek tekintjük.
    if (eps < 1e-15) { // Ahol 1e-15 egy kis, pozitív szám
        return DBL_MAX; 
    }

    double omega = calculate_keplerian_angular_velocity(r, disk_params);

    if (omega <= EPSILON_SIZE) {
        fprintf(stderr, "HIBA [calculate_growth_timescale]: Az omega érték túl kicsi vagy nulla (%.2e). Visszatérés DBL_MAX.\n", omega);
        return DBL_MAX; // Javítás: Végtelen időt ad vissza a nulla helyett, hogy elkerüljük a szokatlan viselkedést.
    }

    // A növekedési idő képlete
    double tau_growth = eps / omega;
    if (isnan(tau_growth) || isinf(tau_growth) || tau_growth < 0.0) {
        return DBL_MAX;
    }
    return tau_growth;
}


double update_particle_size(double particle_radius, double particle_density, double sigma_gas, double sigma_dust, double particle_distance_au, double gas_pressure, double gas_dpdr, double time_step, const disk_t *disk_params) {
    // JAVÍTÁS: A részecskeméret ne legyen túl kicsi vagy negatív a számítások előtt sem.
    if (particle_radius < EPSILON_SIZE) {
        return EPSILON_SIZE;
    }

    // JAVÍTÁS: Hozzáadott robusztus ellenőrzés a bemeneti paraméterekre
    if (sigma_gas < EPSILON_SIZE || sigma_dust < 0.0 || gas_pressure < EPSILON_SIZE || isnan(sigma_gas) || isnan(sigma_dust) || isnan(gas_pressure)) {
        fprintf(stderr, "HIBA [update_particle_size]: Érvénytelen bemeneti paraméterek. (sigma_gas=%.2e, sigma_dust=%.2e, gas_pressure=%.2e).\n", sigma_gas, sigma_dust, gas_pressure);
        return particle_radius; // Visszatérés az eredeti mérethez, ha a bemenet hibás.
    }

    double size_max_turb = calculate_max_size_from_turbulence(sigma_gas, particle_distance_au, particle_density, disk_params);
    double size_max_drift_frag = calculate_max_size_from_drift_fragmentation(sigma_gas, particle_distance_au, gas_pressure, gas_dpdr, particle_density,disk_params);
    double size_max_drift = calculate_max_size_from_drift(sigma_dust, particle_distance_au, gas_pressure, gas_dpdr, particle_density, disk_params);
    
    // JAVÍTÁS: A méretkorlátok helyes meghatározása
    // Csak a pozitív, nem nulla értékeket vesszük figyelembe
    double size_min = DBL_MAX;
    if (size_max_turb > EPSILON_SIZE) {
        if (size_max_turb < size_min) size_min = size_max_turb;
    }
    if (size_max_drift_frag > EPSILON_SIZE) {
        if (size_max_drift_frag < size_min) size_min = size_max_drift_frag;
    }
    if (size_max_drift > EPSILON_SIZE) {
        if (size_max_drift < size_min) size_min = size_max_drift;
    }

    // Ha az összes korlát túl kicsi (vagy nulla), akkor nincs méretkorlát a növekedésre
    if (size_min == DBL_MAX) {
        size_min = DBL_MAX; // Ezt az állapotot a növekedési logika kezeli
    }
    
    // ELLENŐRZÉS: Robusztusabb ellenőrzés a nullával való osztásra.
    double eps = 0.0;
    if (sigma_gas > EPSILON_SIZE) {
        eps = sigma_dust / sigma_gas;
    } else {
        eps = EPSILON_SIZE;
    }
    
    // JAVÍTÁS: Ellenőrizzük, hogy az eps ne legyen negatív.
    if (eps < 0.0) {
        eps = EPSILON_SIZE;
    }

    double tau_growth = calculate_growth_timescale(particle_distance_au, eps, disk_params);
    double new_size = 0.0;
    
    // A méret konvertálása cm-ből AU-ba
    // FONTOS: Csak a size_min-t konvertáljuk, mivel a particle_radius már AU-ban van.
    if (size_min != DBL_MAX) {
        size_min = size_min / AU_TO_CM;
    }
    
    // A növekedési mechanizmus beépítése
    if (tau_growth <= EPSILON_SIZE || tau_growth >= HUGE_VAL) {
        new_size = size_min;
    } else if (particle_radius < size_min) {
        // ELLENŐRZÉS: A logaritmus argumentumának pozitívnak kell lennie.
        double exponent = time_step / tau_growth;
        if (exponent > 100) {
            new_size = size_min;
        } else {
            // A növekedés mértéke
            new_size = particle_radius * exp(exponent);
            // JAVÍTÁS: a new_size nem haladhatja meg a size_min-t
            if (new_size > size_min) {
                new_size = size_min;
            }
        }
    } else {
        new_size = size_min;
    }
    
    // ELLENŐRZÉS: A részecskeméret ne legyen túl kicsi vagy negatív.
    if (new_size < EPSILON_SIZE) {
        new_size = EPSILON_SIZE;
    }
    
    // VISSZATÉRÉSI ÉRTÉK ELLENŐRZÉSE
    if (isnan(new_size) || isinf(new_size) || new_size < 0.0) {
        fprintf(stderr, "HIBA [update_particle_size]: NaN, Inf vagy negatív érték a new_size-ban (%.2e). Visszatérés 0.0. Paraméterek: particle_radius=%.2e, tau_growth=%.2e.\n", new_size, particle_radius, tau_growth);
        return 0.0;
    }

    return new_size;
}

void calculate_dust_density_grid(const ParticleData_t *p_data, disk_t *disk_params, const simulation_options_t *sim_opts) {
    double *sigma_d_temp = (double *)calloc(disk_params->NGRID, sizeof(double));
    if (sigma_d_temp == NULL) {
        fprintf(stderr, "HIBA: Memóriafoglalás sikertelen a sigma_d_temp számára a calculate_dust_density_grid-ben.\n");
        exit(EXIT_FAILURE);
    }
    
    double *sigma_dm_temp = NULL;
    if (sim_opts->twopop == 1.0) {
        sigma_dm_temp = (double *)calloc(disk_params->NGRID, sizeof(double));
        if (sigma_dm_temp == NULL) {
            fprintf(stderr, "HIBA: Memóriafoglalás sikertelen a sigma_dm_temp számára a calculate_dust_density_grid-ben.\n");
            free(sigma_d_temp); // Fontos: Szabadítsuk fel a korábban lefoglalt memóriát
            exit(EXIT_FAILURE);
        }
    }

    if (p_data->num_particles_pop1 > 0) {
        double *rad_pop1_1d = (double *)calloc(p_data->num_particles_pop1, sizeof(double));
        double *mass_pop1_1d = (double *)calloc(p_data->num_particles_pop1, sizeof(double));

        if (rad_pop1_1d == NULL || mass_pop1_1d == NULL) {
            fprintf(stderr, "HIBA: Memóriafoglalás sikertelen a Pop1 ideiglenes tömbökhöz a calculate_dust_density_grid-ben.\n");
            free(sigma_d_temp);
            if (sigma_dm_temp != NULL) free(sigma_dm_temp);
            free(rad_pop1_1d); // A calloc NULL-t ad, ha hibás, így itt a free biztonságos.
            free(mass_pop1_1d);
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < p_data->num_particles_pop1; i++) {
            rad_pop1_1d[i] = p_data->particles_pop1[i].distance_au;
            mass_pop1_1d[i] = p_data->particles_pop1[i].initial_mass_msun;
        }

        calculate_dust_surface_density_profile(sigma_d_temp, disk_params->rvec,
                                                rad_pop1_1d, mass_pop1_1d,
                                                p_data->num_particles_pop1, disk_params->NGRID, disk_params);

        free(rad_pop1_1d);
        free(mass_pop1_1d);
    }

    if (sim_opts->twopop == 1.0 && p_data->num_particles_pop2 > 0) {
        double *rad_pop2_1d = (double *)calloc(p_data->num_particles_pop2, sizeof(double));
        double *mass_pop2_1d = (double *)calloc(p_data->num_particles_pop2, sizeof(double));

        if (rad_pop2_1d == NULL || mass_pop2_1d == NULL) {
            fprintf(stderr, "HIBA: Memóriafoglalás sikertelen a Pop2 ideiglenes tömbökhöz a calculate_dust_density_grid-ben.\n");
            free(sigma_d_temp);
            if (sigma_dm_temp != NULL) free(sigma_dm_temp);
            free(rad_pop2_1d);
            free(mass_pop2_1d);
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < p_data->num_particles_pop2; i++) {
            rad_pop2_1d[i] = p_data->particles_pop2[i].distance_au;
            mass_pop2_1d[i] = p_data->particles_pop2[i].initial_mass_msun;
        }

        calculate_dust_surface_density_profile(sigma_dm_temp, disk_params->rvec,
                                                rad_pop2_1d, mass_pop2_1d,
                                                p_data->num_particles_pop2, disk_params->NGRID, disk_params);

        free(rad_pop2_1d);
        free(mass_pop2_1d);
    }

    for (int i = 0; i < disk_params->NGRID; i++) {
        disk_params->sigmadustvec[i] = sigma_d_temp[i];
        if(sim_opts->twopop == 1.0) {
            disk_params->sigmadustmicrvec[i] = sigma_dm_temp[i];
        }
    }

    free(sigma_d_temp);
    if (sigma_dm_temp != NULL) free(sigma_dm_temp); // Fontos: Csak akkor szabadítsuk fel, ha lefoglaltuk.
}


void update_particle_positions(dust_particle_t *particles_array,
                               int num_particles,
                               double deltat,
                               double t,
                               const simulation_options_t *sim_opts,
                               const disk_t *disk_params) {

    if (disk_params == NULL) {
        fprintf(stderr, "HIBA [update_particle_positions]: a disk_params pointer NULL a belépésnél!\n");
        exit(EXIT_FAILURE);
    }

    if (particles_array == NULL || num_particles <= 0) {
        fprintf(stderr, "HIBA [update_particle_positions]: A függvény korán kilép érvénytelen részecskeadatok miatt.\n");
        return;
    }

    const double *aggregated_sigmad_ptr = disk_params->sigmadustvec;
    const double *aggregated_rdvec_ptr = disk_params->rvec;
    int num_grid_points = disk_params->NGRID;
    double grid_dd = disk_params->DD;
    double grid_rmin = disk_params->RMIN;
    double grid_rmax = disk_params->RMAX;

    #pragma omp parallel for
    for (int i = 0; i < num_particles; ++i) {
        // ELLENŐRZÉS: Ha a részecske már korábban nullázva lett, lépjünk át rajta.
        if (particles_array[i].current_size_au <= 0.0) {
            continue;
        }
        
        double current_r = particles_array[i].distance_au;
        double current_particle_size = particles_array[i].current_size_au;
        
        // JAVÍTÁS: Adjon hozzá egy kis értéket (EPSILON_SIZE) a részecske méretéhez,
        // hogy megakadályozza a nullával való osztást a belső függvényekben.
        if (current_particle_size < EPSILON_SIZE) {
            current_particle_size = EPSILON_SIZE;
        }
        
        double new_r = 0.0;
        double new_size = current_particle_size;

        int_step(t,
                 current_particle_size,
                 current_r,
                 deltat,
                 &new_r,
                 &new_size,
                 aggregated_sigmad_ptr,
                 aggregated_rdvec_ptr,
                 num_grid_points,
                 grid_dd,
                 grid_rmin,
                 disk_params,
                 sim_opts);

        particles_array[i].distance_au = new_r;

        if (sim_opts->growth == 1.0) {
            particles_array[i].current_size_au = new_size;
            // ELLENŐRZÉS: Ellenőrizzük a new_size-t, mielőtt reciprokot számolunk.
            if (particles_array[i].current_size_au > EPSILON_SIZE) {
                particles_array[i].size_reciprocal = 1.0 / particles_array[i].current_size_au;
            } else {
                particles_array[i].size_reciprocal = 0.0;
            }
        }

        // JAVÍTÁS: A részecskék eltávolítása most már akkor történik, ha elhagyják a lemezt.
        if (particles_array[i].distance_au <= grid_rmin || particles_array[i].distance_au >= grid_rmax) {
            fprintf(stderr, "HIBA [update_particle_positions]: A %d. részecske (R=%.4e AU) kívül esett a korongon (RMIN=%.4e, RMAX=%.4e). Eltávolítás.\n",
                            i, particles_array[i].distance_au, grid_rmin, grid_rmax);
            particles_array[i].distance_au = 0.0;
            particles_array[i].current_size_au = 0.0;
            particles_array[i].size_reciprocal = 0.0;
            particles_array[i].drdt = 0.0;
        }
    }
}
