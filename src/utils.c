#include "utils.h"
#include "config.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "particle_data.h"
#include "simulation_types.h" 
#include "dust_physics.h"
#include "gas_physics.h"
#include "boundary_conditions.h"


void linearInterpolation(double *input_value, double *radial_grid, double actual_position, double *interpolated_value, double grid_spacing, const DiskParameters *disk_params) {

	double relative_grid_position, left_cell_radius, linear_slope, interpolated_result;
	int left_cell_index; 

    relative_grid_position = actual_position - disk_params->r_min;
	relative_grid_position = relative_grid_position / grid_spacing;     					
	left_cell_index = (int) floor(relative_grid_position);				
	left_cell_radius = radial_grid[left_cell_index];       		
 	linear_slope = (input_value[left_cell_index + 1] - input_value[left_cell_index]) / grid_spacing; 
	interpolated_result = input_value[left_cell_index] + linear_slope * (actual_position - left_cell_radius);   

	*interpolated_value = interpolated_result;
}

double findMinimumOfAnArray(double value_1, double value_2, double value_3) {

    double minimum = (value_1 < value_2) ? value_1 : value_2; return (minimum < value_3) ? minimum : value_3;

}

double ftcsSecondDerivativeCoefficient(double radial_distance, const DiskParameters *disk_params){					

    double second_derivate_coefficient;
    second_derivate_coefficient = 3.0 * calculateKinematicViscosity(radial_distance, disk_params);
    return second_derivate_coefficient;
}

double ftcsFirstDerivativeCoefficient(double radial_distance, const DiskParameters *disk_params){							

    double first_derivate_coefficient;
    first_derivate_coefficient = 9.0 * calculateKinematicViscosity(radial_distance,disk_params) / (2.0 * radial_distance);
    return first_derivate_coefficient;
}


int identifyPressureTraps(const DiskParameters *disk_params, PressureTrap *traps, int max_traps) {
    PressureTrap temp_list[10];
    int temp_count = 0;
    const double *gradient = disk_params->gas_pressure_gradient_vector;
    const double *grid = disk_params->radial_grid;

    // 1. Összes zérushely megkeresése
    for (int i = 0; i < disk_params->grid_number - 1; i++) {
        if ((gradient[i] * gradient[i+1] <= 0.0) && (gradient[i] > gradient[i+1])) {
            if (temp_count < 10) {
                double slope = (gradient[i+1] - gradient[i]) / (grid[i+1] - grid[i]);
                temp_list[temp_count].radial_position = grid[i] - (gradient[i] / slope);
                temp_list[temp_count].trap_id = -1; 
                temp_count++;
            }
        }
    }

    // Alaphelyzetbe állítás (MINDENKIT nullázunk)
    for(int i = 0; i < max_traps; i++) {
        traps[i].radial_position = 0.0;
        traps[i].trap_id = i;
    }

    if (temp_count == 0) return 0;

    double dze_inner_target = disk_params->r_dze_i;
    double dze_outer_target = disk_params->r_dze_o;
    double tolerance = 5.0; // AU

    // 2. TRAP 0 (Inner DZE) keresése
    int best_inner_idx = -1;
    double min_dist_inner = tolerance;
    for (int i = 0; i < temp_count; i++) {
        double dist = fabs(temp_list[i].radial_position - dze_inner_target);
        if (dist < min_dist_inner) {
            min_dist_inner = dist;
            best_inner_idx = i;
        }
    }
    if (best_inner_idx != -1) {
        traps[0].radial_position = temp_list[best_inner_idx].radial_position;
        temp_list[best_inner_idx].trap_id = 0;
    }

    // 3. TRAP 1 (Outer DZE) keresése
    int best_outer_idx = -1;
    double min_dist_outer = tolerance;
    for (int i = 0; i < temp_count; i++) {
        if (temp_list[i].trap_id != -1) continue; 
        double dist = fabs(temp_list[i].radial_position - dze_outer_target);
        if (dist < min_dist_outer) {
            min_dist_outer = dist;
            best_outer_idx = i;
        }
    }
    if (best_outer_idx != -1) {
        traps[1].radial_position = temp_list[best_outer_idx].radial_position;
        temp_list[best_outer_idx].trap_id = 1;
    }

    // 4. Maradék feltöltése Traps 2+ helyre
    int extra_slot = 2;
    int max_reached_idx = 1; // Legalább a Trap 1-ig nézzük a logban

    for (int i = 0; i < temp_count; i++) {
        if (temp_list[i].trap_id == -1) {
            if (extra_slot < max_traps) {
                traps[extra_slot].radial_position = temp_list[i].radial_position;
                max_reached_idx = extra_slot;
                extra_slot++;
            }
        }
    }

    fprintf(stderr,"DEBUG: Total traps identified = %d\n", max_reached_idx + 1);


    // A legnagyobb használt index + 1-et adjuk vissza, hogy a ciklus mindent lásson
    return (max_reached_idx + 1);
}



/**
 * @brief Kiszámítja a pormasszát egy adott csapdában a strukturált (2D) adatok alapján.
 */
void calculateMassInSpecificTrapStructured(PressureTrap *trap, const StructuredParticleData *data, const SimulationOptions *sim_opts) {
    double primary_mass = 0.0;
    double secondary_mass = 0.0;

    // Végigiterálunk a radiális oszlopokon
    for (size_t i = 0; i < data->n_r; i++) {
        // Mivel minden oszlopban (j) azonos a sugárirányú pozíció, 
        // elég a legalsó (j=0) részecske r_au értékét megnézni.
        double r = data->particles[i][0].r_au;

        // Ellenőrizzük, hogy ez a radiális sáv a csapda határain belül van-e
        if (r >= trap->inner_boundary && r <= trap->outer_boundary) {
            
            // Összeadjuk a teljes vertikális oszlop tömegét (ez a "NewSum")
            double column_mass = 0.0;
            for (size_t j = 0; j < data->n_z; j++) {
                column_mass += data->particles[i][j].mass_g;
            }

            // Pop1 (elsődleges, növekvő por)
            primary_mass += column_mass;

            // Pop2 (szekunder, mikronos por)
            // Megjegyzés: Ha a mikronos port is a Structured-ben tárolod egy másik mezőben, 
            // ide jön az összeadása. Ha külön Structured tömbje van, akkor azt is be kell vonni.
            if (sim_opts->option_for_dust_secondary_population == 1.0) {
                // Példa, ha a DustParticle struktúrában van külön mass_micron_g:
                // secondary_mass += data->particles[i][j].mass_micron_g;
                // Jelenleg feltételezzük, hogy a Pop2-t még fejlesztened kell az új struktúrához.
            }
        }
    }

    trap->primary_dust_mass = primary_mass;
    trap->secondary_dust_mass = secondary_mass;
    trap->total_dust_mass = primary_mass + secondary_mass;



  // DEBUG KIÍRÁS
    fprintf(stderr,"DEBUG: Trap [%.3f - %.3f] AU: primary_mass = %.6e, secondary_mass = %.6e, total_mass = %.6e\n",
           trap->inner_boundary, trap->outer_boundary,
           trap->primary_dust_mass, trap->secondary_dust_mass, trap->total_dust_mass);


}


/**
 * @brief Frissíti a részecskék rácsindexeit.
 * Nincs szükség current_time-ra, mert a mass_g már tartalmazza a fix tömeget.
 */
void updateStructuredParticleGridIndices(StructuredParticleData *sd,  const DiskParameters *disk_params) {
    if (!sd || !sd->particles) return;

    for (size_t i = 0; i < sd->n_r; i++) {
        for (size_t j = 0; j < sd->n_z; j++) {
            DustParticle *p = &sd->particles[i][j];

            // Rácsindex számítása a sugár alapján
            double relative_pos = (p->r_au - disk_params->r_min) / disk_params->delta_r;
            int g_idx = (int)floor(relative_pos + 0.5);

            // Biztonsági korlátok
            if (relative_pos < 0 || isnan(relative_pos)) {
                g_idx = 0;
            } else if (g_idx >= (int)disk_params->grid_number) {
                g_idx = (int)disk_params->grid_number - 1;
            }
            
            p->grid_index = g_idx;
        }
    }
}


void updateDustSurfaceDensityEulerian(StructuredParticleData *data, const DiskParameters *disk_params) {
    size_t n_r = disk_params->grid_number; // Fix gáz-rács mérete
    size_t n_z = data->n_z;

    // 1. Nullázás
    for (size_t i = 0; i < n_r; i++) disk_params->dust_surface_density_euler[i] = 0.0;

    // 2. Binning: A részecskék tömegét a fix gáz-rácshoz adjuk
    for (size_t i = 0; i < data->n_r; i++) { // Végig az összes Lagrange-i oszlopon
        for (size_t j = 0; j < n_z; j++) {
            DustParticle *p = &data->particles[i][j];
            
            if (p->r_au > disk_params->r_min && p->r_au < disk_params->r_max) {
                // Megkeressük a fix Euler-cella indexét
                int idx = (int)((p->r_au - disk_params->r_min) / disk_params->delta_r);
                
                if (idx >= 0 && idx < (int)n_r) {
                    disk_params->dust_surface_density_euler[idx] += p->mass_g;
                }
            }
        }
    }

    // 3. Normalizálás a fix cellák területével
    for (size_t i = 0; i < n_r; i++) {
        double r_cell = disk_params->radial_grid[i];
        double area_cgs = 2.0 * M_PI * r_cell * disk_params->delta_r * (AU_IN_CM * AU_IN_CM);
        if (area_cgs > 0) disk_params->dust_surface_density_euler[i] /= area_cgs;
    }
}



void updateDustSurfaceDensityEulerianCIC(StructuredParticleData *data, const DiskParameters *disk_params) {
    size_t n_r_euler = disk_params->grid_number;
    size_t n_r_lag = data->n_r;
    double dr_e = disk_params->delta_r;
    double r_min = disk_params->r_min;

    // 1. Nullázás
    for (size_t i = 0; i < n_r_euler; i++) disk_params->dust_surface_density_euler[i] = 0.0;

    // 2. Tömeggyűjtés (CSAK a tömeget adjuk össze, nem osztunk semmivel a ciklusban!)
    for (size_t i = 0; i < n_r_lag; i++) {
        for (size_t j = 0; j < data->n_z; j++) {
            DustParticle *p = &data->particles[i][j];
            
            // Hol van a részecske az Euler-rácshoz képest?
            double float_idx = (p->r_au - r_min) / dr_e-0.5;
            int i_low = (int)floor(float_idx);
            int i_high = i_low + 1;
            
            double weight_high = float_idx - (double)i_low;
            double weight_low = 1.0 - weight_high;

            // Alacsonyabb indexű cella
            if (i_low >= 0 && i_low < (int)n_r_euler) {
                disk_params->dust_surface_density_euler[i_low] += (double)(p->mass_g * weight_low);
            }
            // Magasabb indexű cella
            if (i_high >= 0 && i_high < (int)n_r_euler) {
                disk_params->dust_surface_density_euler[i_high] += (double)(p->mass_g * weight_high);
            }
        }
    }

    // 3. Normalizálás: Sigma = Össztömeg_a_cellában / Cella_Területe
    for (size_t i = 0; i < n_r_euler; i++) {
        double r_center = r_min + (i+0.5) * dr_e;
        
        // A cella határai (fél rácsközre a középponttól)
        double r_in = r_center - 0.5 * dr_e;
        double r_out = r_center + 0.5 * dr_e;
        
        // Határvédelem
        if (r_in < r_min) r_in = r_min;
        
        // Terület cm2-ben
        double area = M_PI * (r_out * r_out - r_in * r_in);
        
        if (area > 0) {
            disk_params->dust_surface_density_euler[i] /= area;
        }
    }
}



void updateDustSurfaceDensityEulerianTSC(StructuredParticleData *data, const DiskParameters *disk_params) {
    size_t n_r_euler = disk_params->grid_number;
    size_t n_r_lag = data->n_r;
    double dr_e = disk_params->delta_r;
    double r_min = disk_params->r_min;

    // 1. Nullázás
    for (size_t i = 0; i < n_r_euler; i++) disk_params->dust_surface_density_euler[i] = 0.0;

    // 2. Tömeggyűjtés TSC-vel
    for (size_t i = 0; i < n_r_lag; i++) {
        for (size_t j = 0; j < data->n_z; j++) {
            DustParticle *p = &data->particles[i][j];
            
            // Relatív pozíció a rácson (távolság r_min-től egységnyi dr-ben mérve)
            double u = (p->r_au - r_min) / dr_e;
            int i_nearest = (int)round(u); // A legközelebbi rácspont indexe

            // Távolság a legközelebbi rácsponttól (-0.5 és 0.5 között)
            double d = u - (double)i_nearest;

            // TSC Súlyfüggvények (3 pontra)
            double w_mid  = 0.75 - d * d;
            double w_left  = 0.5 * (0.5 - d) * (0.5 - d);
            double w_right = 0.5 * (0.5 + d) * (0.5 + d);

            int indices[3] = {i_nearest - 1, i_nearest, i_nearest + 1};
            double weights[3] = {w_left, w_mid, w_right};

            for (int k = 0; k < 3; k++) {
                int idx = indices[k];
                if (idx >= 0 && idx < (int)n_r_euler) {
                    disk_params->dust_surface_density_euler[idx] += p->mass_g * weights[k];
                }
            }
        }
    }

    // 3. Normalizálás a pontos körgyűrű területtel
    for (size_t i = 0; i < n_r_euler; i++) {
        double r_in  = r_min + (double)i * dr_e;
        double r_out = r_in + dr_e;
        
        double area = M_PI * (r_out * r_out - r_in * r_in);
        
        if (area > 0) {
            disk_params->dust_surface_density_euler[i] /= area;
        }
    }
}


void updateDustSurfaceDensitySmart(StructuredParticleData *data, const DiskParameters *disk_params) {
    size_t n_r_euler = disk_params->grid_number;
    size_t n_r_lag = data->n_r;
    double dr_e = disk_params->delta_r;
    double r_min = disk_params->r_min;

    // 1. Tömeggyűjtés TSC-vel (ez védi legjobban a csúcsokat)
    for (size_t i = 0; i < n_r_euler; i++) disk_params->dust_surface_density_euler[i] = 0.0;

    for (size_t i = 0; i < n_r_lag; i++) {
        for (size_t j = 0; j < data->n_z; j++) {
            DustParticle *p = &data->particles[i][j];
            double u = (p->r_au - r_min) / dr_e;
            int i_nearest = (int)round(u);
            double d = u - (double)i_nearest;

            double weights[3] = {0.5*(0.5-d)*(0.5-d), 0.75-d*d, 0.5*(0.5+d)*(0.5+d)};
            for (int k = -1; k <= 1; k++) {
                int idx = i_nearest + k;
                if (idx >= 0 && idx < (int)n_r_euler) {
                    disk_params->dust_surface_density_euler[idx] += p->mass_g * weights[k+1];
                }
            }
        }
    }

    // 2. Normalizálás
    for (size_t i = 0; i < n_r_euler; i++) {
        double r_in = r_min + i * dr_e;
        double r_out = r_in + dr_e;
        double area = M_PI * (r_out * r_out - r_in * r_in);
        disk_params->dust_surface_density_euler[i] /= area;
    }

    // 3. "Lyukkitöltés" (Gap filling) - Csak ott avatkozik be, ahol üres cellák vannak
    // Egy másolatot készítünk, hogy ne gerjesszünk instabilitást
    double *temp = malloc(n_r_euler * sizeof(double));
    memcpy(temp, disk_params->dust_surface_density_euler, n_r_euler * sizeof(double));

    for (int i = 1; i < (int)n_r_euler - 1; i++) {
        // Ha a cella értéke gyanúsan alacsony a szomszédokhoz képest (zaj)
        if (temp[i] < 1e-15) { // Vagy egy ésszerű küszöbérték
            // Lineáris interpoláció a szomszédokból, hogy ne legyen nulla
            disk_params->dust_surface_density_euler[i] = 0.5 * (temp[i-1] + temp[i+1]);
        }
    }
    free(temp);
}