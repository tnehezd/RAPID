#include "utils.h"
#include "config.h"
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

    // A legnagyobb használt index + 1-et adjuk vissza, hogy a ciklus mindent lásson
    return (max_reached_idx + 1);
}


/**
 * @brief Calculates dust mass for a specific trap based on Lagrangian particle positions.
 */
void calculateMassInSpecificTrap(PressureTrap *trap, const ParticleData *particle_data, int particle_number, const SimulationOptions *sim_opts) {

    double primary_mass = 0.0;
    double secondary_mass = 0.0;

    for (int i = 0; i < particle_number; i++) {
        // Elsődleges (cm) por
        double radial_distance = particle_data->particle_distance_array[i][0];
        if (radial_distance >= trap->inner_boundary && radial_distance <= trap->outer_boundary) {
            primary_mass += particle_data->dust_particle_mass_array[i][0];
        }

        // Másodlagos (mikron) por
        if (sim_opts->option_for_dust_secondary_population == 1.0) {
            double radial_distance_micron = particle_data->micron_particle_distance_array[i][0];
            if (radial_distance_micron >= trap->inner_boundary && radial_distance_micron <= trap->outer_boundary) {
                secondary_mass += particle_data->micron_dust_particle_mass_array[i][0];
            }
        }
    }

    trap->primary_dust_mass = primary_mass;
    trap->secondary_dust_mass = secondary_mass;
    trap->total_dust_mass = primary_mass + secondary_mass;
}



void updateParticleGridIndices(const ParticleData *particle_data, double current_time,int particle_count,const DiskParameters *disk_params) {

    int particle_index;
    int grid_index;
    double relative_grid_position;

    for (particle_index = 0; particle_index < particle_count; particle_index++) {

        relative_grid_position =
            (particle_data->particle_distance_array[particle_index][0] -
             disk_params->r_min) / disk_params->delta_r;

        grid_index = (int)floor(relative_grid_position + 0.5);

        if (relative_grid_position < 0 || isnan(relative_grid_position))
            grid_index = 0;

        if (particle_count == particle_number)
            particle_data->dust_particle_mass_array[particle_index][0] =
                particle_data->dust_particle_mass_grid[particle_index];

        particle_data->dust_particle_mass_array[particle_index][1] = grid_index;

        if (current_time == 0) {
            particle_data->dust_particle_mass_array[particle_index][2] =
                particle_data->dust_particle_mass_array[particle_index][0];
            particle_data->dust_particle_mass_array[particle_index][3] = 0;
        }
    }
}

void computeParticleRadiusRange(const ParticleData *particle_data,int particle_number,int has_secondary_population,double *min_radius,double *max_radius) {

    double min_primary = HUGE_VAL;
    double max_primary = -HUGE_VAL;

    for (int i = 0; i < particle_number; i++) {
        double particle_radius = particle_data->particle_distance_array[i][0];
        if (particle_radius > 0.0) {
            if (particle_radius < min_primary) min_primary = particle_radius;
            if (particle_radius > max_primary) max_primary = particle_radius;
        }
    }

    double min_secondary = HUGE_VAL;
    double max_secondary = -HUGE_VAL;

    if (has_secondary_population) {
        for (int i = 0; i < particle_number; i++) {
            double particle_radius = particle_data->micron_particle_distance_array[i][0];
            if (particle_radius > 0.0) {
                if (particle_radius < min_secondary) min_secondary = particle_radius;
                if (particle_radius > max_secondary) max_secondary = particle_radius;
            }
        }
    }

    *min_radius = has_secondary_population ? fmin(min_primary, min_secondary) : min_primary;
    *max_radius = has_secondary_population ? fmax(max_primary, max_secondary) : max_primary;
}
