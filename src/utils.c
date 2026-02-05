#include "utils.h"
#include "config.h"
#include <math.h>
#include <stdlib.h>
#include "particle_data.h"
#include "simulation_types.h" 
#include "dust_physics.h"
#include "gas_physics.h"
#include "boundary_conditions.h"


void linearInterpolation(double *invec, double *radial_grid, double pos, double *out, double rd, const DiskParameters *disk_params) {

	double rmid, rindex, coef1, temp;
	int index; 

    rmid = pos - disk_params->r_min;
	rmid = rmid / rd;     					
	index = (int) floor(rmid);				
	rindex = radial_grid[index];       		
 	coef1 = (invec[index + 1] - invec[index]) / rd; 
	temp = invec[index] + coef1 * (pos - rindex);   

	*out = temp;
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

double findMinimumOfAnArray(double s1, double s2, double s3) {

	double minimum = (s1 < s2) ? s1 : s2; return (minimum < s3) ? minimum : s3;

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

void sortAnArray(double *rv,int n) {

	double temp;
	int i, step;

	for(step = 1; step <= (n-1); step++) {

		for(i = 0; i <= (n-2); i++) {

			if(rv[i] > rv[i + 1]) {
				temp = rv[i];
				rv[i] = rv[i + 1];
				rv[i + 1] = temp;
			}
		}
	}
}

void histogram(double r, int *hist, double dd, DiskParameters *disk_params) {

    int index;
    double rmid; 

    if (r < disk_params->r_min) {
        r = disk_params->r_min;
    } else if (r > disk_params->r_max) {
        r = disk_params->r_max;
    }

    rmid = (r - disk_params->r_min) / dd;
    index = (int) floor(rmid);

    if (index < 0) {
        index = 0; 
    }

    if (index >= disk_params->grid_number) { 
        index = disk_params->grid_number - 1;
    }

    hist[index]++;
}

void sortAnArrayarray(double rv[][3],int n) {

	double temp, temp2, temp3;
	int i, step;

	for(step = 1; step <= (n-1); step++) {
		for(i = 0; i <= (n-2); i++) {
			if(rv[i][1] > rv[i + 1][1]) {
				temp = rv[i][1];
				rv[i][1] = rv[i + 1][1];
				rv[i + 1][1] = temp;
				temp2 = rv[i][0];
				rv[i][0] = rv[i + 1][0];
				rv[i + 1][0] = temp2;
				temp3 = rv[i][2];
				rv[i][2] = rv[i + 1][2];
				rv[i + 1][2] = temp3;
			}
		}
	}
}


void roundParticleRadii(double in[][3], int n, const DiskParameters *disk_params) {

	double dd = (disk_params->r_max - disk_params->r_min) / (particle_number-1);
	int dker = (int)(1./dd);//
	dker = dker * ROUNDING_FACTOR;
	double ddker = (double) dker;
	int i;
	int temp;

	for(i = 0; i<n; i++) {
		temp = (int)floor(in[i][1] * ddker+0.5);
		in[i][1] = (double)temp / ddker;
	}
}




void mergeParticlesByRadius(double in[][3], double dd, int n, const DiskParameters *disk_params) {

	int i;
	int j;
	int k;
	double sig = 0, radout[n], sigout[n];

	for(i = 0; i < n; i++) {
		radout[i] = 0;
		sigout[i] = 0;
		in[i][2] = 0;
	}

	i = 0;
	j = 0;
	k = 0;

	roundParticleRadii(in,n,disk_params);
	sortAnArrayarray(in,n);

	do {
		if(in[i][1] != in[i+1][1]) {
			radout[j] = in[i][1];
			sigout[j] = in[i][0];
			sig = 0;
			k = 0;
			j++;
			i++;
		} else {
			do {
				radout[j] = in[i][1];
				sig = sig + in[i+k][0];
				sigout[j] = sig;
				k++;
			} while (in[i][1] == in[i+k][1]);
			i = i+k;
			k = 0;
			j++;
		}

	} while (i < n);

	for(i = 0; i < n; i++) {
		in[i][0] = sigout[i];
		in[i][1] = radout[i];
  		double rmid = (in[i][1] - disk_params->r_min) / dd;     	
		int rindex = (int) floor(rmid);	
		in[i][2] = (double)rindex;
	}
}

void updateParticleGridIndices(const ParticleData *particle_data, double t, int n, const DiskParameters *disk_params) {

    int i, rindex;
    double rmid;  

    for (i = 0; i < n; i++) {   
        rmid = (particle_data->particle_distance_array[i][0] - disk_params->r_min) / disk_params->delta_r; 
        rindex = (int) floor(rmid+0.5);
        if(rmid < 0) rindex = 0;
        if(isnan(rmid)) rindex = 0;
		if(n == particle_number) particle_data->dust_particle_mass_array[i][0] = particle_data->dust_particle_mass_grid[i];	
		particle_data->dust_particle_mass_array[i][1] = rindex;		
		if(t == 0) {
			particle_data->dust_particle_mass_array[i][2] = particle_data->dust_particle_mass_array[i][0];					
			particle_data->dust_particle_mass_array[i][3] = 0;
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
