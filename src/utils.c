#include "utils.h"
#include "config.h"
#include <math.h>
#include <stdlib.h>
#include "particle_data.h"
#include "simulation_types.h" 
#include "dust_physics.h"
#include "gas_physics.h"
#include "boundary_conditions.h"


void linearInterpolation(double *invec, double *radial_grid, double pos, double *out, double rd, int opt, const DiskParameters *disk_params) {

	double rmid, rindex, coef1, temp;
	int index; 

    rmid = pos - disk_params->r_min;
	rmid = rmid / rd;     					
	index = (int) floor(rmid);				
	rindex = radial_grid[index];       		
 	coef1 = (invec[index + 1] - invec[index]) / rd; 
	temp = invec[index] + coef1 * (pos - rindex);   
	if(opt == 1) if(temp < 0) temp = -1.*temp;

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

int countZeroPoints(const DiskParameters *disk_params) {

	int i,count;
	count = 0;

	for(i = 0; i < disk_params->grid_number-1; i++) {
		if(((disk_params->gas_pressure_gradient_vector[i] * disk_params->gas_pressure_gradient_vector[i+1]) <= 0.)  && (disk_params->gas_pressure_gradient_vector[i] > disk_params->gas_pressure_gradient_vector[i+1])) {
			count++;
		} 
	}

	return count;
}

double findZeroPointRadius(double r1, double r2, double dp1, double dp2) {

	double a, b, r_zero;
	a = (dp2 - dp1) / (r2 - r1);
	b = dp1 - a * r1;
	r_zero = - b / a;

	return r_zero;

}

double findZeroPoint(int i, const double *radial_grid, const double *dp) {

	double r;
	
	if(((dp[i] * dp[i+1]) <= 0.) && (dp[i] > dp[i+1])) {	
		r = findZeroPointRadius(radial_grid[i],radial_grid[i+1],dp[i],dp[i+1]);	
	} else {
		r = 0.0;
	}

	return r;

}

void findRAnnulusAroundDZE(double rin, double *ind_ii, double *ind_io,
                            double rout, double *ind_oi, double *ind_oo,
                            const SimulationOptions *sim_opts, const DiskParameters *disk_params) {

    if (disk_params == NULL) {
        fprintf(stderr, "ERROR [findRAnnulusAroundDZE]: disk_params is NULL!\n");
        exit(1);
    }

    int i;
    double rmid, rtemp;
    double roimH, roipH, roomH, roopH;
    double riimH, riipH, riomH, riopH;

    *ind_ii = 0.0;
    *ind_io = 0.0;
    *ind_oi = 0.0;
    *ind_oo = 0.0;

    double h_rin = calculatePressureScaleHeight(rin, disk_params); 
    double h_rout = calculatePressureScaleHeight(rout, disk_params); 
    double rin_minus_h_rin = rin - h_rin;
    double rin_plus_h_rin = rin + h_rin;
    double rout_minus_h_rout = rout - h_rout;
    double rout_plus_h_rout = rout + h_rout;

    riimH = rin_minus_h_rin - disk_params->delta_r / 2.0;
    riipH = rin_minus_h_rin + disk_params->delta_r / 2.0;
    riomH = rin_plus_h_rin - disk_params->delta_r / 2.0;
    riopH = rin_plus_h_rin + disk_params->delta_r / 2.0;
    roimH = rout_minus_h_rout - disk_params->delta_r / 2.0;
    roipH = rout_minus_h_rout + disk_params->delta_r / 2.0;
    roomH = rout_plus_h_rout - disk_params->delta_r / 2.0;
    roopH = rout_plus_h_rout + disk_params->delta_r / 2.0;


    for (i = 0; i < disk_params->grid_number; i++) {

        if (sim_opts->flag_for_deadzone == 1) {
            if (disk_params->radial_grid[i] > riimH && disk_params->radial_grid[i] < riipH) {
                rmid = (disk_params->radial_grid[i] - disk_params->r_min) / disk_params->delta_r;
                rtemp = floor(rmid + 0.5);
                *ind_ii = rtemp;
            }

            if (disk_params->radial_grid[i] > riomH && disk_params->radial_grid[i] < riopH) {
                rmid = (disk_params->radial_grid[i] - disk_params->r_min) / disk_params->delta_r;
                rtemp = floor(rmid + 0.5);
                *ind_io = rtemp;
            }
        }

        if (disk_params->radial_grid[i] > roimH && disk_params->radial_grid[i] < roipH) {
            rmid = (disk_params->radial_grid[i] - disk_params->r_min) / disk_params->delta_r;
            rtemp = floor(rmid + 0.5);
            *ind_oi = rtemp;
        }

        if (disk_params->radial_grid[i] > roomH && disk_params->radial_grid[i] < roopH) {
            rmid = (disk_params->radial_grid[i] - disk_params->r_min) / disk_params->delta_r;
            rtemp = floor(rmid + 0.5);
            *ind_oo = rtemp;
        }

        if (disk_params->radial_grid[i] > roopH) break;
    }
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
