// src/disk_model.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "disk_model.h"   
#include "config.h"       
#include "simulation_types.h"
#include "globals.h"
#include "dust_physics.h" 
#include "io_utils.h"     
#include "utils.h" 


/*	A korong parametereinek beolvasasa	*/
void read_disk_parameters(disk_t *disk_params) {
    // Ellenőrzés, ha a pointer NULL (jó gyakorlat)
    if (disk_params == NULL) {
        fprintf(stderr, "ERROR [read_disk_parameters]: Received NULL disk_params pointer.\n");
        exit(EXIT_FAILURE);
    }

    fprintf(stderr, "Calculating derived disk parameters and writing to output file.\n");

    disk_params->PDENSITYDIMLESS = disk_params->PDENSITY / SUN_MASS_TO_GRAMS * AU_TO_CM * AU_TO_CM * AU_TO_CM;

}


/*	r vektor (gridcellák) inicializálása	*/
void initialize_grid_cells(disk_t *disk_params) {
	int i;
 	for(i = 0; i <= disk_params->NGRID+1; i++) {						
 		disk_params->rvec[i] = disk_params->RMIN + (i-1) * disk_params->DD;
	}
}

/*	a sigmara kezdeti profil betoltese	*/
void initial_gas_surface_density_profile(disk_t *disk_params){		/*	initial profile of sigma		*/

  	int i;
  
  	for(i = 1; i <= disk_params->NGRID; i++) {
    		disk_params->sigmavec[i] = disk_params->SIGMA0 * pow(disk_params->rvec[i],disk_params->SIGMAP_EXP);		/*	sigma0*r^x (x could be eg. -1/2)	*/
    }
  

  	Perem(disk_params->sigmavec,disk_params);

}

void initial_gas_pressure_profile(disk_t *disk_params){		/*	initial profile of pressure		*/

  	int i;
  
  	for(i = 1; i <= disk_params->NGRID; i++) {
    		disk_params->pressvec[i] = calculate_gas_pressure(disk_params->sigmavec[i],disk_params->rvec[i],disk_params);
  	}
  	Perem(disk_params->pressvec,disk_params);


}

void initial_gas_pressure_gradient_profile(disk_t *disk_params){		/*	initial profile of pressure		*/

	calculate_gas_pressure_gradient(disk_params);
   	Perem(disk_params->dpressvec,disk_params);


}

/*	ug vektor feltoltese az u_gas ertekevel	*/
void initial_gas_velocity_profile(disk_t *disk_params){		/*	initial profile of pressure		*/
 	
	u_gas(disk_params);
  	Perem(disk_params->ugvec,disk_params);
}



/*	Lokalis viszkozitas erteke	*/
double calculate_gas_viscosity(double r, const disk_t *disk_params) {
    double nu;
    double cs, H;

    H = calculate_scale_height(r,disk_params);
    cs = calculate_local_sound_speed(r,disk_params);

    nu = calculate_turbulent_alpha(r, disk_params) * cs * H;
    return nu;
}

/*	local scale height	*/
double calculate_scale_height(double r, const disk_t *disk_params) {

    if (disk_params == NULL) {
        fprintf(stderr, "ERROR [calculate_scale_height]: disk_params is NULL!\n");
        return 0.0; // Vagy valamilyen hibakód/NaN
    }

    // Itt van az eredeti számítás
    double calculated_result = pow(r, 1. + disk_params->FLIND) * disk_params->HASP;
    return calculated_result;
}


/*	lokális kepleri sebesség	*/
double calculate_keplerian_velocity(double r, const disk_t *disk_params) {
    return sqrt(G_GRAV_CONST * disk_params->STAR_MASS / r);
}

/*	lokalis kepleri korfrekvencia	*/
double calculate_keplerian_angular_velocity(double r, const disk_t *disk_params) {
    return sqrt(G_GRAV_CONST * disk_params->STAR_MASS / r / r / r);
}


/*	local sound speed		*/
double calculate_local_sound_speed(double r, const disk_t *disk_params) {
    return calculate_keplerian_angular_velocity(r,disk_params) * calculate_scale_height(r,disk_params);
}

/*	Suruseg a midplane-ben	*/
double calculate_midplane_gas_density(double sigma, double r, const disk_t *disk_params) {
    return 1. / sqrt(2.0 * M_PI) * sigma / calculate_scale_height(r,disk_params);
}

/* local pressure of the gas p = rho_gas * cs * cs kepletbol!!	*/
double calculate_gas_pressure(double sigma, double r, const disk_t *disk_params) {
    return calculate_midplane_gas_density(sigma, r, disk_params) * calculate_local_sound_speed(r,disk_params) * calculate_local_sound_speed(r, disk_params);
}

/*	a nyomas derivaltja	*/
void calculate_gas_pressure_gradient(disk_t *disk_params) {
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