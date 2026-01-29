// src/disk_model.c

#include "disk_model.h"   
#include "config.h"       
#include "simulation_types.h"
#include "gas_physics.h"
#include "dust_physics.h" 
#include "io_utils.h"     
#include "utils.h" 
#include "boundary_conditions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>



/*	Reads the parameters of the disk	*/
void readDiskParameters(DiskParameters *disk_params) {
    // Check if the pointer is NULL
    if (disk_params == NULL) {
        fprintf(stderr, "ERROR [readDiskParameters]: Received NULL disk_params pointer.\n");
        exit(EXIT_FAILURE);
    }

    fprintf(stderr, "DEBUG [readDiskParameters]: Calculating derived disk parameters and writing to output file.\n");

    disk_params->particle_density_dimensionless = disk_params->particle_density / SOLAR_MASS_IN_GRAMS * AU_IN_CM * AU_IN_CM * AU_IN_CM;

    fprintf(stderr, "DEBUG [readDiskParameters]: Calculated particle_density = %.2e, particle_density_dimensionless = %.2e.\n",
           disk_params->particle_density, disk_params->particle_density_dimensionless);
}




/*	Initialize radial grid	*/
void createRadialGrid(DiskParameters *disk_params) {
	
	int i;
 	for(i = 0; i <= disk_params->grid_number+1; i++) {						/*	load an array of radii	*/
 		disk_params->radial_grid[i] = disk_params->r_min + (i-1) * disk_params->delta_r;
//        fprintf(stderr, "DEBUG [createRadialGrid]: r: %lg\n", disk_params->radial_grid[i]);
	}
}

/*	Create the initial gas surface density profile	*/
void createInitialGasSurfaceDensity(DiskParameters *disk_params){		/*	initial profile of sigma		*/

  	int i;
  
  	for(i = 1; i <= disk_params->grid_number; i++) {
    		disk_params->sigmavec[i] = disk_params->sigma_0 * pow(disk_params->radial_grid[i],disk_params->sigma_power_law_index);		/*	sigma0*r^x (x could be eg. -1/2)	*/
    }
  

  	applyBoundaryConditions(disk_params->sigmavec,disk_params);

}

void createInitialGasPressure(DiskParameters *disk_params){	

  	int i;
  
  	for(i = 1; i <= disk_params->grid_number; i++) {
    		disk_params->pressvec[i] = calculateGasPressure(disk_params->sigmavec[i],disk_params->radial_grid[i],disk_params);
  	}
  	applyBoundaryConditions(disk_params->pressvec,disk_params);


}

void createInitialGasPressureGradient(DiskParameters *disk_params){

	calculateGasPressureGradient(disk_params);
   	applyBoundaryConditions(disk_params->dpressvec,disk_params);


}

/*	Update radial gas velovity	*/
void createInitialGasVelocity(DiskParameters *disk_params){	
 	
	calculateGasRadialVelocity(disk_params);
  	applyBoundaryConditions(disk_params->ugvec,disk_params);
}



void calculateInitialDustSurfaceDensity(double radin[][2], double *massin, double out[][3], int n, const DiskParameters *disk_params) {

	int i;

	for(i=0;i<n;i++){

/*	Calcualte the surface density of the dust grains	*/
/*	If the dust grain is within the simulated regine (above r_min)
 	the surface density is calculated from the representative mass of the dust grain	*/
		if((radin[i][0] >= disk_params->r_min)) {
			out[i][0] = massin[i] / (2. * (radin[i][0]-disk_params->delta_r/2.) * M_PI * disk_params->delta_r);	// sigma = m /(2 * r * pi * dr) --> dr is the original grid step
			out[i][1] = radin[i][0];																	// Saves the radial distance of the dust grain

  			double rmid = (radin[i][0] - disk_params->r_min) / disk_params->delta_r;     						// 	The integer part of this gives at which index is the body		
			int rindex = (int) floor(rmid);																// 	The "whole part of rmin" --> floor rounds down, +0.5 allows us to solve the rounding correctly
			out[i][2] = (double) rindex;

/*	If the dust grain is drifted below  r_min, the surface density is set to 0*/	
		} else {
			out[i][0] = 0;
			out[i][1] = 0;																				// Set r to 0!
			out[i][2] = 0;
		}
	}
}
