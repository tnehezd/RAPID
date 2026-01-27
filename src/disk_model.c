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
void readDiskParameters(disk_t *disk_params) {
    // Check if the pointer is NULL
    if (disk_params == NULL) {
        fprintf(stderr, "ERROR [readDiskParameters]: Received NULL disk_params pointer.\n");
        exit(EXIT_FAILURE);
    }

    fprintf(stderr, "DEBUG [readDiskParameters]: Calculating derived disk parameters and writing to output file.\n");

    disk_params->PDENSITYDIMLESS = disk_params->PDENSITY / SUN2GR * AU2CM * AU2CM * AU2CM;

    fprintf(stderr, "DEBUG [readDiskParameters]: Calculated PDENSITY = %.2e, PDENSITYDIMLESS = %.2e.\n",
           disk_params->PDENSITY, disk_params->PDENSITYDIMLESS);
}




/*	Initialize radial grid	*/
void createRadialGrid(disk_t *disk_params) {
	
	int i;
 	for(i = 0; i <= disk_params->NGRID+1; i++) {						/*	load an array of radii	*/
 		disk_params->rvec[i] = disk_params->RMIN + (i-1) * disk_params->DD;
//        fprintf(stderr, "DEBUG [createRadialGrid]: r: %lg\n", disk_params->rvec[i]);
	}
}

/*	Create the initial gas surface density profile	*/
void createInitialGasSurfaceDensity(disk_t *disk_params){		/*	initial profile of sigma		*/

  	int i;
  
  	for(i = 1; i <= disk_params->NGRID; i++) {
    		disk_params->sigmavec[i] = disk_params->SIGMA0 * pow(disk_params->rvec[i],disk_params->SIGMAP_EXP);		/*	sigma0*r^x (x could be eg. -1/2)	*/
    }
  

  	Perem(disk_params->sigmavec,disk_params);

}

void createInitialGasPressure(disk_t *disk_params){	

  	int i;
  
  	for(i = 1; i <= disk_params->NGRID; i++) {
    		disk_params->pressvec[i] = press(disk_params->sigmavec[i],disk_params->rvec[i],disk_params);
  	}
  	Perem(disk_params->pressvec,disk_params);


}

void createInitialGasPressureGradient(disk_t *disk_params){

	dpress(disk_params);
   	Perem(disk_params->dpressvec,disk_params);


}

/*	Update radial gas velovity	*/
void createInitialGasVelocity(disk_t *disk_params){	
 	
	u_gas(disk_params);
  	Perem(disk_params->ugvec,disk_params);
}



void calculateInitialDustSurfaceDensity(double radin[][2], double *massin, double out[][3], int n, const disk_t *disk_params) {

	int i;

	for(i=0;i<n;i++){

/*	Calcualte the surface density of the dust grains	*/
/*	If the dust grain is within the simulated regine (above RMIN)
 	the surface density is calculated from the representative mass of the dust grain	*/
		if((radin[i][0] >= disk_params->RMIN)) {
			out[i][0] = massin[i] / (2. * (radin[i][0]-disk_params->DD/2.) * M_PI * disk_params->DD);	// sigma = m /(2 * r * pi * dr) --> dr is the original grid step
			out[i][1] = radin[i][0];																	// Saves the radial distance of the dust grain

  			double rmid = (radin[i][0] - disk_params->RMIN) / disk_params->DD;     						// 	The integer part of this gives at which index is the body		
			int rindex = (int) floor(rmid);																// 	The "whole part of rmin" --> floor rounds down, +0.5 allows us to solve the rounding correctly
			out[i][2] = (double) rindex;

/*	If the dust grain is drifted below  RMIN, the surface density is set to 0*/	
		} else {
			out[i][0] = 0;
			out[i][1] = 0;																				// Set r to 0!
			out[i][2] = 0;
		}
	}
}
