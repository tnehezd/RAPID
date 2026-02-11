// src/disk_model.c

#include "disk_model.h"   
#include "config.h"       
#include "simulation_types.h"
#include "gas_physics.h"
#include "dust_physics.h" 
#include "ascii_output.h"     
#include "utils.h" 
#include "boundary_conditions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


void readDiskParameters(DiskParameters *disk_params) {

    if (disk_params == NULL) {
        fprintf(stderr, "ERROR [readDiskParameters]: Received NULL disk_params pointer.\n");
        exit(EXIT_FAILURE);
    }

    fprintf(stderr, "DEBUG [readDiskParameters]: Calculating derived disk parameters and writing to output file.\n");

    disk_params->particle_density_dimensionless = disk_params->particle_density / SOLAR_MASS_IN_GRAMS * AU_IN_CM * AU_IN_CM * AU_IN_CM;

    fprintf(stderr, "DEBUG [readDiskParameters]: Calculated particle_density = %.2e, particle_density_dimensionless = %.2e.\n",
           disk_params->particle_density, disk_params->particle_density_dimensionless);
}


void createRadialGrid(DiskParameters *disk_params) {
	
	int i;
 	for(i = 0; i <= disk_params->grid_number+1; i++) {						
 		disk_params->radial_grid[i] = disk_params->r_min + (i-1) * disk_params->delta_r;
	}
}

void createInitialGasSurfaceDensity(DiskParameters *disk_params){	

  	int i;
  
  	for(i = 1; i <= disk_params->grid_number; i++) {
    		disk_params->gas_surface_density_vector[i] = disk_params->sigma_0 * pow(disk_params->radial_grid[i],disk_params->sigma_power_law_index);	
    }

  	applyBoundaryConditions(disk_params->gas_surface_density_vector,disk_params);
}

void createInitialGasPressure(DiskParameters *disk_params){	

  	int i;
  
  	for(i = 1; i <= disk_params->grid_number; i++) {
    		disk_params->gas_pressure_vector[i] = calculateGasPressure(disk_params->gas_surface_density_vector[i],disk_params->radial_grid[i],disk_params);
  	}
  	applyBoundaryConditions(disk_params->gas_pressure_vector,disk_params);

}

void createInitialGasPressureGradient(DiskParameters *disk_params){

	calculateGasPressureGradient(disk_params);
   	applyBoundaryConditions(disk_params->gas_pressure_gradient_vector,disk_params);

}

void createInitialGasVelocity(DiskParameters *disk_params){	
 	
	calculateGasRadialVelocity(disk_params);
  	applyBoundaryConditions(disk_params->gas_velocity_vector,disk_params);
}



void calculateDustSurfaceDensityFromRepresentativeMass(double input_dust_radii_array[][2], double *input_mass_array, double output_dust_surfacedensity_array[][3], int particle_number, const DiskParameters *disk_params) {

	int i;

	for(i=0;i<particle_number;i++){

		if((input_dust_radii_array[i][0] >= disk_params->r_min)) {
			output_dust_surfacedensity_array[i][0] = input_mass_array[i] / (2. * (input_dust_radii_array[i][0]-disk_params->delta_r/2.) * M_PI * disk_params->delta_r);	// sigma = m /(2 * r * pi * dr) --> dr is the original grid step
			output_dust_surfacedensity_array[i][1] = input_dust_radii_array[i][0];																	

  			double radial_cell_position = (input_dust_radii_array[i][0] - disk_params->r_min) / disk_params->delta_r;     				
			int radial_index = (int) floor(radial_cell_position);														
			output_dust_surfacedensity_array[i][2] = (double) radial_index;
		} else {
			memset(output_dust_surfacedensity_array[i], 0, 3 * sizeof(double));
		}
	}
}
