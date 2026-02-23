// src/dust_physics.c
#include "gas_physics.h" 
#include "config.h"       
#include "simulation_types.h" 
#include "boundary_conditions.h"
#include "simulation_core.h" 
#include "utils.h"           
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>          

double calculateTurbulentAlpha(double radial_distance, const DiskParameters *disk_params) {
    
    double alpha_r;
    alpha_r = 1.0 - 0.5 * (1.0 - disk_params->alpha_parameter_modification) * (tanh((radial_distance - disk_params->r_dze_i) / disk_params->dr_dze_i) + tanh((disk_params->r_dze_o - radial_distance) / disk_params->dr_dze_o));
    return alpha_r * disk_params->alpha_parameter;
}


double calculateKinematicViscosity(double radial_distance, const DiskParameters *disk_params) {
    
    double gas_viscosity, local_soundspeed, local_pressure_scaleheight;
    local_pressure_scaleheight = calculatePressureScaleHeight(radial_distance,disk_params);
    local_soundspeed = calculateLocalSoundSpeed(radial_distance,disk_params);
    gas_viscosity = calculateTurbulentAlpha(radial_distance, disk_params) * local_soundspeed * local_pressure_scaleheight;
    return gas_viscosity;
}

double calculatePressureScaleHeight(double radial_distance, const DiskParameters *disk_params) {

    double local_pressure_scaleheight = pow(radial_distance, 1. + disk_params->flaring_index) * disk_params->h_aspect_ratio;
    return local_pressure_scaleheight;
}

double calculateKeplerianVelocity(double radial_distance, const DiskParameters *disk_params) {
    
    return sqrt(G_DIMENSIONLESS * disk_params->stellar_mass / radial_distance);
}

double calculateKeplerianFrequency(double radial_distance, const DiskParameters *disk_params) {
    
    return sqrt(G_DIMENSIONLESS * disk_params->stellar_mass / radial_distance / radial_distance / radial_distance);
}

double calculateLocalSoundSpeed(double radial_distance, const DiskParameters *disk_params) {
    
    return calculateKeplerianFrequency(radial_distance,disk_params) * calculatePressureScaleHeight(radial_distance,disk_params);
}

double calcualteMidplaneGasDensity(double gas_surface_density, double radial_distance, const DiskParameters *disk_params) {
   
    return 1. / sqrt(2.0 * M_PI) * gas_surface_density / calculatePressureScaleHeight(radial_distance,disk_params);
}

double calculateGasPressure(double gas_surface_density, double radial_distance, const DiskParameters *disk_params) {
  
    return calcualteMidplaneGasDensity(gas_surface_density, radial_distance, disk_params) * calculateLocalSoundSpeed(radial_distance,disk_params) * calculateLocalSoundSpeed(radial_distance, disk_params);
}

void calculateGasPressureGradient(DiskParameters *disk_params) {

    int i;
    double temporary_gas_pressure, calculated_pressure_gradient[disk_params->grid_number + 2];

    for (i = 1; i <= disk_params->grid_number; i++) {
        temporary_gas_pressure = (disk_params->gas_pressure_vector[i + 1] - disk_params->gas_pressure_vector[i - 1]) / (2.0 * disk_params->delta_r);
        calculated_pressure_gradient[i] = temporary_gas_pressure;

    }

    for (i = 1; i <= disk_params->grid_number; i++) {
        disk_params->gas_pressure_gradient_vector[i] = calculated_pressure_gradient[i];
    }

}

double coefficientForGasRadialVelocity(double gas_surface_density, double radial_distance) {

    return -1.0 * (3.0 / (gas_surface_density * sqrt(radial_distance)));
}

void calculateGasRadialVelocity(DiskParameters *disk_params) {

    double temporary_gas_velocity;
    double gas_velocity_vector[disk_params->grid_number + 2];
    double gas_velocity_vectortemp[disk_params->grid_number + 1];

    int i;

    #pragma omp parallel for private(i)
    for (i = 0; i <= disk_params->grid_number + 1; i++) { 
        gas_velocity_vector[i] = disk_params->gas_surface_density_vector[i] * calculateKinematicViscosity(disk_params->radial_grid[i], disk_params) * sqrt(disk_params->radial_grid[i]);
    }

    #pragma omp parallel for private(i, temporary_gas_velocity)
    for (i = 1; i <= disk_params->grid_number; i++) { 
        temporary_gas_velocity = (gas_velocity_vector[i + 1] - gas_velocity_vector[i - 1]) / (2.0 * disk_params->delta_r); 
        gas_velocity_vectortemp[i] = coefficientForGasRadialVelocity(disk_params->gas_surface_density_vector[i], disk_params->radial_grid[i]) * temporary_gas_velocity;
    }

    for (i = 1; i <= disk_params->grid_number; i++) { 
        disk_params->gas_velocity_vector[i] = gas_velocity_vectortemp[i];
    }
}


void initializeGas2D(DiskParameters *disk_params) {
    int n_r = disk_params->grid_number;
    int n_z = disk_params->vertical_grid_number;

    for (int i_r = 0; i_r < n_r; i_r++) {
        double r = disk_params->radial_grid[i_r];
        double sigma_r = disk_params->gas_surface_density_vector[i_r];
        
        // FIZIKAI SKÁLAMAGASSÁG KISZÁMÍTÁSA (Hg)
        double Hg = calculatePressureScaleHeight(r, disk_params);
        
        // A gáz rácsának határa is legyen a Hg függvénye (pl. 4 * Hg)
        double z_limit = disk_params->vertical_grid_max * Hg; 
        double dz = 2.0 * z_limit / (double)(n_z - 1);

        long double weight_sum = 0.0;
        // 1. menet: Súlyok és koordináták
        for (int iz = 0; iz < n_z; iz++) {
            double z = -z_limit + iz * dz;
            disk_params->vertical_grid[iz] = z;
            // Itt Hg-vel kell osztani!
            weight_sum += exp(-0.5 * (z / Hg) * (z / Hg));
        }

        // 2. menet: Sűrűség és nyomás
        for (int iz = 0; iz < n_z; iz++) {
            double z = -z_limit + iz * dz;
            double f = exp(-0.5 * (z / Hg) * (z / Hg)) / weight_sum;
            
            disk_params->gas_density_2d[i_r][iz] = sigma_r * f;
            disk_params->gas_pressure_2d[i_r][iz] = calculateGasPressure(disk_params->gas_density_2d[i_r][iz], r, disk_params);
        }
    }
}

void refreshGasSurfaceDensityPressurePressureGradient(double delta_t, DiskParameters *disk_params) { 

    double gas_sigma_dot_viscosity, gas_sigma_dot_viscosity_backwards, gas_sigma_dot_viscosity_forward;
    double gas_surface_density_temp[disk_params->grid_number + 2]; 
    double gas_velocity_array[disk_params->grid_number + 2]; 

    int i;

    gas_surface_density_temp[0] = disk_params->gas_surface_density_vector[0];
    gas_surface_density_temp[disk_params->grid_number + 1] = disk_params->gas_surface_density_vector[disk_params->grid_number + 1];

    gas_velocity_array[0] = disk_params->gas_surface_density_vector[0] * calculateKinematicViscosity(disk_params->radial_grid[0], disk_params); 
    gas_velocity_array[disk_params->grid_number + 1] = disk_params->gas_surface_density_vector[disk_params->grid_number + 1] * calculateKinematicViscosity(disk_params->radial_grid[disk_params->grid_number + 1], disk_params); 

    #pragma omp parallel for
    for(i = 1; i <= disk_params->grid_number; i++) { 
        gas_velocity_array[i] = disk_params->gas_surface_density_vector[i] * calculateKinematicViscosity(disk_params->radial_grid[i], disk_params); 
    }

    for (i = 1; i <= disk_params->grid_number; i++) {
        gas_sigma_dot_viscosity = gas_velocity_array[i];
        gas_sigma_dot_viscosity_backwards = gas_velocity_array[i - 1];
        gas_sigma_dot_viscosity_forward = gas_velocity_array[i + 1];

        double gas_sigma_dot_viscosity_temporary = ftcsSecondDerivativeCoefficient(disk_params->radial_grid[i], disk_params) * (gas_sigma_dot_viscosity_forward - 2.0 * gas_sigma_dot_viscosity + gas_sigma_dot_viscosity_backwards) / (disk_params->delta_r * disk_params->delta_r) +
                      ftcsFirstDerivativeCoefficient(disk_params->radial_grid[i], disk_params) * (gas_sigma_dot_viscosity_forward - gas_sigma_dot_viscosity_backwards) / (2.0 * disk_params->delta_r);
        
        gas_surface_density_temp[i] = gas_velocity_array[i] + delta_t * gas_sigma_dot_viscosity_temporary; 
    }

    #pragma omp parallel for
    for (i = 1; i <= disk_params->grid_number; i++) { 
        disk_params->gas_surface_density_vector[i] = gas_surface_density_temp[i] / calculateKinematicViscosity(disk_params->radial_grid[i], disk_params);
        disk_params->gas_pressure_vector[i] = calculateGasPressure(disk_params->gas_surface_density_vector[i], disk_params->radial_grid[i], disk_params);
    }

    calculateGasPressureGradient(disk_params);
    applyBoundaryConditions(disk_params->gas_surface_density_vector, disk_params); 
    applyBoundaryConditions(disk_params->gas_pressure_vector, disk_params);
    applyBoundaryConditions(disk_params->gas_pressure_gradient_vector, disk_params);

}