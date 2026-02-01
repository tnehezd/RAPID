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

/*	alpha turbulens paraméter kiszámolása --> alfa csökkentése alpha_r-rel	*/
double calculateTurbulentAlpha(double radial_distance, const DiskParameters *disk_params) {
    
    double alpha_r;
    alpha_r = 1.0 - 0.5 * (1.0 - disk_params->alpha_parameter_modification) * (tanh((radial_distance - disk_params->r_dze_i) / disk_params->dr_dze_i) + tanh((disk_params->r_dze_o - radial_distance) / disk_params->dr_dze_o));
    return alpha_r * disk_params->alpha_parameter;
}


/*  Lokalis viszkozitas erteke  */
double calculateKinematicViscosity(double radial_distance, const DiskParameters *disk_params) {
    
    double gas_viscosity, local_soundspeed, local_pressure_scaleheight;
    local_pressure_scaleheight = calculatePressureScaleHeight(radial_distance,disk_params);
    local_soundspeed = calculateLocalSoundSpeed(radial_distance,disk_params);
    gas_viscosity = calculateTurbulentAlpha(radial_distance, disk_params) * local_soundspeed * local_pressure_scaleheight;
    return gas_viscosity;
}

/*  local scale height  */
double calculatePressureScaleHeight(double radial_distance, const DiskParameters *disk_params) {

    double local_pressure_scaleheight = pow(radial_distance, 1. + disk_params->flaring_index) * disk_params->h_aspect_ratio;
    return local_pressure_scaleheight;
}

/*  lokális kepleri sebesség    */
double calculateKeplerianVelocity(double radial_distance, const DiskParameters *disk_params) {
    
    return sqrt(G_DIMENSIONLESS * disk_params->stellar_mass / radial_distance);
}

/*  lokalis kepleri korfrekvencia   */
double calculateKeplerianFrequency(double radial_distance, const DiskParameters *disk_params) {
    
    return sqrt(G_DIMENSIONLESS * disk_params->stellar_mass / radial_distance / radial_distance / radial_distance);
}

/*  local sound speed       */
double calculateLocalSoundSpeed(double radial_distance, const DiskParameters *disk_params) {
    
    return calculateKeplerianFrequency(radial_distance,disk_params) * calculatePressureScaleHeight(radial_distance,disk_params);
}

/*  Suruseg a midplane-ben  */
double calcualteMidplaneGasDensity(double gas_surface_density, double radial_distance, const DiskParameters *disk_params) {
   
    return 1. / sqrt(2.0 * M_PI) * gas_surface_density / calculatePressureScaleHeight(radial_distance,disk_params);
}

/* local pressure of the gas p = rho_gas * cs * cs kepletbol!!  */
double calculateGasPressure(double gas_surface_density, double radial_distance, const DiskParameters *disk_params) {
  
    return calcualteMidplaneGasDensity(gas_surface_density, radial_distance, disk_params) * calculateLocalSoundSpeed(radial_distance,disk_params) * calculateLocalSoundSpeed(radial_distance, disk_params);
}

/*  a nyomas derivaltja */
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

/*  calculateGasRadialVelocity kiszamolasahoz eltarolt koefficiens   */
double coefficientForGasRadialVelocity(double gas_surface_density, double radial_distance) {

    return -1.0 * (3.0 / (gas_surface_density * sqrt(radial_distance)));
}

/*  calculateGasRadialVelocity = -3/(Sigma*R^0.5)*(d/dR)(nu*Sigma*R^0.5) kiszamolasa */
void calculateGasRadialVelocity(DiskParameters *disk_params) {

    double temporary_gas_velocity;
    // Lokális tömbök, méret grid_number-hez igazítva disk_params-ból
    double gas_velocity_vector[disk_params->grid_number + 2];
    double gas_velocity_vectortemp[disk_params->grid_number + 1]; // Eredeti kód grid_number+1-et használt

    int i;

    // Első ciklus: feltölti a lokális gas_velocity_vector tömböt
    #pragma omp parallel for private(i)
    for (i = 0; i <= disk_params->grid_number + 1; i++) { // Használd a disk_params->grid_number-et
        // Hozzáférés a disk_params tagjaihoz
        gas_velocity_vector[i] = disk_params->gas_surface_density_vector[i] * calculateKinematicViscosity(disk_params->radial_grid[i], disk_params) * sqrt(disk_params->radial_grid[i]);
        // Megjegyzés: A sqrt() függvénynek általában csak egy double paramétere van.
        // Ha valami komplexebb számítást akarsz, akkor lehet, hogy egy saját
        // függvényt hívsz, amihez disk_params is kell. Ellenőrizd a sqrt prototípusát!
    }

    // Második ciklus: feltölti a lokális gas_velocity_vectortemp tömböt
    #pragma omp parallel for private(i, temporary_gas_velocity)
    for (i = 1; i <= disk_params->grid_number; i++) { // Használd a disk_params->grid_number-et
        temporary_gas_velocity = (gas_velocity_vector[i + 1] - gas_velocity_vector[i - 1]) / (2.0 * disk_params->delta_r); // Használd a disk_params->delta_r-t
        // coefficientForGasRadialVelocity hívása, ha szükséges, átadva neki a disk_params-ot
        gas_velocity_vectortemp[i] = coefficientForGasRadialVelocity(disk_params->gas_surface_density_vector[i], disk_params->radial_grid[i]) * temporary_gas_velocity;
    }

    // Harmadik ciklus: Az eredményt bemásolja a disk_params->gas_velocity_vector-be
    for (i = 1; i <= disk_params->grid_number; i++) { // Használd a disk_params->grid_number-et
        disk_params->gas_velocity_vector[i] = gas_velocity_vectortemp[i]; // Így éri el a struktúrán belüli gas_velocity_vector-et
    }
}

/*  Fuggveny a sigma, p, dp kiszamolasara   */
void refreshGasSurfaceDensityPressurePressureGradient(const SimulationOptions *sim_opts, DiskParameters *disk_params) { // Added sim_opts

    double gas_sigma_dot_viscosity, gas_sigma_dot_viscosity_backwards, gas_sigma_dot_viscosity_forward;
    double gas_surface_density_temp[disk_params->grid_number + 2]; // Use disk_params->grid_number
    double gas_velocity_array[disk_params->grid_number + 2];     // Use disk_params->grid_number

    int i;

    // Boundary conditions - access via disk_params
    gas_surface_density_temp[0] = disk_params->gas_surface_density_vector[0];
    gas_surface_density_temp[disk_params->grid_number + 1] = disk_params->gas_surface_density_vector[disk_params->grid_number + 1];

    // gas_velocity_array temporary array initialization
    gas_velocity_array[0] = disk_params->gas_surface_density_vector[0] * calculateKinematicViscosity(disk_params->radial_grid[0], disk_params); // Use disk_params->radial_grid
    gas_velocity_array[disk_params->grid_number + 1] = disk_params->gas_surface_density_vector[disk_params->grid_number + 1] * calculateKinematicViscosity(disk_params->radial_grid[disk_params->grid_number + 1], disk_params); // Use disk_params->radial_grid

    #pragma omp parallel for
    for(i = 1; i <= disk_params->grid_number; i++) { // Use disk_params->grid_number
        gas_velocity_array[i] = disk_params->gas_surface_density_vector[i] * calculateKinematicViscosity(disk_params->radial_grid[i], disk_params); // Use disk_params->gas_surface_density_vector and disk_params->radial_grid
    }

    // This loop is critical due to data dependencies. Keep it sequential for correctness
    for (i = 1; i <= disk_params->grid_number; i++) { // Use disk_params->grid_number
        gas_sigma_dot_viscosity = gas_velocity_array[i];
        gas_sigma_dot_viscosity_backwards = gas_velocity_array[i - 1];
        gas_sigma_dot_viscosity_forward = gas_velocity_array[i + 1];

        double gas_sigma_dot_viscosity_temporary = ftcsSecondDerivativeCoefficient(disk_params->radial_grid[i], disk_params) * (gas_sigma_dot_viscosity_forward - 2.0 * gas_sigma_dot_viscosity + gas_sigma_dot_viscosity_backwards) / (disk_params->delta_r * disk_params->delta_r) +
                      ftcsFirstDerivativeCoefficient(disk_params->radial_grid[i], disk_params) * (gas_sigma_dot_viscosity_forward - gas_sigma_dot_viscosity_backwards) / (2.0 * disk_params->delta_r);
        
        gas_surface_density_temp[i] = gas_velocity_array[i] + sim_opts->user_defined_time_step * gas_sigma_dot_viscosity_temporary; // Use sim_opts->user_defined_time_step for deltat
    }

    // This loop is parallelizable
    #pragma omp parallel for
    for (i = 1; i <= disk_params->grid_number; i++) { // Use disk_params->grid_number
        // Update disk_params' own arrays
        disk_params->gas_surface_density_vector[i] = gas_surface_density_temp[i] / calculateKinematicViscosity(disk_params->radial_grid[i], disk_params);
        disk_params->gas_pressure_vector[i] = calculateGasPressure(disk_params->gas_surface_density_vector[i], disk_params->radial_grid[i], disk_params); // Assuming calculateGasPressure takes disk_params
    }

    calculateGasPressureGradient(disk_params); // Assuming dpress takes arrays and disk_params
    applyBoundaryConditions(disk_params->gas_surface_density_vector, disk_params); // First argument is the array, second is the DiskParameters pointer
    applyBoundaryConditions(disk_params->gas_pressure_vector, disk_params);
    applyBoundaryConditions(disk_params->gas_pressure_gradient_vector, disk_params);
}