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
#include <omp.h>             // OpenMP támogatáshoz

// Globális változó deklarációk, ha nem lennének meg máshol (pl. config.h)
// Fontos: ezeknek a típusoknak egyezniük kell a config.h-ban deklaráltakkal!
// Ha már szerepelnek a config.h-ban, akkor ezeket innen törölni kell,
// vagy csak az extern kulcsszót meghagyni!




/*	alpha turbulens paraméter kiszámolása --> alfa csökkentése alpha_r-rel	*/
double calculateTurbulentAlpha(double r, const DiskParameters *disk_params) {
    double alpha_r;
    alpha_r = 1.0 - 0.5 * (1.0 - disk_params->alpha_parameter_modification) * (tanh((r - disk_params->r_dze_i) / disk_params->dr_dze_i) + tanh((disk_params->r_dze_o - r) / disk_params->dr_dze_o));
    return alpha_r * disk_params->alpha_parameter;
}



/*  Lokalis viszkozitas erteke  */
double calculateKinematicViscosity(double r, const DiskParameters *disk_params) {
    double nu;
    double cs, H;

    H = calculatePressureScaleHeight(r,disk_params);
    cs = calculateLocalSoundSpeed(r,disk_params);

    nu = calculateTurbulentAlpha(r, disk_params) * cs * H;
    return nu;
}

/*  local scale height  */
double calculatePressureScaleHeight(double r, const DiskParameters *disk_params) {

    if (disk_params == NULL) {
        fprintf(stderr, "ERROR [scale_height]: disk_params is NULL!\n");
        return 0.0; // Vagy valamilyen hibakód/NaN
    }

    // Itt van az eredeti számítás
    double calculated_result = pow(r, 1. + disk_params->flaring_index) * disk_params->h_aspect_ratio;
    return calculated_result;
}

/*  lokális kepleri sebesség    */
double calculateKeplerianVelocity(double r, const DiskParameters *disk_params) {
    return sqrt(G_DIMENSIONLESS * disk_params->stellar_mass / r);
}

/*  lokalis kepleri korfrekvencia   */
double calculateKeplerianFrequency(double r, const DiskParameters *disk_params) {
    return sqrt(G_DIMENSIONLESS * disk_params->stellar_mass / r / r / r);
}

/*  local sound speed       */
double calculateLocalSoundSpeed(double r, const DiskParameters *disk_params) {
    return calculateKeplerianFrequency(r,disk_params) * calculatePressureScaleHeight(r,disk_params);
}

/*  Suruseg a midplane-ben  */
double calcualteMidplaneGasDensity(double sigma, double r, const DiskParameters *disk_params) {
    return 1. / sqrt(2.0 * M_PI) * sigma / calculatePressureScaleHeight(r,disk_params);
}

/* local pressure of the gas p = rho_gas * cs * cs kepletbol!!  */
double calculateGasPressure(double sigma, double r, const DiskParameters *disk_params) {
    return calcualteMidplaneGasDensity(sigma, r, disk_params) * calculateLocalSoundSpeed(r,disk_params) * calculateLocalSoundSpeed(r, disk_params);
}

/*  a nyomas derivaltja */
void calculateGasPressureGradient(DiskParameters *disk_params) {
    int i;
    double ptemp, pvec[disk_params->grid_number + 2];

    for (i = 1; i <= disk_params->grid_number; i++) {
        ptemp = (disk_params->gas_pressure_vector[i + 1] - disk_params->gas_pressure_vector[i - 1]) / (2.0 * disk_params->delta_r);
        pvec[i] = ptemp;

    }
    for (i = 1; i <= disk_params->grid_number; i++) {
        disk_params->gas_pressure_gradient_vector[i] = pvec[i];
    }


}

/*  calculateGasRadialVelocity kiszamolasahoz eltarolt koefficiens   */
double coefficientForGasRadialVelocity(double sigma, double r) {
    return -1.0 * (3.0 / (sigma * sqrt(r)));
}

/*  calculateGasRadialVelocity = -3/(Sigma*R^0.5)*(d/dR)(nu*Sigma*R^0.5) kiszamolasa */
void calculateGasRadialVelocity(DiskParameters *disk_params) {

    double tempug;
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
    #pragma omp parallel for private(i, tempug)
    for (i = 1; i <= disk_params->grid_number; i++) { // Használd a disk_params->grid_number-et
        tempug = (gas_velocity_vector[i + 1] - gas_velocity_vector[i - 1]) / (2.0 * disk_params->delta_r); // Használd a disk_params->delta_r-t
        // coefficientForGasRadialVelocity hívása, ha szükséges, átadva neki a disk_params-ot
        gas_velocity_vectortemp[i] = coefficientForGasRadialVelocity(disk_params->gas_surface_density_vector[i], disk_params->radial_grid[i]) * tempug;
    }

    // Harmadik ciklus: Az eredményt bemásolja a disk_params->gas_velocity_vector-be
    for (i = 1; i <= disk_params->grid_number; i++) { // Használd a disk_params->grid_number-et
        disk_params->gas_velocity_vector[i] = gas_velocity_vectortemp[i]; // Így éri el a struktúrán belüli gas_velocity_vector-et
    }
}

/*  Fuggveny a sigma, p, dp kiszamolasara   */
void refreshGasSurfaceDensityPressurePressureGradient(const SimulationOptions *sim_opts, DiskParameters *disk_params) { // Added sim_opts

    double u, u_bi, u_fi;
    double sigma_temp[disk_params->grid_number + 2]; // Use disk_params->grid_number
    double uvec[disk_params->grid_number + 2];     // Use disk_params->grid_number

    int i;

    // Boundary conditions - access via disk_params
    sigma_temp[0] = disk_params->gas_surface_density_vector[0];
    sigma_temp[disk_params->grid_number + 1] = disk_params->gas_surface_density_vector[disk_params->grid_number + 1];

    // uvec temporary array initialization
    uvec[0] = disk_params->gas_surface_density_vector[0] * calculateKinematicViscosity(disk_params->radial_grid[0], disk_params); // Use disk_params->radial_grid
    uvec[disk_params->grid_number + 1] = disk_params->gas_surface_density_vector[disk_params->grid_number + 1] * calculateKinematicViscosity(disk_params->radial_grid[disk_params->grid_number + 1], disk_params); // Use disk_params->radial_grid

    #pragma omp parallel for
    for(i = 1; i <= disk_params->grid_number; i++) { // Use disk_params->grid_number
        uvec[i] = disk_params->gas_surface_density_vector[i] * calculateKinematicViscosity(disk_params->radial_grid[i], disk_params); // Use disk_params->gas_surface_density_vector and disk_params->radial_grid
    }

    // This loop is critical due to data dependencies. Keep it sequential for correctness
    for (i = 1; i <= disk_params->grid_number; i++) { // Use disk_params->grid_number
        u = uvec[i];
        u_bi = uvec[i - 1];
        u_fi = uvec[i + 1];

        // Access delta_r and deltat through the appropriate structs
        // Assuming ftcsSecondDerivativeCoefficient and ftcsFirstDerivativeCoefficient also take disk_params (and sim_opts if they need it)
        double temp = ftcsSecondDerivativeCoefficient(disk_params->radial_grid[i], disk_params) * (u_fi - 2.0 * u + u_bi) / (disk_params->delta_r * disk_params->delta_r) +
                      ftcsFirstDerivativeCoefficient(disk_params->radial_grid[i], disk_params) * (u_fi - u_bi) / (2.0 * disk_params->delta_r);
        
        sigma_temp[i] = uvec[i] + sim_opts->user_defined_time_step * temp; // Use sim_opts->user_defined_time_step for deltat
    }

    // This loop is parallelizable
    #pragma omp parallel for
    for (i = 1; i <= disk_params->grid_number; i++) { // Use disk_params->grid_number
        // Update disk_params' own arrays
        disk_params->gas_surface_density_vector[i] = sigma_temp[i] / calculateKinematicViscosity(disk_params->radial_grid[i], disk_params);
        disk_params->gas_pressure_vector[i] = calculateGasPressure(disk_params->gas_surface_density_vector[i], disk_params->radial_grid[i], disk_params); // Assuming calculateGasPressure takes disk_params
    }

    // These calls likely remain sequential or require their own internal OpenMP if large
    // If applyBoundaryConditions, dpress also update members of disk_params, they should take disk_params as a parameter.
    // And if they are modifying the *content* of the arrays within disk_params, then disk_params should NOT be const in *their* parameter list.
    // However, since refreshGasSurfaceDensityPressurePressureGradient is modifying them, disk_params *here* cannot be const.
    // Let's remove 'const' from disk_params in refreshGasSurfaceDensityPressurePressureGradient signature if it modifies them.
    // void refreshGasSurfaceDensityPressurePressureGradient(DiskParameters *disk_params, const SimulationOptions *sim_opts) { ... }
    
    // Assuming these helper functions need disk_params to access *its* internal arrays
    calculateGasPressureGradient(disk_params); // Assuming dpress takes arrays and disk_params
    applyBoundaryConditions(disk_params->gas_surface_density_vector, disk_params); // First argument is the array, second is the DiskParameters pointer
    applyBoundaryConditions(disk_params->gas_pressure_vector, disk_params);
    applyBoundaryConditions(disk_params->gas_pressure_gradient_vector, disk_params);
}