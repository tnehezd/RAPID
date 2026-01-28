// src/dust_physics.c
#include "gas_physics.h" // A saját headerjét mindig includolni kell
#include "config.h"       // Szükséges lehet a globális konstansokhoz (pl. PARTICLE_NUMBER, AU2CM, RMIN, RMAX, NGRID, G_GRAV_CONST, STAR, SDCONV, CMPSECTOAUPYRP2PI, uFrag, fFrag, PDENSITYDIMLESS, HASP, M_PI, DD, sim_opts->dzone, sim_opts->twopop, RMIN, RMAX, FLIND, alpha_visc, a_mod, r_dze_i, r_dze_o, Dr_dze_i, Dr_dze_o)
#include "simulation_types.h" // Például output_files_t, disk_t struktúrákhoz
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
double calculateTurbulentAlpha(double r, const disk_t *disk_params) {
    double alpha_r;
    alpha_r = 1.0 - 0.5 * (1.0 - disk_params->a_mod) * (tanh((r - disk_params->r_dze_i) / disk_params->Dr_dze_i) + tanh((disk_params->r_dze_o - r) / disk_params->Dr_dze_o));
    return alpha_r * disk_params->alpha_visc;
}



/*  Lokalis viszkozitas erteke  */
double calculateKinematicViscosity(double r, const disk_t *disk_params) {
    double nu;
    double cs, H;

    H = calculatePressureScaleHeight(r,disk_params);
    cs = calculateLocalSoundSpeed(r,disk_params);

    nu = calculateTurbulentAlpha(r, disk_params) * cs * H;
    return nu;
}

/*  local scale height  */
double calculatePressureScaleHeight(double r, const disk_t *disk_params) {

    if (disk_params == NULL) {
        fprintf(stderr, "ERROR [scale_height]: disk_params is NULL!\n");
        return 0.0; // Vagy valamilyen hibakód/NaN
    }

    // Itt van az eredeti számítás
    double calculated_result = pow(r, 1. + disk_params->FLIND) * disk_params->HASP;
    return calculated_result;
}

/*  lokális kepleri sebesség    */
double calculateKeplerianVelocity(double r, const disk_t *disk_params) {
    return sqrt(G_GRAV_CONST * disk_params->STAR_MASS / r);
}

/*  lokalis kepleri korfrekvencia   */
double calculateKeplerianFrequency(double r, const disk_t *disk_params) {
    return sqrt(G_GRAV_CONST * disk_params->STAR_MASS / r / r / r);
}

/*  local sound speed       */
double calculateLocalSoundSpeed(double r, const disk_t *disk_params) {
    return calculateKeplerianFrequency(r,disk_params) * calculatePressureScaleHeight(r,disk_params);
}

/*  Suruseg a midplane-ben  */
double calcualteMidplaneGasDensity(double sigma, double r, const disk_t *disk_params) {
    return 1. / sqrt(2.0 * M_PI) * sigma / calculatePressureScaleHeight(r,disk_params);
}

/* local pressure of the gas p = rho_gas * cs * cs kepletbol!!  */
double calculateGasPressure(double sigma, double r, const disk_t *disk_params) {
    return calcualteMidplaneGasDensity(sigma, r, disk_params) * calculateLocalSoundSpeed(r,disk_params) * calculateLocalSoundSpeed(r, disk_params);
}

/*  a nyomas derivaltja */
void calculateGasPressureGradient(disk_t *disk_params) {
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

/*  calculateGasRadialVelocity kiszamolasahoz eltarolt koefficiens   */
double coefficientForGasRadialVelocity(double sigma, double r) {
    return -1.0 * (3.0 / (sigma * sqrt(r)));
}

/*  calculateGasRadialVelocity = -3/(Sigma*R^0.5)*(d/dR)(nu*Sigma*R^0.5) kiszamolasa */
void calculateGasRadialVelocity(disk_t *disk_params) {

    double tempug;
    // Lokális tömbök, méret NGRID-hez igazítva disk_params-ból
    double ugvec[disk_params->NGRID + 2];
    double ugvectemp[disk_params->NGRID + 1]; // Eredeti kód NGRID+1-et használt

    int i;

    // Első ciklus: feltölti a lokális ugvec tömböt
    #pragma omp parallel for private(i)
    for (i = 0; i <= disk_params->NGRID + 1; i++) { // Használd a disk_params->NGRID-et
        // Hozzáférés a disk_params tagjaihoz
        ugvec[i] = disk_params->sigmavec[i] * calculateKinematicViscosity(disk_params->rvec[i], disk_params) * sqrt(disk_params->rvec[i]);
        // Megjegyzés: A sqrt() függvénynek általában csak egy double paramétere van.
        // Ha valami komplexebb számítást akarsz, akkor lehet, hogy egy saját
        // függvényt hívsz, amihez disk_params is kell. Ellenőrizd a sqrt prototípusát!
    }

    // Második ciklus: feltölti a lokális ugvectemp tömböt
    #pragma omp parallel for private(i, tempug)
    for (i = 1; i <= disk_params->NGRID; i++) { // Használd a disk_params->NGRID-et
        tempug = (ugvec[i + 1] - ugvec[i - 1]) / (2.0 * disk_params->DD); // Használd a disk_params->DD-t
        // coefficientForGasRadialVelocity hívása, ha szükséges, átadva neki a disk_params-ot
        ugvectemp[i] = coefficientForGasRadialVelocity(disk_params->sigmavec[i], disk_params->rvec[i]) * tempug;
    }

    // Harmadik ciklus: Az eredményt bemásolja a disk_params->ugvec-be
    for (i = 1; i <= disk_params->NGRID; i++) { // Használd a disk_params->NGRID-et
        disk_params->ugvec[i] = ugvectemp[i]; // Így éri el a struktúrán belüli ugvec-et
    }
}

/*  Fuggveny a sigma, p, dp kiszamolasara   */
void refreshGasSurfaceDensityPressurePressureGradient(const simulation_options_t *sim_opts, disk_t *disk_params) { // Added sim_opts

    double u, u_bi, u_fi;
    double sigma_temp[disk_params->NGRID + 2]; // Use disk_params->NGRID
    double uvec[disk_params->NGRID + 2];     // Use disk_params->NGRID

    int i;

    // Boundary conditions - access via disk_params
    sigma_temp[0] = disk_params->sigmavec[0];
    sigma_temp[disk_params->NGRID + 1] = disk_params->sigmavec[disk_params->NGRID + 1];

    // uvec temporary array initialization
    uvec[0] = disk_params->sigmavec[0] * calculateKinematicViscosity(disk_params->rvec[0], disk_params); // Use disk_params->rvec
    uvec[disk_params->NGRID + 1] = disk_params->sigmavec[disk_params->NGRID + 1] * calculateKinematicViscosity(disk_params->rvec[disk_params->NGRID + 1], disk_params); // Use disk_params->rvec

    #pragma omp parallel for
    for(i = 1; i <= disk_params->NGRID; i++) { // Use disk_params->NGRID
        uvec[i] = disk_params->sigmavec[i] * calculateKinematicViscosity(disk_params->rvec[i], disk_params); // Use disk_params->sigmavec and disk_params->rvec
    }

    // This loop is critical due to data dependencies. Keep it sequential for correctness
    for (i = 1; i <= disk_params->NGRID; i++) { // Use disk_params->NGRID
        u = uvec[i];
        u_bi = uvec[i - 1];
        u_fi = uvec[i + 1];

        // Access DD and deltat through the appropriate structs
        // Assuming ftcsSecondDerivativeCoefficient and ftcsFirstDerivativeCoefficient also take disk_params (and sim_opts if they need it)
        double temp = ftcsSecondDerivativeCoefficient(disk_params->rvec[i], disk_params) * (u_fi - 2.0 * u + u_bi) / (disk_params->DD * disk_params->DD) +
                      ftcsFirstDerivativeCoefficient(disk_params->rvec[i], disk_params) * (u_fi - u_bi) / (2.0 * disk_params->DD);
        
        sigma_temp[i] = uvec[i] + sim_opts->DT * temp; // Use sim_opts->DT for deltat
    }

    // This loop is parallelizable
    #pragma omp parallel for
    for (i = 1; i <= disk_params->NGRID; i++) { // Use disk_params->NGRID
        // Update disk_params' own arrays
        disk_params->sigmavec[i] = sigma_temp[i] / calculateKinematicViscosity(disk_params->rvec[i], disk_params);
        disk_params->pressvec[i] = calculateGasPressure(disk_params->sigmavec[i], disk_params->rvec[i], disk_params); // Assuming calculateGasPressure takes disk_params
    }

    // These calls likely remain sequential or require their own internal OpenMP if large
    // If applyBoundaryConditions, dpress also update members of disk_params, they should take disk_params as a parameter.
    // And if they are modifying the *content* of the arrays within disk_params, then disk_params should NOT be const in *their* parameter list.
    // However, since refreshGasSurfaceDensityPressurePressureGradient is modifying them, disk_params *here* cannot be const.
    // Let's remove 'const' from disk_params in refreshGasSurfaceDensityPressurePressureGradient signature if it modifies them.
    // void refreshGasSurfaceDensityPressurePressureGradient(disk_t *disk_params, const simulation_options_t *sim_opts) { ... }
    
    // Assuming these helper functions need disk_params to access *its* internal arrays
    calculateGasPressureGradient(disk_params); // Assuming dpress takes arrays and disk_params
    applyBoundaryConditions(disk_params->sigmavec, disk_params); // First argument is the array, second is the disk_t pointer
    applyBoundaryConditions(disk_params->pressvec, disk_params);
    applyBoundaryConditions(disk_params->dpressvec, disk_params);
}