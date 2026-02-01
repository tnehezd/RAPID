#ifndef GAS_PHYSICS_H
#define GAS_PHYSICS_H

#include "simulation_types.h" 

/*	Calculates turbulent alpha parameter	*/
double calculateTurbulentAlpha(double radial_distance, const DiskParameters *disk_params);

/*	Lokalis viszkozitas erteke	*/
double calculateKinematicViscosity(double radial_distance, const DiskParameters *disk_params);

/*	local scale height	*/
double calculatePressureScaleHeight(double radial_distance, const DiskParameters *disk_params);

/*	lokális kepleri sebesség	*/
double calculateKeplerianVelocity(double radial_distance, const DiskParameters *disk_params);

/*	lokalis kepleri korfrekvencia	*/
double calculateKeplerianFrequency(double radial_distance, const DiskParameters *disk_params);

/*	local sound speed		*/
double calculateLocalSoundSpeed(double radial_distance, const DiskParameters *disk_params);

/*	Suruseg a midplane-ben	*/
double calcualteMidplaneGasDensity(double gas_surface_density, double radial_distance, const DiskParameters *disk_params);

/* local pressure of the gas p = rho_gas * cs * cs kepletbol!!	*/
double calculateGasPressure(double gas_surface_density, double radial_distance, const DiskParameters *disk_params);

/*	a nyomas derivaltja	*/
void calculateGasPressureGradient(DiskParameters *disk_params);

/*	calculateGasRadialVelocity kiszamolasahoz eltarolt koefficiens	*/
double coefficientForGasRadialVelocity(double gas_surface_density, double radial_distance); // Feltételezve, hogy coefficientForGasRadialVelocity a dust_physics.c-ben van

/*	calculateGasRadialVelocity = -3/(Sigma*R^0.5)*(d/dR)(nu*Sigma*R^0.5) kiszamolasa	*/
void calculateGasRadialVelocity(DiskParameters *disk_params); // disk_params nem const, mert módosítva van az gas_velocity_vector tagja

/*	Fuggveny a sigma, p, dp kiszamolasara	*/
// disk_params nem const, mert módosítva van a gas_surface_density_vector, gas_pressure_vector, gas_pressure_gradient_vector tagjai
void refreshGasSurfaceDensityPressurePressureGradient(const SimulationOptions *sim_opts, DiskParameters *disk_params);

#endif // GAS_PHYSICS_H
