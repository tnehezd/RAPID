/**
 * @file gas_physics.h
 * @brief Physical models for gas dynamics in the protoplanetary disk.
 *
 * This module provides the core routines that describe the physical state and
 * evolution of the gas component in a 1D protoplanetary disk. It includes
 * calculations of turbulent viscosity (alpha prescription), kinematic viscosity,
 * pressure scale height, Keplerian quantities, sound speed, midplane density,
 * gas pressure, and the radial pressure gradient.
 *
 * In addition, the module contains routines for computing the viscous gas
 * radial velocity and for updating the gas surface density, pressure, and
 * pressure gradient during the simulation. These functions operate on the
 * DiskParameters structure and supply the gas‑related quantities required by
 * the dust physics and the global time‑integration loop.
 */

#ifndef GAS_PHYSICS_H
#define GAS_PHYSICS_H

#include "simulation_types.h" 

/**
 * @brief Calculates the turbulent alpha parameter \f$\alpha(r)\f$.
 *
 * The turbulent alpha parameter controls the strength of angular momentum
 * transport in the disk. It may depend on radius through a prescribed model
 * (e.g., dead zone).
 *
 * @param radial_distance   Radial position where α is evaluated.
 * @param disk_params       Pointer to disk parameter structure.
 * @return Turbulent alpha parameter at the given radius.
 */
double calculateTurbulentAlpha(double radial_distance, const DiskParameters *disk_params);

/**
 * @brief Computes the kinematic viscosity \f$\nu(r)\f$.
 *
 * The kinematic viscosity is given by \f$\nu = \alpha c_s H\f$, where \f$c_s\f$ is the sound speed
 * and \f$H\f$ is the pressure scale height. This function evaluates ν at a given
 * radial position.
 *
 * @param radial_distance   Radial position.
 * @param disk_params       Pointer to disk parameter structure.
 * @return Kinematic viscosity at the given radius.
 */
double calculateKinematicViscosity(double radial_distance, const DiskParameters *disk_params);

/**
 * @brief Computes the pressure scale height \f$H(r)\f$.
 *
 * The pressure scale height is defined as: \f$ H(r) = c_s(r)/\Omega_K(r)\f$
 * where \f$ c_s \f$ is the local sound speed and \f$ \Omega_K \f$ is the Keplerian frequency.
 * It characterizes the vertical thickness of the gas disk.
 *
 * @param radial_distance   Radial position.
 * @param disk_params       Pointer to disk parameter structure.
 * @return Pressure scale height at the given radius.
 */
double calculatePressureScaleHeight(double radial_distance, const DiskParameters *disk_params);

/**
 * @brief Computes the Keplerian orbital velocity: \f$v_K(r)\f$.
 *
 * The Keplerian velocity is given by \f$v_K = \sqrt{(GM_\star/r)}\f$. It is used in gas and
 * dust dynamical calculations.
 *
 * @param radial_distance   Radial position.
 * @param disk_params       Pointer to disk parameter structure.
 * @return Keplerian velocity at the given radius.
 */
double calculateKeplerianVelocity(double radial_distance, const DiskParameters *disk_params);

/**
 * @brief Computes the Keplerian orbital frequency \f$\Omega_K(r)\f$.
 *
 * The Keplerian frequency is defined as \f$\Omega_K = \sqrt(GM_star/r^3)\f$. It is used in
 * scale height, viscosity, and dynamical calculations.
 *
 * @param radial_distance   Radial position.
 * @param disk_params       Pointer to disk parameter structure.
 * @return Keplerian frequency at the given radius.
 */
double calculateKeplerianFrequency(double radial_distance, const DiskParameters *disk_params);

/**
 * @brief Computes the local sound speed \f$c_s(r)\f$.
 *
 * The sound speed depends on the disk temperature or scale height. It is used
 * in viscosity, pressure, and dust–gas coupling calculations.
 *
 * @param radial_distance   Radial position.
 * @param disk_params       Pointer to disk parameter structure.
 * @return Local sound speed at the given radius.
 */
double calculateLocalSoundSpeed(double radial_distance, const DiskParameters *disk_params);

/**
 * @brief Computes the midplane gas density \f$\rho_g(r)\f$.
 *
 * The midplane density is obtained from the surface density \f$\Sigma\f$ and the scale
 * height \f$H\f$ via \f$\rho_g = \Sigma / (sqrt(2\pi) H)\f$.
 *
 * @param gas_surface_density   Gas surface density at the given radius.
 * @param radial_distance       Radial position.
 * @param disk_params           Pointer to disk parameter structure.
 * @return Midplane gas density at the given radius.
 */
double calcualteMidplaneGasDensity(double gas_surface_density, double radial_distance, const DiskParameters *disk_params);

/**
 * @brief Computes the gas pressure \f$P(r)\f$.
 *
 * The gas pressure is computed using \f$P = \rho_g c_s^2\f$, where \f$\rho_g\f$ is the midplane
 * gas density and \f$c_s\f$ is the local sound speed.
 *
 * @param gas_surface_density   Gas surface density.
 * @param radial_distance       Radial position.
 * @param disk_params           Pointer to disk parameter structure.
 * @return Gas pressure at the given radius.
 */
double calculateGasPressure(double gas_surface_density, double radial_distance, const DiskParameters *disk_params);

/**
 * @brief Computes the radial pressure gradient \f$dP/dr\f$ across the grid.
 *
 * This function updates the gas_pressure_gradient_vector field in the
 * DiskParameters structure using finite differences.
 *
 * @param disk_params Pointer to disk parameter structure.
 */
void calculateGasPressureGradient(DiskParameters *disk_params);

/**
 * @brief Computes a coefficient used in the viscous gas radial velocity formula.
 *
 * This helper function evaluates the coefficient appearing in the expression
 * for the viscous radial velocity:
 *   \f$v_r = -3 / (\Sigma r^{1/2}) \cdot d/dr (\nu \Sigma r^{1/2})\f$
 *
 * @param gas_surface_density   Gas surface density.
 * @param radial_distance       Radial position.
 * @return Coefficient used in the radial velocity calculation.
 */
 double coefficientForGasRadialVelocity(double gas_surface_density, double radial_distance); 

/**
 * @brief Computes the viscous radial velocity of the gas.
 *
 * Updates the gas_velocity_vector field in the DiskParameters structure using
 * the standard viscous accretion formula:
 *   \f$v_r = -3 / (\Sigma r^{1/2}) \cdot d/dr (\nu \Sigma r^{1/2})\f$
 *
 * @param disk_params Pointer to disk parameter structure.
 */
 void calculateGasRadialVelocity(DiskParameters *disk_params);

/**
 * @brief Updates gas surface density, pressure, and pressure gradient.
 *
 * This function recomputes \f$\Sigma(r)\f$, \f$P(r)\f$, and \f$dP/dr\f$ after the gas surface density
 * has been evolved by the simulation. It ensures that all dependent quantities
 * remain consistent.
 *
 * @param sim_opts     Pointer to simulation options.
 * @param disk_params  Pointer to disk parameter structure.
 */
void refreshGasSurfaceDensityPressurePressureGradient(SimulationOptions *sim_opts, DiskParameters *disk_params);

#endif // GAS_PHYSICS_H
