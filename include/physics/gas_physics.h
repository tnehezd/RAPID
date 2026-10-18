/**
 * @file gas_physics.h
 * @brief Physical models for gas dynamics in the protoplanetary disk.
 *
 * This module provides the core routines that describe the physical state and
 * evolution of the gas component in a 1D protoplanetary disk. It includes
 * calculations of turbulent viscosity (alpha prescription), kinematic viscosity,
 * Keplerian quantities and the radial pressure gradient. Thermal closure,
 * sound speed, scale height, midplane density and pressure are provided by eos.h.
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
#include "eos.h"

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

/**
 * @brief Updates the disk self-gravity contribution to the rotation curve.
 *
 * Sums the softened radial gravitational pull of all gas rings (complete elliptic
 * integral kernel, softening 0.6 H) and stores \f$g_{sg}/r\f$ in
 * DiskParameters::selfgravity_omega2_extra, so that
 * \f$\Omega^2 = \Omega_K^2 + g_{sg}/r\f$ in calculateKeplerianFrequency().
 * Does nothing if self-gravity is disabled.
 *
 * @param disk_params Pointer to disk parameter structure.
 */
void updateSelfGravity(DiskParameters *disk_params);

#endif // GAS_PHYSICS_H
