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
 * @brief Calculates the turbulent alpha parameter α(r).
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
 * @brief Computes the kinematic viscosity ν(r).
 *
 * The kinematic viscosity is given by ν = α c_s H, where c_s is the sound speed
 * and H is the pressure scale height. This function evaluates ν at a given
 * radial position.
 *
 * @param radial_distance   Radial position.
 * @param disk_params       Pointer to disk parameter structure.
 * @return Kinematic viscosity at the given radius.
 */
double calculateKinematicViscosity(double radial_distance, const DiskParameters *disk_params);

/**
 * @brief Computes the pressure scale height H(r).
 *
 * The pressure scale height is defined as H = c_s / Ω_K, where c_s is the local
 * sound speed and Ω_K is the Keplerian frequency. It characterizes the vertical
 * thickness of the gas disk.
 *
 * @param radial_distance   Radial position.
 * @param disk_params       Pointer to disk parameter structure.
 * @return Pressure scale height at the given radius.
 */
double calculatePressureScaleHeight(double radial_distance, const DiskParameters *disk_params);

/**
 * @brief Computes the Keplerian orbital velocity v_K(r).
 *
 * The Keplerian velocity is given by v_K = sqrt(GM_star/r). It is used in gas and
 * dust dynamical calculations.
 *
 * @param radial_distance   Radial position.
 * @param disk_params       Pointer to disk parameter structure.
 * @return Keplerian velocity at the given radius.
 */
double calculateKeplerianVelocity(double radial_distance, const DiskParameters *disk_params);

/**
 * @brief Computes the Keplerian orbital frequency Ω_K(r).
 *
 * The Keplerian frequency is defined as Ω_K = sqrt(GM_star/r^3). It is used in
 * scale height, viscosity, and dynamical calculations.
 *
 * @param radial_distance   Radial position.
 * @param disk_params       Pointer to disk parameter structure.
 * @return Keplerian frequency at the given radius.
 */
double calculateKeplerianFrequency(double radial_distance, const DiskParameters *disk_params);

/**
 * @brief Computes the local sound speed c_s(r).
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
 * @brief Computes the midplane gas density ρ_g(r).
 *
 * The midplane density is obtained from the surface density Σ and the scale
 * height H via ρ_g = Σ / (sqrt(2π) H).
 *
 * @param gas_surface_density   Gas surface density at the given radius.
 * @param radial_distance       Radial position.
 * @param disk_params           Pointer to disk parameter structure.
 * @return Midplane gas density at the given radius.
 */
double calcualteMidplaneGasDensity(double gas_surface_density, double radial_distance, const DiskParameters *disk_params);

/**
 * @brief Computes the gas pressure P(r).
 *
 * The gas pressure is computed using P = ρ_g c_s^2, where ρ_g is the midplane
 * gas density and c_s is the local sound speed.
 *
 * @param gas_surface_density   Gas surface density.
 * @param radial_distance       Radial position.
 * @param disk_params           Pointer to disk parameter structure.
 * @return Gas pressure at the given radius.
 */
double calculateGasPressure(double gas_surface_density, double radial_distance, const DiskParameters *disk_params);

/**
 * @brief Computes the radial pressure gradient dP/dr across the grid.
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
 *   v_r = -3 / (Σ r^{1/2}) * d/dr (ν Σ r^{1/2})
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
 *   v_r = -3 / (Σ r^{1/2}) * d/dr (ν Σ r^{1/2})
 *
 * @param disk_params Pointer to disk parameter structure.
 */
 void calculateGasRadialVelocity(DiskParameters *disk_params);

/**
 * @brief Updates gas surface density, pressure, and pressure gradient.
 *
 * This function recomputes Σ(r), P(r), and dP/dr after the gas surface density
 * has been evolved by the simulation. It ensures that all dependent quantities
 * remain consistent.
 *
 * @param sim_opts     Pointer to simulation options.
 * @param disk_params  Pointer to disk parameter structure.
 */
void refreshGasSurfaceDensityPressurePressureGradient(const SimulationOptions *sim_opts, DiskParameters *disk_params);

#endif // GAS_PHYSICS_H
