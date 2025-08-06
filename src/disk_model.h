C
/**
 * @file disk_model.h
 * @brief Declarations for functions handling gas disk properties and evolution.
 *
 * This header file provides the declarations for functions related to
 * initializing the protoplanetary gas disk, setting up its initial profiles,
 * calculating various physical properties, and evolving the disk over time.
 * The corresponding implementations are found in `disk_model.c`.
 */

#ifndef DISK_MODEL_H
#define DISK_MODEL_H

#include "config.h"
#include "simulation_types.h" // Assumed to contain the definition of disk_t
#include "globals.h"

/**
 * @brief Initializes the radial grid cell positions.
 *
 * This function populates the `rvec` array within the `disk_t` structure
 * with the radial positions of the grid cell centers. It includes the
 * setup for boundary (ghost) cells.
 *
 * @param disk_params Pointer to the `disk_t` structure containing grid parameters (`RMIN`, `DD`, `NGRID`).
 */
void initialize_grid_cells(disk_t *disk_params);

/**
 * @brief Sets up the initial gas surface density profile.
 *
 * This function calculates the initial gas surface density (`sigmavec`)
 * for each grid cell, typically based on a power-law or other predefined profile.
 * It also applies necessary boundary conditions.
 *
 * @param disk_params Pointer to the `disk_t` structure containing initial gas parameters
 * and grid data (`rvec`, `NGRID`).
 */
void initial_gas_surface_density_profile(disk_t *disk_params);

/**
 * @brief Sets up the initial gas pressure profile.
 *
 * This function calculates the initial gas pressure (`pressvec`) for each
 * grid cell, typically derived from the gas surface density and local sound speed.
 * Boundary conditions are applied after calculation.
 *
 * @param disk_params Pointer to the `disk_t` structure containing gas and grid data.
 */
void initial_gas_pressure_profile(disk_t *disk_params);

/**
 * @brief Sets up the initial gas pressure gradient profile.
 *
 * This function calculates the initial radial gas pressure gradient (`dpressvec`)
 * across the disk, which is crucial for gas dynamics. Boundary conditions are applied.
 *
 * @param disk_params Pointer to the `disk_t` structure containing pressure data.
 */
void initial_gas_pressure_gradient_profile(disk_t *disk_params);

/**
 * @brief Sets up the initial gas radial velocity profile.
 *
 * This function calculates the initial radial velocity (`ugvec`) of the gas
 * across the disk, consistent with the initial surface density and pressure profiles.
 * Boundary conditions are applied.
 *
 * @param disk_params Pointer to the `disk_t` structure containing relevant disk properties.
 */
void initial_gas_velocity_profile(disk_t *disk_params);

/**
 * @brief Calculates the turbulent alpha parameter, including the dead zone.
 *
 * This function computes the dimensionless turbulent alpha viscosity parameter
 * ($\alpha$) at a given radial position `r`. The calculation accounts for the
 * presence of a dead zone, where the turbulence is suppressed.
 *
 * @param r The radial distance at which to calculate alpha.
 * @param disk_params Pointer to the `disk_t` structure containing disk parameters.
 * @return The calculated turbulent alpha value.
 */
double calculate_turbulent_alpha(double r, const disk_t *disk_params);


/**
 * @brief Calculates the local kinematic viscosity of the gas.
 *
 * This function computes the turbulent kinematic viscosity (`nu`) at a given
 * radial position `r` within the disk. It typically depends on the local sound speed,
 * scale height, and an alpha viscosity parameter.
 *
 * @param r The radial distance at which to calculate viscosity.
 * @param disk_params Pointer to the `disk_t` structure containing disk parameters.
 * @return The calculated kinematic viscosity at radius `r`.
 */
double calculate_gas_viscosity(double r, const disk_t *disk_params);

/**
 * @brief Calculates the local gas disk scale height.
 *
 * This function computes the vertical scale height `H` of the gas disk
 * at a given radial position `r`. This quantity represents the vertical extent of the disk.
 *
 * @param r The radial distance at which to calculate the scale height.
 * @param disk_params Pointer to the `disk_t` structure containing relevant disk properties.
 * @return The calculated scale height at radius `r`.
 */
double calculate_scale_height(double r, const disk_t *disk_params);

/**
 * @brief Calculates the local Keplerian orbital velocity.
 *
 * This function computes the orbital speed (`v_K`) that a test particle
 * would have if it were orbiting in a circular path under the sole influence
 * of the central star's gravity (Keplerian motion).
 * The formula used is: $v_K = \sqrt{G M_{star} / r}$.
 *
 * @param r The radial distance at which to calculate the velocity.
 * @param disk_params Pointer to a structure containing disk parameters, specifically `STAR_MASS`.
 * @return The Keplerian orbital velocity in consistent units.
 */
double calculate_keplerian_velocity(double r, const disk_t *disk_params);

/**
 * @brief Calculates the local Keplerian angular velocity.
 *
 * This function computes the angular velocity ($\omega_K$) that a test particle
 * would have if it were orbiting in a circular path under the sole influence
 * of the central star's gravity (Keplerian motion).
 * The formula used is: $\omega_K = \sqrt{G M_{star} / r^3}$.
 *
 * @param r The radial distance at which to calculate the angular velocity.
 * @param disk_params Pointer to a structure containing disk parameters, specifically `STAR_MASS`.
 * @return The Keplerian angular velocity in consistent units (e.g., radians per second).
 */
double calculate_keplerian_angular_velocity(double r, const disk_t *disk_params);

/**
 * @brief Calculates the local isothermal sound speed of the gas.
 *
 * This function computes the sound speed (`c_s`) at a given radial position `r`.
 * In an isothermal disk, the sound speed is directly related to the Keplerian
 * angular velocity and the disk scale height.
 *
 * @param r The radial distance at which to calculate the sound speed.
 * @param disk_params Pointer to the `disk_t` structure containing relevant disk properties.
 * @return The calculated local sound speed at radius `r`.
 */
double calculate_local_sound_speed(double r, const disk_t *disk_params);

/**
 * @brief Calculates the gas density in the midplane of the disk.
 *
 * This function computes the midplane gas density (`rho_gas`) assuming a
 * vertically isothermal and Gaussian density distribution. It is derived from
 * the surface density and the local scale height.
 * The formula used is: $\rho_{midplane} = \frac{1}{\sqrt{2\pi}} \frac{\Sigma}{H}$.
 *
 * @param sigma The gas surface density at radius `r`.
 * @param r The radial distance.
 * @param disk_params Pointer to the `disk_t` structure containing relevant disk properties.
 * @return The calculated midplane gas density.
 */
double calculate_midplane_gas_density(double sigma, double r, const disk_t *disk_params);

/**
 * @brief Calculates the local gas pressure.
 *
 * This function computes the local gas pressure `p` using the ideal gas law,
 * relating midplane density and sound speed: $p = \rho_{gas} c_s^2$.
 *
 * @param sigma The gas surface density at radius `r`.
 * @param r The radial distance.
 * @param disk_params Pointer to the `disk_t` structure containing relevant disk properties.
 * @return The calculated local gas pressure.
 */
double calculate_gas_pressure(double sigma, double r, const disk_t *disk_params);

/**
 * @brief Calculates the radial gas pressure gradient.
 *
 * This function computes the radial pressure gradient (`dpressvec`) for each
 * grid cell, typically using a finite difference approximation.
 * The result is stored internally in `disk_params->dpressvec`.
 *
 * @param disk_params Pointer to the `disk_t` structure containing the `pressvec` and grid parameters.
 */
void calculate_gas_pressure_gradient(disk_t *disk_params);

/**
 * @brief Calculates a coefficient used in the gas radial velocity calculation.
 *
 * This function computes a specific coefficient, typically of the form $-3 / (\Sigma \sqrt{R})$,
 * which is a part of the expression for the gas radial velocity in viscous accretion disks.
 *
 * @param sigma The local gas surface density.
 * @param r The radial distance.
 * @return The calculated coefficient. Returns 0.0 if `sigma` is zero or `r` is non-positive.
 */
double coefficient_for_gas_velocity(double sigma, double r);

/**
 * @brief Calculates the radial gas velocity profile.
 *
 * This function computes the radial velocity (`ugvec`) of the gas for each grid cell
 * based on the diffusion equation for viscous accretion disks. It uses numerical
 * differentiation to approximate the terms.
 * The result is stored internally in `disk_params->ugvec`.
 *
 * @param disk_params Pointer to the `disk_t` structure containing gas properties
 * (`sigmavec`, `rvec`, `ugvec`) and grid parameters (`NGRID`, `DD`).
 */
void calculate_gas_velocity(disk_t *disk_params);

/**
 * @brief Updates the gas surface density, pressure, and pressure gradient for the next time step.
 *
 * This is a core time-stepping function responsible for evolving the gas disk's
 * properties. It calculates the new surface density based on a numerical scheme
 * (e.g., a finite difference solution to the diffusion equation), and then updates
 * the pressure and pressure gradient accordingly. Boundary conditions are applied
 * after the updates.
 *
 * @param sim_opts Pointer to the `simulation_options_t` structure containing simulation time step (`DT`).
 * @param disk_params Pointer to the `disk_t` structure containing all gas disk properties
 * (e.g., `sigmavec`, `pressvec`, `dpressvec`, `rvec`, `NGRID`) which will be modified.
 */
void get_gas_surface_density_pressure_pressure_gradient(const simulation_options_t *sim_opts, disk_t *disk_params);

#endif // DISK_MODEL_H