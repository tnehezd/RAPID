#ifndef DUST_PHYSICS_H
#define DUST_PHYSICS_H

#include <stdio.h> 
#include "simulation_types.h"
#include "particle_data.h"

/**
 * @brief Computes the Stokes number of a dust particle.
 *
 * The Stokes number characterizes the aerodynamic coupling between a dust grain
 * and the surrounding gas. It depends on the particle size, local gas surface
 * density, and disk physical parameters.
 *
 * @param particle_radius       Radius of the dust particle (cm or AU, depending on code units).
 * @param gas_surfacedensity    Local gas surface density at the particle position.
 * @param disk_params           Pointer to disk parameter structure.
 * @return The Stokes number at the given location.
 */
double calculateStokesNumber(double particle_radius, double gas_surfacedensity, DiskParameters *disk_params);

/**
 * @brief Computes particle masses for the two-population dust model.
 *
 * This function fills the output arrays with the representative masses
 * of the inner/outer populations based on the particle distribution.
 *
 * @param number_of_particles Total number of dust particles.
 * @param partmassind         2D array storing intermediate mass values.
 * @param indii, indio, indoi, indoo  Index boundaries for the four mass bins.
 * @param massiout            Output pointer for inner population mass.
 * @param massoout            Output pointer for outer population mass.
 * @param sim_opts            Pointer to simulation options.
 */
void calculateParticleMass(int number_of_particles, double (*partmassind)[5], int indii, int indio, int indoi, int indoo, double *massiout, double *massoout, const SimulationOptions *sim_opts);

/**
 * @brief Computes the radial drift barrier size.
 *
 * The drift barrier is the maximum grain size limited by radial drift
 * before significant growth can occur.
 *
 * @param dust_surfacedensity Local dust surface density.
 * @param radial_distance     Radial position of the particle.
 * @param p                   Gas pressure at the particle location.
 * @param dp                  Radial pressure gradient.
 * @param rho_p               Intrinsic particle density.
 * @param disk_params         Pointer to disk parameter structure.
 * @return Maximum particle size allowed by the drift barrier.
 */
double calculateRadialDriftBarrier(double dust_surfacedensity, double radial_distance, double p, double dp, double rho_p, const DiskParameters *disk_params);

/**
 * @brief Computes the turbulent fragmentation barrier size.
 *
 * The fragmentation barrier is the maximum grain size set by destructive
 * collisions driven by turbulent relative velocities.
 *
 * @param gas_surfacedensity  Local gas surface density.
 * @param r                   Radial position.
 * @param rho_p               Intrinsic particle density.
 * @param disk_params         Pointer to disk parameter structure.
 * @return Maximum particle size allowed by turbulent fragmentation.
 */
double calculateTurbulentFragmentationBarrier(double gas_surfacedensity, double r, double rho_p, const DiskParameters *disk_params);

/**
 * @brief Computes the drift-induced fragmentation barrier size.
 *
 * This barrier arises when drift-driven relative velocities exceed the
 * fragmentation threshold.
 *
 * @param gas_surfacedensity  Local gas surface density.
 * @param r                   Radial position.
 * @param p                   Gas pressure.
 * @param dp                  Pressure gradient.
 * @param rho_p               Intrinsic particle density.
 * @param disk_params         Pointer to disk parameter structure.
 * @return Maximum particle size allowed by drift-induced fragmentation.
 */
double calculateDriftInducedFragmentationBarrier(double gas_surfacedensity, double r, double p, double dp, double rho_p, const DiskParameters *disk_params);

/**
 * @brief Computes the dust growth timescale.
 *
 * The growth timescale determines how quickly dust grains grow due to
 * coagulation processes at a given radial location.
 *
 * @param r       Radial position.
 * @param eps     Dust-to-gas ratio.
 * @param disk_params Pointer to disk parameter structure.
 * @return Growth timescale at the given location.
 */
double calculateGrowthTimescale(double r, double eps, const DiskParameters *disk_params);

/**
 * @brief Computes the dust particle size using the Birnstiel et al. model.
 *
 * This function evaluates the maximum allowed particle size based on
 * fragmentation, drift, and growth limits.
 *
 * @param prad                  Current particle radius.
 * @param pdens                 Intrinsic particle density.
 * @param gas_surfacedensity    Local gas surface density.
 * @param dust_surfacedensity   Local dust surface density.
 * @param y                     Radial position of the particle.
 * @param p                     Gas pressure.
 * @param dpress_val            Pressure gradient.
 * @param actual_timestep       Time step.
 * @param disk_params           Pointer to disk parameter structure.
 * @return Updated particle radius after growth.
 */
double calculateDustParticleSize(double prad, double pdens, double gas_surfacedensity, double dust_surfacedensity, double y, double p, double dpress_val, double actual_timestep, const DiskParameters *disk_params);

/**
 * @brief Computes the dust surface density profile from particle data.
 *
 * This function reconstructs the dust surface density on the gas grid
 * based on the positions and masses of dust particles.
 *
 * @param max_param     Upper radial bound.
 * @param min_param     Lower radial bound.
 * @param particle_data Pointer to particle data structure.
 * @param sim_opts      Pointer to simulation options.
 * @param disk_params   Pointer to disk parameter structure.
 */
void calculateDustSurfaceDensity(double max_param, double min_param, const ParticleData *particle_data, const SimulationOptions *sim_opts, const DiskParameters *disk_params);

/**
 * @brief Updates and stores the new radial positions of dust particles.
 *
 * This function integrates the dust particle motion and writes the updated
 * positions to the particle data structure.
 *
 * @param nev                 Output filename or identifier.
 * @param opt       Output mode selector.
 * @param particle_data Pointer to particle data structure.
 * @param actual_timestep     Time step.
 * @param actual_time         Current simulation time.
 * @param n         Number of particles.
 * @param sim_opts  Pointer to simulation options.
 * @param disk_params Pointer to disk parameter structure.
 */
void calculateDustDistance(const char *nev, int opt, ParticleData *particle_data, double actual_timestep, double actual_time, int n, const SimulationOptions *sim_opts, const DiskParameters *disk_params);



#endif // DUST_PHYSICS_H