/**
 * @file steady_state_test.h
 * @brief Declaration for the steady-state accretion disk benchmark test.
 */

#ifndef STEADY_STATE_TEST_H
#define STEADY_STATE_TEST_H

#include "simulation_types.h"

/**
 * @brief Initializes the disk parameters to a steady-state accretion profile.
 * 
 * @param dp Pointer to DiskParameters structure to be initialized.
 * @param opt Pointer to SimulationOptions structure for simulation settings.
 * @param mdot0 Target steady-state mass accretion rate [M_sun/yr].
 */
void initializeSteadyStateProfile(DiskParameters *dp,
                                  SimulationOptions *opt,
                                  double mdot0);    


/**
 * @brief Refreshes the steady-state profile without evolving the surface density.
 * 
 * @param sim_opts Pointer to SimulationOptions structure.
 * @param dp Pointer to DiskParameters structure.
 */
void refreshSteadyState(SimulationOptions *sim_opts, DiskParameters *dp);

/**
 * @brief Runs the steady-state accretion disk test (constant mass accretion rate).
 * 
 * @param disk_params Pointer to DiskParameters.
 * @param sim_opts Pointer to SimulationOptions.
 */
void runSteadyStateTest(DiskParameters *disk_params, SimulationOptions *sim_opts);

#endif // STEADY_STATE_TEST_H