/**
 * @file ring_viscosity.h
 * @brief Declarations for the viscous ring spreading test routine.
 */

#ifndef RING_VISCOSITY_H
#define RING_VISCOSITY_H

#include "simulation_types.h"

/**
 * @brief Runs the viscous ring spreading test (Lynden-Bell & Pringle benchmark).
 * 
 * @param disk_params Pointer to DiskParameters.
 * @param sim_opts Pointer to SimulationOptions.
 */
void runRingViscosityTest(DiskParameters *disk_params, SimulationOptions *sim_opts);

#endif // RING_VISCOSITY_H