#ifndef SIMULATION_CORE_H
#define SIMULATION_CORE_H

#include "disk_model.h"        
#include "simulation_types.h"  

void calculate1DDustDrift(double particle_radius, double pressure_gradient, double gas_surface_density, double gas_velocity, double radial_distance, double *drift_velocity, const DiskParameters *disk_params);
double calculateTimeStep(const DiskParameters *disk_params);
void timeIntegrationForTheSystem(SnapshotMode mode, DiskParameters *disk_params, const SimulationOptions *sim_opts, OutputFiles *output_files);

#endif // SIMULATION_CORE_H