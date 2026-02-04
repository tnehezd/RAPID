#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "disk_model.h"      
#include "simulation_types.h"

void integrateParticleRungeKutta4(double time, double prad, const double *sigmad, const double *particle_distance_grid, double step, double y, double *ynew, double *pradnew, const DiskParameters *disk_params, const SimulationOptions *sim_opts);

#endif // INTEGRATOR_H