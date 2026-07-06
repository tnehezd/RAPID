// include/simulation_io.h
#ifndef SIMULATION_IO_H
#define SIMULATION_IO_H

#include "simulation_types.h"

void printStatus(int step,
                      double deltat,
                      double current_time_years,
                      double internal_time,
                      double output_time,
                      const char *mode,
                      int was_snapshot,
                      double current_mass,
                      double target_mass,
                      double initial_mass, 
                      double last_snapshot_time, 
                      double interval,
                      SimulationOptions *sim_opts);

#endif // SIMULATION_IO_H