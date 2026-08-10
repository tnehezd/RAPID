#ifndef RING_TEST_H
#define RING_TEST_H

#include "simulation_types.h"

void setupViscousRingTestCustom(
    SimulationOptions *sim,
    DiskParameters *disk,
    double ring_center,
    double ring_width,
    double ring_mass
);

#endif
