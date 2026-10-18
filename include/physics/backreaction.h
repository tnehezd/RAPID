#ifndef BACKREACTION_H
#define BACKREACTION_H

#include "simulation_types.h"
#include "particle_data.h"

typedef struct BackreactionState {
    int cells;
    double *loading_a, *loading_b;
    double *radial, *azimuthal, *flux, *viscous;
    double values[];
} BackreactionState;

/* Local stationary Keplerian drag balance; azimuthal velocity is relative to Keplerian. */
void solveStationaryDrag(double a, double b, double viscous_velocity,
                         double pressure_velocity, double *radial, double *azimuthal);
double backreactionCellArea(int i, const DiskParameters *disk);
void refreshBackreaction(const ParticleData *particles, int count,
                         DiskParameters *disk, SimulationOptions *options);
double backreactionDustVelocity(double radius, double size, double sigma,
                                const DiskParameters *disk);
double backreactionTimeStep(const ParticleData *particles, int count,
                            const DiskParameters *disk, const SimulationOptions *options);
void advanceBackreactionGas(double dt, DiskParameters *disk, SimulationOptions *options);

#endif
