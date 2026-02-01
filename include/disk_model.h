// src/disk_model.h

#ifndef DISK_MODEL_H
#define DISK_MODEL_H

#include "config.h"
#include "simulation_types.h" 

void readDiskParameters(DiskParameters *disk_params);
void createRadialGrid(DiskParameters *disk_params);
void createInitialGasSurfaceDensity(DiskParameters *disk_params);
void createInitialGasPressure(DiskParameters *disk_params);
void createInitialGasPressureGradient(DiskParameters *disk_params);
void createInitialGasVelocity(DiskParameters *disk_params);
void calculateDustSurfaceDensityFromRepresentativeMass(double input_dust_radii_array[][2], double *input_mass_array, double output_dust_surfacedensity_array[][3], int particle_number, const DiskParameters *disk_params);

#endif // DISK_MODEL_H
