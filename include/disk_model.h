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
void calculateDustSurfaceDensityFromRepresentativeMass(double radin[][2], double *massin, double out[][3], int n, const DiskParameters *disk_params);

#endif // DISK_MODEL_H
