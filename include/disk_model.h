// src/disk_model.h

#ifndef DISK_MODEL_H
#define DISK_MODEL_H



#include "config.h"
#include "simulation_types.h" 

void readDiskParameters(disk_t *disk_params);
void createRadialGrid(disk_t *disk_params);
void createInitialGasSurfaceDensity(disk_t *disk_params);
void createInitialGasPressure(disk_t *disk_params);
void createInitialGasPressureGradient(disk_t *disk_params);
void createInitialGasVelocity(disk_t *disk_params);
void calculateInitialDustSurfaceDensity(double radin[][2], double *massin, double out[][3], int n, const disk_t *disk_params);

#endif // DISK_MODEL_H
