#ifndef DISK_THERMAL_H
#define DISK_THERMAL_H

#include "simulation_types.h"

/* Internal time unit: one year / (2*pi), in seconds. */
double thermalTimeUnitSeconds(void);
double interpolateThermalSoundSpeed(double radius, const DiskParameters *disk);
double interpolateThermalScaleHeight(double radius, const DiskParameters *disk);
void refreshRadiativeDisk(DiskParameters *disk, SimulationOptions *options);
void evolveRadiativeEnergy(DiskParameters *disk, SimulationOptions *options);
/* Outward r*Sigma*v_r at the face between cells j and j+1, in code units. */
double radiativeMassFlux(int j, const DiskParameters *disk, const SimulationOptions *options);
double totalDiskInternalEnergy(const DiskParameters *disk);
void writeEnergyBudget(const DiskParameters *disk, const SimulationOptions *options,
                       double time_years, double output_time);
void writeThermalSnapshot(const DiskParameters *disk, const SimulationOptions *options,
                          double time_years, double output_time);

#endif
