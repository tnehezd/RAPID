#ifndef UTILS_H
#define UTILS_H

#include "particle_data.h"
#include "simulation_types.h" 

void linearInterpolation(double *invec, double *radial_grid, double pos, double *out, double rd, const DiskParameters *disk_params);
double findMinimumOfAnArray(double s1, double s2, double s3);
double ftcsSecondDerivativeCoefficient(double radial_distance, const DiskParameters *disk_params);
double ftcsFirstDerivativeCoefficient(double radial_distance, const DiskParameters *disk_params);
int identifyPressureTraps(const DiskParameters *disk_params, PressureTrap *traps, int max_traps);
int countZeroPoints(const DiskParameters *disk_params);
double findZeroPointRadius(double r1, double r2, double dp1, double dp2);
double findZeroPoint(int i, const double *radial_grid, const double *dp);
void findRAnnulusAroundDZE(double rin_val, double *ind_ii, double *ind_io, double rout_val, double *ind_oi, double *ind_oo, const SimulationOptions *sim_opts, const DiskParameters *disk_params);
void sortAnArray(double *rv,int n);
void histogram(double r, int *hist, double dd, DiskParameters *disk_params);
void sortAnArrayarray(double rv[][3],int n);
void roundParticleRadii(double in[][3], int n, const DiskParameters *disk_params);
void mergeParticlesByRadius(double in[][3], double dd, int n, const DiskParameters *disk_params); 
void updateParticleGridIndices(const ParticleData *particle_data, double t, int n, const DiskParameters *disk_params);
void computeParticleRadiusRange(
    const ParticleData *particle_data,
    int particle_number,
    int has_secondary_population,
    double *min_radius,
    double *max_radius
);

#endif // UTILS_H