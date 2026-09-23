#ifndef PHOTOEVAPORATION_H
#define PHOTOEVAPORATION_H

#include "simulation_types.h"

/* Main dispatcher called from dust_physics.c */
void computePhotoevaporationSink(DiskParameters *disk);

/* Hole / gap detection */
void PhotoevapSearchHole(const double *radius_array,
                          const double *sigma_array,
                          DiskParameters *disk);

/* Owen 2012 */
double PhotoevapNormOwen2012(const double *radius_array,
                               int hole,
                               double r_hole,
                               const double *dr_array,
                               const DiskParameters *disk);

double PhotoevapFuncOwen2012(double r,
                               double norm,
                               int hole,
                               double r_hole,
                               const DiskParameters *disk);

void PhotoevaporationOwen2012(double *evap_array,
                               const double *radius_array,
                               double norm,
                               int hole,
                               double r_hole,
                               const double *dr_array,
                               const DiskParameters *disk);

/* Picogna 2019 */
double PhotoevapFuncPicogna2019(double r, double L_x, int hole, double r_hole);

void PhotoevaporationPicogna2019(double *evap_array,
                                  const double *radius_array,
                                  double L_x,
                                  const double *dr_array,
                                  const DiskParameters *disk);

#endif
